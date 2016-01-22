__author__ = 'etseng@pacificbiosciences.com'

import os, sys, subprocess, time
import logging
from cPickle import *
from multiprocessing.pool import ThreadPool
from pbcore.io.FastaIO import FastaReader, FastaWriter
from pbcore.io.FastqIO import FastqReader
from pbcore.util.Process import backticks

class DazzIDHandler:
    """
    For any kind of input fasta, convert & maintain ID mapping

    <input>.fasta --> <input>.dazz.fasta, <input>.dazz.fasta.pickle

    If the input is fastq, automatically output as fasta since dalign only takes fasta

    By default, converted is False, so go ahead with converting.
    Otherwise, just read the mapping pickle.

    NOTE: Consistent with dazz's indexing, the IDs will be 1-based!!!!
    """
    def __init__(self, input_filename, converted=False):

        i = input_filename.rfind('.')

        self.input_filename = input_filename
        self.file_prefix = input_filename[:i]
        self.file_suffix = input_filename[i+1:].lower()

        self.dazz_movie_name = 'prolog'
        self.dazz_filename = self.file_prefix + '.dazz.fasta'
        self.dazz_mapping = {} # index --> original sequence ID   ex: 1 --> movie/zmw/start_end_CCS


        if self.file_suffix in ('fa', 'fasta'):
            self.filetype = 'fasta'
        elif self.file_suffix in ('fq', 'fastq'):
            self.filetype = 'fastq'
        else:
            raise Exception, "Unrecognized file suffix! Must be .fa, .fasta, .fq or .fastq!"

        if not converted:
            self.convert_to_dazz_fasta()
        else:
            self.read_dazz_pickle()


    def convert_to_dazz_fasta(self):
        """
        Convert input fasta/fastq file to daligner-compatibe fasta with ids:
        <prefix>/<index>/0_<seqlen>

        Also write out mappings to pickle
        """
        i = 1
        reader = FastaReader(self.input_filename) if self.filetype == 'fasta' else \
            FastqReader(self.input_filename)

        f = FastaWriter(self.dazz_filename)

        for r in reader:
            f.writeRecord("{p}/{i}/0_{len}".format(p=self.dazz_movie_name, i=i, len=len(r.sequence)), r.sequence)
            self.dazz_mapping[i] = r.id
            i += 1

        f.close()

        with open(self.dazz_filename + '.pickle', 'w') as f:
            dump(self.dazz_mapping, f)


    def read_dazz_pickle(self):
        filename = self.dazz_filename + '.pickle'
        with open(filename) as f:
            self.dazz_mapping = load(f)


    def keys(self):
        return self.dazz_mapping.keys()

    def __getitem__(self, key):
        """
        key should be a single integer from the id
        ex: prolog/1/0_1324 --> key should be 1
        """
        return self.dazz_mapping[key]


class DalignerRunner:
    def __init__(self, query_filename, db_filename, is_FL, same_strand_only, query_converted=False, db_converted=False, query_made=False, db_made=False, use_sge=True, sge_opts=None, script_dir="scripts/", cpus=24):
        self.query_filename = query_filename
        self.db_filename = db_filename
        self.is_FL = is_FL
        self.same_strand_only = same_strand_only

        self.cpus = cpus

        # always assume they are converted
        self.query_dazz_handler = DazzIDHandler(query_filename, converted=query_converted)
        # DB may have already been converted (if shared)
        self.db_dazz_handler = DazzIDHandler(db_filename, converted=db_converted)
        if not db_made:
            DalignerRunner.make_db(self.db_dazz_handler.dazz_filename)
        if not query_made:
            DalignerRunner.make_db(self.query_dazz_handler.dazz_filename)
        self.db_blocks = self.get_num_blocks(self.db_dazz_handler.dazz_filename)
        self.query_blocks = self.get_num_blocks(self.query_dazz_handler.dazz_filename)

        self.use_sge = use_sge
        self.sge_opts = sge_opts
        self.script_dir = os.path.realpath(script_dir)

        if not os.path.exists(self.script_dir):
            os.makedirs(self.script_dir)

    @staticmethod
    def make_db(filename):
        """
        1. fasta2DB
        2. DBsplit
        3. get & store number of blocks
        """
        db_filename = filename + '.db'
        if os.path.exists(db_filename):
            subprocess.check_call("DBrm " + filename, shell=True)

        cmd = "fasta2DB {0} {0}".format(filename)
        subprocess.check_call(cmd, shell=True)

        cmd = "DBsplit -s200 {0}".format(filename)
        subprocess.check_call(cmd, shell=True)

    def get_num_blocks(self, filename):
        db_filename = filename + '.db'
        with open(db_filename) as f:
            f.readline()
            f.readline()
            x = f.readline().strip()
            assert x.startswith('blocks =')
            num_blocks = int(x.split('=')[1])
        return num_blocks

    def runHPC(self, min_match_len=300, output_dir='.', sensitive_mode=False):
        """
        ex: daligner -h35 -e.80 -l500 -s100 <query> <db>

        if self.use_sge --- writes to <scripts>/daligner_job_#.sh
        else --- run locally, dividing into self.cpus/4 tasks (capped max at 4)

        NOTE 1: when using SGE, be careful that multiple calls to this might end up writing to the SAME job.sh files,
        this should be avoided by changing <scripts> directory

        NOTE 2: more commonly this should be invoked locally (since ice_partial.py i/one be qsub-ed),
        in that case it is more recommended to keep self.cpus = 4 so that each daligner job is run consecutively
        and that the original qsub job should have been called with qsub -pe smp 4 (set by --blasr_nproc 4)
        In this way, the daligner jobs are called consecutively, but LA4Ice (LAshow) is parallelized 4X
        """
        old_dir = os.path.realpath(os.path.curdir)
        os.chdir(output_dir)
        cmds_daligner = []
        cmds_show = []
        las_filenames = []
        las_out_filenames = []

        for i in xrange(self.query_blocks):
            for j in xrange(self.db_blocks):
                # avoid unnecessary redundant calls to daligner when doing self hits
                if self.query_filename == self.db_filename and i > j:
                    continue

                # old DALIGNER param is not sensitive enough for > 5 kb CCS reads
                # but the more sensitive param seems to do badly on 1 - 2 kb reads. WTH =_=
                if not sensitive_mode:
                    cmd = "timeout 600 daligner -w12 -h24 -k24 -e.80 -l{m} -s100 -t10 {q}.{i} {db}.{j}"
                    #cmd = "timeout 600 daligner -h35 -k16 -e.80 -l{m} -s100 -t10 {q}.{i} {db}.{j}"
                else:
                    cmd = "timeout 1200 daligner -w12 -h16 -k12 -e.70 -l{m} -s100 -t10 {q}.{i} {db}.{j}"
                    #cmd = "timeout 600 daligner -w12 -h24 -k24 -e.70 -l{m} -s100 -t10 {q}.{i} {db}.{j}"

                cmd = cmd.format(q=self.query_dazz_handler.dazz_filename, i=i+1, \
                    db=self.db_dazz_handler.dazz_filename, j=j+1, \
                    m=min_match_len)
                cmds_daligner.append(cmd)

                cmd = "LA4Ice"
                if self.is_FL: cmd += " -E"
                for k in xrange(4): # NTHREADS --- hardcoded!!!
                    p = "{q}.{i}.{db}.{j}.N{k}.las".format(\
                        q=self.query_dazz_handler.dazz_filename, i=i+1, \
                        db=os.path.basename(self.db_dazz_handler.dazz_filename), j=j+1, k=k)
                    cmd_ice = cmd + " -a -m -i0 -w100000 -b0 {dbq} {db} {p} > {p}.out".format(\
                        p=p, db=self.db_dazz_handler.dazz_filename,\
                        dbq=self.query_dazz_handler.dazz_filename)
                    cmds_show.append(cmd_ice)
                    las_filenames.append(p)
                    las_out_filenames.append(p + '.out')

                    if not self.same_strand_only:
                        p = "{q}.{i}.{db}.{j}.C{k}.las".format(\
                            q=self.query_dazz_handler.dazz_filename, i=i+1, \
                            db=os.path.basename(self.db_dazz_handler.dazz_filename), j=j+1, k=k)
                        cmd_ice = cmd + " -a -m -i0 -w100000 -b0 {dbq} {db} {p} > {p}.out".format(\
                            p=p, db=self.db_dazz_handler.dazz_filename,\
                            dbq=self.query_dazz_handler.dazz_filename)
                        cmds_show.append(cmd_ice)
                        las_filenames.append(p)
                        las_out_filenames.append(p + '.out')

        assert len(cmds_show) > 0

        logging.info("CMD: " + "\n".join(cmds_daligner) + '\n')

        # (a) run all daligner jobs
        start_t = time.time()
        done_script = os.path.join(self.script_dir, "daligner_job_wait.sh")
        with open(done_script, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("touch DALIGN.DONE\n")
        if self.use_sge:
            qsub_job_runner(cmds_daligner, os.path.join(self.script_dir, "daligner_job_{i}.sh"), done_script, self.sge_opts)
        else:
            local_job_runner(cmds_daligner, num_processes=max(1, min(self.cpus/4, 4))) # max 4 at a time to avoid running out of mem..
        logging.info("daligner jobs took {0} sec.".format(time.time()-start_t))


        #print >> sys.stderr, "\n".join(cmds_show)

        # (b) run all LA4Ice jobs
        start_t = time.time()
        #print cmds_show
        if self.use_sge:
            qsub_job_runner(cmds_show, os.path.join(self.script_dir, "LA4Ice_job_{i}.sh"), done_script, self.sge_opts)
        else:
            local_job_runner(cmds_show, num_processes=max(1, min(self.cpus, 4))) # max 4 at a time to avoid running out of memory...
        #print >> sys.stderr, "LA4Ice jobs took {0} sec".format(time.time()-start_t)
        os.chdir(old_dir)
        return las_filenames, las_out_filenames


def local_job_runner(cmds_list, num_processes):
    run_cmd_in_shell = lambda x: subprocess.check_call(x, shell=True)
    pool = ThreadPool(processes=num_processes)
    rets = pool.map(run_cmd_in_shell, cmds_list)
    pool.close()
    pool.join()
    for i, cmd in enumerate(cmds_list):
        if rets[i] != 0:
            raise RuntimeError, "CMD failed:", cmd


def qsub_job_runner(cmds_list, sh_file_format, done_script, sge_opts, qsub_retry=3, run_local_if_qsub_fail=True):
    """
    cmds_list -- list of commands to run (each in a separate file)
    sh_file_format ---- ex: test_script.{i}.sh

    ToDo:
    (1) add in ways to gracefully fail if SGE submits fail -- resubmit? wait? run local?
    (2) add in ways to monitor if certain qsub jobs died or hung --- resubmit? kill? run local?
    """
    jids = []
    for i, cmd in enumerate(cmds_list):
        f = open(sh_file_format.format(i=i), 'w')
        f.write("#!/bin/bash\n")
        f.write(cmd + '\n')
        f.close()

        # hard-coded to 4 CPUS because hard-coded in daligner!
        qsub_cmd = "qsub"
        if sge_opts.queue_name is not None:
            qsub_cmd += " -q " + sge_opts.queue_name
        qsub_cmd += " -cwd -V -S /bin/bash -pe {env} 4 -e {out}.elog -o {out}.olog {out}".format(\
            env=sge_opts.sge_env_name, out=f.name)
        try_times = 1
        while try_times <= qsub_retry:
            _out, _code, _msg = backticks(qsub_cmd)
            if _code == 0: # succeeded, break
                break
            else:
                # failed, sleep for a little, try again
                time.sleep(10)
                try_times += 1
        if try_times > qsub_retry:
            if run_local_if_qsub_fail:
                raise NotImplementedError, "Not yet implemented to not use SGE!"
            else:
                raise RuntimeError, "Unable to qsub. Abort!:", qsub_cmd
        # ex: # Your job 596028 ("a.sh") has been submitted
        jids.append(str(_out).split()[2])

    # use a qsub job to wait for the commands to finish
    # ToDo: this is NOT bullet proof! watch for cases where the job may have died or been killed or hung
    wait_cmd = "qsub "
    if sge_opts.queue_name is not None:
        wait_cmd += " -q " + sge_opts.queue_name
    wait_cmd += " -sync y -pe {2} 1 -cwd -S /bin/bash -V -e /dev/null -o /dev/null -hold_jid {0} {1}".format(",".join(jids), done_script, sge_opts.sge_env_name)
    _out, _code, _msg = backticks(wait_cmd)
    if _code != 0: # failed, just wait manually then
        active_jids = [x.split()[0] for x in os.popen("qstat").read().strip().split('\n')[2:]]
        while True:
            if any(x in jids for x in active_jids): # some jobs are still running
                time.sleep(10)
            else:
                break


