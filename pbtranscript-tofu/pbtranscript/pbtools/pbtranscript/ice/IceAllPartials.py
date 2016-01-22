#!/usr/bin/env python
"""
IceAllPartials procedures:

(1) For each input fasta file (i.e., each splitted nfl read file), map its
    reads to unpolished consensus isoforms in ref_fasta
    (i.e., final.consensus.fa), then create a pickle file
    (e.g. *.partial_uc.pickle).

(2) Wait for all pickle files to be created

(3) Merge all pickle files and dump to a big pickle.

"""

import os.path as op
import logging
import time
from pbtools.pbtranscript.PBTranscriptOptions import \
    add_sge_arguments, add_fofn_arguments
from pbtools.pbtranscript.Utils import realpath, mkdir, real_upath, ln
from pbtools.pbtranscript.ice.IceFiles import IceFiles
from pbtools.pbtranscript.ice.IceUtils import combine_nfl_pickles
from pbtools.pbtranscript.ice.IceUtils import ice_fq2fa


class IceAllPartials(IceFiles):

    """Align non-full-length reads to unpolished consensus
    isoforms and assign them to isoforms based on blasr alignments,
    (build partial uc) and save results to a pickle file.
    """

    # Define description of IceAllPartials
    desc = "Align and assign all non-full-length reads to unpolished " + \
           "consensus isoforms, and save results to a pickle file."

    prog = "ice_partial.py all " # used by ice_partial.py

    def __init__(self, root_dir, fastq_filenames, ref_fasta,
                 out_pickle, sge_opts, sa_file=None, ccs_fofn=None):
        """
        fastq_filenames --- a list of splitted nfl fastq files.

        ref_fasta --- (unpolished) consensus isoforms

        out_pickle --- a pickle file with all nfl fasta reads

        ccs_fofn --- should be reads_of_insert.fofn or None

        root_dir --- ICE root output directory

        sge_opts --- params for SGE environment, including
            use_sge    : use SGE or not
            max_sge_jobs: maximum number of sub-jobs submitted
            unique_id  : unique qsub job id, important that this
                        DOES NOT CONFLICT!
            blasr_nproc: blasr -nproc param, number of threads per cpu.
        """
        self.prog_name = "IceAllPartials"
        IceFiles.__init__(self, prog_name=self.prog_name, root_dir=root_dir)

        self.add_log("DEBUG: in IceAllPartials, ccs_fofn is {0}.".format(ccs_fofn), level=logging.INFO)

        self.fastq_filenames, self.ref_fasta, self.ccs_fofn, self.sa_file = \
            self._validate_inputs(fastq_filenames=fastq_filenames,
                                  ref_fasta=ref_fasta,
                                  ccs_fofn=ccs_fofn,
                                  sa_file=sa_file)

        self.out_pickle = out_pickle

        self.sge_opts = sge_opts

        self.add_log("Making dir for mapping noFL reads: " + self.nfl_dir)
        mkdir(self.nfl_dir)

        self.add_log("input fasta files are: " +
                     ", ".join(self.fastq_filenames))
        self.add_log("temp pickle files are: " +
                     ", ".join(self.pickle_filenames))
        self.add_log("out pickle file is: " + self.out_pickle)

    def _validate_inputs(self, fastq_filenames, ref_fasta, ccs_fofn, sa_file):
        """Validate input files."""
        for f in fastq_filenames:
            if not op.exists(f):
                raise IOError("Input fastq {f} does not exist.".format(f=f))
        if ref_fasta is None or not op.exists(ref_fasta):
            raise IOError("Reference {r} does not exist.".format(r=ref_fasta))
        if ccs_fofn is not None and not op.exists(ccs_fofn):
            raise IOError("ccs_fofn file {ccs_fofn} does not exist.".format(
                ccs_fofn=ccs_fofn))
        if sa_file is not None and not op.exists(sa_file):
            raise IOError("suffix array {s} does not exist.".format(s=sa_file))
        return ([realpath(f) for f in fastq_filenames],
                realpath(ref_fasta),
                realpath(ccs_fofn),
                realpath(sa_file) if sa_file is not None else None)

    def cmd_str(self):
        """Return a cmd string."""
        return self._cmd_str(fastq_filenames=self.fastq_filenames,
                             ref_fasta=self.ref_fasta,
                             out_pickle=self.out_pickle,
                             root_dir=self.root_dir,
                             ccs_fofn=self.ccs_fofn,
                             sa_file=self.sa_file,
                             sge_opts=self.sge_opts)

    def _cmd_str(self, fastq_filenames, ref_fasta, out_pickle,
                 root_dir, ccs_fofn, sa_file, sge_opts):
        """Return a cmd string."""
        cmd = self.prog + \
              "{fns} ".format(fns=",".join(fastq_filenames)) + \
              "{ref} ".format(ref=ref_fasta) + \
              "{out} ".format(out=out_pickle)
        if root_dir is not None:
            cmd += "--root_dir={d} ".format(d=root_dir)
        if ccs_fofn is not None:
            cmd += "--ccs_fofn={f} ".format(f=ccs_fofn)
        if sa_file is not None:
            cmd += "--sa={sa} ".format(sa=sa_file)
        cmd += sge_opts.cmd_str(show_blasr_nproc=True)

        return cmd

    @property
    def pickle_filenames(self):
        """pickle files for each fasta file."""
        return [op.join(self.nfl_dir, op.basename(f) + ".partial_uc.pickle")
                for f in self.fastq_filenames]

    @property
    def done_filenames(self):
        """done files to indicate that pickles are done."""
        return [op.join(self.nfl_dir, op.basename(f) + ".DONE")
                for f in self.pickle_filenames]

    @property
    def script_filenames(self):
        """scripts to generate pickles from fasta files."""
        return [op.join(self.script_dir, op.basename(f) + ".partial_uc.sh")
                for f in self.fastq_filenames]

    def createPickles(self):
        """For each file in fastq_filenames, call 'ice_partial.py one' to
        build clusters and to save results to a pickle file. When all pickles
        are done, union all pickles.
        """
        self.add_log("Mapping non-full-length reads to consensus isoforms.")
        self.add_log("Creating pickles...", level=logging.INFO)

        # using --blasr_nproc=4 because DALIGNER uses only 4 cores
        for idx, fq in enumerate(self.fastq_filenames):
            # for each splitted non-full-length reads fasta file, build #
            # partial_uc.pickle
            cmd = "ice_partial.py one {i} ".format(i=real_upath(fq)) + \
                  "{r} ".format(r=real_upath(self.ref_fasta)) + \
                  "{o} ".format(o=real_upath(self.pickle_filenames[idx])) + \
                  "--blasr_nproc={n} ".format(n=4) + \
                  "--done={d} ".format(d=real_upath(self.done_filenames[idx]))
            if self.ccs_fofn is not None:
                cmd += "--ccs_fofn={f} ".format(f=real_upath(self.ccs_fofn))
            if self.sa_file is not None:
                cmd += "--sa={sa} ".format(sa=real_upath(self.sa_file))

            self.add_log("Writing command to script {fsh}".
                         format(fsh=self.script_filenames[idx]))
            self.add_log("CMD: {0}".format(cmd))
            with open(self.script_filenames[idx], 'w') as fsh:
                fsh.write(cmd + "\n")

            # determine elog & olog
            partial_log_fn = op.join(self.log_dir,
                                     'IcePartial.{idx}'.format(idx=idx))
            elog = partial_log_fn + ".elog"
            olog = partial_log_fn + ".olog"
            jid = "ice_partial_{unique_id}_{name}".format(
                unique_id=self.sge_opts.unique_id,
                name=op.basename(fq))

            qsub_cmd = "qsub"
            if self.sge_opts.sge_queue is not None:
                qsub_cmd += " -q " + self.sge_opts.sge_queue
            qsub_cmd += " -pe {env} {n} ".format(env=self.sge_opts.sge_env_name, n=4) + \
                       "-cwd -S /bin/bash -V " + \
                       "-e {elog} ".format(elog=real_upath(elog)) + \
                       "-o {olog} ".format(olog=real_upath(olog)) + \
                       "-N {jid} ".format(jid=jid) + \
                       "{sh}".format(sh=real_upath(self.script_filenames[idx]))

            self.add_log("Creating a pickle for {f}".format(f=fq))

            if self.sge_opts.use_sge is True:
                self.qsub_cmd_and_log(qsub_cmd)
            else:
                cmd += " 1>{olog} 2>{elog}".format(olog=real_upath(olog),
                                                   elog=real_upath(elog))
                self.run_cmd_and_log(cmd=cmd, olog=olog, elog=elog)

    def waitForPickles(self, pickle_filenames, done_filenames):
        """Wait for *.pickle and *.pickle.DONE to be created."""
        self.add_log("Waiting for pickles {ps} to be created.".
                     format(ps=", ".join(pickle_filenames)),
                     level=logging.INFO)
        stop = False
        sleep_time = 10
        while stop is not True:
            stop = all(op.exists(p) for p in pickle_filenames) and \
                all(op.exists(d) for d in done_filenames)
            sleep_time = min(60, sleep_time + 10) # wait in increments of 10 sec, up to 1 min
            time.sleep(sleep_time)
            self.add_log("Waiting for pickles to be created: {ps}".
                         format(ps=", ".join([p for p in pickle_filenames
                                              if op.exists(p)])))

    def combinePickles(self, pickle_filenames, out_pickle):
        """Combine all *.pickle files to one and dump to self.out_pickle."""
        combine_nfl_pickles(pickle_filenames, out_pickle)

    def run(self):
        """Assigning nfl reads to consensus isoforms and merge."""
        # Call ice_partial.py to create a pickle for each splitted nfl fasta
        self.createPickles()
        # Wait for pickles to be created, if SGE is used.
        self.waitForPickles(pickle_filenames=self.pickle_filenames,
                            done_filenames=self.done_filenames)
        # Combine all pickles to a big pickle file: nfl_all_pickle_fn.
        self.combinePickles(pickle_filenames=self.pickle_filenames,
                            out_pickle=self.nfl_all_pickle_fn)
        # Create symbolic link if necessary
        ln(self.nfl_all_pickle_fn, self.out_pickle)

        # Close log
        self.close_log()


def add_ice_all_partials_arguments(parser):
    """Add IceAllPartials argument parser."""
    parser.add_argument("fastq_filenames",
                        type=str,
                        help="comma delimited fasta files of " +
                             "splitted non-full-length reads")
    parser.add_argument("ref_fasta",
                        type=str,
                        help="Reference fasta file, most likely " +
                             "ref_consensus.fa from ICE output")
    parser.add_argument("out_pickle",
                        type=str,
                        help="Output pickle file.")
    parser.add_argument("--root_dir",
                        dest="root_dir",
                        default="",
                        help="A directory for saving intermediate files.")

    parser = add_sge_arguments(parser, blasr_nproc=True)
    parser = add_fofn_arguments(parser, ccs_fofn=True)

    parser.add_argument("--sa",
                        dest="sa_file",
                        default=None,
                        help="saffix array of ref_fasta")
    parser.add_argument("--done", dest="done_filename", type=str,
                        help="An empty file generated to indicate that " +
                             "out_pickle is done.")
    return parser

#
# import sys
# import os
# from pbcore.util.ToolRunner import PBToolRunner
# from pbtools.pbtranscript.__init__ import get_version
# from pbtools.pbtranscript.ClusterOptions import SgeOptions
#
# class IceAllPartialsRunner(PBToolRunner):
#
#     """IceAllPartials Runner"""
#
#     def __init__(self):
#         PBToolRunner.__init__(self, IceAllPartials.desc)
#         self.parser = add_ice_all_partials_arguments(self.parser)
#
#     def getVersion(self):
#         """Get version string."""
#         return get_version()
#
#     def run(self):
#         """Run"""
#         logging.info("Running {f} v{v}.".format(f=op.basename(__file__),
#                                                 v=self.getVersion()))
#         cmd_str = ""
#         try:
#             args = self.args
#             sge_opts = SgeOptions(os.getpid(),
#                                   blasr_nproc=args.blasr_nproc,
#                                   use_sge=args.use_sge,
#                                   max_sge_jobs=args.max_sge_jobs)
#
#             icep = IceAllPartials(
#                 root_dir=args.root_dir,
#                 fasta_filenames=args.fasta_filenames.split(','),
#                 ref_fasta=args.ref_fasta,
#                 out_pickle=args.out_pickle,
#                 sge_opts=sge_opts,
#                 sa_file=args.sa_file,
#                 ccs_fofn=args.ccs_fofn)
#
#             cmd_str = obj.cmd_str()
#             icep.run()
#         except:
#             logging.exception("Exiting {cmd_str} with return code 1.".
#                               format(cmd_str=cmd_str))
#             return 1
#         return 0
#
#
# def main():
#     """Main function"""
#     runner = IceAllPartialsRunner()
#     return runner.start()
#
# if __name__ == "__main__":
#     sys.exit(main())
