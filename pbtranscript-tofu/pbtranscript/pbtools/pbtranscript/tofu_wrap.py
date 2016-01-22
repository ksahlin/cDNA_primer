#!/usr/bin/env python
import os, sys, argparse, subprocess, binascii, time
from cPickle import *
from pbcore.io.FastaIO import FastaReader
from pbcore.io.FastqIO import FastqReader, FastqWriter
from pbtools.pbtranscript.Utils import realpath
from pbtools.pbtranscript.Cluster import Cluster
from pbtools.pbtranscript.ClusterOptions import IceOptions, SgeOptions, \
            IceQuiverHQLQOptions
from pbtools.pbtranscript.PBTranscriptOptions import add_cluster_arguments
from pbtools.pbtranscript.ice.IceUtils import convert_fofn_to_fasta
from pbtools.pbtranscript.counting import get_read_count_from_collapsed as sp
from bisect import bisect_right
from collections import defaultdict
from pbtools.pbtranscript.__init__ import get_version

def sep_flnc_by_primer(flnc_filename, root_dir, output_filename='isoseq_flnc.fastq'):
    """
    Separate flnc fastq by primer. Useful for targeted sequencing.
    ex: make <root_dir>/primer0/isoseq_flnc.fastq ... etc ...
    """
    def get_primer(r):
        for x in r.name.split(';'):
            if x.startswith('primer='):
                return x.split('=')[1]
        return None

    primers = set()
    for r in FastqReader(flnc_filename):
        p = get_primer(r)
        if p is None:
            raise Exception, "ERROR: Unable to find primer information from sequence ID for {0}! Abort!".format(r.name)
        primers.add(p)

    handles = {}
    for p in primers:
        dirname = os.path.join(root_dir, "primer{0}".format(p))
        if os.path.exists(dirname):
            print >> sys.stderr, "WARNING: {0} already exists.".format(dirname)
        else:
            os.makedirs(dirname)
        handles[p] = FastqWriter(os.path.join(dirname, output_filename))
        #handles[p] = open(os.path.join(dirname, output_filename), 'w')

    for r in FastqReader(flnc_filename):
        p = get_primer(r)
        handles[p].writeRecord(r.name, r.sequence, r.quality)
        #handles[p].write(">{0}\n{1}\n".format(r.name, r.sequence))

    for f in handles.itervalues(): f.close()
    primers = list(primers)
    primers.sort(key=lambda x: int(x))
    names = [handles[x].file.name for x in primers]
    return filter(lambda x: os.stat(x).st_size > 0, names)


def sep_flnc_by_size(flnc_filename, root_dir, bin_size_kb=1, bin_manual=None, max_base_limit_MB=600, output_filename='isoseq_flnc.fastq'):
    """
    Separate flnc fastq into different size bins
    ex: make <root_dir>/0to2k/isoseq_flnc.fastq ... etc ...

    If <bin_manual> (ex: (0, 2, 4, 12)) is given, <bin_size_kb> is ignored.

    If max_base_limit_MB is not None, it caps the per-partition # of bases (in Mb)
    So could have 0to1k_part1, 0to2k_part2...etc...
    """
    # first check min - max size range
    min_size = 0
    max_size = 0
    base_in_each_size = defaultdict(lambda: 0) # size(in kb)--> number of bases
    for r in FastqReader(flnc_filename):
        seqlen = len(r.sequence)
        min_size = min(min_size, seqlen)
        max_size = max(max_size, seqlen)
        base_in_each_size[seqlen/1000] += len(r.sequence)
    
    min_size_kb = min_size/1000
    max_size_kb = max_size/1000 + (1 if max_size%1000 > 1 else 0)
    if bin_manual is not None:
        if bin_manual[0] > min_size_kb:
            raise Exception, "Min sequence length is {0} kb, below the bin!".format(min_size)
        if bin_manual[-1] < max_size_kb:
            raise Exception, "Max sequence length is {0} kb, above the bin!".format(max_size)
        bins = bin_manual
    else:
        bins = range(min_size_kb, max_size_kb+1, bin_size_kb)

    print >> sys.stderr, bins

    max_bin = len(bins)-1
    # for each bin, see if further partitioning is needed
    base_in_each_bin = defaultdict(lambda: 0)
    num_parts_in_each_bin = defaultdict(lambda: 1)
    part_counter_in_each_bin = defaultdict(lambda: 0)
    if max_base_limit_MB is not None:
        for kb_size, num_bases in base_in_each_size.iteritems():
            i = min(max_bin, max(0, bisect_right(bins, kb_size)-1))
            base_in_each_bin[i] += num_bases
        for bin_i, num_bases in base_in_each_bin.iteritems():
            num_parts_in_each_bin[bin_i] = int((num_bases*1. / 10**6) / max_base_limit_MB) + \
                    (1 if (num_bases*1. / 10**6) % max_base_limit_MB > 0 else 0)

    handles = {}
    for i in xrange(len(bins)-1):
        for part in xrange(num_parts_in_each_bin[i]):
            dirname = os.path.join(root_dir, "{0}to{1}kb_part{2}".format(bins[i], bins[i+1], part))
            if os.path.exists(dirname):
                print >> sys.stderr, "WARNING: {0} already exists.".format(dirname)
            else:
                os.makedirs(dirname)
            handles[(i,part)] = FastqWriter(os.path.join(dirname, output_filename))


    for r in FastqReader(flnc_filename):
        kb_size = len(r.sequence)/1000
        i = min(max_bin, max(0, bisect_right(bins, kb_size)-1))
        part = part_counter_in_each_bin[i] % num_parts_in_each_bin[i]
        part_counter_in_each_bin[i] += 1
        handles[(i,part)].writeRecord(r.name, r.sequence, r.quality)
        print >> sys.stderr, "putting {0} in {1}".format(len(r.sequence), handles[(i,part)].file.name)
    for h in handles.itervalues(): h.close()

    names = []
    for i in xrange(len(bins)-1):
        for part in xrange(num_parts_in_each_bin[i]):
            names.append(handles[(i,part)].file.name)
    return filter(lambda x: os.stat(x).st_size > 0, names)

def combine_quiver_results(split_dirs, output_dir, hq_filename, lq_filename, tofu_prefix=''):
    """
    For each size bin result, ex: clusterOut/0to2k/all.quiveredXXX.fastq
    combine it together, remember to add a prefix (ex: i0|c12, i1|c13/....)
    """
    prefix_dict_hq  = {}
    prefix_dict_lq = {}
    fout_hq = FastqWriter(os.path.join(output_dir, 'all_sizes.quivered_hq.fastq'))
    fout_lq = FastqWriter(os.path.join(output_dir, 'all_sizes.quivered_lq.fastq'))
    for i,d in enumerate(split_dirs):
        file_hq = os.path.join(d, hq_filename) #'all_quivered_hq.100_30_0.99.fastq')
        file_lq = os.path.join(d, lq_filename) #'all_quivered_lq.fastq')
        print >> sys.stderr, "Adding prefix i{0}| to {1},{2}...".format(i, file_hq, file_lq)
        prefix_dict_hq["i{i}HQ_{p}".format(i=i,p=tofu_prefix)] = os.path.abspath(d)
        prefix_dict_lq["i{i}LQ_{p}".format(i=i,p=tofu_prefix)] = os.path.abspath(d)
        for r in FastqReader(file_hq):
            _name_ = "i{i}HQ_{p}|{n}".format(p=tofu_prefix, i=i, n=r.name)
            fout_hq.writeRecord(_name_, r.sequence, r.quality)
        for r in FastqReader(file_lq):
            _name_ = "i{i}LQ_{p}|{n}".format(p=tofu_prefix, i=i, n=r.name)
            fout_lq.writeRecord(_name_, r.sequence, r.quality)
    fout_hq.close()
    fout_lq.close()
    print >> sys.stderr, "HQ quivered output combined to:", fout_hq.file.name
    print >> sys.stderr, "LQ quivered output combined to:", fout_lq.file.name
    return fout_hq.file.name,fout_lq.file.name,prefix_dict_hq,prefix_dict_lq

def run_collapse_sam(fastq_filename, gmap_db_dir, gmap_db_name, cpus=24, min_coverage=0.99, min_identity=0.95, dun_merge_5_shorter=False, max_fuzzy_junction=5):
    """
    Wrapper for running collapse script
    (a) run GMAP
    (b) sort GMAP sam
    (c) run collapse_isoforms_by_sam
    """
    if not os.path.exists(fastq_filename + '.sam'):
        cmd = "gmap -D {d} -d {name} -n 0 -t {cpus} -z sense_force --cross-species -f samse {fq} > {fq}.sam 2> {fq}.sam.log".format(\
                d=gmap_db_dir, name=gmap_db_name, cpus=cpus, fq=fastq_filename)
        print >> sys.stderr, "CMD:", cmd
        subprocess.check_call(cmd, shell=True)
    cmd = "sort -k 3,3 -k 4,4n {fq}.sam > {fq}.sorted.sam".format(fq=fastq_filename)
    print >> sys.stderr, "CMD:", cmd
    subprocess.check_call(cmd, shell=True) 
    cmd = "collapse_isoforms_by_sam.py --input {fq} --fq -s {fq}.sorted.sam --max_fuzzy_junction {j} -c {c} -i {i}".format(\
            c=min_coverage, i=min_identity, 
            fq=fastq_filename, j=max_fuzzy_junction)
    if dun_merge_5_shorter:
        cmd += " --dun-merge-5-shorter -o {fq}.no5merge".format(fq=fastq_filename)
    else:
        cmd += " -o {fq}.5merge".format(fq=fastq_filename)
    print >> sys.stderr, "CMD:", cmd
    subprocess.check_call(cmd, shell=True)

    if dun_merge_5_shorter:
        return fastq_filename + '.no5merge.collapsed'
    else:
        return fastq_filename + '.5merge.collapsed'


def get_abundance(collapse_prefix, prefix_dict, output_prefix, restricted_movies=None):
    cid_info = sp.read_group_filename(collapse_prefix + ".group.txt", is_cid=True,\
            sample_prefixes=prefix_dict.keys())

    fl_pickles = []
    nfl_pickles = []
    for i, d in prefix_dict.iteritems():
        file = os.path.join(d, 'output/final.pickle')
        if not os.path.exists(file):
            raise Exception, "Expected FL pickle file {0} but not found!".format(file)
        fl_pickles.append((i, file))
        file = os.path.join(d, 'output/map_noFL/nfl.all.partial_uc.pickle')
        if not os.path.exists(file):
            raise Exception, "Expected nFL pickle file {0} but not found!".format(file)
        nfl_pickles.append((i, file))

    sp.output_read_count_FL(cid_info, fl_pickles, output_prefix + '.read_stat.txt', 'w', restricted_movies=restricted_movies)
    sp.output_read_count_nFL(cid_info, nfl_pickles, output_prefix + '.read_stat.txt', 'a', restricted_movies=restricted_movies)
    sp.make_abundance_file(output_prefix + '.read_stat.txt', output_prefix + '.abundance.txt', restricted_movies=restricted_movies)
    print >> sys.stderr, "Abundance file written to", output_prefix + '.abundance.txt'

def run_filtering_by_count(input_prefix, output_prefix, min_count):
    cmd = "filter_by_count.py {i} {o} --min_count {c}".format(i=input_prefix,
                                                              o=output_prefix,
                                                              c=min_count)
    print >> sys.stderr, "CMD:", cmd
    subprocess.check_call(cmd, shell=True)

def run_filtering_away_subsets(input_prefix, output_prefix, fuzzy_junction):
    cmd = "filter_away_subset_in_no5merge.py {i} {o} --fuzzy_junction {j}".format(i=input_prefix,
                                                                                  o=output_prefix,
                                                                                  j=fuzzy_junction)
    print >> sys.stderr, "CMD:", cmd
    subprocess.check_call(cmd, shell=True)

def tofu_wrap_main():
    parser = argparse.ArgumentParser(prog='tofu_wrap')
    add_cluster_arguments(parser, show_sge_env_name=True, show_sge_queue=True)

    parser.add_argument("--bin_size_kb", default=1, type=int, help="Bin size by kb (default: 1)")
    parser.add_argument("--bin_manual", default=None, help="Bin manual (ex: (1,2,3,5)), overwrites bin_size_kb")
    parser.add_argument("--bin_by_primer", default=False, action="store_true", help="Instead of binning by size, bin by primer (overwrites --bin_size_kb and --bin_manual)")
    parser.add_argument("--max_base_limit_MB", default=600, type=int, help="Maximum number of bases per partitioned bin, in MB (default: 600)")
    parser.add_argument("--gmap_name", default="hg19", help="GMAP DB name (default: hg19)")
    parser.add_argument("--gmap_db", default="/home/UNIXHOME/etseng/share/gmap_db_new/", help="GMAP DB location (default: /home/UNIXHOME/etseng/share/gmap_db_new/)")
    parser.add_argument("--output_seqid_prefix", type=str, default=None, help="Output seqid prefix. If not given, a random ID is generated")
    parser.add_argument("--mem_debug", default=False, action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--max_fuzzy_junction", default=5, type=int, help="Max fuzzy junction (default: 5 bp)")
    parser.add_argument("--version", action='version', version='%(prog)s ' + str(get_version()))
    args = parser.parse_args()

    # PRINT VERSION AND EXIT
#    if args.version:
#        print >> sys.stderr, get_version()
#        sys.exit(0)
    # DEBUG
    if args.mem_debug:
        from memory_profiler import memory_usage
    
    # #################################################################
    # SANITY CHECKS
    if not args.quiver:
        print >> sys.stderr, "--quiver must be turned on for tofu_wrap. Quit."
        sys.exit(-1)
    if args.nfl_fa is None:
        print >> sys.stderr, "--nfl_fa must be provided for tofu_wrap. Quit."
        sys.exit(-1)
    if not os.path.exists(args.gmap_db):
        print >> sys.stderr, "GMAP DB location not valid: {0}. Quit.".format(args.gmap_db)
        sys.exit(-1)
    if not os.path.exists(os.path.join(args.gmap_db, args.gmap_name)):
        print >> sys.stderr, "GMAP name not valid: {0}. Quit.".format(args.gmap_name)
        sys.exit(-1)
    # #################################################################

    tofu_prefix = binascii.b2a_hex(os.urandom(3)) if args.output_seqid_prefix is None else args.output_seqid_prefix

    ice_opts = IceOptions(quiver=args.quiver,
            use_finer_qv=args.use_finer_qv,
            targeted_isoseq=args.targeted_isoseq,
            ece_penalty=args.ece_penalty,
            ece_min_len=args.ece_min_len,
    )
    sge_opts = SgeOptions(unique_id=args.unique_id,
            use_sge=args.use_sge,
            max_sge_jobs=args.max_sge_jobs,
            blasr_nproc=args.blasr_nproc,
            quiver_nproc=args.quiver_nproc,
            gcon_nproc=args.gcon_nproc,
            sge_env_name=args.sge_env_name,
            sge_queue=args.sge_queue)
    ipq_opts = IceQuiverHQLQOptions(qv_trim_5=args.qv_trim_5,
            qv_trim_3=args.qv_trim_3,
            hq_quiver_min_accuracy=args.hq_quiver_min_accuracy,
            hq_isoforms_fa=args.hq_isoforms_fa,
            hq_isoforms_fq=args.hq_isoforms_fq,
            lq_isoforms_fa=args.lq_isoforms_fa,
            lq_isoforms_fq=args.lq_isoforms_fq)

    # ex: all_quivered_hq.100_30_0.99.fastq
    quiver_hq_filename = "all_quivered_hq.{0}_{1}_{2:.2f}.fastq".format(\
            args.qv_trim_5,args.qv_trim_3,args.hq_quiver_min_accuracy)
    quiver_lq_filename = "all_quivered_lq.fastq"

    # (1) separate input flnc into size bins or primers
    if args.bin_by_primer:
        split_files = sep_flnc_by_primer(args.flnc_fa, os.path.abspath(args.root_dir))
    else:
        bin_manual = eval(args.bin_manual) if args.bin_manual is not None else None
        split_files = sep_flnc_by_size(args.flnc_fa, args.root_dir, bin_size_kb=args.bin_size_kb, bin_manual=bin_manual, max_base_limit_MB=args.max_base_limit_MB)
    print >> sys.stderr, "split input {0} into {1} bins".format(args.flnc_fa, len(split_files))

    # (2) if fasta_fofn already is there, use it; otherwise make it first
    if args.quiver and args.fasta_fofn is None:
        print >> sys.stderr, "Making fasta_fofn now"
        nfl_dir = os.path.abspath(os.path.join(args.root_dir, "fasta_fofn_files"))
        if not os.path.exists(nfl_dir):
            os.makedirs(nfl_dir)
        args.fasta_fofn = os.path.join(nfl_dir, 'input.fasta.fofn')
        print >> sys.stderr, "fasta_fofn", args.fasta_fofn
        print >> sys.stderr, "nfl_dir", nfl_dir
        convert_fofn_to_fasta(fofn_filename=args.bas_fofn,
                            out_filename=args.fasta_fofn,
                            fasta_out_dir=nfl_dir,
                            cpus=args.blasr_nproc)
    else:
        if not os.path.exists(args.fasta_fofn):
            raise Exception, "fasta_fofn {0} does not exist!".format(args.fasta_fofn)
        for line in open(args.fasta_fofn):
            file = line.strip()
            if len(file) > 0 and not os.path.exists(file):
                raise Exception, "File {0} does not exists in {1}".format(file, args.fasta_fofn)

    # (3) run ICE/Quiver (the whole thing), providing the fasta_fofn
    split_dirs = []
    for cur_file in split_files:
        cur_dir = os.path.abspath(os.path.dirname(cur_file))
        split_dirs.append(cur_dir)
        cur_out_cons = os.path.join(cur_dir, args.consensusFa)
        
        hq_quiver = os.path.join(cur_dir, quiver_hq_filename)
        if os.path.exists(hq_quiver):
            print >> sys.stderr, "{0} already exists. SKIP!".format(hq_quiver)
            continue
        print >> sys.stderr, "running ICE/Quiver on", cur_dir
        start_t = time.time()

        obj = Cluster(root_dir=cur_dir,
                flnc_fa=cur_file,
                nfl_fa=realpath(args.nfl_fa),
                bas_fofn=realpath(args.bas_fofn),
                ccs_fofn=realpath(args.ccs_fofn),
                fasta_fofn=realpath(args.fasta_fofn),
                out_fa=cur_out_cons,
                sge_opts=sge_opts,
                ice_opts=ice_opts,
                ipq_opts=ipq_opts,
                report_fn=args.report_fn,
                summary_fn=args.summary_fn,
                nfl_reads_per_split=args.nfl_reads_per_split)
        
        # DEBUG
        if args.mem_debug: 
            mem_usage = memory_usage(obj.run, interval=60)
            end_t = time.time()
            with open('mem_debug.log', 'a') as f:
                f.write("Running ICE/Quiver on {0} took {1} secs.\n".format(cur_dir, end_t-start_t))
                f.write("Maximum memory usage: {0}\n".format(max(mem_usage)))
                f.write("Memory usage: {0}\n".format(mem_usage))
        else:
            obj.run()

    combined_dir = os.path.join(args.root_dir, 'combined')
    if not os.path.exists(combined_dir):
        os.makedirs(combined_dir)
    # (4) combine quivered HQ/LQ results
    hq_filename, lq_filename, hq_pre_dict, lq_pre_dict = \
            combine_quiver_results(split_dirs, combined_dir, quiver_hq_filename, quiver_lq_filename,\
            tofu_prefix)
    with open(os.path.join(args.root_dir, 'combined', 'combined.hq_lq_pre_dict.pickle'), 'w') as f:
        dump({'HQ': hq_pre_dict, 'LQ': lq_pre_dict}, f)
    # (5) collapse quivered HQ results
    collapse_prefix_hq = run_collapse_sam(hq_filename, args.gmap_db, args.gmap_name, cpus=args.blasr_nproc, max_fuzzy_junction=args.max_fuzzy_junction, dun_merge_5_shorter=True)
    # (6) make abundance 
    get_abundance(collapse_prefix_hq, hq_pre_dict, collapse_prefix_hq)
    # (7) run filtering & removing subsets in no5merge
    if args.targeted_isoseq:
        run_filtering_by_count(collapse_prefix_hq, collapse_prefix_hq+'.min_fl_5', min_count=5)
        run_filtering_away_subsets(collapse_prefix_hq+'.min_fl_5', collapse_prefix_hq+'.min_fl_5.filtered', args.max_fuzzy_junction)
    else:
        run_filtering_by_count(collapse_prefix_hq, collapse_prefix_hq+'.min_fl_2', min_count=2)
        run_filtering_away_subsets(collapse_prefix_hq+'.min_fl_2', collapse_prefix_hq+'.min_fl_2.filtered', args.max_fuzzy_junction)
    # Now do it for LQ (turned OFF for now)
    #collapse_prefix_lq = run_collapse_sam(lq_filename, args.gmap_db, args.gmap_name, cpus=args.blasr_nproc)
    #get_abundance(collapse_prefix_lq, lq_pre_dict, collapse_prefix_lq)

if __name__ == "__main__":
    tofu_wrap_main()
