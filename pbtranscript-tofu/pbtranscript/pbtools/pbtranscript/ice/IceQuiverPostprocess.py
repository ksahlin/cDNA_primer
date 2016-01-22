#!/usr/bin/env python
###############################################################################
# Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY
# THIS LICENSE.  THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR
# ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
###############################################################################

"""
Overview:
    pbtranscript cluster contains two main components:
    * (1) ICE (iterative clustering and error correction) to predict
      unpolished consensus isoforms.
    * (2) Polish, to use nfl reads and quiver to polish those predicted
      unpolished isoforms. Polish contains three steps:
      + (2.1) IceAllPartials (ice_partial.py all)
              Align and assign nfl reads to unpolished isoforms, and
              save results to a pickle file.
      + (2.2) IceQuiver (ice_quiver.py i and ice_quiver.py merge)
              Call quiver to polish each isoform based on alignments
              created by mapping its associated fl and nfl reads to
              this isoform.
      + (2.3) IceQuiverPostprocess (ice_quiver.py postprocess)
              Collect and post process quiver results, and classify
              HQ/LQ isoforms.

    In order to handle subtasks by SMRTPipe instead of pbtranscript
    itself, we will refactor the polish phase including
    (2.1) (2.2) and (2.3). The refactor of (2.1) is described in
    ice_partial.py.

    (2.2) IceQuiver will be refactored to
       + (2.2.1) IceQuiverI (ice_quiver.py i)
                 Split all unpolished isoforms into N chunks and
                 call Quiver to polish isoforms of the i-th chunk
                 at a time
       + (2.2.2) IceQuiverMerge (ice_quiver.py merge)
                 When all splitted quiver jobs are done,
                 collect all submitted jobs and save to
                 root_dir/log/submitted_quiver_jobs.txt
    (2.3) IceQuiverPostProcess's entry will be renamed from
          ice_post_quiver.py to:
       + (2.3.1) ice_quiver.py postprocess

   *** Here we are focusing on (2.3.1) 'ice_quiver.py postprocess' ***

Description:
    (2.3.1) ice_quiver.py postprocess

    Assumption:
     * Phase (1) ICE is done, unpolished isoforms are created, fl reads
       are assigned to isoforms, and saved to a pickle (i.e., final.pickle)
     * Step (2.1) IceAllPartials is done, all nfl reads are assigned
       to unpolished isoforms, saved to a pickle (i.e., nfl_all_pickle_fn)
     * Step (2.2.1) and (2.2.2) are done. All Quiver jobs finished and
       results merged.

    Process:
       Given root_dir, collect and output high QV and low QV isoforms, write
       a report and a summary.

    Input:
        Positional:
            root_dir, an output directory for running pbtranscript cluster.

    Output:
        Polished high QV and low QV isoforms in fasta/fastq, a cluster report
        and a summary.

    Hierarchy:
        pbtranscript = IceIterative

        pbtranscript --quiver = IceIterative + \
                                ice_polish.py

        ice_polish.py =  ice_make_fasta_fofn.py + \
                         ice_partial.py all + \
                         ice_quiver.py all

        ice_partial.py all = ice_partial.py split + \
                             ice_partial.py i + \
                             ice_partial.py merge

        (ice_partial.py one --> only apply ice_partial on a given input fasta)

        ice_quiver.py all = ice_quiver.py i + \
                            ice_quiver.py merge + \
                            ice_quiver.py postprocess

    Example:
        ice_quiver.py postprocess root_dir
"""

"""Post quiver, pick up the high QV and low QV conesnsus isoforms."""
import os
import re
import logging
import os.path as op
from collections import defaultdict
from cPickle import load
from time import sleep
from pbtools.pbtranscript.ClusterOptions import IceQuiverHQLQOptions
from pbtools.pbtranscript.PBTranscriptOptions import \
    add_cluster_root_dir_as_positional_argument, \
    add_ice_post_quiver_hq_lq_arguments, \
    add_cluster_summary_report_arguments
from pbtools.pbtranscript.Utils import phred_to_qv, \
    get_all_files_in_dir, ln, nfs_exists
from pbtools.pbtranscript.ice.IceFiles import IceFiles
from pbtools.pbtranscript.ice.IceUtils import cid_with_annotation, locally_run_failed_quiver_jobs
from pbcore.io import FastaWriter, FastqReader, FastqWriter


class IceQuiverPostprocess(IceFiles):

    """check if quiver jobs are finished and quiver results are compeleted.
       If quiver jobs are completed, pick up high QV consensus isoforms.
       If use_sge is True and quiver jobs are still running,
           * If quit_if_not_done is True, exit.
           * If quit_if_not_done is False, wait till quiver jobs are finished.

    """

    desc = "Post-quiver process, pick up high QV and low QV consensus isoforms."
    prog = "ice_quiver.py postprocess "

    def __init__(self, root_dir, ipq_opts,
                 use_sge=False, quit_if_not_done=True,
                 summary_fn=None, report_fn=None):
        super(IceQuiverPostprocess, self).__init__(
                prog_name="IceQuiverPostprocess",
                root_dir=root_dir)
        self.use_sge = use_sge
        self.quit_if_not_done = quit_if_not_done

        assert(type(ipq_opts) is IceQuiverHQLQOptions)
        self.ipq_opts = ipq_opts
        self.hq_isoforms_fa = ipq_opts.hq_isoforms_fa
        self.hq_isoforms_fq = ipq_opts.hq_isoforms_fq
        self.lq_isoforms_fa = ipq_opts.lq_isoforms_fa
        self.lq_isoforms_fq = ipq_opts.lq_isoforms_fq
        # Quiver usually can't call accurate QV on both ends
        # because of very well + less coverage
        # Ignore QV of the first 100 bp on 5' end
        self.qv_trim_5 = ipq_opts.qv_trim_5
        # Ignore QV of the last 30 bp on the 3' end"""
        self.qv_trim_3 = ipq_opts.qv_trim_3
        # The max allowed average base error rate within
        # seq[qv_trim_5:-qv_trim_3]
        self.hq_quiver_min_accuracy = ipq_opts.hq_quiver_min_accuracy

        self.fq_filenames = []

        self.report_fn = report_fn
        self.summary_fn = summary_fn

        self.validate_inputs()

    def get_existing_binned_quivered_fq(self):
        """Return all existing quivered fq files for binned clusters."""
        pattern = r"c(\d+)to(\d+)"  # e.g. c0to214
        fs = get_all_files_in_dir(self.quivered_dir,
                                  extension="quivered.fq")
        return [f for f in fs if re.search(pattern, f) is not None]

    def validate_inputs(self):
        """Validate if logs and pickle for non-full-length reads exist."""
        errMsg = ""

        if not nfs_exists(self.nfl_all_pickle_fn):
            errMsg = "Pickle file {f} ".format(f=self.nfl_all_pickle_fn) + \
                     "which assigns all non-full-length reads to isoforms " + \
                     "does not exist. Please check 'ice_partial.py *' are " + \
                     "all done."
        elif not nfs_exists(self.final_pickle_fn):
            errMsg = "Pickle file {f} ".format(f=self.final_pickle_fn) + \
                     "which assigns full-length non-chimeric reads to " + \
                     "isoforms does not exist."
        elif not nfs_exists(self.submitted_quiver_jobs_log):
            errMsg = "Log file {f}".format(f=self.submitted_quiver_jobs_log) + \
                     " of all submitted quiver jobs {f} does not exist."

        if errMsg != "":
            self.add_log(errMsg, level=logging.ERROR)
            raise IOError(errMsg)

    def check_quiver_jobs_completion(self):
        """Check whether quiver jobs are completed.
        submitted_quiver_jobs.txt should have format like:
        <job_id> \t ./quivered/<range>.sh

        (1) if all jobs are done and files are there return True
        (2) if all jobs are done but some files incomplete ask if to resubmit
        (3) if not all jobs are done, just quit
        fq_filenames contains all the finished fastq files.
        """
        self.add_log("Checking if quiver jobs are completed.")
        done_flag = True
        bad_sh = []
        self.fq_filenames = []
        submitted = {}
        self.add_log("Submitted quiver jobs are at {f}:".
                     format(f=self.submitted_quiver_jobs_log))

        sge_used = False
        with open(self.submitted_quiver_jobs_log, 'r') as f:
            for line in f:
                a, b = line.strip().split('\t')
                if a == 'local':
                    submitted[b] = b
                else:
                    sge_used = True
                    submitted[a] = b

        if sge_used is True and self.use_sge is True:
            stuff = os.popen("qstat").read().strip().split('\n')
            # first two lines are header
            running_jids = []
            for x in stuff[2:]:
                job_id = x.split()[0]
                running_jids.append(job_id)
                if job_id in submitted:
                    self.add_log("job {0} is still running.".format(job_id))
                    done_flag = False

        for job_id, sh_name in submitted.iteritems():
            fq_filename = op.join(self.quivered_dir,
                                  op.basename(sh_name).replace('.sh', '.quivered.fq'))

            if not nfs_exists(fq_filename) or \
                    os.stat(fq_filename).st_size == 0:
                if job_id in running_jids:  # still running, pass
                    done_flag = False
                else:
                    self.add_log("job {0} is completed but {1} is still empty!".
                                 format(job_id, fq_filename))
                    bad_sh.append(submitted[job_id])
            else:
                self.add_log("job {0} is done".format(job_id))
                self.fq_filenames.append(fq_filename)

        if not done_flag:
            if len(bad_sh) == 0:
                return "RUNNING"
            else:
                self.add_log("Some Quiver jobs failed. Attempt to rerun locally.\n")
                still_bad_sh = locally_run_failed_quiver_jobs(bad_sh)
                if len(still_bad_sh) > 0:
                    self.add_log("The following jobs were completed but " +
                             "no output file. Please check and resubmit: " +
                             "\n{0}\n".format('\n'.join(still_bad_sh)))
                    return "FAILED"
                else:
                    return "DONE"
        else:
            return "DONE"

    @property
    def quivered_good_fa(self):
        """Return $root_dir/all_quivered.hq.a_b_c.fasta"""
        return op.join(self.root_dir,
                       "all_quivered_hq.{a}_{b}_{c}.fasta".format(
                           a=self.qv_trim_5,
                           b=self.qv_trim_3,
                           c=self.hq_quiver_min_accuracy))

    @property
    def quivered_good_fq(self):
        """Return $root_dir/all_quivered_hq.a_b_c.fq"""
        return op.join(self.root_dir,
                       "all_quivered_hq.{a}_{b}_{c}.fastq".format(
                           a=self.qv_trim_5,
                           b=self.qv_trim_3,
                           c=self.hq_quiver_min_accuracy))

    @property
    def quivered_bad_fa(self):
        """Return $root_dir/all_quivered_lq.fa"""
        return op.join(self.root_dir, "all_quivered_lq.fasta")

    @property
    def quivered_bad_fq(self):
        """Return $root_dir/all_quivered_lq.fq"""
        return op.join(self.root_dir, "all_quivered_lq.fastq")

    def pickup_best_clusters(self, fq_filenames):
        """Pick up hiqh QV clusters."""
        self.add_log("Picking up the best clusters according to QVs from {fs}.".
                     format(fs=", ".join(fq_filenames)))
        a = load(open(self.final_pickle_fn))
        uc = a['uc']
        quivered = {}

        for fq in fq_filenames:
            self.add_log("Looking at quivered fq {f}".format(f=fq))
            for r in FastqReader(fq):
                # possible ID: c0/0_1611|quiver
                cid = r.name.split('|')[0]
                if cid.endswith('_ref'):
                    cid = cid[:-4]
                i = cid.find('/')
                if i > 0:
                    cid = cid[:i]
                cid = int(cid[1:])
                quivered[cid] = r

        good = []

        for cid, r in quivered.iteritems():
            qv_len = max(0, len(r.quality) - self.qv_trim_5 - self.qv_trim_3)
            if qv_len != 0:
                q = [phred_to_qv(x) for x in r.quality]
                err_sum = sum(q[self.qv_trim_5: -self.qv_trim_3])
                # LIZ HACK: definitely of HQ must include # of FL >= 2 !!!
                if 1.0 - (err_sum / float(qv_len)) >= self.hq_quiver_min_accuracy and len(uc[cid]) >= 2:
                    good.append(cid)

        partial_uc = load(open(self.nfl_all_pickle_fn))['partial_uc']
        partial_uc2 = defaultdict(lambda: [])
        partial_uc2.update(partial_uc)

        if self.report_fn is not None:
            self.write_report(report_fn=self.report_fn,
                              uc=uc, partial_uc=partial_uc2)

        self.add_log("Writing hiqh-quality isoforms to {f}|fq".
                     format(f=self.quivered_good_fa))
        self.add_log("Writing low-quality isoforms to {f}|fq".
                     format(f=self.quivered_bad_fa))
        with FastaWriter(self.quivered_good_fa) as good_fa_writer, \
                FastaWriter(self.quivered_bad_fa) as bad_fa_writer, \
                FastqWriter(self.quivered_good_fq) as good_fq_writer, \
                FastqWriter(self.quivered_bad_fq) as bad_fq_writer:
            for cid in quivered:
                r = quivered[cid]
                newname = "c{cid}/f{flnc_num}p{nfl_num}/{read_len}".\
                    format(cid=cid,
                           flnc_num=len(uc[cid]),
                           nfl_num=len(partial_uc2[cid]),
                           read_len=len(r.sequence))
                newname = cid_with_annotation(newname)

                if cid in good:
                    self.add_log("processing quivered cluster {c} --> good.".
                                 format(c=cid))
                    good_fa_writer.writeRecord(newname, r.sequence)
                    good_fq_writer.writeRecord(newname, r.sequence, r.quality)
                else:
                    self.add_log("processing quivered cluster {c} --> bad.".
                                 format(c=cid))
                    bad_fa_writer.writeRecord(newname, r.sequence)
                    bad_fq_writer.writeRecord(newname, r.sequence, r.quality)

        self.add_log("-" * 60, level=logging.INFO)
        self.add_log("High-quality Quivered consensus written " +
                     "to:\n{0}\n{1}".format(self.quivered_good_fa,
                                            self.quivered_good_fq))
        self.add_log("Low-qulality Quivered consensus written " +
                     "to:\n{0}\n{1}".format(self.quivered_bad_fa,
                                            self.quivered_bad_fq))
        self.add_log("-" * 60, level=logging.INFO)

    def cmd_str(self):
        """Return a cmd string (ice_quiver.py postprocess)."""
        return self._cmd_str(root_dir=self.root_dir,
                             ipq_opts=self.ipq_opts,
                             use_sge=self.use_sge,
                             quit_if_not_done=self.quit_if_not_done,
                             summary_fn=self.summary_fn,
                             report_fn=self.report_fn)

    def _cmd_str(self, root_dir, ipq_opts, use_sge, quit_if_not_done,
                 summary_fn, report_fn):
        """Return a cmd string (ice_quiver.py postprocess)."""
        cmd = self.prog + \
              "{d} ".format(d=root_dir) + \
              ipq_opts.cmd_str()
        if use_sge is True:
            cmd += "--use_sge "
        if quit_if_not_done is True:
            cmd += "--quit_if_not_done "
        if summary_fn is not None:
            cmd += "--summary={f} ".format(f=summary_fn)
        if report_fn is not None:
            cmd += "--report={f} ".format(f=report_fn)
        return cmd

    def run(self):
        """Check all quiver jobs are running, failed or done. Write high-quality
        consensus and low-quality consensus to all_quivered.good|bad.fa|fq.
        """
        job_stats = self.check_quiver_jobs_completion()
        self.add_log("quiver job status: {s}".format(s=job_stats))

        if self.use_sge is not True and job_stats != "DONE":
            self.add_log("quiver jobs were not submitted via sge, " +
                         "however are still incomplete. Please check.",
                         level=logging.ERROR)
            return -1
        elif self.use_sge is True:
            while job_stats != "DONE":
                self.add_log("Sleeping for 180 seconds.")
                sleep(180)
                job_stats = self.check_quiver_jobs_completion()
                if job_stats == "DONE":
                    break
                elif job_stats == "FAILED":
                    self.add_log("There are some failed jobs. Please check.",
                                 level=logging.ERROR)
                    return 1
                elif job_stats == "RUNNING":
                    self.add_log("There are jobs still running, waiting...",
                                 level=logging.INFO)
                    if self.quit_if_not_done is True:
                        return 0
                else:
                    msg = "Unable to recognize job_stats {s}".format(job_stats)
                    self.add_log(msg, logging.ERROR)
                    raise ValueError(msg)

        self.pickup_best_clusters(self.fq_filenames)

        self.add_log("Creating polished high quality consensus isoforms.")
        if self.hq_isoforms_fa is not None:
            ln(self.quivered_good_fa, self.hq_isoforms_fa)
        if self.hq_isoforms_fq is not None:
            ln(self.quivered_good_fq, self.hq_isoforms_fq)

        self.add_log("Creating polished low quality consensus isoforms.")
        if self.lq_isoforms_fa is not None:
            ln(self.quivered_bad_fa, self.lq_isoforms_fa)
        if self.lq_isoforms_fq is not None:
            ln(self.quivered_bad_fq, self.lq_isoforms_fq)

        if self.summary_fn is not None:
            self.write_summary(summary_fn=self.summary_fn,
                               isoforms_fa=self.final_consensus_fa,
                               hq_fa=self.hq_isoforms_fa,
                               lq_fa=self.lq_isoforms_fa)

        self.close_log()


def add_ice_quiver_postprocess_arguments(parser):
    """Add arugments for IceQuiverPostprocess (ice_quiver.py postprocess)."""
    parser = add_cluster_root_dir_as_positional_argument(parser)
    parser = add_ice_post_quiver_hq_lq_arguments(parser)
    parser = add_cluster_summary_report_arguments(parser)

    parser.add_argument("--use_sge",
                        default=False,
                        dest="use_sge",
                        action="store_true",
                        help="quiver jobs have been submitted to sge."
                             "Check qstat")
    parser.add_argument("--quit_if_not_done",
                        default=False,
                        dest="quit_if_not_done",
                        action="store_true",
                        help="Quit if quiver jobs haven't been completed.")
    return parser


# import sys
# from pbcore.util.ToolRunner import PBToolRunner
# from pbtools.pbtranscript.__init__ import get_version
#
#
# class IceQuiverPostprocessRunner(PBToolRunner):
#
#     """IceQuiverPostprocess runner"""
#
#     def __init__(self):
#         PBToolRunner.__init__(self, IceQuiverPostprocess.desc)
#         add_ice_quiver_postprocess_arguments(self.parser)
#
#     def getVersion(self):
#         """Get version string"""
#         return get_version()
#
#     def run(self):
#         """Run"""
#         logging.info("Running {f} v{v}.".format(f=op.basename(__file__),
#                                                 v=self.getVersion()))
#         args = self.args
#         cmd_str = ""
#         try:
#             ipq_opts = IceQuiverHQLQOptions(
#                 hq_isoforms_fa=args.hq_isoforms_fa,
#                 hq_isoforms_fq=args.hq_isoforms_fq,
#                 lq_isoforms_fa=args.lq_isoforms_fa,
#                 lq_isoforms_fq=args.lq_isoforms_fq,
#                 qv_trim_5=args.qv_trim_5,
#                 qv_trim_3=args.qv_trim_3,
#                 hq_quiver_min_accuracy=args.hq_quiver_min_accuracy)
#
#             obj = IceQuiverPostprocess(root_dir=args.root_dir,
#                 use_sge=args.use_sge,
#                 ipq_opts=ipq_opts,
#                 quit_if_not_done=args.quit_if_not_done,
#                 summary_fn=args.summary_fn,
#                 report_fn=args.report_fn)
#             cmd_str = obj.cmd_str()
#             obj.run()
#         except:
#             logging.exception("Exiting {cmd_str} with return code 1.".
#                               format(cmd_str=cmd_str))
#             return 1
#         return 0
#
#
# def main():
#     """Main function."""
#     runner = IceQuiverPostprocessRunner()
#     return runner.start()
#
# if __name__ == "__main__":
#     sys.exit(main())
