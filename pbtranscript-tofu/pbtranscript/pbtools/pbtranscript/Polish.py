#!/usr/bin/env python

"""Polish clustered isoforms using quiver."""

import sys
import os.path as op
import logging
#from pbtools.pbtranscript.io.FastaSplitter import splitFasta
from pbtools.pbtranscript.io.FastqSplitter import splitFastq
from pbtools.pbtranscript.__init__ import get_version
from pbtools.pbtranscript.Utils import realpath
from pbtools.pbtranscript.ClusterOptions import IceOptions, SgeOptions, \
    IceQuiverHQLQOptions
from pbtools.pbtranscript.ice.IceUtils import build_sa, \
    convert_fofn_to_fasta
from pbtools.pbtranscript.PBTranscriptOptions import add_fofn_arguments, \
    add_ice_post_quiver_hq_lq_arguments, add_sge_arguments, \
    add_nfl_fa_argument, add_cluster_root_dir_as_positional_argument
from pbtools.pbtranscript.ice.IceFiles import IceFiles
from pbtools.pbtranscript.ice.IceAllPartials import IceAllPartials
from pbtools.pbtranscript.ice.IceQuiver import IceQuiver
from pbtools.pbtranscript.ice.IceQuiverPostprocess import IceQuiverPostprocess
from pbtools.pbtranscript.icedalign.IceDalignUtils import DazzIDHandler, DalignerRunner


class Polish(IceFiles):

    """Polish isoforms clusters using Quiver."""

    def __init__(self, root_dir, nfl_fa, bas_fofn, ccs_fofn,
            ice_opts, sge_opts, ipq_opts, fasta_fofn=None, nfl_reads_per_split=30000):
        """
        root_dir --- IceFiles.root_dir, usually data/clusterOutDir
        nfl_fa    --- non-full-length reads in fasta, e.g., isoseq_nfl.fa
        bas_fofn --- e.g. input.fofn of bas|bax.h5 files
        ccs_fofn --- e.g. reads_of_insert.fofn of ccs files.

        ipq_opts --- IceQuiverHQLQOptions
                     qv_trim_5: ignore QV of n bases in the 5' end
                     qv_trim_3: ignore QV of n bases in the 3' end
                     hq_quiver_min_accuracy: minimum allowed quiver accuracy
                                      to mark an isoform as high quality
                     hq_isoforms_fa|fq: polished, hiqh quality consensus
                                        isoforms in fasta|q
                     lq_isoforms_fa|fq: polished, low quality consensus
                                        isoforms in fasta|q
        """
        IceFiles.__init__(self, prog_name="IcePolish", root_dir=root_dir,
                          bas_fofn=bas_fofn, ccs_fofn=ccs_fofn,
                          fasta_fofn=fasta_fofn)

        #self.add_log("DEBUG: in Polish ccs_fofn is {0}".format(self.ccs_fofn))
        #self.add_log("DEBUG: in Polish fasta_fofn is {0}".format(self.fasta_fofn))
        #self.add_log("DEBUG: in Polish bas_fofn is {0}".format(self.bas_fofn))
        self.nfl_fa = realpath(nfl_fa)
        self.nfl_reads_per_split = nfl_reads_per_split
        self.ice_opts = ice_opts
        self.sge_opts = sge_opts
        self.ipq_opts = ipq_opts

        self.add_log("ece_penalty: {0}, ece_min_len: {1}".format(self.ice_opts.ece_penalty, self.ice_opts.ece_min_len))

        self.icep = None   # IceAllPartials.
        self.iceq = None   # IceQuiver
        self.icepq = None  # IceQuiverPostprocess
        self._nfl_splitted_fas = None

        self.validate_inputs()

    def validate_inputs(self):
        """
        Validate input directories: root_dir, and
        files: nfl_fa, bas_fofn, ccs_fofn.
        """
        self.add_log("Validating inputs.")
        errMsg = ""
        if not op.exists(self.root_dir):
            errMsg = "Root dir {d} is not an existing directory!".\
                format(d=self.root_dir)
        if not op.exists(self.nfl_fa):
            errMsg = "Failed to find non-full-length reads {f}!".\
                format(f=self.nfl_fa)
        if not op.exists(self.bas_fofn):
            errMsg = "Failed to find bas fofn {f}!".format(f=self.bas_fofn)
        if errMsg != "":
            self.add_log(errMsg, level=logging.ERROR)
            raise ValueError(errMsg)

    # def get_sa_file(self):
    #     """Return saffix array of final_consensus_fa."""
    #     ret = True
    #     if not op.exists(self.final_consensus_sa):
    #         ret = build_sa(self.final_consensus_fa, self.final_consensus_sa)
    #
    #     if ret is True:
    #         self.add_log("Successfully created a suffix array for {f} at {sa}.".
    #                      format(f=self.final_consensus_fa,
    #                             sa=self.final_consensus_sa))
    #     else:
    #         self.add_log("Failed to generate suffix array for {f}".format(
    #             f=self.final_consensus_fa), level=logging.WARNING)
    #     if op.exists(self.final_consensus_sa):
    #         return self.final_consensus_sa
    #     else:
    #         return None

    def run(self):
        """
        First, split non-full-length (nfl) fasta files into smaller
        chunks, assign nfl reads in each splitted fasta file
        into unpolished isoform clusters and then merge all pickles
        into self.nfl_all_pickle_fn.
        Second, bin every 100 clusters, for each bin, call blasr,
        samto5h, loadPulses, cmph5tools to create cmp.h5 files and
        call quiver to polish each isoforms within each bin.
        Finally, pick up good isoform clusters whose QV errors is less
        than a threshold.
        Save all high quality isoforms to hq_isoforms_fa|fq if they are not None
        Save all low quality isoforms to lq_isoforms_fa|fq if they are not None
        """
        # Create final.consensus.fa.sa
        #self.add_log("Generating suffix array for {f}".format(
        #             f=self.final_consensus_sa), level=logging.INFO)
        #sa_file = self.get_sa_file()

        # Create input.fasta.fofn from bas_fofn
        self.add_log("Creating fasta fofn from bas/bax.h5 fofn",
                     level=logging.INFO)
        if self.fasta_fofn is None:
            self.fasta_fofn = op.join(self.nfl_dir, "input.fasta.fofn")
        self.add_log("bas fofn={f}".format(f=self.bas_fofn))
        self.add_log("fasta fofn={f}".format(f=self.fasta_fofn))
        if op.exists(self.fasta_fofn):
            self.add_log("No need to run convert_fofn_to_fasta.")
        else:
            convert_fofn_to_fasta(fofn_filename=self.bas_fofn,
                                out_filename=self.fasta_fofn,
                                fasta_out_dir=self.nfl_dir,
                                cpus=self.sge_opts.blasr_nproc)

        # Split non-full-length reads into smaller fasta files
        # and save files to root_dir/nfl_00.fa, ..., .
        self.add_log("Splitting {nfl} into ".format(nfl=self.nfl_fa) +
                     "smaller files each containing {n} reads.".format(
                     n=self.nfl_reads_per_split),
                     level=logging.INFO)
        self._nfl_splitted_fqs = splitFastq(input_fastq=self.nfl_fa,
                                            reads_per_split=self.nfl_reads_per_split,
                                            out_dir=self.nfl_dir,
                                            out_prefix="input.split")
        msg = "Splitted files are: " + "\n".join(self._nfl_splitted_fqs)
        self.add_log(msg, level=logging.INFO)

        # Generating dazz DB for final.consensus.fasta
        ref_obj = DazzIDHandler(self.final_consensus_fa, False)
        DalignerRunner.make_db(ref_obj.dazz_filename)
        msg = "Dazz DB made for: " + ref_obj.dazz_filename
        self.add_log(msg, level=logging.INFO)

        # Process nfl reads in each splitted fasta.
        self.add_log("Initializing IceAllPartials.", level=logging.INFO)
        #sa_file = self.final_consensus_sa \
        #    if op.exists(self.final_consensus_fa) else None

        self.icep = IceAllPartials(
            root_dir=self.root_dir,
            fastq_filenames=self._nfl_splitted_fqs,
            ref_fasta=self.final_consensus_fa,
            out_pickle=self.nfl_all_pickle_fn,
            sge_opts=self.sge_opts,
            sa_file=None,  # since we are switching to daligner, just give it as None now; remove sa_file completely later when daligner is mature (ToDo)
            ccs_fofn=self.ccs_fofn)
        self.add_log("IceAllPartials log: {f}.".format(f=self.icep.log_fn),
                     level=logging.INFO)
        self.icep.run()
        self.add_log("IceAllPartials completed.", level=logging.INFO)

        self.add_log("Initializing IceQuiver.", level=logging.INFO)
        self.iceq = IceQuiver(root_dir=self.root_dir,
                              bas_fofn=self.bas_fofn,
                              fasta_fofn=self.fasta_fofn,
                              sge_opts=self.sge_opts)
        self.add_log("IceQuiver log: {f}.".format(f=self.iceq.log_fn),
                     level=logging.INFO)
        self.iceq.run()
        self.add_log("IceQuiver finished.", level=logging.INFO)

        self.add_log("Initializing IceQuiverPostprocess.", level=logging.INFO)
        self.icepq = IceQuiverPostprocess(root_dir=self.root_dir,
                                          use_sge=self.sge_opts.use_sge,
                                          quit_if_not_done=False,
                                          ipq_opts=self.ipq_opts)
        self.add_log("IceQuiverPostprocess log: {f}.".
                     format(f=self.icepq.log_fn), level=logging.INFO)
        self.icepq.run()
        self.add_log("IceQuiverPostprocess finished.", level=logging.INFO)


def add_ice_polish_arguments(parser):
    """Set up argument parser."""
    parser = add_cluster_root_dir_as_positional_argument(parser)
    parser = add_nfl_fa_argument(parser, positional=True)
    parser = add_fofn_arguments(parser, ccs_fofn=True, bas_fofn=True, fasta_fofn=True)
    parser = add_ice_post_quiver_hq_lq_arguments(parser)
    parser = add_sge_arguments(parser, quiver_nproc=True, blasr_nproc=True, sge_env_name=True, sge_queue=True)
    return parser


from pbcore.util.ToolRunner import PBToolRunner


class PolishRunner(PBToolRunner):

    """Polish Runner"""

    def __init__(self):
        desc = "Call quiver to polish consensus isoforms with " + \
               "non-full-length reads."
        PBToolRunner.__init__(self, desc)
        add_ice_polish_arguments(self.parser)

    def getVersion(self):
        """Get version string."""
        return get_version()

    def run(self):
        """Run"""
        logging.info("Running {f} v{v}.".format(f=op.basename(__file__),
                                                v=self.getVersion()))
        args = self.args

        sge_opts = SgeOptions(unique_id=args.unique_id,
                              use_sge=args.use_sge,
                              max_sge_jobs=args.max_sge_jobs,
                              quiver_nproc=args.quiver_nproc,
                              blasr_nproc=args.blasr_nproc,
                              sge_env_name=args.sge_env_name,
                              sge_queue=args.sge_queue)
        ipq_opts = IceQuiverHQLQOptions(
            hq_isoforms_fa=args.hq_isoforms_fa,
            hq_isoforms_fq=args.hq_isoforms_fq,
            lq_isoforms_fa=args.lq_isoforms_fa,
            lq_isoforms_fq=args.lq_isoforms_fq,
            qv_trim_5=args.qv_trim_5,
            qv_trim_3=args.qv_trim_3,
            hq_quiver_min_accuracy=args.hq_quiver_min_accuracy)
        try:
            obj = Polish(root_dir=args.root_dir,
                         nfl_fa=args.nfl_fa,
                         bas_fofn=args.bas_fofn,
                         ccs_fofn=args.ccs_fofn,
                         fasta_fofn=args.fasta_fofn,
                         sge_opts=sge_opts,
                         ice_opts=IceOptions(),
                         ipq_opts=ipq_opts,
                         nfl_reads_per_split=args.nfl_reads_per_split)
            obj.run()
        except Exception as e:
            logging.error(str(e))
            import traceback
            traceback.print_exc()
            return 1
        return 0


def main():
    """Main function."""
    runner = PolishRunner()
    return runner.start()

if __name__ == "__main__":
    sys.exit(main())
