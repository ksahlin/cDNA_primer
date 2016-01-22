"""Define class `Cluster` and `ClusterException`."""
#!/usr/bin/env python
import os
import os.path as op
import logging
import cPickle
import time
from pbtools.pbtranscript.PBTranscriptException import PBTranscriptException
from pbtools.pbtranscript.io.FastqSplitter import splitFastq
from pbtools.pbtranscript.Utils import realpath, ln
from pbtools.pbtranscript.Polish import Polish
from pbtools.pbtranscript.ice.ProbModel import ProbFromModel, ProbFromQV, ProbFromFastq
from pbtools.pbtranscript.ice.IceFiles import IceFiles
from pbtools.pbtranscript.ice.IceInit import IceInit
from pbtools.pbtranscript.ice.IceIterative import IceIterative
from pbtools.pbtranscript.ice.IceUtils import ice_fa2fq, ice_fq2fa
from pbtools.pbtranscript.__init__ import get_version


class ClusterException(PBTranscriptException):

    """
    Exception class for Classifier.
    """

    def __init__(self, msg):
        PBTranscriptException.__init__(self, "cluster", msg)


class Cluster(IceFiles):

    """
    An object of `Cluster` calls the ICE algorithm to
    generate consensus isoforms.
    """

    def __init__(self, root_dir, flnc_fa, nfl_fa,
                 bas_fofn, ccs_fofn, out_fa,
                 sge_opts, ice_opts, ipq_opts,
                 report_fn=None, summary_fn=None, fasta_fofn=None, nfl_reads_per_split=30000):
        super(Cluster, self).__init__(prog_name="Cluster",
                                      root_dir=root_dir,
                                      bas_fofn=bas_fofn,
                                      ccs_fofn=ccs_fofn,
                                      fasta_fofn=fasta_fofn)

        self.sge_opts = sge_opts  # SGE, CPU arguments and etc
        self.ice_opts = ice_opts  # ICE clustering algorithm arguments
        self.ipq_opts = ipq_opts  # IceQuiver HQ/LQ isoform arguments
        self.nfl_reads_per_split = int(nfl_reads_per_split) # ToDo: setsanity check here?

        self.flnc_fa, self.nfl_fa, self.ccs_fofn, self.fasta_fofn = self._validate_inputs(
            flnc_fa, nfl_fa, ccs_fofn, fasta_fofn=fasta_fofn, quiver=self.ice_opts.quiver)

        self.root_dir, self.out_fa = self._validate_outputs(
            root_dir, out_fa)

        self.sanity_check()

        self._probqv = None     # probability & quality value

        self._flnc_splitted_fqs = []  # split flnc_fa into smaller files.
        #self._nflncSplittedFqs = []  # split nfl_fa into smaller files.  # (Liz) this is not used!
        self._logConfigs()      # Log configurations

        self.iceinit = None
        self.icec = None
        self.iceq = None
        self.pol = None

        self.add_log("Setting ece_penalty: {0} ece_min_len: {1}".format(ice_opts.ece_penalty, ice_opts.ece_min_len),\
                     level=logging.INFO)


        self.report_fn = realpath(report_fn) if report_fn is not None \
            else op.join(self.root_dir, "cluster_report.csv")
        self.summary_fn = realpath(summary_fn) if summary_fn is not None \
            else op.join(self.root_dir, "cluster_summary.txt")

        self.add_log("A Cluster Object created.", level=logging.INFO)

    def _validate_inputs(self, _flnc_fa, _nfl_fa, _ccs_fofn, fasta_fofn=None, quiver=False):
        """Validate input files and return absolute expaneded paths."""
        flnc_fa, nfl_fa, ccs_fofn = _flnc_fa, _nfl_fa, _ccs_fofn
        self.add_log("Checking input files.", level=logging.INFO)
        if flnc_fa is None:
            raise ClusterException("Input full-length non-chimeric reads " +
                    "files (i.e., flnc_fa) needs to be specified.")
        else:
            flnc_fa = realpath(flnc_fa)
            if not op.exists(flnc_fa):
                raise ClusterException("Unable to find full-length " +
                    "non-chimeric reads: {fn}".format(fn=flnc_fa))

        if nfl_fa is None:
            if quiver is True:
                raise ClusterException("Input non-full-length reads file " +
                    "(i.e., nfl_fa) needs to be specified for isoform polish.")
        else:
            nfl_fa = realpath(nfl_fa)
            if not op.exists(nfl_fa):
                raise ClusterException("Unable to find non-full-length " +
                    "non-chimeric reads: {fn}".format(fn=nfl_fa))

        if ccs_fofn is not None and not op.exists(ccs_fofn):
            raise ClusterException("Unable to find CCS FOFN file: " +
                                   "{fn}".format(fn=ccs_fofn))

        if fasta_fofn is not None and quiver:
            if not os.path.exists(fasta_fofn):
                raise ClusterException("Unable to find FASTA_FOFN file: {0}".format(fasta_fofn))
            for line in open(fasta_fofn):
                if not os.path.exists(line.strip()):
                    raise ClusterException("Unable to locate {0} in {1}".format(\
                            line.strip(), fasta_fofn))
        return (flnc_fa, nfl_fa, ccs_fofn, fasta_fofn)

    def _validate_outputs(self, _root_dir, _out_fa):
        """Validate outputs, create root_dir if it does not exist."""
        self.add_log("Checking outputs.", level=logging.INFO)
        root_dir, out_fa = _root_dir, _out_fa
        if root_dir is None:
            self.add_log("Output directory needs to be specified.",
                         level=logging.ERROR)
        if out_fa is None:
            self.add_log("Output consensus fasta needs to be specified.",
                         level=logging.ERROR)

        root_dir = realpath(root_dir)
        out_fa = realpath(out_fa)

        if op.exists(root_dir):
            self.add_log("Output directory {d} already exists.".
                         format(d=root_dir))
        else:
            self.add_log("Creating output directory {d}.".format(d=root_dir))
            os.mkdir(root_dir)
        if op.exists(out_fa):
            raise ClusterException("Consensus FASTA file {f} already exists.".
                                   format(f=out_fa))
        return root_dir, out_fa

    def sanity_check(self):
        """Do sanity check before stat to run."""
        errMsg = ""
        if self.ice_opts.quiver is True:
            if self.bas_fofn is None:
                errMsg = "A fofn of bas/bax.h5 files, e.g., input.fofn, " + \
                         "is required in order to polish consensus " + \
                         "isoforms using quiver."
            if self.nfl_fa is None:
                errMsg = "Non-full-length reads are required for polishing " + \
                         "consensus isoforms using quiver."
        if errMsg != "":
            self.add_log(errMsg, level=logging.ERROR)
            raise ValueError(errMsg)

    @property
    def configFN(self):
        """Return configuration file of the current run."""
        return op.join(self.root_dir, "run_ice_config.txt")

    def _logConfigs(self):
        """Log configuration."""
        with open(self.configFN, 'w', 0) as f:
            f.write('pbtranscript ' + get_version() + "\n")
            f.write(str(self.ice_opts) + "\n")
            f.write(str(self.sge_opts) + "\n")

    @property
    def initPickleFN(self):
        """Return path to pickle file with initial clusters."""
        return op.join(self.root_dir, "init.uc.pickle")

    def _setProbQV_ccs(self, ccs_fofn, firstSplitFa):
        """Set ProbQV from .ccs.h5"""
        start_t = time.time()
        self._probqv = ProbFromQV(ccs_fofn, firstSplitFa)
        self.add_log("Loading QVs from {f} + {c} took {t} sec.".
                format(f=firstSplitFa, c=ccs_fofn, t=time.time()-start_t),
                level=logging.INFO)

    def _setProbQV_fq(self, firstSplitFq=None):
        """Set probability and quality values from FASTQ"""
        start_t = time.time()
        self._probqv = ProbFromFastq(firstSplitFq)
        self.add_log("Loading QVs from {f} took {t} sec.".
                format(f=firstSplitFq, t=time.time()-start_t),
                level=logging.INFO)

    def run(self):
        """Call ICE to cluster consensus isoforms."""
        self.add_log("Start to run cluster.", level=logging.INFO)

        if self.ice_opts.targeted_isoseq:
            first_split = 1000
            self.ice_opts.flnc_reads_per_split = 10000
            self.add_log("targeted_isoseq: further splitting JUST first split to 1000. Changing flnc_reads_per_split=10000.")
        else:
            first_split = None

        # Split flnc_fa into smaller files and save files to _flnc_splitted_fas.
        self.add_log("Splitting {flnc} into ".format(flnc=self.flnc_fa) +
                     "smaller files each containing {n} reads.".format(
                     n=self.ice_opts.flnc_reads_per_split),
                     level=logging.INFO)

        self._flnc_splitted_fqs = splitFastq(
            input_fastq=self.flnc_fa,
            reads_per_split=self.ice_opts.flnc_reads_per_split,
            out_dir=self.root_dir,
            out_prefix="input.split",
            first_split=first_split)


        self.add_log("Splitted files are: " +
                     "\n".join(self._flnc_splitted_fqs),
                     level=logging.INFO)

        firstSplit_fq = self._flnc_splitted_fqs[0]
        firstSplit_fa = firstSplit_fq[:firstSplit_fq.rfind('.')] + '.fasta'
        ice_fq2fa(firstSplit_fq, firstSplit_fa)


        # -------- commenting out below because we already have FASTQ!
        #firstSplit = self._flnc_splitted_fas[0]
        #firstSplit_fq = firstSplit[:firstSplit.rfind('.')] + '.fastq'
        #self.add_log("Converting first split file {0} + {1} into fastq\n".format(\
        #        firstSplit, self.ccs_fofn), level=logging.INFO)
        ## Convert this into FASTQ
        #ice_fa2fq(firstSplit, self.ccs_fofn, firstSplit_fq)
        #
        ## Set up probabbility and quality value model
        #if self.ice_opts.use_finer_qv:
        #    self._setProbQV_ccs(self.ccs_fofn, firstSplit)
        #else:
        #    self._setProbQV_fq(firstSplitFq=firstSplit_fq)
        self._setProbQV_fq(firstSplitFq=firstSplit_fq)

        # Initialize cluster by clique
        if os.path.exists(self.initPickleFN):
            self.add_log("Reading existing uc pickle: {0}".format(self.initPickleFN), level=logging.INFO)
            with open(self.initPickleFN) as f:
                uc = cPickle.load(f)
        else:
            self.add_log("Finding maximal cliques: initializing IceInit.",
                         level=logging.INFO)
            self.iceinit = IceInit(readsFa=firstSplit_fa,
                                   qver_get_func=self._probqv.get_smoothed,
                                   ice_opts=self.ice_opts,
                                   sge_opts=self.sge_opts,
                                   qvmean_get_func=self._probqv.get_mean)
            uc = self.iceinit.uc

            # Dump uc to a file
            self.add_log("Dumping initial clusters to {f}".format(
                         f=self.initPickleFN), level=logging.INFO)
            with open(self.initPickleFN, 'w') as f:
                cPickle.dump(uc, f)

        # Run IceIterative.
        self.add_log("Iterative clustering: initializing IceIterative.",
                     level=logging.INFO)
        #self.add_log("In Cluster. DEBUG: Calling Cluster with {0} {1} {2} ".format(self.bas_fofn, self.ccs_fofn, self.fasta_fofn))


        self.icec = IceIterative(
            fasta_filename=firstSplit_fa,
            fastq_filename=firstSplit_fq,
            fastq_filenames_to_add=self._flnc_splitted_fqs[1:],
            all_fasta_filename=self.flnc_fa,
            ccs_fofn=self.ccs_fofn,
            root_dir=self.root_dir,
            ice_opts=self.ice_opts,
            sge_opts=self.sge_opts,
            uc=uc,
            probQV=self._probqv,
            use_ccs_qv=self.ice_opts.use_finer_qv)
        self.add_log("IceIterative log: {f}.".format(f=self.icec.log_fn))
        self.icec.run()
        self.add_log("IceIterative completed.", level=logging.INFO)

        # IceIterative done, write predicted (unplished) consensus isoforms
        # to an output fasta
        self.add_log("Creating a link to unpolished consensus isoforms.")
        ln(self.icec.final_consensus_fa, self.out_fa)

        #self.add_log("In Cluster. DEBUG: End Cluster with {0} {1} {2} ".format(self.bas_fofn, self.ccs_fofn, self.fasta_fofn))


        # Call quiver to polish predicted consensus isoforms.
        if self.ice_opts.quiver is not True:
            self.add_log("Creating a link to cluster report.",
                         level=logging.INFO)
            ln(src=self.icec.report_fn, dst=self.report_fn)

            # Summarize cluster and write to summary_fn.
            self.write_summary(summary_fn=self.summary_fn,
                               isoforms_fa=self.out_fa)
        else:  # self.ice_opts.quiver is True
            self.add_log("Polishing clusters: initializing IcePolish.",
                         level=logging.INFO)
            #self.add_log("In Cluster. DEBUG: Calling Polish with {0} {1} {2} ".format(self.bas_fofn, self.ccs_fofn, self.fasta_fofn))
            self.pol = Polish(root_dir=self.root_dir,
                              nfl_fa=self.nfl_fa,
                              bas_fofn=self.bas_fofn,
                              ccs_fofn=self.ccs_fofn,
                              fasta_fofn=self.fasta_fofn,
                              ice_opts=self.ice_opts,
                              sge_opts=self.sge_opts,
                              ipq_opts=self.ipq_opts,
                              nfl_reads_per_split=self.nfl_reads_per_split)
            self.add_log("IcePolish log: {f}.".format(f=self.pol.log_fn),
                         level=logging.INFO)
            self.pol.run()
            self.add_log("IcePolish completed.", level=logging.INFO)

            # cluster report
            self.add_log("Creating a link to cluster report.",
                         level=logging.INFO)
            ln(src=self.pol.iceq.report_fn, dst=self.report_fn)

            # Summarize cluster & polish and write to summary_fn.
            self.write_summary(summary_fn=self.summary_fn,
                               isoforms_fa=self.out_fa,
                               hq_fa=self.pol.icepq.quivered_good_fa,
                               lq_fa=self.pol.icepq.quivered_bad_fa)

        # Create log file.
        self.close_log()
        return 0
