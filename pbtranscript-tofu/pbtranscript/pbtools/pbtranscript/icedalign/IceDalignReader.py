__author__ = 'etseng@pacificbiosciences.com'

import pdb
import os, sys, subprocess
import numpy
from pbtools.pbtranscript.io.BLASRRecord import BLASRRecord
from pbtools.pbtranscript.ice.IceUtils import HitItem, eval_blasr_alignment, alignment_has_large_nonmatch
from pbtools.pbtranscript.ice.c_IceAlign import get_ece_arr_from_alignment
class LAshowAlignReader:
    """
    Reader for using LAshow with -a option that gives the alignments
    """
    def __init__(self, las_out_filename):
        self.las_out_filename = las_out_filename
        self.f = open(las_out_filename)

    def __iter__(self):
        return self

    def next(self):
        """
        Should be the printed out of running:
        LA4Ice -m -i0 -w100000 -b0 -a:{db} {las}
        """
        # first line is BLASR-like output
        # ex: 000000002 000002845 -1192 85.12 0 1855 3082 3082 0 2324 3516 3517 overlap
        # if is - strand, then strand=1, start=S, end=E means the sequence is
        # seq[S:E].reverse_complement()

        raw = self.f.readline().strip().split()
        if raw[0] == '+' and raw[1] == '+': # FALCON-added EOF signature
                raise StopIteration
        qID = int(raw[0]) + 1 # convert to 1-based
        sID = int(raw[1]) + 1 # convert to 1-based
        score = int(raw[2])
        iden = float(raw[3])
        qStrand = int(raw[4])
        qStart = int(raw[5]) # 0-based
        qEnd = int(raw[6])
        qLen = int(raw[7])
        sStrand = int(raw[8])
        sStart = int(raw[9])
        sEnd = int(raw[10])
        sLen = int(raw[11])

        self.f.readline() # blank line
        _qStart, qAln = self.f.readline().strip().split()
        assert (qStrand == 0 and int(_qStart)-1 == qStart) or (qStrand == 1 and int(_qStart)-1 == qLen-qEnd)
        alnStr = self.f.readline().strip()
        _sStart, sAln = self.f.readline().strip().split()[:2]
        assert (sStrand == 0 and int(_sStart)-1 == sStart) or (sStrand == 1 and int(_sStart)-1 == sLen-sEnd)
        return BLASRRecord(qID, qLen, qStart, qEnd, qStrand, sID, sLen, sStart, sEnd, sStrand, score, None, \
                    qAln=qAln, alnStr=alnStr, sAln=sAln, identity=iden, strand='+' if qStrand==sStrand else '-')


def dalign_against_ref(dazz_query_obj, dazz_db_obj, las_out_filename, is_FL, sID_starts_with_c,
                      qver_get_func, qvmean_get_func, qv_prob_threshold=.03,
                      ece_penalty=1, ece_min_len=20, same_strand_only=True, no_qv_or_aln_checking=False,
                      max_missed_start=50, max_missed_end=50):
    """
    Excluding criteria:
    (1) self hit
    (2) opposite strand hit  (should already be in the same orientation;
        can override with <same_strand_only> set to False)
    (3) less than 90% aligned or more than 50 bp missed

    qver_get_func --- should be basQV.basQVcacher.get() or
                      .get_smoothed(), or can just pass in
                      lambda (x, y): 1. to ignore QV
    """
    for r in LAshowAlignReader(las_out_filename):
        #pdb.set_trace()
        missed_q = r.qStart + r.qLength - r.qEnd
        missed_t = r.sStart + r.sLength - r.sEnd

        r.qID = dazz_query_obj[r.qID]
        r.sID = dazz_db_obj[r.sID]
        if sID_starts_with_c:
            # because all consensus should start with
            # c<cluster_index>
            assert r.sID.startswith('c')
            if r.sID.find('/') > 0:
                r.sID = r.sID.split('/')[0]
            if r.sID.endswith('_ref'):
                # probably c<cid>_ref
                cID = int(r.sID[1:-4])
            else:
                cID = int(r.sID[1:])
        else:
            cID = r.sID


        # self hit, useless!
        # (identity is removed here, NOT trustworthy using Jason's code calculations)
        # opposite strand not allowed!
        if (cID == r.qID or
                (r.strand == '-' and same_strand_only)):
            yield HitItem(qID=r.qID, cID=cID)
            continue

        # (Liz) this is used for partial_uc/nFL reads only
        # simply accepts hits from daligner for the nFL partial hits
        # testing shows that it does not affect much the Quiver consensus calling
        if no_qv_or_aln_checking:
            yield HitItem(qID=r.qID, cID=cID,
                              qStart=r.qStart, qEnd=r.qEnd,
                              missed_q=missed_q * 1. / r.qLength,
                              missed_t=missed_t * 1. / r.sLength,
                              fakecigar=1,
                              ece_arr=1)
            continue


        # full-length case: allow up to <max_missed_start> bp of 5' not aligned
        # and 50bp of 3' not aligned
        if (is_FL and (r.sStart > max_missed_start or r.qStart > max_missed_start or
                      (r.sLength - r.sEnd > max_missed_end) or
                      (r.qLength - r.qEnd > max_missed_end))):
            yield HitItem(qID=r.qID, cID=cID)
        else:
            cigar_str, ece_arr = eval_blasr_alignment(
                    record=r,
                    qver_get_func=qver_get_func,
                    sID_starts_with_c=sID_starts_with_c,
                    qv_prob_threshold=qv_prob_threshold,
                    qvmean_get_func=qvmean_get_func)
            #else: # don't use QV, just look at alignment


            if alignment_has_large_nonmatch(ece_arr,
                                            ece_penalty, ece_min_len):
                yield HitItem(qID=r.qID, cID=cID)
            else:
                yield HitItem(qID=r.qID, cID=cID,
                              qStart=r.qStart, qEnd=r.qEnd,
                              missed_q=missed_q * 1. / r.qLength,
                              missed_t=missed_t * 1. / r.sLength,
                              fakecigar=cigar_str,
                              ece_arr=ece_arr)


