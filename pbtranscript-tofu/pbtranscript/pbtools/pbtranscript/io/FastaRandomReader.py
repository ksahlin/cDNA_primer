"""
Define miscBio.FastaRandomReader.
It will be repaced by pbcore.io.FastaIO.FastaRandomReader.
"""
from pbcore.io.FastaIO import FastaRecord
from collections import namedtuple


Interval = namedtuple('Interval', ['start', 'end'])


class FastaRandomReader:

    """
        This is meant to substitute for the Bio.SeqIO.to_dict method since some fasta files
        are too big to fit entirely to memory. The only requirement is that every id line
        begins with the symbol >. It is ok for the sequences to stretch multiple lines.
        The sequences, when read, are returned as Bio.SeqRecord objects.

        Example:
            r = FastaRandomReader('output/test.fna')
            r['6C_49273_NC_008578/2259031-2259297'] ==> this shows a FastaRecord
    """

    def __init__(self, fasta_filename):
        self.f = open(fasta_filename)
        self.d = {}
        self.locations = []
        self.locations_d_key_map = {}

        while 1:
            line = self.f.readline()
            if len(line) == 0:
                break
            if line.startswith('>'):
                sid = line.strip()[1:].split(None, 1)[0]
                # the header MUST be just 1 line
                # if id in self.d:
                #    print "duplicate id {0}!!".format(id)
                self.d[sid] = self.f.tell()

    def __getitem__(self, k):
        if k not in self.d:
            errMsg = "key {k} not in {f}!".format(k=k, f=self.f.name)
            raise ValueError(errMsg)
        self.f.seek(self.d[k])
        content = ''
        for line in self.f:
            if line.startswith('>'):
                break
            content += line.strip()
        # return SeqRecord(Seq(content), id=k)
        return FastaRecord(k, sequence=content)

    def __len__(self):
        return len(self.d)

    def __delitem__(self, key):
        errMsg = "FastaRandomReader.__delitem__ not defined."
        raise NotImplementedError(errMsg)

    def __setitem__(self, key):
        errMsg = "FastaRandomReader.__setitem__ not defined."
        raise NotImplementedError(errMsg)

    def keys(self):
        """Return d.keys."""
        return self.d.keys()


class MetaSubreadFastaReader(object):

    """Reader for reading PabBio subreads in a list of fasta files."""

    def __init__(self, fasta_filenames):
        self.meta_f = {}
        self.meta_zmw_d = {}
        # record just the zmw and let the subread reader handle it
        for fn in fasta_filenames:
            self.meta_f[fn] = SubreadFastaReader(fn)
            # combine all the keys
            for k in self.meta_f[fn].zmw_d.keys():
                self.meta_zmw_d[k] = fn

    def __getitem__(self, k):
        """
        k -- could be zmw or subread id
        """
        if k.count('/') == 2:
            zmw = k[:k.rfind('/')]
        else:
            zmw = k
        return self.meta_f[self.meta_zmw_d[zmw]][k]

    def __len__(self):
        return len(sum([len(d) for d in self.meta_zmw_d]))

    def __delitem__(self, key):
        errMsg = "MetaSubreadFastaReader.__delitem__ not defined."
        raise NotImplementedError(errMsg)

    def __setitem__(self, key):
        errMsg = "MetaSubreadFastaReader.__setitem__ not defined."
        raise NotImplementedError(errMsg)


class SubreadFastaReader(object):

    """Reader for reading PabBio subreads in a fasta file."""

    def __init__(self, fasta_filename):
        self.f = open(fasta_filename)
        self.d = {}
        self.zmw_d = {}
        while 1:
            line = self.f.readline()
            if len(line) == 0:
                break
            if line.startswith('>'):
                # ex: m140..._s1_p0/155/0_1673 RQ=0.845
                rid = line.strip()[1:].split(None, 1)[0]
                # the header MUST be just 1 line
                zmw = rid[:rid.rfind('/')]
                self.d[rid] = (rid, self.f.tell())
                if zmw not in self.zmw_d:
                    self.zmw_d[zmw] = []
                self.zmw_d[zmw].append((rid, self.f.tell()))

    def __getitem__(self, k):
        """
        k --- should be <movie>/<zmw> or <movie>/<zmw>/<start_end>
        If former, return a list of records associated with that ZMW
        If latter, return just that record but still in a list
        """
        if k.count('/') == 2:  # is a subread
            if k not in self.d:
                raise ValueError("key {0} not in dictionary!".format(k))
            locations = [self.d[k]]
        else:  # is a ZMW
            if k not in self.zmw_d:
                raise ValueError("key {0} not in dictionary!".format(k))
            locations = self.zmw_d[k]
        output = []
        for seqid, loc in locations:
            self.f.seek(loc)
            content = ''
            for line in self.f:
                if line.startswith('>'):
                    break
                content += line.strip()
            output.append(FastaRecord(seqid, sequence=content))
        return output

    def keys(self):
        """return keys (subreads)."""
        return self.d.keys()

    def __len__(self):
        return len(self.d)

    def __delitem__(self, key):
        errMsg = "SubreadFastaReader.__delitem__ not defined."
        raise NotImplementedError(errMsg)

    def __setitem__(self, key):
        errMsg = "SubreadFastaReader.__setitem__ not defined."
        raise NotImplementedError(errMsg)
