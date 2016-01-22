"""Define util functions."""
import os.path as op
import os
import shutil
import logging
import sys
from time import sleep
from pbcore.util.Process import backticks
from pbcore.io.FastaIO import FastaReader
from pbcore.io.FastqIO import FastqReader


def check_ids_unique(fa_or_fq_filename, is_fq=False):
    """
    Confirm that a FASTA/FASTQ file has all unique IDs
    (used probably by collapse or fusion finding script)
    """
    if is_fq:
        reader = FastqReader(fa_or_fq_filename)
    else:
        reader = FastaReader(fa_or_fq_filename)
    seen = set()
    for r in reader:
        if r.id in seen:
            raise Exception, "Duplicate id {0} detected. Abort!".format(r.id)
        seen.add(r.id)


def revcmp(seq):
    """Given a sequence return its reverse complement sequence."""
    NTMAP = {'a': 't', 'c': 'g', 't': 'a', 'g': 'c',
             'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C',
             '*': '*', 'n': 'n', 'N': 'N'}
    return "".join([NTMAP[x] for x in seq.rstrip()])[::-1]


def realpath(f):
    """Return absolute, user expanded path."""
    if f is None:
        return None
    return op.abspath(op.expanduser(f))


def real_ppath(fn):
    """Return real 'python-style' path of a file.
    Consider files with white spaces in their paths, such as
    'res\ with\ space/out.sam' or 'res with space/out.sam',
    'res\ with\ space/out.sam' is unix-style file path.
    'res with space/out.sam' is python style file path.

    We need to convert all '\_' in path to ' ' so that python
    can handle files with space correctly, which means that
    'res\ with\ space/out.sam' will be converted to
    'res with space/out.sam'.
    """
    if fn is None:
        return None
    return op.abspath(op.expanduser(fn)).replace(r'\ ', ' ')


def real_upath(fn):
    """Return real 'unix-style' path of a file.
    Consider files with white spaces in their paths, such as
    'res\ with\ space/out.sam' or 'res with space/out.sam',
    'res\ with\ space/out.sam' is unix-style file path.
    'res with space/out.sam' is python style file path.

    We need to convert all ' ' to '\ ' so that unix can handle
    files with space correctly, which means that
    'res with space/out.sam' will be converted to
    'res\ with\ space/out.sam'.
    """
    if fn is None:
        return None
    return real_ppath(fn).replace(' ', r'\ ')

def nfs_exists(fn):
    """Detect whether a NFS file or a directory exists or not.
    In rare cases, a file f, which is created by a node X,
    may not be detected by node Y immediately due to NFS latency.
    In even more rare cases, node Y's shell script may be able
    to detect f's existence while Y's python can not, because of
    a python file system cache latency. So in order to eliminate
    this problem, first call 'ls' to trigger mount, then detect
    existence of fn using python.open(fn, 'r'), and try twice
    before finally give up.

    This script should return True, if fn is either
    an existing file or an existing directory created before
    nfs_exists is called. However, this script may return an
    arbitrary value, if fn is created or deleted at the same
    time or after nfs_exists is called.
    """
    # Call ls just to trigger mount, don't trust the return value.
    _o, _c, _m = backticks("ls {f}".format(f=fn))

    ERROR_NO_SUCH_FILE_OR_DIRECTORY = 2
    ERROR_IS_DIRECTORY = 21
    # Try to open fn with read-only mode
    try:
        with open(fn, 'r') as _reader:
            pass
        return True # fn is an existing file
    except IOError as err:
        if err.errno == ERROR_NO_SUCH_FILE_OR_DIRECTORY:
            # Wait 15 seconds and try it again
            sleep(15)
            try:
                with open(fn, 'r') as _reader:
                    pass
                return True # An newly detected existing file
            except IOError as err2:
                if err2.errno == ERROR_NO_SUCH_FILE_OR_DIRECTORY:
                    return False # fn does not exist
                elif err2.errno == ERROR_IS_DIRECTORY:
                    return True # fn is an existing directory
        elif err.errno == ERROR_IS_DIRECTORY:
            return True # fn is an existing directory
        else:
            return False # other IOErrors
    return False # other errors


def mkdir(path):
    """Create a directory if it does not pre-exist,
    otherwise, pass."""
    if not op.exists(path) and not op.lexists(path):
        try:
            os.makedirs(path)
        except OSError as e:
            # "File exists error" can happen when python
            # fails to syncronize with NFS or multiple
            # processes are trying to make the same dir.
            if e.errno == 17:
                pass
            else:
                raise OSError(e)


def mknewdir(path):
    """Create a new directory if it does not pre-exist,
    otherwise, delete it and then re-create it."""
    if op.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)


def touch(path):
    """touch a file."""
    if op.exists(path):
        os.utime(path, None)
    else:
        open(path, 'a').close()


def generateChunkedFN(out_dir, prefix, num_chunks):
    """Generate n chunked file names, e.g.
    outDir/$prefix.0, outDir/$prefix.1, ..., outDir/$prefix.num_chunks-1
    """
    return [op.join(out_dir, prefix + "." + str(i))
            for i in xrange(0, num_chunks)]


def get_files_from_fofn(fofn_filename):
    """Return a list of file names within a fofn file."""
    fns = []
    try:
        with open(fofn_filename, 'r') as fofn:
            for line in fofn:
                fns.append(realpath(line.strip()))
    except (IOError, OSError) as e:
        raise IOError("Failed to read from fofn file {fofn}.\n".
                      format(fofn=fofn_filename) + str(e))
    return fns


def write_files_to_fofn(file_names, fofn_filename):
    """Write files in list file_names to fofn_filename."""
    try:
        with open(fofn_filename, 'w') as fofn:
            for fn in file_names:
                fofn.write(str(fn) + "\n")
    except (IOError, OSError) as e:
        raise IOError("Failed to files to fofn file {fofn}.\n".
                      format(fofn=fofn_filename) + str(e))


def validate_fofn(fofn_filename):
    """Validate existence of FOFN and files within the FOFN.

    :param fofn: (str) Path to File of file names or None.
    :raises: IOError if any file is not found.
    :return: (str) input fofn or None
    """
    if fofn_filename is None:
        return None

    if nfs_exists(fofn_filename):
        fns = get_files_from_fofn(fofn_filename)
        for fn in fns:
            if not nfs_exists(fn):
                raise IOError("Unable to find {f} in FOFN {fofn}.".
                              format(f=fn, fofn=fofn_filename))
        return fofn_filename
    else:
        raise IOError("Unable to find FOFN {fofn}.".
                      format(fofn=fofn_filename))


def setup_log(alog, file_name=None, level=logging.DEBUG, str_formatter=None):
    """
    Copied from mkocher's pbreports/utils.py.
    Util function for setting up logging.

    Due to how smrtpipe logs, the default behavior is that the stdout
    is where the logging is redirected. If a file name is given the log
    will be written to that file.

    :param log: (log instance) Log instance that handlers and filters will
    be added.
    :param file_name: (str, None), Path to file. If None, stdout will be used.
    :param level: (int) logging level
    """
    if file_name is None:
        handler = logging.StreamHandler(sys.stdout)
    else:
        handler = logging.FileHandler(file_name)

    if str_formatter is None:
        str_formatter = '[%(levelname)s] %(asctime)-15s ' + \
                        '[%(name)s %(funcName)s %(lineno)d] %(message)s'

    formatter = logging.Formatter(str_formatter)
    handler.setFormatter(formatter)
    alog.addHandler(handler)
    alog.setLevel(level)


def now_str():
    """Return string of current time."""
    import datetime
    return str(datetime.datetime.now()).split(".")[0]


def phred_to_qv(phred):
    """Phred value to quality value."""
    return 10 ** -(phred / 10.0)


def cat_files(src, dst):
    """Concatenate files in src and save to dst.
       src --- source file names in a list
       dst --- destinate file name
    """
    if src is None or len(src) == 0:
        raise ValueError("src should contain at least one file.")
    if dst in src:
        raise IOError("Unable to cat a file and save to itself.")

    with open(real_ppath(dst), 'w') as writer:
        for src_f in src:
            with open(real_ppath(src_f), 'r') as reader:
                for line in reader:
                    writer.write(line.rstrip() + '\n')


def get_all_files_in_dir(dir_path, extension=None):
    """return all files in a directory."""
    fs = []
    for f in os.listdir(dir_path):
        if extension is None:
            fs.append(f)
        else:
            if f.endswith(extension):
                fs.append(f)
    return fs


class CIGAR(object):

    """Cigar string."""

    def __init__(self, cigar_str):
        self.cigar_str = cigar_str
        self.num_match = 0
        self.num_mismatch = 0
        self.num_insert = 0
        self.num_deletion = 0
        self.num_hardclip = 0
        self.num_softclip = 0
        self.num_padding = 0
        self.num_unknown = 0
        self.parse(self.cigar_str)

    def parse(self, cigar_str):
        """Parse cigar string."""
        s = cigar_str
        if s == "*" or s == "=":
            return
        while(len(s) > 0):
            i = 0
            while s[i].isdigit() and i < len(s):
                i += 1
            num = int(s[:i])
            action = s[i]
            s = s[i + 1:]
            if action == "M":
                self.num_match += num
            elif action == "X":
                self.num_mismatch += num
            elif action == "I":
                self.num_insert += num
            elif action == "D":
                self.num_deletion += num
            elif action == "H":
                self.num_hardclip += num
            elif action == "S":
                self.num_softclip += num
            elif action == "P":
                self.num_padding += num
            elif action == "N":
                self.num_unknown += num
            else:
                raise ValueError("Can not parse CIGAR string " +
                                 "{s}".format(s=cigar_str))
# using regular expression is 20% slower than naive way
#        pattern = r"^(\d+)(M|I|D|X|S|H|P|N)(.*)$"
#        while(len(s) > 0):
#            m = re.search(pattern, s)
#            if m:
#                num = int(m.groups()[0])
#                action = m.groups()[1]

    def __str__(self):
        return "Match = {m}, ".format(m=self.num_match) + \
               "Mismatch = {m}, ".format(m=self.num_mismatch) + \
               "Insert = {m}, ".format(m=self.num_insert) + \
               "Deletion = {m}, ".format(m=self.num_deletion) + \
               "HardClipping = {m}, ".format(m=self.num_hardclip) + \
               "SoftClipping = {m}, ".format(m=self.num_softclip) + \
               "Padding = {m}, ".format(m=self.num_padding)

    def match_seq(self, seq):
        """Return if this cigar string matches the sequence."""
        return self.cigar_str == "*" or self.cigar_str == "=" or \
            (self.num_match + self.num_insert + self.num_softclip == len(seq))


def cigar_match_seq(sam_str):
    """Return True if cigar length match sequence length, otherwise, False"""
    fields = sam_str.split('\t')
    cigar, seq = CIGAR(fields[5]), fields[9]
    if cigar.match_seq(seq):
        return True
    else:
        return False


def filter_sam(in_sam, out_sam):
    """Filter sam alignments with bad cigar string."""
    if not op.exists(in_sam):
        raise IOError("Unable to find input sam {f}".format(f=in_sam))
    if realpath(in_sam) == realpath(out_sam):
        raise IOError("in_sam and out_sam can not be identical.")
    with open(in_sam, 'r') as reader, \
            open(out_sam, 'w') as writer:
        for l, in_line in enumerate(reader):
            logging.info("Processing {l}".format(l=in_line))
            if in_line.startswith("#") or \
               in_line.startswith("@") or \
               cigar_match_seq(in_line):
                writer.write(in_line)
            else:
                logging.warn("line {l}, cigar does not match sequence.".
                             format(l=l + 1))


def ln(src, dst):
    """if src and dst are identical, pass. Otherwise, create dst, a soft
    symbolic link pointing to src."""
    if realpath(src) != realpath(dst):
        if op.exists(dst) or op.lexists(dst):
            os.remove(dst)
        logging.debug("Creating a symbolic link {dst} pointing to {src}".
                      format(dst=dst, src=src))
        os.symlink(src, dst)
