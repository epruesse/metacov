# cython: c_string_type=unicode, c_string_encoding=ascii
cimport cython

from cpython.array cimport array, clone, resize
from libc.stdio cimport fopen, fdopen, fclose, fwrite, getline, FILE, perror, setbuf
from libc.stdlib cimport malloc, free, atoi
from libc.stdint cimport uint8_t, uint64_t
from libc.errno cimport errno

from subprocess import Popen, PIPE, check_output
from fcntl import fcntl, F_SETFL, F_GETFL
import os, re, sys

from pysam.libchtslib cimport BAM_FREAD1, BAM_FREAD2, BAM_FPAIRED
from pysam.libcalignmentfile cimport AlignmentFile, IteratorRowRegion, IteratorRowAll
from pysam.libcalignedsegment cimport pysam_bam_get_cigar, pysam_bam_get_seq, \
    pysam_bam_get_qual, pysam_bam_get_qname, PileupColumn

from metacov.compat import FileNotFoundError

# A->0
# C->1
# G->2
# T->3
# rest -> 4


cdef uint8_t *iupac_to_nt4_table = [
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, # 0
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, # 16
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, # 32
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, # 48

    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, # 64
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, # 80
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, # 96
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, # 112

    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, # 128
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, #
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, #
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, #

    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, #
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, #
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, #
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4  #
]
cdef inline uint8_t iupac_to_nt4(unsigned char n) nogil:
    return iupac_to_nt4_table[n]


cdef long gzip_get_size(filename):
    """Get the (uncompressed) size of a gzip'ed file"""
    cdef:
       long size = os.path.getsize(filename)
       long guess = 0
    with open(filename, "rb") as f:
       f.seek(size-4)
       guess = int.from_bytes(f.read(4), 'little')
    while guess < size * 2:
        guess = guess + 2**32
    return guess


cdef class FastQFile:
    """
    FastQ File Class (read only)
    """
    def __cinit__(self):
        cdef int i
        self.proc = None # for unzip process
        self.fstream = NULL
        for i in range(4):
            self.buf[i] = NULL
            self.buf_len[i] = 0
            self.len[i] = 0
        self._file_size = 0
        self._pos = 0

        self.rseq_arr = clone(array('B'), 100, False)
        self.rseq = self.rseq_arr

    def __init__(self, str filename, int max_linelen=1000):
        """
        Arguments:
            filename:    name of file to open
            max_linelen: maximum length of each line/read
        """
        self.set_max_linelen(max_linelen)
        self.filename = filename

    def __dealloc__(self):
        cdef int i
        for i in range(4):
            if self.buf[i]:
                free(self.buf[i])
                self.buf[i] = NULL

    def __enter__(self):
        if self.filename.endswith(".gz"):
            try:
                self.proc = Popen(["unpigz", "-c", self.filename],
                                  stdout=PIPE, shell=False)
            except FileNotFoundError:
                self.proc = Popen(["gunzip", "-c", self.filename],
                                      stdout=PIPE, shell=False)
            self.fstream = fdopen(self.proc.stdout.fileno(), "r")
        else:
            self.fstream = fopen(self.filename.encode("UTF-8"), "r")
        if not self.fstream:
            raise OSError(errno, "Failed to open", self.filename)
        setbuf(self.fstream, <char*>malloc(8*1024))
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        if self.fstream:
            fclose(self.fstream)
        if self.proc:
            self.proc.wait(5)

    def __iter__(self):
        return self

    def __next__(self):
        if self.cnext() <= 0:
            raise StopIteration()
        return self

    @property
    def rlen(self):
        """Length of current read"""
        return self.len[1]

    @property
    def seq(self):
        """List of bases of current read (ACGTN->01234)"""
        cdef int i
        self.get_seq()
        return [self.rseq[i] for i in range(self.rlen)]

    @property
    def char_seq(self):
        return self.buf[1][:self.len[1]]

    @property
    def pos(self):
        """Position in file in kb"""
        return long(self._pos / 1024)

    @property
    def size(self):
        """Size of file in kb"""
        if not self._file_size:
            if self.filename.endswith(".gz"):
                self._file_size = gzip_get_size(self.filename)
            else:
                self._file_size = os.path.getsize(self.filename)
        return long(self._file_size / 1024)

    cdef void set_max_linelen(self, int max_linelen):
        """Change the maximum read/line length"""
        self.max_linelen = max_linelen
        resize(self.rseq_arr, max_linelen)
        self.rseq = self.rseq_arr

    cdef int cnext(self) nogil:
        """Parse the next read"""
        cdef int i, l
        if not self.fstream:
            return -2
        for i in range(4):
            l = getline(&self.buf[i], &self.buf_len[i], self.fstream)
            if l == -1: return -1
            self._pos += l
            self.len[i] = l
        return 1

    @cython.boundscheck(False)
    cdef public uint8_t[:] get_seq(self) nogil:
        """Convert current read from text to numeric (ACGTN->01234)"""
        cdef int i
        for i in range(self.len[1]):
            self.rseq[i] = iupac_to_nt4(self.buf[1][i])
        return self.rseq


    cdef public int get_len(self) nogil:
        """Get the length of the current read"""
        return self.len[1]

    cdef public int get_flags(self) nogil:
        return 0


cdef class FastQFilePair(FastQFile):
    """
    FastQ File Pair Class

    Essentially `zip(FastQFile, FastQFile)`
    """
    def __cinit__(self):
        pass

    def __init__(self, str read1, str read2, int max_linelen=1000):
        self.read1 = FastQFile(read1, max_linelen)
        self.read2 = FastQFile(read2, max_linelen)
        self.flags1 = BAM_FPAIRED | BAM_FREAD1
        self.flags2 = BAM_FPAIRED | BAM_FREAD2

    def __enter__(self):
        self.read1.__enter__()
        self.read2.__enter__()
        self.cur = 0
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.read1.__exit__(exception_type, exception_value, traceback)
        self.read2.__exit__(exception_type, exception_value, traceback)

    def __iter__(self):
        return self

    def __next__(self):
        if self.cnext() <= 0:
            raise StopIteration()
        return self

    @property
    def rlen(self):
        if self.cur > 0:
            return self.read2.rlen
        else:
            return self.read1.rlen

    @property
    def seq(self):
        if self.cur > 0:
            return self.read2.seq
        else:
            return self.read1.seq
    @property
    def char_seq(self):
        if self.cur > 0:
            return self.read2.char_seq
        else:
            return self.read1.char_seq

    @property
    def pos(self):
        return long(self.read1.pos + self.read2.pos)

    @property
    def size(self):
        return long(self.read1.size + self.read2.size)

    cdef int cnext(self) nogil:
        self.cur = self.cur ^ 1
        if self.cur > 0:
            return self.read2.cnext()
        else:
            return self.read1.cnext()

    @cython.boundscheck(False)
    cdef public uint8_t[:] get_seq(self) nogil:
        if self.cur > 0:
            return self.read2.get_seq()
        else:
            return self.read1.get_seq()

    cdef public int get_len(self) nogil:
        if self.cur > 0:
            return self.read2.get_len()
        else:
            return self.read1.get_len()

    cdef public int get_flags(self) nogil:
        if self.cur > 0:
            return self.flags2
        else:
            return self.flags1

cdef class FastQWriter:
    """
    FastQ File Writer
    """
    def __cinit__(self):
        self.proc = None  # for gzip process
        self.file = None
        self.fstream = NULL
        self.filename = None

    def __init__(self, object fileobj):
        """
        Arguments:
            filename:  name of file to write to
        """
        if hasattr(fileobj, 'fileno'):
            self.file = fileobj
        else:
            self.filename = fileobj

    def __enter__(self):
        if self.filename is not None:
            self.file = open(self.filename, "wb")
        if self.filename.endswith(".gz"):
            try:
                self.proc = Popen(["pigz", "-c"],
                                  stdin=PIPE,
                                  stdout=self.file,
                                  shell=False)
            except FileNotFoundError:
                self.proc = Popen(["gzip", "-c"],
                                  stdin=PIPE,
                                  stdout=self.file,
                                  shell=False)
            self.fstream = fdopen(self.proc.stdin.fileno(), "wb")
        else:
            self.fstream = fdopen(self.file.fileno(), "wb")
        if not self.fstream:
            raise OSError(errno, "Failed to open", self.filename)

        return self

    def __exit__(self, exception_type, exception_value, traceback):
        if self.fstream:
            fclose(self.fstream)
        if self.proc:
            self.proc.wait(5)
        if self.filename is not None:
            try:
                self.file.close()
            except OSError:
                pass

    cdef public int write_from_FastQFile(self, FastQFile f) nogil:
        cdef size_t res
        for i in range(4):
            res = fwrite(f.buf[i], 1, f.len[i], self.fstream)
            if res != f.len[i]:
                return False
        return True

    cdef public int write_from_AlignmentFile(self, AlignmentFile bam) nogil:
        cdef int l = bam.b.core.l_qseq
        #cdef char* name = bam_get_qname(bam.b)

    def write(self, read):
        if isinstance(read, FastQFile):
            return self.write_from_FastQFile(read)
        else:
            raise Exception("writing {} not implemented".format(type(read)))
