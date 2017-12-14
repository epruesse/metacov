import numpy as np
cimport numpy as np

from scan cimport *
from pyfq cimport FastQFile

cimport cython
from cpython.array cimport array, clone, resize
from cython.parallel import parallel, prange
from pysam.libchtslib cimport *

from pysam.libcalignmentfile cimport AlignmentFile, IteratorRowRegion, IteratorRowAll
from pysam.libcalignedsegment cimport pysam_bam_get_cigar, pysam_bam_get_seq, \
    pysam_bam_get_qual, PileupColumn
from pysam.libcfaidx cimport FastaFile

from copy import copy
from itertools import chain, islice

## Fast:

cdef char* NT4_ACGTN_table = 'ACGTN'
cdef inline char nt4_to_ascii(uint8_t n) nogil:
    "Convert NT4 coded base to ASCII"
    return NT4_ACGTN_table[n]

cdef uint8_t *nt16_to_nt4_table = [ 4, 0, 1, 4,  2, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4 ]
cdef inline uint8_t nt16_to_nt4(unsigned char n) nogil:
    "Convert NT16 (bits) coded base to NT4"
    return nt16_to_nt4_table[n]

cdef char* NT16_ACGTN_table = 'NACNGNNNTNNNNNNN'
cdef inline char nt16_to_ascii(uint8_t n) nogil:
    return NT16_ACGTN_table[n]

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


cdef inline uint8_t nt4_comp(uint8_t n) nogil:
    return 3-n + ((3-n & 4)>>2)*5

## Not Fast:

cdef str bamseq2str(char* seq, int lseq):
    "Convert BAM nt16 encoded base string to str"
    return "".join([chr(nt16_to_ascii(bam_seqi(seq,i))) for i in range(lseq)])

cpdef str kmer_base2_to_ascii(uint32_t kmer, int l):
    cdef int n
    return "".join([chr(nt4_to_ascii(kmer>>n & 3)) for n in range(0, 2*l, 2) ])

## Flags

cdef class Flag:
    def __cinit__(self, int flag, str name_true, str name_false, str name_col):
        self.flag = flag
        self.name_true = name_true
        self.name_false = name_false
        self.name_col = name_col


cdef:
    # 0x01 template having multiple segments
    Flag _FLAG_PAIRED = Flag(BAM_FPAIRED, "Paired", "Unpaired", "Paired")
    # 0x02 each segment properly aligned
    Flag _FLAG_PROPER_PAIR = Flag(BAM_FPROPER_PAIR, "Paired", "Unpaired", "ProperPair")
    # 0x04 segment unmapped
    Flag _FLAG_MAPPED = Flag(BAM_FUNMAP, "Unmapped", "Mapped", "Mapped")
    # 0x08 next segment unmapped
    Flag _FLAG_MMAPPED = Flag(BAM_FMUNMAP, "Unmapped", "Mapped", "MateMapped")
    # 0x10 read mapped to reverse strand
    Flag _FLAG_REVERSE = Flag(BAM_FREVERSE, "Reverse", "Forward", "Readdir")
    # 0x20 next segment mapped to reverse strand
    Flag _FLAG_MREVERSE = Flag(BAM_FMREVERSE, "Reverse", "Forward", "MateReaddir")
    # 0x40 first segment in template (read1)
    Flag _FLAG_READ1 = Flag(BAM_FREAD1, "R1", "", "IsRead1")
    # 0x80 last segment in template (read2)
    Flag _FLAG_READDIR = Flag(BAM_FREAD2, "R2", "R1", "R")
    # 0x100 secondary alignment
    Flag _FLAG_SECONDARY = Flag(BAM_FSECONDARY, "Secondary", "Primary", "Alignment")
    # 0x200 qc failure
    Flag _FLAG_QCFAIL = Flag(BAM_FQCFAIL, "Fail", "Pass", "QC")
    # 0x400 optical or PCR duplicate
    Flag _FLAG_DUP = Flag(BAM_FDUP, "Duplicate", "Singleton", "Duplicate")


FLAG_PAIRED = _FLAG_PAIRED
FLAG_PROPER_PAIR = _FLAG_PROPER_PAIR
FLAG_MAPPED = _FLAG_MAPPED
FLAG_MMAPPED = _FLAG_MMAPPED
FLAG_REVERSE = _FLAG_REVERSE
FLAG_MREVERSE = _FLAG_MREVERSE
FLAG_READ1 = _FLAG_READ1
FLAG_READDIR = _FLAG_READDIR
FLAG_SECONDARY = _FLAG_SECONDARY
FLAG_QCFAIL = _FLAG_QCFAIL
FLAG_DUP = _FLAG_DUP


## Iterating over reads

cdef class ReadIterator:
    """
    Abstract Base Class for iterating over reads
    """
    def __cinit__(self):
        self.max_readlen = 50
        self.rseq_arr = clone(array('B'), self.max_readlen, False)
        self.rseq = self.rseq_arr

    cdef void set_max_readlen(self, int max_readlen):
        self.max_readlen = max_readlen
        resize(self.rseq_arr, max_readlen)
        self.rseq = self.rseq_arr

    cdef int cnext(self) nogil:
        return -1

    cdef uint8_t[:] get_seq(self) nogil:
        return None

    cdef uint8_t[:] get_ref(self) nogil:
        return None

    cdef int get_len(self) nogil:
        return 0

    cdef int get_flags(self) nogil:
        return 0

    cdef int get_pos(self) nogil:
        return 0

    cdef int get_isize(self) nogil:
        return 0

    cdef char* get_name(self) nogil:
        return ""

    cdef int get_tid(self) nogil:
        return 0

    cdef char* get_rname(self) nogil:
        return ""

    cdef int is_reverse(self) nogil:
        return 0

    cdef int is_paired(self) nogil:
        return 0


cdef class AlignmentFileIterator(ReadIterator):
    """
    Iterates over all reads of a BAM/SAM/CRAM file
    """
    cdef:
       AlignmentFile bam
       IteratorRowAll row
       FastaFile fasta
       int tid
       array curseq_data
       uint8_t[:] curseq
       int curseq_len

    def __cinit__(self, AlignmentFile bam, FastaFile fasta):
        # super __cinit__ called automatically
        self.bam = bam
        self.row = IteratorRowAll(bam)
        self.fasta = fasta
        self.tid = -1
        self.curseq_data = clone(array('B'), 1, False)
        self.curseq = self.curseq_data
        self.curseq_len = 0

    @cython.boundscheck(False)
    cdef int cnext(self) nogil:
        cdef int res, length
        cdef char* seq
        with gil:
            res = self.row.cnext()

        # if we have fasta and are switching reference templates,
        # load the new current sequence
        if self.fasta is not None and self.tid != self.get_tid():
            self.tid = self.get_tid()
            length = faidx_seq_len(self.fasta.fastafile, self.get_rname())
            if length > 0:
                seq = faidx_fetch_seq(self.fasta.fastafile, self.get_rname(),
                                      0, length, &length)
                with gil:
                    resize(self.curseq_data, length)
                self.curseq = self.curseq_data
                for i in range(length):
                    self.curseq[i] = iupac_to_nt4(seq[i])
                self.curseq_len = length
            else:
                self.curseq_len = 0

        return res

    cdef uint8_t[:] get_ref(self) nogil:
        return self.curseq

    @cython.boundscheck(False)
    cdef uint8_t[:] get_seq(self) nogil:
        cdef:
           char* rqseq = bam_get_seq(self.row.b)
           int rlen = self.row.b.core.l_qseq
           int i

        if rlen > self.max_readlen:
            with gil:
                self.set_max_readlen(rlen)

        if self.is_reverse():
            # undo reverse complementing of read by mapper
            for i in range(rlen):
                self.rseq[rlen-i-1] = nt4_comp(nt16_to_nt4(bam_seqi(rqseq, i)))
        else:
            for i in range(rlen):
                self.rseq[i] = nt16_to_nt4(bam_seqi(rqseq, i))

        return self.rseq

    cdef int get_len(self) nogil:
        return self.row.b.core.l_qseq

    cdef int get_flags(self) nogil:
        return self.row.b.core.flag

    cdef int get_isize(self) nogil:
        if self.is_paired():
            return self.row.b.core.isize
        else:
            return 0

    cdef int get_pos(self) nogil:
        if self.get_flags() & _FLAG_REVERSE.flag == 0:
            return self.row.b.core.pos
        else:
            return self.row.b.core.pos + self.row.b.core.l_qseq

    cdef char* get_name(self) nogil:
        return <char*>self.row.b.data

    cdef int get_tid(self) nogil:
        return self.row.b.core.tid

    cdef char* get_rname(self) nogil:
        if self.get_tid() < 0:
            return ""
        return self.bam.header.target_name[self.get_tid()]

    cdef int is_reverse(self) nogil:
        return self.row.b.core.flag & BAM_FREVERSE

    cdef int is_paired(self) nogil:
        return self.row.b.core.flag & BAM_FPROPER_PAIR

cdef class FastQFileIterator(ReadIterator):
    """
    Iterates over all reads of a fq/fg.gz file
    """
    cdef:
        FastQFile fq

    def __cinit__(self, FastQFile fq):
        self.fq = fq

    cdef int cnext(self) nogil:
        return self.fq.cnext()

    cdef uint8_t[:] get_seq(self) nogil:
        return self.fq.get_seq()

    cdef uint8_t[:] get_ref(self) nogil:
        return None

    cdef int get_len(self) nogil:
        return self.fq.get_len()

    cdef int get_flags(self) nogil:
        return self.fq.get_flags()

    cdef int get_isize(self) nogil:
        return -1

    cdef int get_pos(self) nogil:
        return -1

    cdef char* get_name(self) nogil:
        return ""

    cdef int get_tid(self) nogil:
        return -1

    cdef char* get_rname(self) nogil:
        return ""

    cdef int is_reverse(self) nogil:
        return 0 # FIXME

    cdef int is_paired(self) nogil:
        return 0 # FIXME


## Counter Classes

cdef class ReadProcessor:
   """Base class for read stats accumulators"""
   cpdef public void set_max_readlen(self, int rlen):
       pass
   
   cdef public void process_read(self, int rlen, uint8_t[:] read, int flags,
                                 ReadIterator it) nogil:
       with gil: print("ArgHHH")


cdef class ReadProcessorList(ReadProcessor):
   """Container class multiplexing over several read stats accumulators"""
   def __init__(self, list processors):
       self.processors = processors
       for processor in processors:
           assert isinstance(processor, ReadProcessor)

   def __copy__(self):
       return ReadProcessorList([copy(processor) for processor in self.processors])

   cpdef public void set_max_readlen(self, int rlen):
       for processor in self.processors:
           (<ReadProcessor>processor).set_max_readlen(rlen)

   cdef public void process_read(self, int rlen, uint8_t[:] read, int flags,
                                 ReadIterator it) nogil:
       with gil:
           for processor in self.processors:
               (<ReadProcessor>processor).process_read(rlen, read, flags, it)

   def get_rows(self, int i):
       return self.processors[i].get_rows()



cdef class ByFlag(ReadProcessorList):
    def __init__(self, object processor, list flags):
        if isinstance(processor, list):
            processor = ReadProcessorList(processor)
        assert isinstance(processor, ReadProcessor)
        
        self.nflags = len(flags)
        self.flags = flags
        #self.processors = [ copy.deepcopy(processor) for _ in range(2 ** self.nflags) ]
        self.processors = [ copy(processor) for _ in range(2 ** self.nflags) ]
        super().__init__(self.processors)

    def get_rows(self, int i):
        tag_head = [flag.name_col for flag in reversed(self.flags)]
        yield next(self.processors[0].get_rows(i)) + tag_head
        for n, processor in enumerate(self.processors):
            tag = []
            for m, flag in enumerate(reversed(self.flags)):
                if 1<<m & n:
                    tag.append(flag.name_true)
                else:
                    tag.append(flag.name_false)
            for l in islice(processor.get_rows(i), 1, None):
                yield l + tag
            

    @cython.boundscheck(False)
    cdef public void process_read(self, int rlen, uint8_t[:] read, int flags,
                                  ReadIterator it) nogil:
        cdef:
            int flag
            int i
            int n = 0
        with gil:
            for i in range(self.nflags):
                #n = (n<<1) | (<Flag>self.flags[i]).flag & flags > 0
                n <<= 1
                if (<Flag>self.flags[i]).flag & flags:
                    n += 1
        
            (<ReadProcessor>self.processors[n]).process_read(rlen, read, flags, it)
            

cdef class BaseHist(ReadProcessor):
    def __cinit__(self, int start_pos):
        """
        Computes a profile of base counts along the length of the reads

        Params:
          start_pos:  include this many bases from the reference ("left" of read)
        """
        self.start_pos = start_pos
        self._counts_data = np.zeros((10, 5), dtype=np.uint32)
        self._counts = self._counts_data

    def __copy__(self):
        return BaseHist(self.start_pos)

    cpdef public void set_max_readlen(self, int rlen):
        self._counts_data.resize((rlen+self.start_pos, 5), refcheck=False)
        self._counts = self._counts_data

    @cython.boundscheck(False)
    cdef public void process_read(self, int rlen, uint8_t[:] read, int flags,
                                  ReadIterator it) nogil:
        cdef int i, x = 0
        cdef int pos = it.get_pos()
        cdef uint8_t[:] ref = it.get_ref()
        cdef int mismatch = 0

        if pos < self.start_pos:
            return

        if it.is_reverse():
            # verify match
            for i in range(rlen):
                if read[i] != nt4_comp(ref[pos - i - 1]):
                    mismatch += 1
            if mismatch * 32 > rlen:
                return

            # count in reference
            for i in range(self.start_pos):
                self._counts[i, nt4_comp(ref[pos - i - 1 + self.start_pos])] += 1
        else:
            # verify match
            for i in range(rlen):
                if read[i] != ref[pos + i]:
                    mismatch += 1
            if mismatch * 32 > rlen:
                return

            # count in reference
            for i in range(self.start_pos):
                self._counts[i, ref[pos + i - self.start_pos]] += 1

        # count in read
        for i in range(rlen):
            self._counts[i+self.start_pos, read[i]] += 1

    @property
    def counts(self):
        return self._counts_data

    def get_rows(self):
        cdef int i
        yield ["Pos", "A", "G", "C", "T", "N"]
        for i in range(self._counts.shape[0]):
            yield [i-self.start_pos] + list(self._counts[i])


cdef class KmerHist(ReadProcessor):
    def __cinit__(self, int K, int NK, int STEP, int OFFSET):
        self.K = K
        self.NK = NK
        self.STEP = STEP
        self.OFFSET = OFFSET
        self._counts_data = np.zeros((4 ** K + 1, NK), dtype=np.uint32)
        self._counts = self._counts_data

    def __copy__(self):
        return KmerHist(self.K, self.NK, self.STEP, self.OFFSET)

    @cython.boundscheck(False)
    cdef public void process_read(self, int rlen, uint8_t[:] read, int flags,
                                  ReadIterator it) nogil:
        cdef:
            int i, j, k
            uint8_t c

        if rlen < self.OFFSET + self.STEP * self.NK:
            return

        for i in range(self.NK): # count self.NK kmers
            k = 0
            for j in range(self.K): # for each base of kmer
                c = read[self.OFFSET + i * self.STEP + j]
                if c > 3:
                    k = 4 ** self.K
                    break
                k = k | c << (2*j)
            
            self._counts[k, i] += 1

    @property
    def counts(self):
        return self._counts_data

    def get_rows(self):
        cdef int i
        yield ["kmer"] + ["n{}".format(i) for i in range(self.NK)]
        yield ['N'*self.K ] + list(self.counts[4**self.K])
        for i in range(4**self.K):
            yield [kmer_base2_to_ascii(i, self.K)] + list(self._counts_data[i])


cdef class MirrorHist(ReadProcessor):
    def __cinit__(self, int OFFSET=4, int N=10):
        self.OFFSET = OFFSET
        self.N = N
        self._counts_data = np.zeros((N+1,2), dtype=np.uint32)
        self._counts = self._counts_data

    def __copy__(self):
        return MirrorHist(self.OFFSET, self.N)

    @cython.boundscheck(False)
    cdef public void process_read(self, int rlen, uint8_t[:] read, int flags,
                                  ReadIterator it) nogil:
        cdef:
            int pos = it.get_pos() + self.OFFSET
            uint8_t[:] ref = it.get_ref()
            int mismatch_plain = 0
            int mismatch_comp = 0
            
            int i = 0

        if pos < self.N - self.OFFSET:
            return

        for i in range(self.N):
            if ref[pos + i + 1] != ref[pos - i - 1]:
                mismatch_plain += 1
            if ref[pos + i + 1] != nt4_comp(ref[pos - i - 1]):
                mismatch_comp += 1

        with gil:
            self._counts[mismatch_plain,0] += 1
            self._counts[mismatch_comp,1] += 1

    @property
    def counts(self):
        return self._counts_data

    def get_rows(self):
        cdef int i
        yield ["n","plain","comp"]
        for i in range(self.N):
            yield [i, self._counts_data[i,0], self._counts_data[i,1]]


cdef class IsizeHist(ReadProcessor):
    def __cinit__(self):
        self._counts_data = np.zeros(128, dtype=np.uint32)
        self._counts = self._counts_data
        self.max_isize = 0

    def __copy__(self):
        return IsizeHist()

    @cython.boundscheck(False)
    cdef public void process_read(self, int rlen, uint8_t[:] read, int flags,
                                  ReadIterator it) nogil:
        cdef:
           int isize = it.get_isize()
           int asize = self._counts.shape[0]

        if isize < 0:
            isize = -isize

        if isize > self.max_isize:
            self.max_isize = isize

        if asize <= isize:
            while asize <= isize:
                asize *= 2
            with gil:
                self._counts_data.resize(asize, refcheck=False)
                self._counts = self._counts_data

        self._counts[isize] += 1

    @property
    def counts(self):
        return self._counts_data

    def get_rows(self):
        cdef int i
        yield ["n","count"]
        for i in range(self.max_isize + 1):
            yield [i, self._counts_data[i]]


cpdef long scan_reads(
    object infile, object fasta, object counters,
    int progress_interval=10000000, object progress_cb=None,
    maxreads=0):
    cdef:
        ReadIterator it
        ReadProcessor processor

        long readno = 0
        int rlen
        uint8_t[:] rseq
        int rflags
        int max_readlen = 50

    if isinstance(infile, AlignmentFile):
        it = AlignmentFileIterator(infile, fasta)
    elif isinstance(infile, FastQFile):
        it = FastQFileIterator(infile)
    else:
        raise Exception("meh")

    if isinstance(counters, list):
        processor = ReadProcessorList(counters)
    elif isinstance(counters, ReadProcessor):
        processor = counters
    else:
        raise Exception("mah")
        
    processor.set_max_readlen(max_readlen)

    while it.cnext() > 0:
        readno += 1
        if readno % progress_interval == 0:
            if progress_cb:
                progress_cb()
        
        rlen = it.get_len()
        rseq = it.get_seq()
        rflags = it.get_flags()

        if rlen > max_readlen:
            max_readlen = rlen
            processor.set_max_readlen(max_readlen)

        processor.process_read(rlen, rseq, rflags, it)

        if maxreads and readno >= maxreads:
            break
        
    return readno
        

    
    



