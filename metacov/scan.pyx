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


## Not Fast:

cdef str bamseq2str(char* seq, int lseq):
    "Convert BAM nt16 encoded base string to str"
    return "".join([chr(nt16_to_ascii(bam_seqi(seq,i))) for i in range(lseq)])

cpdef str kmer_base2_to_ascii(uint32_t kmer, int l):
    cdef int n
    return "".join([chr(nt4_to_ascii(kmer>>n & 3)) for n in range(0, 2*l, 2) ])


## Counter Classes

cdef class ReadProcessor:
   """Base class for read stats accumulators"""
   cpdef public void set_max_readlen(self, int rlen):
       pass
   
   cdef public void process_read(self, int rlen, uint8_t[:] read, int flags) nogil:
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

   cdef public void process_read(self, int rlen, uint8_t[:] read, int flags) nogil:
       with gil:
           for processor in self.processors:
               (<ReadProcessor>processor).process_read(rlen, read, flags)

   def get_rows(self, int i):
       return self.processors[i].get_rows()

cdef class Flag:
    def __init__(self, int flag, str name_true, str name_false, str name_col):
        self.flag = flag
        self.name_true = name_true
        self.name_false = name_false
        self.name_col = name_col

FLAG_PAIRED = Flag(BAM_FPAIRED, "Paired", "Unpaired", "Paired")
FLAG_PROPER_PAIR = Flag(BAM_FPROPER_PAIR, "Paired", "Unpaired", "ProperPair")
FLAG_MAPPED = Flag(BAM_FUNMAP, "Unmapped", "Mapped", "Mapped")
FLAG_READDIR = Flag(BAM_FREAD2, "R2", "R1", "R")
FLAG_SECONDARY = Flag(BAM_FSECONDARY, "Secondary", "Primary", "Alignment")
FLAG_DUP = Flag(BAM_FDUP, "Duplicate", "Singleton", "Duplicate")
               
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
    cdef public void process_read(self, int rlen, uint8_t[:] read, int flags) nogil:
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
        
            (<ReadProcessor>self.processors[n]).process_read(rlen, read, flags)
            
            
               
cdef class BaseHist(ReadProcessor):
    def __cinit__(self):
        self._counts_data = np.zeros((10, 5), dtype=np.uint32)
        self._counts = self._counts_data

    def __copy__(self):
        return BaseHist()

    cpdef public void set_max_readlen(self, int rlen):
        self._counts_data.resize((rlen, 5), refcheck=False)
        self._counts = self._counts_data

    @cython.boundscheck(False)
    cdef public void process_read(self, int rlen, uint8_t[:] read, int flags) nogil:
        cdef int i
        for i in range(rlen):
            self._counts[i, read[i]] += 1

    @property
    def counts(self):
        return self._counts_data

    def get_rows(self):
        cdef int i
        yield ["Pos", "A", "G", "C", "T", "N"]
        for i in range(self._counts.shape[0]):
            yield [i+1] + list(self._counts[i])


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
    cdef public void process_read(self, int rlen, uint8_t[:] read, int flags) nogil:
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
            yield [kmer_base2_to_ascii(i, self.K)] + list(self.counts[i])



            
## scanning algos

cdef class ReadIterator:
    cdef:
        int max_readlen
        array rseq_arr
        uint8_t[:] rseq
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
    
    cdef int get_len(self) nogil:
        return 0

    cdef int get_flags(self) nogil:
        return 0


cdef class AlignmentFileIterator(ReadIterator):
    cdef:
       AlignmentFile bam
       IteratorRowAll row

    def __cinit__(self, AlignmentFile bam):
        self.bam = bam
        self.row = IteratorRowAll(bam)

    cdef int cnext(self) nogil:
        with gil:
            return self.row.cnext()

    @cython.boundscheck(False)
    cdef uint8_t[:] get_seq(self) nogil:
        cdef:
           char* rqseq = bam_get_seq(self.row.b)
           int rlen = self.row.b.core.l_qseq
           int i

        if rlen > self.max_readlen:
            with gil:
                self.set_max_readlen(rlen)
           
        for i in range(rlen):
            self.rseq[i] = nt16_to_nt4(bam_seqi(rqseq, i))

        return self.rseq

    cdef int get_len(self) nogil:
        return self.row.b.core.l_qseq

    cdef int get_flags(self) nogil:
        return self.row.b.core.flag


cdef class FastQFileIterator(ReadIterator):
    cdef:
        FastQFile fq
    def __cinit__(self, FastQFile fq):
        self.fq = fq

    cdef int cnext(self) nogil:
        return self.fq.cnext()

    cdef uint8_t[:] get_seq(self) nogil:
        return self.fq.get_seq()

    cdef int get_len(self) nogil:
        return self.fq.get_len()

    cdef int get_flags(self) nogil:
        return self.fq.get_flags()


cpdef long scan_reads(
    object infile, object counters,
    int progress_interval=10000000, object progress_cb=None,
    maxreads=0):
    cdef:
        ReadIterator it
        ReadProcessor proc

    if isinstance(infile, AlignmentFile):
        it = AlignmentFileIterator(infile)
    elif isinstance(infile, FastQFile):
        it = FastQFileIterator(infile)
    else:
        raise Exception("meh")

    if isinstance(counters, list):
        proc = ReadProcessorList(counters)
    elif isinstance(counters, ReadProcessor):
        proc = counters
    else:
        raise Exception("mah")
        

    return scan(it, proc, progress_interval, progress_cb, maxreads)


cpdef long scan(
    ReadIterator it, ReadProcessor processor,
    int progress_interval=10000000, object progress_cb=None,
    maxreads=0):
    cdef:
        long readno = 0
        int rlen
        uint8_t[:] rseq
        int rflags
        int max_readlen = 50

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

        processor.process_read(rlen, rseq, rflags)

        if maxreads and readno >= maxreads:
            break
        
    return readno
        

cpdef scan_reads_bam(
    AlignmentFile bam, counters,
    int progress_interval=10000000, object progress_cb=None,
    maxreads=0):

    cdef:
        IteratorRowAll row = IteratorRowAll(bam)
        int readno = 0
        int rlen, i
        int max_rlen = 50
        char *rqseq
        np.ndarray rseq_arr = np.zeros((200), dtype=np.uint8)
        uint8_t[:] rseq = rseq_arr

    for counter in counters:
            (<ReadProcessor>counter).set_max_readlen(max_rlen)


    while row.cnext() > 0:
        readno += 1
        if readno % progress_interval == 0:
            if progress_cb:
                progress_cb()

        rlen = row.b.core.l_qseq
        rqseq = bam_get_seq(row.b)

        if rlen > max_rlen:
            max_rlen = rlen
            for counter in counters:
                (<ReadProcessor>counter).set_max_readlen(rlen)
            
        for i in range(rlen):
            rseq[i] = nt16_to_nt4(bam_seqi(rqseq, i))

        for counter in counters:
            (<ReadProcessor>counter).process_read(rlen, rseq, 0)

        if maxreads and readno >= maxreads:
            break

        
cpdef scan_reads_fq(
    FastQFile fq, counters,
    int progress_interval=10000000, object progress_cb=None,
    maxreads=0):
    cdef:
        int readno = 0
        int rlen
        int max_rlen = 50
        uint8_t[:] rseq

    for counter in counters:
        (<ReadProcessor>counter).set_max_readlen(max_rlen)

    while fq.cnext() > 0:
        readno += 1
        if readno % progress_interval == 0:
            if progress_cb:
                progress_cb()

        rlen = fq.get_len()
        rseq = fq.get_seq()
        
        if rlen > max_rlen:
            max_rlen = rlen
            for counter in counters:
                (<ReadProcessor>counter).set_max_readlen(rlen)

        for counter in counters:
            (<ReadProcessor>counter).process_read(rlen, rseq, 0)

        if maxreads and readno >= maxreads:
            break

    
    


