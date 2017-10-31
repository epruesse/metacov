from pysam.libchtslib cimport *
import numpy as np
cimport numpy as np
from cpython.array cimport array


cdef class ReadProcessor(object):
    """Base class for read stats accumulators"""

    cpdef void set_max_readlen(self, int rlen)
    cdef:
        public void process_read(self, int rlen, uint8_t[:] read, int flags) nogil


cdef class ReadProcessorList(ReadProcessor):
    """Base class for read stats accumulators"""
    cdef:
        list processors

    cpdef void set_max_readlen(self, int rlen)
    cdef:
        public void process_read(self, int rlen, uint8_t[:] read, int flags) nogil


cdef class Flag:
    cdef public:
        int flag
        str name_true
        str name_false
        str name_col


cdef class ByFlag(ReadProcessorList):
    cdef:
        int nflags
        list flags
        array flag_ints_arr
        uint32_t[:] flag_ints

    cpdef void set_max_readlen(self, int rlen)
    cdef:
        public void process_read(self, int rlen, uint8_t[:] read, int flags) nogil


cdef class BaseHist(ReadProcessor):
    cdef:
        np.ndarray _counts_data
        uint32_t[:,:] _counts
        public void process_read(self, int rlen, uint8_t[:] read, int flags) nogil

    cpdef void set_max_readlen(self, int rlen)


cdef class KmerHist(ReadProcessor):
    cdef:
        np.ndarray _counts_data
        uint32_t[:,:] _counts
        int K
        int NK
        int OFFSET
        int STEP

        public void process_read(self, int rlen, uint8_t[:] read, int flags) nogil

    cpdef void set_max_readlen(self, int rlen)


"""
  - insert size
  - length
  - fwd/rev?
  - quality per base
  - mapq

"""
