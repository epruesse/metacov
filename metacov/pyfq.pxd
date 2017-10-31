from libc.stdio cimport FILE
from libc.stdint cimport uint8_t
from cpython.array cimport array

cdef class FastQFile(object):
   cdef:
      str filename
      object proc
      FILE *fstream
      int max_linelen
      char *buf[4]
      size_t buf_len[4]
      ssize_t len[4]
      long _file_size
      long _pos

      array rseq_arr
      uint8_t[:] rseq

      public void set_max_linelen(self, int max_linelen)
      public int cnext(self) nogil
      public uint8_t[:] get_seq(self) nogil
      public int get_len(self) nogil
      public int get_flags(self) nogil

cdef class FastQFilePair(FastQFile):
   cdef:
      FastQFile  read1, read2
      int flags1, flags2
      int cur

      public int cnext(self) nogil
      public uint8_t[:] get_seq(self) nogil
      public int get_len(self) nogil
      public int get_flags(self) nogil
