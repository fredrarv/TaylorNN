      subroutine gauss_odf_wrap(orifile, ckofile, status)
      implicit none
      character*255 orifile, ckofile
      integer status
C-----------------------------------------------------------------------
C     A wrapper function for gauss_odf() with no unit number argument.
C     On return, 'status' contains one of the following error codes:
C       0     success
C       1     cannot open orientation file
C-----------------------------------------------------------------------
      integer uin, n, m
      character*1 nul

      status = 0

C     Accept both NULL-terminated strings and fortran strings
      nul = char(0)
      n = index(orifile, nul) - 1
      m = index(ckofile, nul) - 1
      if (n.lt.0) n = len_trim(orifile)
      if (m.lt.0) m = len_trim(ckofile)

      uin = 26
      open(unit=uin, file=orifile(1:n), status='old', err=1)
      call gauss_odf(uin, orifile(1:n), ckofile(1:m))

      return

 1    continue
      status = 1
      return

      end
      
