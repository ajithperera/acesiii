
c Single-line parser of Z-matrix elements ONLY.

c INPUT
c    zline(80) : the Z-matrix line to parse

c OUTPUT
c    izl(2,7) : the beginning (1,*) and ending (2,*) indices of each
c               element (non-existent elements are at zline(0:0))
c               EXAMPLE:
c                  zline = '  X  1  R3  2  A2  3  D1   '
c                           123456789+123456789+1234567
c                  izl(1:2,1:7) = 3  6   9  13  16  20  23
c                                 3  6  10  13  17  20  24




      subroutine parsez(zline,izl)
      implicit none

c     Maximum string length of terminal lines
      INTEGER LINELEN
      PARAMETER (LINELEN=80)

      character*(linelen) zline
      integer izl(2,7)

      character*1 czTmp
      integer ndx, count
      logical find_char, not_done

c "parameters"
      character*1 achar, czTab, czSpace, czComment
      integer    max_cols
      parameter (max_cols = 80)

c ----------------------------------------------------------------------

      czTab     = achar(9)
      czSpace   = achar(32)
      czComment = achar(35)

      izl(1,1) = 0
      izl(2,1) = 0
      izl(1,2) = 0
      izl(2,2) = 0
      izl(1,3) = 0
      izl(2,3) = 0
      izl(1,4) = 0
      izl(2,4) = 0
      izl(1,5) = 0
      izl(2,5) = 0
      izl(1,6) = 0
      izl(2,6) = 0
      izl(1,7) = 0
      izl(2,7) = 0

c   o the zline pointer
      ndx = 1

c   o the izl pointer
      count = 1

c   o start looking for a char
      find_char = .true.

      not_done = .true.
      do while (not_done)

         czTmp = zline(ndx:ndx)
         if (find_char) then
            if ( (czTmp.ne.czSpace) .and.
     &           (czTmp.ne.czTab  )       ) then
               if (czTmp.eq.czComment) then
                  not_done = .false.
               else
                  izl(1,count) = ndx
                  find_char = .false.
               end if
            end if
         else
            if ( (czTmp.eq.czSpace).or.
     &           (czTmp.eq.czTab  )    ) then
               izl(2,count) = ndx - 1
               if (count.eq.7) then
                  return
               else
                  find_char = .true.
                  count = count + 1
               end if
            else
               if (czTmp.eq.czComment) then
                  izl(2,count) = ndx - 1
                  not_done = .false.
               end if
            end if
         end if

         if (not_done.and.(ndx.eq.max_cols)) then
            not_done = .false.
            if (.not.find_char) izl(2,count) = ndx
         else
            ndx = ndx + 1
         end if

c     end do while (not_done)
      end do

c      write(*,*) (izl(1,count),count=1,7)
c      write(*,*) (izl(2,count),count=1,7)

      return
      end

