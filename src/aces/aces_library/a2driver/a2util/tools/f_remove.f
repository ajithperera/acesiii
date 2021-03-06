
c This routine removes a file (or directory) without having to open it and
c then close it with status='DELETE'.

c WARNING:
c Fortran blank-pads all strings. This routine will simply append the
c null character to szFile and pass it to f_remove_core. If there are trailing
c blanks in the file name, then f_remove_core will not remove it and the process
c will die.



      subroutine f_remove(szFile)
      implicit none

c ARGUMENTS
      character*(*) szFile

c EXTERNAL FUNCTIONS
      integer f_remove_core
      character*1 achar

c INTERNAL VARIABLES
      integer iLength, iTmp
      character*(256) sz

c ----------------------------------------------------------------------


      iTmp = 0
c   o assert szFile fits in sz
      if (len(szFile).ge.256) then
         print *, '@F_REMOVE: Assertion failed.'
         print *, '   szFile  = "',szFile,'"'
         print *, '   len(sz) = ',256
         iTmp = 1
      end if
      if (iTmp.ne.0) call c_exit(iTmp)


      iLength = 1
      do while (szFile(iLength:iLength).ne.' '.and.
     &          iLength.le.len(szFile))
         iLength = iLength + 1
      end do
      iLength = iLength - 1
      if (iLength.eq.0) return

c ----------------------------------------------------------------------

      if (iLength.lt.256) then



         sz   = szFile(1:iLength)//achar(0)
         iTmp = f_remove_core(sz)

         if (iTmp.eq.0) return
         print *, '@F_REMOVE: The file "',szFile,
     &            '" could not be removed.'
         print *, '           error code = ',iTmp
      else
         print *, '@F_REMOVE: The sz buffer is too small ',
     &            'to contain the input string.'
         print *, '           Recompile with at least ',iLength+1,
     &            ' characters in the buffer.'
         print *, '           (Currently ',256,' characters.)'
      end if

      call c_exit(1)
      end

