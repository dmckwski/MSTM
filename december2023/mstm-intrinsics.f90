      module intrinsics
!
! compiler-dependent intrinsic functions.
!
!
!  last revised: 15 January 2011
!
      implicit none
      contains
!
!   system clock
!
         real function mytime()
         implicit none
         real :: etime
         real(4), parameter :: x=0.
         real :: t(2)
!         mytime=secnds(x)
!         mytime=etime(t)
         call cpu_time(t(1))
         mytime=t(1)
         end function mytime

!   flush is one of the best named functions in fortran, although it does not do what I think it should.   This function flushes
!   the buffer to print unit i, so when the unit is flushed, all data written to the open file will appear in the file.
!   Not all compilers have this function.

         subroutine flush(i)
         implicit none
         integer :: i
         flush(i)
         end subroutine flush
!
!   number of command-line arguments.
!
         integer function mstm_nargs()
         implicit none
         integer nargs
!         mstm_nargs=nargs()
!         mstm_nargs=iargc()
         mstm_nargs=command_argument_count()
         end function mstm_nargs
!
!   command line argument retrieval
!
         subroutine mstm_getarg(char)
         implicit none
         integer :: istat
         character(*) :: char
!         call getarg(1,char,istat)
!         call getarg(1,char)
         call get_command_argument(1,char)
         end subroutine mstm_getarg
!
!  file positioning subroutine: uses gfortran intinsic fseek
!
         subroutine mstm_fseek(unit,position,whence,ierr)
         implicit none
         integer :: unit,position,ierr,whence,i
         character :: c*127

!         call fseek(unit,position,whence,ierr)
         rewind(unit)
         do i=1,position
            read(unit,'(a)',iostat=ierr) c
            if(ierr.ne.0) exit
         enddo
         end subroutine mstm_fseek
      end module intrinsics
