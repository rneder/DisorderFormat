module error_mod
!
!  Handle error messages
!
integer, parameter :: NMSG = 3   ! Number of lines in ier_msg
!
!*******************************************************************************
!
contains
!
!*******************************************************************************
!
subroutine error_message(ier_num, ier_msg)
!-
! Print error messages
!+
integer                              , intent(out) :: ier_num
character(len=*)  , dimension(NMSG)  , intent(out) :: ier_msg
!
integer, parameter :: IU = -4
integer, parameter :: IO =  0
!
integer :: i,j
character(len=45), dimension(IU: IO) :: error
!
data error(IU:IO) / &
 'Data set has wrong dimensions               ', & ! -4
 'Data set has wrong rank                     ', & ! -3
 'Data set does not exist in HDF5 file        ', & ! -2
 'File not found                              ', & ! -1
 'No error                                    '  & !  0
 /
!
write(*,'(a5, a45, a5,i4)') ' *** ',error(ier_num),' *** ', ier_num
do i  = 1, NMSG
  if(ier_msg(i)/= ' ') then
     j = min(len(ier_msg), 45)
     write(*,'(a5, a45, a5,i4)') ' *** ',ier_msg(i)(1:j), ' *** ', ier_num
  endif
enddo
!
end subroutine error_message
!
!*******************************************************************************
!
end module error_mod
