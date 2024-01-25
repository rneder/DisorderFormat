module top_data_mod
!-
!  Summary of all(?) relevant variables
!+
!
use precision_mod         ! Data precision for real, integer, character
!
implicit none
!
integer                            :: value      ! The data are PDF or intensities
LOGICAL :: laver
character(LEN=200)                 :: outfile    ! Output file name
integer           , dimension(3)   :: out_inc
real(kind=PREC_DP), dimension(3,4) :: eck ! (3,4) Corners (HKL, number)
real(kind=PREC_DP), dimension(3,3) :: vi  ! Increment vectors (HKL, number)
integer                            :: extr_abs
integer                            :: extr_ord
integer                            :: extr_top
real(kind=PREC_DP), dimension(3)   :: cr_a0      ! Direct space lattice parameters
real(kind=PREC_DP), dimension(3)   :: cr_win     ! Direct space angles
real(kind=PREC_DP), dimension(:,:,:), allocatable :: qvalues   ! Actual data values
integer, parameter      :: VAL_PDF    = 14
integer, parameter      :: VAL_3DPDF  = 15 
real(kind=PREC_DP)      :: valmax
integer                 :: ier_num     ! An error number
integer                 :: ier_typ     ! An error type
integer, parameter      :: ER_IO   = 2 ! Input/output error
integer, parameter      :: ER_APPL = 6 ! Error in the user program
!
end module top_data_mod
