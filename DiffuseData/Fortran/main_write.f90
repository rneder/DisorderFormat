program main_write
!-
! Writes/reads an example of 3D-Diffuse scattering into an HDF5 file format
!
! The data consist of a simple 3D grid with corners:
!                       h     k     l
! lower left  bottom:  -3.0, -3.0, -3.0   in  eck(:,1)
! lower right bottom:   3.0, -3.0, -3.0   in  eck(:,2)
! upper left  bottom:  -3.0,  3.0, -3.0   in  eck(:,3)
! upper right bottom:   3.0,  3.0, -3.0
! lower left  top   :  -3.0, -3.0, -3.0   in  eck(:,4)
! lower right top   :   3.0, -3.0,  3.0
! upper left  top   :  -3.0,  3.0,  3.0
! upper right top   :   3.0,  3.0,  3.0
! Increment from  lower left  bottom to lower right bottom:  0.10, 0.00, 0.00
! Increment from  lower left  bottom to upper right bottom:  0.00, 0.20, 0.00
! Increment from  lower left  bottom to lower left  top   :  0.00, 0.00, 0.25
! the array size is thus: abscissa:  61
!                         ordinate:  31
!                         top axis:  25
! The asymmetry is intended to make the directions sitinguishable
!
! The data has values of
!                   a 3DGaussian with sigma = 1.0 along each axis
!                   10.0 at each integer hkl
! 
!+
use top_data_mod
use gen_hdf_write_mod
use lib_hdf5_read_mod
!
use precision_mod         ! Data precision for real, integer, character
use prompt_mod
!
implicit none
!
integer :: i,j,k    ! Dummy indices
real(kind=PREC_DP), dimension(3) :: h   ! Dummy vector
real(kind=PREC_DP) :: sigma
!
outfile = 'example.h5'
eck(:,1) = (/ -3.0_PREC_DP, -3.0_PREC_DP, -3.0_PREC_DP /)    ! Lower left bottom
eck(:,2) = (/  3.0_PREC_DP, -3.0_PREC_DP, -3.0_PREC_DP /)    ! lower right bottom
eck(:,3) = (/ -3.0_PREC_DP,  3.0_PREC_DP, -3.0_PREC_DP /)    ! Upper left bottom
eck(:,4) = (/ -3.0_PREC_DP, -3.0_PREC_DP,  3.0_PREC_DP /)    ! Lower left top
!
vi(:,1)  = (/  0.10_PREC_DP, 0.00_PREC_DP, 0.00_PREC_DP /)   ! Steps along abscissa
vi(:,2)  = (/  0.00_PREC_DP, 0.20_PREC_DP, 0.00_PREC_DP /)   ! Steps along ordinate
vi(:,3)  = (/  0.00_PREC_DP, 0.00_PREC_DP, 0.25_PREC_DP /)   ! Steps along top
!
out_inc  = (/  61, 31, 25 /)
extr_abs = 1    ! Allows for exchange of axes, standard should be as typed
extr_ord = 2
extr_top = 3
!
value = VAL_3DPDF     ! Direct space data
!
cr_a0  = (/ 90.00_PREC_DP,90.00_PREC_DP,90.00_PREC_DP /)
cr_win = (/ 90.00_PREC_DP,90.00_PREC_DP,90.00_PREC_DP /)
!
valmax = 0.0_PREC_DP  ! Allows to scale data to this maximum value if > 0.0
!
sigma = 1.0_PREC_DP
!
allocate(qvalues(out_inc(1), out_inc(2), out_inc(3)))
qvalues = 0.0_PREC_DP
!
do k=1, out_inc(3)
   do j=1, out_inc(2)
      do i=1, out_inc(1)
         h(1) = eck(1,1) + (i-1)*vi(1,1) + (j-1)*vi(1,2) + (k-1)*vi(1,3)
         h(2) = eck(2,1) + (i-1)*vi(2,1) + (j-1)*vi(2,2) + (k-1)*vi(2,3)
         h(3) = eck(3,1) + (i-1)*vi(3,1) + (j-1)*vi(3,2) + (k-1)*vi(3,3)
         qvalues(i,j,k) = exp(-0.5_PREC_DP*(h(1)/sigma)**2) *  &
                          exp(-0.5_PREC_DP*(h(2)/sigma)**2) *  &
                          exp(-0.5_PREC_DP*(h(3)/sigma)**2)
      enddo
   enddo
enddo
!
do k=1, out_inc(3), 4
   do j=1, out_inc(2), 5 
      do i=1, out_inc(1), 10
         qvalues(i,j,k) = 10.0_PREC_DP
      enddo
   enddo
enddo
!
call gen_hdf5_write (value, laver, outfile, out_inc, eck, vi, &
                       extr_abs, extr_ord, extr_top,                       &
                       cr_a0, cr_win, qvalues, VAL_PDF, VAL_3DPDF, valmax, &
                       ier_num, ier_typ, ER_IO, ER_APPL)
!
write(output_io, '(a,a)') ' Wrote data to file: ', outfile(1:len_trim(outfile))
write(*,'(a,3i5)')   ' Aray size ', out_inc
write(*,'(a,3f9.3)') ' LeftLowerBottom  ', eck(:,1)
write(*,'(a,3f9.3)') 'RightLowerBottom  ', eck(:,2)
write(*,'(a,3f9.3)') ' LeftUpperBottom  ', eck(:,3)
write(*,'(a,3f9.3)') ' LeftLowerTop     ', eck(:,4)
write(*,'(a,3f9.3)') ' Steps abscissa   ', vi (:,1)
write(*,'(a,3f9.3)') ' Steps ordinate   ', vi (:,2)
write(*,'(a,3f9.3)') ' Steps top        ', vi (:,3)
write(*,'(a,i5)')    ' Direct=1, rec=0  ', value
write(*,'(a,3f9.3)') ' Data min/max     ', minval(qvalues), maxval(qvalues)
write(*,*)
!
! Clear all arrays and read file 
!
eck = 0.0_PREC_DP
vi  = 0.0_PREC_DP
out_inc = 0
cr_a0  = 0.0_PREC_DP
cr_win = 0.0_PREC_DP
qvalues = 0.0_PREC_DP
!
write(*,*) ' READ DATA '
call hdf5_read(outfile, out_inc, eck, vi, cr_a0, cr_win, value, qvalues)
write(*,*) ' GOT  DATA '
write(*,'(a,3i5)')   ' Aray size ', out_inc
write(*,'(a,3f9.3)') ' LeftLowerBottom  ', eck(:,1)
write(*,'(a,3f9.3)') 'RightLowerBottom  ', eck(:,2)
write(*,'(a,3f9.3)') ' LeftUpperBottom  ', eck(:,3)
write(*,'(a,3f9.3)') ' LeftLowerTop     ', eck(:,4)
write(*,'(a,3f9.3)') ' Steps abscissa   ', vi (:,1)
write(*,'(a,3f9.3)') ' Steps ordinate   ', vi (:,2)
write(*,'(a,3f9.3)') ' Steps top        ', vi (:,3)
write(*,'(a,i5)')    ' Direct=1, rec=0  ', value
write(*,'(a,3f9.3)') ' Data min/max     ', minval(qvalues), maxval(qvalues)
!
deallocate(qvalues)
!
end program main_write
