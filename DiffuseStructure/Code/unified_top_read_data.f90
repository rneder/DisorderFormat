program unified_top_read_data
!
!  Dummy top program that defines a crystal structure and uses the h5fortran library
!  to write a unified data file
!  
! Manually set up an example structure with all essential input. 
! Fixed example in analogy to test_write.py and the DISCUS macro
!
use error_mod
use lib_unified_chars_mod
use lib_nx_transfer_mod
use unified_read_mod
!
use precision_mod         ! Defines KIND of real numbers, strings etc
!
implicit none
!
!real(kind=prec_DP), parameter :: twopi = 2.0D0 * 3.1415926535897932384626433832795028841971693993751D0
character(len=PREC_STRING) ::  infile           !  Input file name
!character(len=PREC_STRING) :: program_version   ! Which main program wrote the structure
!character(len=PREC_STRING) :: author_name       ! Authors
!
real(kind=PREC_DP), dimension(3)   :: unit_cell_lengths    ! (a, b, c)
real(kind=PREC_DP), dimension(3)   :: unit_cell_angles     ! (alpha, beta, gamma)
!real(kind=PREC_DP), dimension(3,3) :: metric_tensor        ! Direct space metric tensor
!
character(len=32)          :: symmetry_H_M      ! Hermann-Mauguin Symbol
integer                    :: symmetry_origin   ! Origin choice 1 or 2
character(len=3)           :: symmetry_abc      ! For orthorhombic space groups abc, cba etc
integer                    :: symmetry_n_mat    ! Number of symmetry metrices
real(kind=PREC_DP), dimension(:,:,:), allocatable :: symmetry_mat         ! Direct space symmetry matrices
!
character(len=32)                                        :: data_type_experiment ! 'experimental';  'calculated'
character(len=32)                                        :: data_type_style      ! 'powder_diffraction', 'powder_pdf', 'single_diffraction', 'single_pdf' ...
character(len=32)                                        :: data_type_axes       ! 'hkl', 'Q', 'theta', 'theta', 'uvw', 'xyz', 'r' ...
character(len=32)                                        :: data_type_content    ! 'intensity', '3d-delta-pdf', 'amplitide', 'density' ...
character(len=32)                                        :: data_type_reciprocal ! 'reciprocal', 'patterson', 'direct'
character(len=32)                                        :: data_type_with_bragg ! 'bragg_present'; 'bragg_subtracted'
character(len=32)                                        :: data_type_symmetrized! 'none'; 'point'; 'laue'; 'space', 'sym_elem'
character(len=32)                                        :: data_type_number     ! 'real';  'complex'
character(len=32)                                        :: data_rad_radiation   ! 'xray'; 'neutron'; 'electron'
character(len=32)                                        :: data_rad_symbol      ! CU, CUA1, CU12
real(kind=PREC_DP)        , dimension(3)                 :: data_rad_length      ! Numerical value
integer                   , dimension(3)                 :: data_dimension   ! Data have dimensions along (HKL) / (UVW)
integer                                                  :: data_abs_is_hkl      ! Abscissa is 1=h 2=k 3=l
integer                                                  :: data_ord_is_hkl      ! Ordinate is 1=h 2=k 3=l
integer                                                  :: data_top_is_hkl      ! top-axis is 1=h 2=k 3=l
real(kind=PREC_DP)        , dimension(3   )              :: data_corner          ! Lower left bottom corner in fractional coordinates
real(kind=PREC_DP)        , dimension(3, 3)              :: data_vector          ! Increment vectors abs: (:,1); ord: (:,2); top: (:,3)
real(kind=PREC_DP)        , dimension(:, :, :), allocatable :: data_values  ! Actual data array
real(kind=PREC_DP)        , dimension(:, :, :), allocatable :: data_imag    ! Optional imaginary part of complex data values
!
character(len=PREC_STRING), dimension(  5)      :: crystal_meta     ! Dictionary type, version, date, creation_program, author
!
integer                                         :: ier_num = 0
character(len=PREC_STRING), dimension(NMSG)     :: ier_msg = ' '
!
!character(len=8)        :: date
!
!integer :: h, k,l  ! Dummy index
integer :: i
!real(kind=PREC_DP), dimension(3) :: hkl
!
!  Set an example structure
!
 infile           = 'example_fortran_data.hdf5'   ! Arbitrary name
!
write(*,*) ' CALLING read_data', allocated(symmetry_mat)

l_dump = .false.
call unified_read_data(  infile, &
                             unit_cell_lengths, unit_cell_angles,               &
                             symmetry_H_M, symmetry_origin, symmetry_abc, symmetry_n_mat, &
                             symmetry_mat,                                                &
                             data_type_experiment, &
                             data_type_style     , &
                             data_type_axes      , &
                             data_type_content   , &
                             data_type_reciprocal, &
                             data_type_with_bragg, &
                             data_type_symmetrized, &
                             data_type_number     , &
                             data_rad_radiation , data_rad_symbol, data_rad_length, &
                             data_dimension  , &
                             data_abs_is_hkl     , &
                             data_ord_is_hkl     , &
                             data_top_is_hkl     , &
                             data_corner     , &
                             data_vector     , &
                             data_values     , &
                             NMSG, ier_num, ier_msg,                                      &
                             crystal_meta,                  &
                             data_imag             &
                             )
!
if(ier_num/=0) call error_message(ier_num, ier_msg)
!
write(*,*)
write(*,'(a)') '      Partial dump of unified data '
write(*,*)
write(*,'(a)') '   Crystal Meta information'
do i= 1, N_META
   write(*,'(a26, 1x, a)') c_meta(i), crystal_meta(i)(1:len_trim(crystal_meta(i)))
enddo
write(*,*)
!
write(*,'(a,3f12.6)')  ' Unit cell params ', unit_cell_lengths
write(*,'(a,3f12.6)')  ' Unit cell angles ', unit_cell_angles
!
write(*,'(a, a     )') ' Symmetry H_M     ', symmetry_H_M
write(*,'(a, i5    )') ' Symmetry Origin  ', symmetry_origin
write(*,'(a, a     )') ' Symmetry abc     ', symmetry_abc
write(*,'(a, i5    )') ' Symmetry n MAT   ', symmetry_n_mat
write(*,'(a, 4f11.6)') ' Symmetry 1 fst   ', symmetry_mat(1,:,1)
write(*,'(a, 4f11.6)') ' Symmetry 1 scnd  ', symmetry_mat(2,:,1)
write(*,'(a, 4f11.6)') ' Symmetry 1 lst   ', symmetry_mat(3,:,1)
if(symmetry_n_mat>10) then
write(*,'(a, i3, a5, 4f11.6)') ' Symmetry ', symmetry_n_mat/2, 'fst  ',symmetry_mat(1,:,symmetry_n_mat/2)
write(*,'(a, i3, a5, 4f11.6)') ' Symmetry ', symmetry_n_mat/2, 'scnd ',symmetry_mat(2,:,symmetry_n_mat/2)
write(*,'(a, i3, a5, 4f11.6)') ' Symmetry ', symmetry_n_mat/2, 'lst  ',symmetry_mat(3,:,symmetry_n_mat/2)
endif
write(*,'(a, 4f11.6)') ' Symmetry N fst   ', symmetry_mat(1,:,symmetry_n_mat)
write(*,'(a, 4f11.6)') ' Symmetry N scnd  ', symmetry_mat(2,:,symmetry_n_mat)
write(*,'(a, 4f11.6)') ' Symmetry N lst   ', symmetry_mat(3,:,symmetry_n_mat)
write(*,*)
write(*,'(a, a     )') ' Data experiment  ', data_type_experiment
write(*,'(a, a     )') ' Data with  bragg ', data_type_with_bragg
write(*,'(a, a     )') ' Data symmetrized ', data_type_symmetrized
write(*,'(a, a     )') ' Data number      ', data_type_number
write(*,'(a, a     )') ' Data axes        ', data_type_axes
write(*,'(a, a     )') ' Data radiation   ', data_rad_radiation
write(*,'(a, 3i5   )') ' Data dimensions  ', data_dimension
write(*,'(a, 3i5   )') ' Data abs,ord,top ', data_abs_is_hkl, data_ord_is_hkl, data_top_is_hkl
write(*,'(a, a     )') ' Data style       ', data_type_style
write(*,'(a, a     )') ' Data content     ', data_type_content
write(*,'(a, a     )') ' Data symbol      ', data_rad_symbol
write(*,'(a, 3f11.6)') ' Data length      ', data_rad_length
write(*,*)
write(*,'(a, 3f11.6)') ' Data corner      ', data_corner(:)
write(*,*)
write(*,'(a, 3f11.6)') ' Data abscissa    ', data_vector(:,1)
write(*,'(a, 3f11.6)') ' Data ordinate    ', data_vector(:,2)
write(*,'(a, 3f11.6)') ' Data top_axis    ', data_vector(:,3)
write(*,*)
write(*,*) ' DATA_VALUES', data_values(1:2,1,1)
write(*,*) ' DATA_MINMAX', minval(data_values), maxval(data_values)
if(data_type_number == 'complex' .and. allocated(data_imag)) then
write(*,*)
write(*,*) ' DATA_IMAG  ', data_imag  (1:2,1,1)
write(*,*) ' DATA_MINMAX', minval(data_imag  ), maxval(data_imag  )
endif
write(*,*)
write(*,*)

!
end program unified_top_read_data
