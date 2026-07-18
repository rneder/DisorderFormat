program unified_top_write_data
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
use unified_write_mod
!
use precision_mod         ! Defines KIND of real numbers, strings etc
!
implicit none
!
real(kind=prec_DP), parameter :: twopi = 2.0D0 * 3.1415926535897932384626433832795028841971693993751D0
character(len=PREC_STRING) :: outfile           ! Output file name
character(len=PREC_STRING) :: program_version   ! Which main program wrote the structure
character(len=PREC_STRING) :: author_name       ! Authors
!
real(kind=PREC_DP), dimension(3)   :: unit_cell_lengths    ! (a, b, c)
real(kind=PREC_DP), dimension(3)   :: unit_cell_angles     ! (alpha, beta, gamma)
real(kind=PREC_DP), dimension(3,3) :: metric_tensor        ! Direct space metric tensor
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
character(len=32)                                        :: coordinate_unit      !'basecell_fractional', 'supercell_fractional', ...
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
character(len=8)        :: date
!
integer :: h, k,l  ! Dummy index
real(kind=PREC_DP), dimension(3) :: hkl
!
!  Set an example structure
!
outfile           = 'example_fortran_data.hdf5'   ! Arbitrary name
!
program_version   = 'unified_top_write_data.f90'  ! This program
author_name       = 'R.B.Neder'
!
call date_and_time(date=date)
crystal_meta(1) = 'Disorder unified data'
crystal_meta(2) = '0.0.0'
crystal_meta(3) = date(1:4) // '-' // date(5:6) // '-' // date(7:8)
crystal_meta(4) = program_version
crystal_meta(5) = author_name
!
unit_cell_lengths    =  10.0_PREC_DP
unit_cell_angles     =  90.0_PREC_DP
metric_tensor        =   0.0_PREC_DP
metric_tensor(1,1)   = 100.0_PREC_DP
metric_tensor(2,2)   = 100.0_PREC_DP
metric_tensor(3,3)   = 100.0_PREC_DP
!
symmetry_H_M         = 'P m m 2'
symmetry_origin      = 1
symmetry_abc         = 'abc'
symmetry_n_mat       = 4
allocate(symmetry_mat(3,4,symmetry_n_mat))
symmetry_mat         = 0.0_PREC_DP
symmetry_mat(1,1,1)  = 1.0_PREC_DP     ! 1
symmetry_mat(2,2,1)  = 1.0_PREC_DP
symmetry_mat(3,3,1)  = 1.0_PREC_DP
!
symmetry_mat(1,1,2) =-1.0_PREC_DP      ! 2 (00z)
symmetry_mat(2,2,2) =-1.0_PREC_DP
symmetry_mat(3,3,2) = 1.0_PREC_DP
!
symmetry_mat(1,1,3) = 1.0_PREC_DP      ! m (x0z)
symmetry_mat(2,2,3) =-1.0_PREC_DP
symmetry_mat(3,3,3) = 1.0_PREC_DP
!
symmetry_mat(1,1,4) =-1.0_PREC_DP      ! m (0yz)
symmetry_mat(2,2,4) = 1.0_PREC_DP
symmetry_mat(3,3,4) = 1.0_PREC_DP
!
data_type_experiment = 'calculated'
data_type_style      = 'single_diffraction'
data_type_axes       = 'hkl'
data_type_content    = 'intensity'
data_type_reciprocal = 'reciprocal'
data_type_with_bragg = 'bragg_subtracted'
data_type_symmetrized = 'none'
data_type_number      = 'complex'
data_rad_radiation    = 'xray'
data_rad_symbol       = 'MOA1'
data_rad_length       = 0.709260D0
data_dimension        = (/41, 31, 21 /)
data_abs_is_hkl       = 1
data_ord_is_hkl       = 2
data_top_is_hkl       = 3
coordinate_unit       = 'basecell_fractional'
data_corner(:  )      = (/ -2.0D0, -1.5D0, -1.0D0 /)
data_vector(:,1)      = (/  0.1D0,  0.0D0,  0.0D0 /)
data_vector(:,2)      = (/  0.0D0,  0.1D0,  0.0D0 /)
data_vector(:,3)      = (/  0.0D0,  0.0D0,  0.1D0 /)
!
allocate(data_values(data_dimension(1), data_dimension(2), data_dimension(3)))
allocate(data_imag  (data_dimension(1), data_dimension(2), data_dimension(3)))
!
do l=1, data_dimension(3)
   hkl(3) = data_corner(3) + (l-1)*data_vector(3,3)
   do k=1, data_dimension(2)
      hkl(2) = data_corner(2) + (k-1)*data_vector(2,2)
      do h=1, data_dimension(1)
         hkl(1) = data_corner(1) + (h-1)*data_vector(1,1)
         data_values(h,k,l) = cos(twopi *hkl(1)) * &
                              cos(twopi *hkl(2)) * &
                              cos(twopi *hkl(3))
         data_imag  (h,k,l) = sin(twopi *hkl(1)) * &
                              sin(twopi *hkl(2)) * &
                              cos(twopi *hkl(3))
      enddo
   enddo
enddo
!write(*,*) ' ======================================'
!write(*,*)
!write(*,'(a, a     )') ' Data experiment  ', data_type_experiment
!write(*,'(a, a     )') ' Data with  bragg ', data_type_with_bragg
!write(*,'(a, a     )') ' Data symmetrized ', data_type_symmetrized
!write(*,'(a, a     )') ' Data number      ', data_type_number
!write(*,'(a, a     )') ' Data axes        ', data_type_axes
!write(*,*) ' ======================================'

!l_dump = .true.
call unified_write_data( outfile, unit_cell_lengths, unit_cell_angles,               &
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
                             coordinate_unit , &
                             data_corner     , &
                             data_vector     , &
                             data_values     , &
                             NMSG, ier_num, ier_msg,                                      &
                                            crystal_meta,                  &
                             data_imag = data_imag &
                             )
!
if(ier_num/=0) call error_message(ier_num, ier_msg)
!
end program unified_top_write_data
