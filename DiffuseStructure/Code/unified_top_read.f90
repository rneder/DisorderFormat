program unified_top_read
!
!  Dummy top program that uses the DISCUS routines and FORPY to read 
!  a unified structure file from an existing HDF5 file
!  
!
use error_mod
use lib_nx_transfer_mod
use lib_unified_chars_mod
!
use unified_read_mod
use precision_mod          ! Defines KIND of real numbers, strings etc
!
implicit none
!
!
character(len=PREC_STRING) ::  infile           !  Input file name
!
real(kind=PREC_DP), dimension(3)   :: unit_cell_lengths    ! (a, b, c)
real(kind=PREC_DP), dimension(3)   :: unit_cell_angles     ! (alpha, beta, gamma)
real(kind=PREC_DP), dimension(3,3) :: metric_tensor        ! Direct space metric tensor
!
character(len=32)          :: symmetry_H_M                 ! Hermann-Mauguin Symbol
integer                    :: symmetry_origin              ! Origin choice 1 or 2
character(len=3)           :: symmetry_abc                 ! For orthorhombic space groups abc, cba etc
integer                    :: symmetry_n_mat               ! Number of symmetry metrices
real(kind=PREC_DP), dimension(:,:,:), allocatable :: symmetry_mat         ! Direct space symmetry matrices
!
integer           , dimension(3  )              :: unit_cells   ! Number of unit cells along a, b, c  !DISCUSSION
!
integer                                         :: number_of_types     ! Number of atom types
character(len=4)  , dimension(:),   allocatable :: types_names
integer           , dimension(:),   allocatable :: types_ordinal
integer           , dimension(:),   allocatable :: types_charge
integer           , dimension(:),   allocatable :: types_isotope
real(kind=PREC_DP), dimension(:),   allocatable :: types_occupancy
!
integer                                         :: number_of_atoms  ! Number of atoms
integer           , dimension(:)  , allocatable :: atom_id
real(kind=PREC_DP), dimension(:,:), allocatable :: atom_pos
integer           , dimension(:)  , allocatable :: atom_type
integer           , dimension(:,:), allocatable :: atom_unit_cell
integer           , dimension(:)  , allocatable :: atom_site
integer           , dimension(:)  , allocatable :: atom_property
!
real(kind=PREC_DP), dimension(:,:), allocatable :: magnetic_spins
!
logical :: l_property
logical :: l_anisotropic_adp
logical :: l_molecules
logical :: l_average_struc
logical :: l_magnetic_spins
logical :: l_types_occupancy
!
character(len=PREC_STRING), dimension(  5)      :: crystal_meta
integer                                         :: ier_num = 0
character(len=PREC_STRING), dimension(NMSG)     :: ier_msg = ' '
!
type(anis_adp_type)     :: anisotropic_adp
type(molecule_data)     :: molecules
type(average_structure) :: average_struc
logical, dimension(2,6) :: status_flags
!
integer :: i   ! Dummy index
!
!
infile           = 'example_fortran.hdf5'   ! Arbitrary name
!
!
call unified_read_structure( infile, unit_cell_lengths, unit_cell_angles,                 &
                             symmetry_H_M, symmetry_origin, symmetry_abc, symmetry_n_mat, &
                             symmetry_mat, unit_cells, number_of_types, types_names,      &
                             types_ordinal, types_charge, types_isotope, number_of_atoms, &
                             atom_type, atom_pos, atom_unit_cell, atom_site,              &
                             atom_property, status_flags, crystal_meta,                   &
                             anisotropic_adp, molecules, average_struc, magnetic_spins,   &
                             types_occupancy,                                             &
                             l_property, l_anisotropic_adp, l_molecules, l_average_struc, &
                             l_magnetic_spins, l_types_occupancy,                         &
                             NMSG, ier_num, ier_msg)
if(ier_num/=0) then
   call error_message(ier_num, ier_msg)
   stop
endif
allocate(atom_id(number_of_atoms))
!
write(*,'(a,3f12.6)') ' UNIT_CELL_LENGTHS   ', unit_cell_lengths
write(*,'(a,3f12.6)') ' UNIT_CELL_ANGLES    ', unit_cell_angles
metric_tensor(1,1) = unit_cell_lengths(1)**2
metric_tensor(2,2) = unit_cell_lengths(2)**2
metric_tensor(3,3) = unit_cell_lengths(3)**2
metric_tensor(1,2) = unit_cell_lengths(1)*unit_cell_lengths(2)*cosd(unit_cell_angles(3))
metric_tensor(2,1) = unit_cell_lengths(1)*unit_cell_lengths(2)*cosd(unit_cell_angles(3))
metric_tensor(1,3) = unit_cell_lengths(1)*unit_cell_lengths(3)*cosd(unit_cell_angles(2))
metric_tensor(3,1) = unit_cell_lengths(1)*unit_cell_lengths(3)*cosd(unit_cell_angles(2))
metric_tensor(2,3) = unit_cell_lengths(2)*unit_cell_lengths(3)*cosd(unit_cell_angles(1))
metric_tensor(3,2) = unit_cell_lengths(2)*unit_cell_lengths(3)*cosd(unit_cell_angles(1))
write(*,'(a,3f12.6)') ' Metric Tensor (1,:) ', metric_tensor(1,:)
write(*,'(a,3f12.6)') ' Metric Tensor (1,:) ', metric_tensor(2,:)
write(*,'(a,3f12.6)') ' Metric Tensor (1,:) ', metric_tensor(3,:)
write(*,'(a,a     )') ' Symmetry H_M        ', symmetry_H_M(1:len_trim(symmetry_H_M))
write(*,'(a,i8)'    ) ' Space group origin  ', symmetry_origin
write(*,'(a,a     )') ' Symmetry abc        ', symmetry_abc(1:len_trim(symmetry_abc))
write(*,'(a,i8)'    ) ' Number of Symmetry  ', symmetry_n_mat
do i=1,symmetry_n_mat
write(*,'(a,i8    )') ' Symmetry matrix No. ', i
write(*,'(a,4f12.6)') ' Symmetry matrix     ', symmetry_mat(1,:,i)
write(*,'(a,4f12.6)') ' Symmetry matrix     ', symmetry_mat(2,:,i)
write(*,'(a,4f12.6)') ' Symmetry matrix     ', symmetry_mat(3,:,i)
enddo
write(*,'(a,3i8   )') ' Unit cells          ', unit_cells(1)
write(*,'(a,3i8   )') ' Unit cells          ', unit_cells(2)
write(*,'(a,3i8   )') ' Unit cells          ', unit_cells(3)
write(*,'(a,i8)'    ) ' Number of types     ', number_of_types
write(*,'(a,20(a4,a1))') ' Atom names        ', (types_names(i), ' ',i=1, number_of_types)
write(*,'(a,20(i4,a1))') ' Atom ordinal      ', (types_ordinal(i), ' ',i=1, number_of_types)
write(*,'(a,20(i4,a1))') ' Atom charge       ', (types_charge (i), ' ',i=1, number_of_types)
write(*,'(a,20(i4,a1))') ' Atom isotope      ', (types_isotope(i), ' ',i=1, number_of_types)
if(.not.allocated(types_occupancy)) then
   write(*,'()') ' Atom occupancy    ', ' not present'
else
   write(*,'(a,20(f8.4,a1 ))') ' Atom occupancy    ', (types_occupancy(i), ' ',i=1, number_of_types)
endif
!
if(.not.l_property) then                  ! Properties were missing
   allocate(atom_property(number_of_atoms))    ! Silently allocate property
   atom_property = 1                           ! Silently assume normal atoms, no further properties
endif
write(*,'(a,i8)'    ) ' Number of atoms     ', number_of_atoms
write(*,'(a)') ' Atom:       No.      Id.   Type     x           y           z          Prop       ____Unit_cell____    Site       Spin_vector'
do i=1, number_of_atoms
   atom_id(i) = i
   if(l_magnetic_spins) then                   ! Magnetic spins were present
      write(*,'(a,3i8,3f12.6,5i8, 3f8.4)') 'Atom:   ', i, atom_id(i), atom_type(i), &
               atom_pos(:,i), atom_property(i), atom_unit_cell(:,i), atom_site(i),      &
               magnetic_spins(:,i)
   else
      write(*,'(a,3i8,3f12.6,5i8)') 'Atom:   ', i, atom_id(i), atom_type(i), &
               atom_pos(:,i), atom_property(i), atom_unit_cell(:,i), atom_site(i)
   endif
enddo
do i=1, 6
   write(*,'(2a,2l2    )') ' Crystal Flags ', c_flags(i), status_flags(:,i)
enddo
do i=1, 5
   write(*,'(a, a26, 2a )') ' Meta data     ', c_meta(i)(1:len_trim(c_meta(i))), ': ',&
      crystal_meta(i)(1:len_trim(crystal_meta(i)))
enddo
if(l_anisotropic_adp) then                     ! Anisotropic ADP were present
   write(*,'(a,i8)'    ) ' Number of ADP type  ', anisotropic_adp%anisotropic_n_type
   do i=1, anisotropic_adp%anisotropic_n_type
      write(*,'(a,i3, 6f10.6,2x, f9.6)') ' ADP Uij       ', anisotropic_adp%anisotropic_is_iso(i), anisotropic_adp%anisotropic_adp(:,i)
   enddo
   write(*,'(a,  10i5)') ' Atoms have ADP type : ', anisotropic_adp%atom_index( 1:10)
   write(*,'(23x,10i5)')                            anisotropic_adp%atom_index(11:20)
   write(*,'(23x,10i5)')                            anisotropic_adp%atom_index(21:30)
   write(*,'(a,2i5)')   ' Mole: number, types   ', molecules%number_moles, molecules%number_types
else
   write(*,'(a)') ' Anisotropic ADPs were absent'
endif
if(l_molecules      ) then                     ! Anisotropic ADP were present
   do i=1, molecules%number_moles
      write(*,'(a,4i5)'  )    ' Mole: Nr; Ty, CH, LEN ', i, molecules%mole_int(:,i)
      write(*,'(a,5x,3f9.6)') '           Ueqv, CL, CQ',    molecules%mole_real(:,molecules%mole_int(1,i))
      write(*,'(a,5x,20i5)' ) '           content     ',    molecules%atom_index(:,i)
   enddo
else
   write(*,'(a)') ' Molecules        were absent'
endif
!
end program unified_top_read
