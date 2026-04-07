program unified_top_write
!
!  Dummy top program that defines a crystal structure and uses the h5fortran library
!  to write a unified structure file
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
integer           , dimension(3)                :: unit_cells   ! Number of unit cells along a, b, c  !DISCUSSION
!
integer                                         :: number_of_types     ! Number of atom types
character(len=4)  , dimension(:),   allocatable :: types_names         ! Atom type names as "O"; "O2-"
integer           , dimension(:),   allocatable :: types_ordinal       ! Ordinal numbers in Periodic table
integer           , dimension(:),   allocatable :: types_charge        ! Atomic charge as +2 or -2 ...
integer           , dimension(:),   allocatable :: types_isotope       ! Zero or specific isotpoe number
real(kind=PREC_DP), dimension(:),   allocatable :: types_occupancy     ! Occupancy for the atom type, not for an individual atom
!
integer                                         :: number_of_atoms  ! Number of atoms
integer           , dimension(:)  , allocatable :: atom_id          ! Sequential number   (might be omitted ???)
real(kind=PREC_DP), dimension(:,:), allocatable :: atom_pos         ! Fractional coordinates 1.0 is a subcell unit cell length
integer           , dimension(:)  , allocatable :: atom_type        ! This atom corresponds to types_names(i)
integer           , dimension(:,:), allocatable :: atom_unit_cell   ! Atom n is in unit cell (1:3, n)
integer           , dimension(:)  , allocatable :: atom_site        ! Atom n is on this site in its unit cell
integer           , dimension(:)  , allocatable :: atom_property    ! The atom has these properties (bits in  atom_property)
!
real(kind=PREC_DP), dimension(:,:), allocatable :: magnetic_spins   ! The atom has these magnetic spins
!
character(len=PREC_STRING), dimension(  5)      :: crystal_meta     ! Dictionary type, version, date, creation_program, author
!
integer                                         :: anis_n
integer                                         :: ier_num = 0
character(len=PREC_STRING), dimension(NMSG)     :: ier_msg = ' '
!
character(len=8)        :: date
!
logical           , dimension(8)                :: optional_intended  ! TRUE if the respective optional component is transfered
!
type(anis_adp_type)     :: anisotropic_adp      ! Details in lib_nx_transfer_mod.f90
type(molecule_data)     :: molecules            ! Details in lib_nx_transfer_mod.f90
type(average_structure) :: average_struc        ! Details in lib_nx_transfer_mod.f90
logical, dimension(2,6) :: crystal_flags        ! is_super_structure, is_asymmetric_unit, is periodic_x,y,z, is_homogeneous
!
integer :: i   ! Dummy index
!
!  Set an example structure
!
outfile           = 'example_fortran.hdf5'   ! Arbitrary name
!
program_version   = 'unified_top_write.f90'  ! This program
author_name       = 'R.B.Neder'
!
call date_and_time(date=date)
crystal_meta(1) = 'Disorder structure'
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
unit_cells           =   0
unit_cells(1)        =   3
unit_cells(2)        =   1
unit_cells(3)        =   1
!
number_of_types         = 4
allocate(types_names    (number_of_types))
allocate(types_ordinal  (number_of_types))
allocate(types_charge   (number_of_types))
allocate(types_isotope  (number_of_types))
allocate(types_occupancy(number_of_types))
!
types_names     = (/'O', 'H', 'N', 'H'/)
types_ordinal   = (/ 8 ,  1 ,  7 ,  1 /)
types_charge    = 0
types_isotope   = (/16 ,  1 , 14 ,  1 /)
types_occupancy = 1.0_PREC_DP
!
number_of_atoms = 30
allocate(atom_id       (  number_of_atoms))
allocate(atom_type     (  number_of_atoms))
allocate(atom_pos      (3,number_of_atoms))
allocate(atom_unit_cell(3,number_of_atoms))
allocate(atom_site     (  number_of_atoms))
allocate(atom_property (  number_of_atoms))
!
do i=1, number_of_atoms
   atom_id(i) = i
enddo
atom_type     = (/1, 2, 2, 2, 2, 3, 4, 4, 4, 4,  &
                  1, 2, 2, 2, 2, 3, 4, 4, 4, 4,  &
                  1, 2, 2, 2, 2, 3, 4, 4, 4, 4   &
                /)
atom_pos      (:, 1)    = 0.0_PREC_DP     ! ( 0, 0, 0)
atom_pos      (:, 2)    = (/ -0.87_PREC_DP,  0.17_PREC_DP,  0.00_PREC_DP  /)
atom_pos      (:, 3)    = (/  1.13_PREC_DP, -0.17_PREC_DP,  0.00_PREC_DP  /)
atom_pos      (:, 4)    = (/ -0.87_PREC_DP, -0.17_PREC_DP,  0.00_PREC_DP  /)
atom_pos      (:, 5)    = (/ -1.13_PREC_DP,  0.17_PREC_DP,  0.00_PREC_DP  /)
atom_pos      (:, 6)    = 0.5_PREC_DP     ! ( 0.5, 0.5, 0.5)
atom_pos      (:, 7)    = (/ -0.33_PREC_DP,  0.67_PREC_DP,  0.00_PREC_DP  /)
atom_pos      (:, 8)    = (/ -0.67_PREC_DP,  0.33_PREC_DP,  0.00_PREC_DP  /)
atom_pos      (:, 9)    = (/ -0.33_PREC_DP,  0.33_PREC_DP,  0.00_PREC_DP  /)
atom_pos      (:,10)    = (/ -0.67_PREC_DP,  0.63_PREC_DP,  0.00_PREC_DP  /)
do i=1, 10
   atom_pos(1, 10+i) = atom_pos(1, i) + 1.0_PREC_DP
   atom_pos(2, 10+i) = atom_pos(2, i)
   atom_pos(3, 10+i) = atom_pos(3, i)
   atom_pos(1, 20+i) = atom_pos(1, i) + 2.0_PREC_DP
   atom_pos(2, 20+i) = atom_pos(2, i)
   atom_pos(3, 20+i) = atom_pos(3, i)
enddo
atom_unit_cell(1, 1:10) = 1   ! First ten atoms are in unit cell (1,1,1)
atom_unit_cell(2, 1:10) = 1   ! First ten atoms are in unit cell (1,1,1)
atom_unit_cell(3, 1:10) = 1   ! First ten atoms are in unit cell (1,1,1)
atom_unit_cell(1,11:20) = 2   ! Secnd ten atoms are in unit cell (1,1,2)
atom_unit_cell(2,11:20) = 1   ! Secnd ten atoms are in unit cell (1,1,2)
atom_unit_cell(3,11:20) = 1   ! Secnd ten atoms are in unit cell (1,1,2)
atom_unit_cell(1,21:30) = 3   ! Third ten atoms are in unit cell (1,1,3)
atom_unit_cell(2,21:30) = 1   ! Third ten atoms are in unit cell (1,1,3)
atom_unit_cell(3,21:30) = 1   ! Third ten atoms are in unit cell (1,1,3)
atom_site     = (/   1,  2,  3,  4,  5,   6,  7,  8,  9, 10, &
                     1,  2,  3,  4,  5,   6,  7,  8,  9, 10, &
                     1,  2,  3,  4,  5,   6,  7,  8,  9, 10  &
                 /)
atom_property = 3     ! DISCUS: normal atom + in_molecule
!
crystal_flags = .true.
crystal_flags(2,2) = .false.
!
! Anisotropic Atomic displacment parameters
!
anis_n       = 4
allocate(anisotropic_adp%anisotropic_adp(7, anis_n))
allocate(anisotropic_adp%anisotropic_is_iso(anis_n))
allocate(anisotropic_adp%atom_index(number_of_atoms))
anisotropic_adp%anisotropic_n_type = anis_n
anisotropic_adp%anisotropic_n_atom = number_of_atoms
anisotropic_adp%anisotropic_is_iso = 0
anisotropic_adp%anisotropic_adp    = 0.0_PREC_DP
anisotropic_adp%anisotropic_adp(1,1) = 0.1/(8.*3.141519**2)
anisotropic_adp%anisotropic_adp(1,2) = 0.2/(8.*3.141519**2)
anisotropic_adp%anisotropic_adp(1,3) = 0.3/(8.*3.141519**2)
anisotropic_adp%anisotropic_adp(1,4) = 0.4/(8.*3.141519**2)
anisotropic_adp%anisotropic_adp(7,:) = anisotropic_adp%anisotropic_adp(1,:)
anisotropic_adp%atom_index    = (/1, 2, 2, 2, 2, 3, 4, 4, 4, 4,  &
                                  1, 2, 2, 2, 2, 3, 4, 4, 4, 4,  &
                                  1, 2, 2, 2, 2, 3, 4, 4, 4, 4   &
                                /)
!
! Molecules
!
molecules%number_moles = 6
molecules%number_types = 2
allocate(molecules%mole_int  (3, molecules%number_moles))
allocate(molecules%mole_real (3, molecules%number_types))
allocate(molecules%atom_index(5, molecules%number_moles))
molecules%mole_real       = 0.0_PREC_DP
molecules%mole_real(1,1)  = 0.1_PREC_DP/(8.*3.141519**2)
!
molecules%mole_int (1,1)  = 1 !Molecule types for each molecule
molecules%mole_int (1,2)  = 2
molecules%mole_int (1,3)  = 1
molecules%mole_int (1,4)  = 2
molecules%mole_int (1,5)  = 1
molecules%mole_int (1,6)  = 2
molecules%mole_int (3,:)  = 5 ! All molecules contain 5 atoms
molecules%atom_index(:,1) = (/  1,  2,  3,  4,  5/)  ! Atoms 1 to 5 are in molecule # 1
molecules%atom_index(:,2) = (/  6,  7,  8,  9, 10/)
molecules%atom_index(:,3) = (/ 11, 12, 13, 14, 15/)
molecules%atom_index(:,4) = (/ 16, 17, 18, 19, 20/)
molecules%atom_index(:,5) = (/ 21, 22, 23, 24, 25/)
molecules%atom_index(:,6) = (/ 26, 27, 28, 29, 30/)
!
! Average structure
!
average_struc%aver_n_atoms = 10
allocate(average_struc%atom_type  (   average_struc%aver_n_atoms))
allocate(average_struc%position   (3, average_struc%aver_n_atoms))
allocate(average_struc%occupancy  (   average_struc%aver_n_atoms))
allocate(average_struc%anis_adp   (7, average_struc%aver_n_atoms))
allocate(average_struc%site_number(   average_struc%aver_n_atoms))
!
average_struc%atom_type = (/ 1, 2, 2, 2, 2, 3, 4, 4, 4, 4 /)
average_struc%position  = atom_pos(:,1:10) ! 0.0_PREC_DP
average_struc%occupancy = 1.0_PREC_DP
average_struc%anis_adp(1,1) = 0.1/(8.*3.141519**2)
average_struc%anis_adp(1,2) = 0.1/(8.*3.141519**2)
average_struc%anis_adp(1,3) = 0.1/(8.*3.141519**2)
average_struc%anis_adp(1,4) = 0.1/(8.*3.141519**2)
average_struc%anis_adp(7,:) = average_struc%anis_adp(1,:)
average_struc%site_number   = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 /)
!
allocate(magnetic_spins(3,number_of_atoms))
magnetic_spins = 1.0_PREC_DP
magnetic_spins(1,:) = -1.0_PREC_DP
!
optional_intended =  .true.    ! All optional parts are present
!
call unified_write_structure(outfile, unit_cell_lengths, unit_cell_angles,                &
                             symmetry_H_M, symmetry_origin, symmetry_abc, symmetry_n_mat, &
                             symmetry_mat, unit_cells, number_of_types, types_names,      &
                             types_ordinal, types_charge, types_isotope, number_of_atoms, &
                             atom_type, atom_pos, atom_unit_cell, atom_site,              &
                             NMSG, ier_num, ier_msg,                                      &
                             optional_intended,                                           &
                             atom_property=atom_property,                                 &
                             crystal_flags=crystal_flags,                                 &
                             crystal_meta=crystal_meta,                                   &
                             anisotropic_adp=anisotropic_adp,                             &
                             molecules=molecules,                                         &
                             average_struc=average_struc,                                 &
                             magnetic_spins=magnetic_spins,                               &
                             types_occupancy=types_occupancy) 
!
if(ier_num/=0) call error_message(ier_num, ier_msg)
!
end program unified_top_write
