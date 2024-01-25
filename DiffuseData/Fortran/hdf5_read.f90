module lib_hdf5_read_mod
!
! Load an HDF5 file into the lib_f90 data structure
!
!use errlist_mod
!use lib_data_struc_h5
use lib_hdf5_params_mod
use hdf5_def_mod
use precision_mod
!
private
public hdf5_read
!
CHARACTER(LEN=PREC_STRING)                            :: h5_infile         ! input file
CHARACTER(LEN=PREC_STRING), DIMENSION(:), ALLOCATABLE :: h5_datasets       ! Names of the data set in file
INTEGER                                               :: ndims             ! Number of dimensions
INTEGER                                               :: one_ndims         ! Number of dimensions
INTEGER                                               :: H5_MAX_DATASETS   ! Current MAX data sets
INTEGER                                               :: h5_data_type      ! Data type 
INTEGER                                               :: h5_n_datasets     ! Current actual data sets
INTEGER                                               :: h5_layer=1        ! Current layer in data set
INTEGER(KIND=LIB_HSIZE_T), DIMENSION(3)               :: h5_dims           ! Actual dimensions
INTEGER                  , DIMENSION(3)               :: d5_dims           ! Actual dimensionsA in transposed sequence
!NTEGER(KIND=LIB_HSIZE_T), DIMENSION(3)               :: maxdims           ! Maximum dimensions
INTEGER(KIND=LIB_HSIZE_T), DIMENSION(3)               :: one_dims          ! Actual dimensions
INTEGER(KIND=LIB_HSIZE_T), DIMENSION(3)               :: one_maxdims       ! Maximum dimensions
LOGICAL                                               :: h5_direct         ! Direct space == TRUE
LOGICAL                                               :: h5_is_grid=.true. ! Data on periodic grid
LOGICAL                                               :: h5_has_dxyz=.false. ! Data on periodic grid
LOGICAL                                               :: h5_has_dval=.false. ! Data on periodic grid
REAL(KIND=PREC_DP)   , DIMENSION(3,4)                 :: h5_corners        ! steps in H, K, L
REAL(KIND=PREC_DP)   , DIMENSION(3,3)                 :: h5_vectors        ! steps in H, K, L
REAL(KIND=PREC_DP)   , DIMENSION(6)                   :: h5_unit           ! Lattice parameters
REAL(KIND=PREC_DP)   , DIMENSION(:)    , ALLOCATABLE  :: h5_x              ! Actual x-coordinates
REAL(KIND=PREC_DP)   , DIMENSION(:)    , ALLOCATABLE  :: h5_y              ! Actual y-coordinates
REAL(KIND=PREC_DP)   , DIMENSION(:)    , ALLOCATABLE  :: h5_z              ! Actual z-coordinates
REAL(KIND=PREC_DP)   , DIMENSION(:)    , ALLOCATABLE  :: h5_dx             ! Actual x-coordinates
REAL(KIND=PREC_DP)   , DIMENSION(:)    , ALLOCATABLE  :: h5_dy             ! Actual y-coordinates
REAL(KIND=PREC_DP)   , DIMENSION(:)    , ALLOCATABLE  :: h5_dz             ! Actual z-coordinates
REAL(KIND=PREC_SP)   , DIMENSION(:,:,:), ALLOCATABLE  :: h5_data           ! Actual diffraction data
REAL(KIND=PREC_SP)   , DIMENSION(:,:,:), ALLOCATABLE  :: h5_sigma          ! Actual diffraction data
REAL(KIND=PREC_DP)   , DIMENSION(:,:,:), ALLOCATABLE  :: d5_data           ! Actual diffraction data  Double precision
REAL(KIND=PREC_DP)   , DIMENSION(:,:,:), ALLOCATABLE  :: d5_sigma          ! Actual diffraction sigma Double precision
REAL(KIND=PREC_DP)   , DIMENSION(3)                   :: h5_llims          ! Lower limits
REAL(KIND=PREC_DP)   , DIMENSION(3)                   :: h5_steps          ! steps in H, K, L
REAL(KIND=PREC_DP)   , DIMENSION(3,3)                 :: h5_steps_full     ! steps in H, K, L
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE hdf5_read(infile, out_inc, out_eck, out_vi, cr_a0, cr_win, value, qvalues)
!
!use top_data_mod
use hdf5
use iso_c_binding
!
use precision_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: infile
integer           , dimension(3)  , intent(out) :: out_inc
real(kind=PREC_DP), dimension(3,4), intent(out) :: out_eck 
real(kind=PREC_DP), dimension(3,3), intent(out) :: out_vi
real(kind=PREC_DP), dimension(3)  , intent(out) :: cr_a0
real(kind=PREC_DP), dimension(3)  , intent(out) :: cr_win
integer                           , intent(out) :: value
real(kind=PREC_DP), dimension(:,:,:), allocatable, intent(out) :: qvalues 
!
!integer, dimension(3), intent(out) :: dims
integer                            :: extr_abs
integer                            :: extr_ord
integer                            :: extr_top
integer                            :: ier_num
integer                            :: ier_typ
integer, parameter      :: ER_IO   = 2 ! Input/output error
integer, parameter      :: ER_APPL = 6 ! Error in the user program
!
!
CHARACTER(LEN=14)   :: dataname    ! Dummy name for HDF5 datasets
INTEGER(KIND=HID_T) :: file_id     ! File identifier
INTEGER(KIND=HID_T) :: dset_id     ! dataset identifier
INTEGER(KIND=HID_T) :: space_id    ! space identifier
INTEGER             :: hdferr      ! Error number
INTEGER(KIND=HSIZE_T) :: idx
INTEGER             :: ret_value   ! Error number
TYPE(C_FUNPTR)      :: funptr      ! Pointer to display function
TYPE(C_PTR), DIMENSION(:), ALLOCATABLE, TARGET :: rdata   ! READ buffer
TYPE(C_PTR)         :: f_ptr        ! Pointer to data elements
TYPE(C_PTR) :: ptr
CHARACTER(LEN=1024, KIND=c_char),  POINTER ::  rstring
INTEGER(KIND=2), PARAMETER :: SHORT_TWO = 2
INTEGER(KIND=2), TARGET :: r_is_direct
!
REAL(KIND=PREC_DP)        , DIMENSION(3)    :: steps     ! dummy steps in H, K, L
INTEGER :: i,j,k                    ! Dummy loop indices
!
integer, parameter :: idims=5
character(len=PREC_STRING), dimension(idims) :: ier_msg
!
!
!
real(kind=PREC_DP), dimension(3,3) :: temp_vi    ! Temporary copy
!
! Make sure all old data are gone
!
if(allocated(h5_datasets)) deallocate(h5_datasets)
if(allocated(h5_data)) deallocate(h5_data)
if(allocated(d5_data)) deallocate(d5_data)
if(allocated(h5_sigma)) deallocate(h5_sigma)
if(allocated(d5_sigma)) deallocate(d5_sigma)
if(allocated(h5_x)) deallocate(h5_x)
if(allocated(h5_y)) deallocate(h5_y)
if(allocated(h5_z)) deallocate(h5_z)
if(allocated(h5_dx)) deallocate(h5_dx)
if(allocated(h5_dy)) deallocate(h5_dy)
if(allocated(h5_dz)) deallocate(h5_dz)
!
h5_infile = infile
dataname = ' '
h5_steps_full = 0.0D0 
!
h5_steps      = 0.0D0
!
!
H5_MAX_DATASETS = 10                                        ! Initial estimate of dataset number
ALLOCATE(h5_datasets(H5_MAX_DATASETS))
h5_n_datasets = 0                                           ! Currently no datasets found
!
h5_dims    = 1
one_maxdims = 1
CALL H5open_f(hdferr)                                       ! Open access to HDF5 stream
CALL H5Eset_auto_f(1, hdferr)                              ! Turn Error messages off
!
CALL H5Fopen_f(h5_infile, H5F_ACC_RDWR_F, file_id, hdferr)     ! Open existing file
IF(hdferr/=0) THEN
   CALL hdf5_error(h5_infile, file_id, dataname, dset_id, ier_num, ier_typ, idims,ier_msg,&
                   -2, ER_IO, 'Could not open ''H5'' file')
   RETURN
ENDIF
!
idx = 0
funptr = C_FUNLOC(op_func) ! call back function
ptr    = C_NULL_PTR
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Iterate across file to obtain format, and list of datasets
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
CALL H5Literate_f(file_id, H5_INDEX_NAME_F, H5_ITER_NATIVE_F, idx, funptr, ptr, ret_value, hdferr)
IF(hdferr/=0) THEN
   dataname = ' '
   dset_id  = 0
   CALL hdf5_error(h5_infile, file_id, dataname, dset_id, ier_num, ier_typ, idims,ier_msg,&
                   -70, ER_APPL, 'Initial iteration failed')
   RETURN
ENDIF
!
IF(h5_n_datasets<1) THEN
   dataname = ' '
   dset_id  = 0
   CALL hdf5_error(h5_infile, file_id, dataname, dset_id, ier_num, ier_typ, idims,ier_msg,&
                   -69, ER_APPL, 'No datasets in H5 file')
   RETURN
ENDIF
!write(*,*) 'NUMBER of data sets found ', h5_n_datasets
yd_present = .false.              ! Assume all data set missing
DO i=1,h5_n_datasets
   loop_ydnd:do j=1, YD_ND
      if(h5_datasets(i)==yd_datasets(j)) then
         yd_present(j) = .true.
         exit loop_ydnd
      endif
   enddo loop_ydnd
!   WRITE(*,*) ' DATASETS ', h5_datasets(i)(1:LEN_TRIM(h5_datasets(i)))
ENDDO
!do j=1, YD_ND
!   write(*,*) ' DATASETS ', yd_present(j),yd_datasets(j)
!enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get the format identifier
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dataname = 'format'
CALL H5Dopen_f(file_id, dataname, dset_id, hdferr)           ! Open the 'format'
IF(hdferr/=0) THEN
   CALL hdf5_error(h5_infile, file_id, dataname, dset_id, ier_num, ier_typ, idims,ier_msg,&
                   -69, ER_APPL, 'Dataset does not exist')
   RETURN
ENDIF
!CALL H5Dget_storage_size_f(dset_id, nummer, hdferr)          ! Get rough storage size
!
ALLOCATE(rdata(1:1024))
f_ptr = C_LOC(rdata(1))
CALL H5Dread_F(dset_id, H5T_STRING, f_ptr, hdferr)
CALL C_F_POINTER(rdata(1), rstring)
j = 0
DO
   IF(j>8) EXIT
   IF(rstring(j+1:j+1)==C_NULL_CHAR) EXIT
   j = j+1
ENDDO
!write(*,*) ' RSTRING ', rstring(1:j), ' J ', j, LEN_TRIM(rstring)
CALL H5Dclose_f(dset_id , hdferr)
!
DEALLOCATE(rdata)
!
IF(rstring(1:j)/='Yell 1.0') THEN
   CALL hdf5_error(h5_infile, file_id, rstring, dset_id, ier_num, ier_typ, idims,ier_msg,&
                   -69, ER_APPL, 'format dataset has wrong value')
   RETURN
ENDIF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get the data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
one_maxdims = 1
dataname='data'
CALL H5Dopen_f(file_id, dataname, dset_id, hdferr)           ! Open the dataset
CALL H5Dget_space_f(dset_id, space_id, hdferr)
CALL H5Sget_simple_extent_ndims_f(space_id, ndims, hdferr)   ! Get the number of dimensions in data set
CALL H5Sget_simple_extent_dims_f(space_id, h5_dims, one_maxdims, hdferr)   ! Get the dimensions in data set
IF(ALLOCATED(h5_data)) DEALLOCATE(h5_data)
ALLOCATE(h5_data(h5_dims(1), h5_dims(2), h5_dims(3)))
CALL H5Dread_f(dset_id, H5T_NATIVE_REAL, h5_data, h5_dims, hdferr)
CALL H5Dclose_f(dset_id, hdferr)                             ! Close the dataset file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get the limits
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
dataname = 'lower_limits'
one_dims    = 1
one_maxdims = 1
CALL H5Dopen_f(file_id, dataname, dset_id, hdferr)           ! Open the dataset
CALL H5Dget_space_f(dset_id, space_id, hdferr)
CALL H5Sget_simple_extent_ndims_f(space_id, one_ndims, hdferr)   ! Get the number of dimensions in data set
CALL H5Sget_simple_extent_dims_f(space_id, one_dims, one_maxdims, hdferr)   ! Get the dimensions in data set
CALL H5Dread_f(dset_id, H5T_NATIVE_DOUBLE, h5_llims, one_dims, hdferr)
CALL H5Dclose_f(dset_id, hdferr)                             ! Close the dataset file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get the is_direct dataset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
dataname = 'is_direct'
CALL H5Dopen_f(file_id, dataname, dset_id, hdferr)           ! Open the dataset
CALL H5Dget_space_f(dset_id, space_id, hdferr)
f_ptr = C_LOC(r_is_direct)
CALL H5Dread_F(dset_id, H5T_STD_I8BE , f_ptr, hdferr)
CALL H5Dclose_f(dset_id, hdferr)
!h5_direct = 1 == (r_is_direct-8192)
h5_direct = (abs(mod(r_is_direct,SHORT_TWO))==1)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get the unit cell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
dataname = 'unit_cell'
one_dims    = 1
one_maxdims = 1
CALL H5Dopen_f(file_id, dataname, dset_id, hdferr)           ! Open the dataset
CALL H5Dget_space_f(dset_id, space_id, hdferr)
CALL H5Sget_simple_extent_ndims_f(space_id, one_ndims, hdferr)   ! Get the number of dimensions in data set
CALL H5Sget_simple_extent_dims_f(space_id, one_dims, one_maxdims, hdferr)   ! Get the dimensions in data set
CALL H5Dread_f(dset_id, H5T_NATIVE_DOUBLE, h5_unit , one_dims, hdferr)
CALL H5Dclose_f(dset_id, hdferr)                             ! Close the dataset file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get the stepss
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
dataname = 'step_sizes'
one_dims    = 1
one_maxdims = 1
CALL H5Dopen_f(file_id, dataname, dset_id, hdferr)           ! Open the dataset
CALL H5Dget_space_f(dset_id, space_id, hdferr)
CALL H5Sget_simple_extent_ndims_f(space_id, one_ndims, hdferr)   ! Get the number of dimensions in data set
CALL H5Sget_simple_extent_dims_f(space_id, one_dims, one_maxdims, hdferr)   ! Get the dimensions in data set
CALL H5Dread_f(dset_id, H5T_NATIVE_DOUBLE, h5_steps, one_dims, hdferr)
CALL H5Dclose_f(dset_id, hdferr)                             ! Close the dataset file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get the steps DISCUS style
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
if(yd_present(YD_step_sizes_abs)) then
   dataname = 'step_sizes_abs'
   one_dims    = 1
   one_maxdims = 1
   CALL H5Dopen_f(file_id, dataname, dset_id, hdferr)           ! Open the dataset
   if(hdferr==0) then
      CALL H5Dget_space_f(dset_id, space_id, hdferr)
      CALL H5Sget_simple_extent_ndims_f(space_id, one_ndims, hdferr)   ! Get the number of dimensions in data set
      CALL H5Sget_simple_extent_dims_f(space_id, one_dims, one_maxdims, hdferr)   ! Get the dimensions in data set
      CALL H5Dread_f(dset_id, H5T_NATIVE_DOUBLE, steps, one_dims, hdferr)
      h5_steps_full(:,1) = steps
      CALL H5Dclose_f(dset_id, hdferr)                             ! Close the dataset file
   endif
else
   h5_steps_full(1,1) = h5_steps(1)
endif
!
if(yd_present(YD_STEP_SIZES_ORD)) then
   dataname = 'step_sizes_ord'
   one_dims    = 1
   one_maxdims = 1
   CALL H5Dopen_f(file_id, dataname, dset_id, hdferr)           ! Open the dataset
   if(hdferr==0) then
      CALL H5Dget_space_f(dset_id, space_id, hdferr)
      CALL H5Sget_simple_extent_ndims_f(space_id, one_ndims, hdferr)   ! Get the number of dimensions in data set
      CALL H5Sget_simple_extent_dims_f(space_id, one_dims, one_maxdims, hdferr)   ! Get the dimensions in data set
      CALL H5Dread_f(dset_id, H5T_NATIVE_DOUBLE, steps, one_dims, hdferr)
      h5_steps_full(:,2) = steps
      CALL H5Dclose_f(dset_id, hdferr)                             ! Close the dataset file
   endif
else
   h5_steps_full(2,2) = h5_steps(2)
endif
!
if(yd_present(YD_STEP_SIZES_TOP)) then
   dataname = 'step_sizes_top'
   one_dims    = 1
   one_maxdims = 1
   CALL H5Dopen_f(file_id, dataname, dset_id, hdferr)           ! Open the dataset
   if(hdferr==0) then
      CALL H5Dget_space_f(dset_id, space_id, hdferr)
      CALL H5Sget_simple_extent_ndims_f(space_id, one_ndims, hdferr)   ! Get the number of dimensions in data set
      CALL H5Sget_simple_extent_dims_f(space_id, one_dims, one_maxdims, hdferr)   ! Get the dimensions in data set
      CALL H5Dread_f(dset_id, H5T_NATIVE_DOUBLE, steps, one_dims, hdferr)
      h5_steps_full(:,3) = steps
      CALL H5Dclose_f(dset_id, hdferr)                             ! Close the dataset file
   endif
else
   h5_steps_full(3,3) = h5_steps(3)
endif
!
CALL h5fclose_f(file_id, hdferr)                             ! Close the input file
!
CALL H5close_f(hdferr)                                    ! Close HDF interface
!
extr_abs = 1
extr_ord = 2
extr_top = 3
h5_steps(1) = h5_steps_full(extr_abs,1)
h5_steps(2) = h5_steps_full(extr_ord,2)
h5_steps(3) = h5_steps_full(extr_top,3)

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copy into H5 storage   Do transpose from C-style to Fortran style
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
d5_dims(1) = int(h5_dims(3))
d5_dims(2) = int(h5_dims(2))
d5_dims(3) = int(h5_dims(1))
h5_is_grid  = .true.
h5_has_dxyz = .false.
h5_has_dval = .false.
!
allocate(d5_data(d5_dims(1), d5_dims(2), d5_dims(3)))
do i=1, d5_dims(1)
   do j=1, d5_dims(2)
      do k=1, d5_dims(3)
         d5_data(i,j,k) = real(h5_data(k,j,i),kind=PREC_DP)
!        d5_data(i,j,k) = real(h5_data(i,j,k),kind=PREC_DP)
      enddo
   enddo
enddo
if(allocated(h5_sigma)) then
   allocate(d5_sigma(d5_dims(1), d5_dims(2), d5_dims(3)))
   do i=1, d5_dims(1)
      do j=1, d5_dims(2)
         do k=1, d5_dims(3)
            d5_sigma(i,j,k) = h5_sigma(k,j,i)
         enddo
      enddo
   enddo
endif
h5_vectors      = h5_steps_full
h5_corners(:,1) = h5_llims                                          ! Lower left
h5_corners(:,2) = h5_corners(:,1) + (d5_dims(1)-1)* h5_vectors(:,1)   ! Lower right
h5_corners(:,3) = h5_corners(:,1) + (d5_dims(2)-1)* h5_vectors(:,2)   ! Upper left
h5_corners(:,4) = h5_corners(:,1) + (d5_dims(3)-1)* h5_vectors(:,3)   ! Top left
!
!
! DUMMY SECTION copy into output
!
out_inc =  d5_dims
out_eck = h5_corners
out_vi  = h5_vectors
cr_a0   = h5_unit(1:3)
cr_win  = h5_unit(4:6)
if(       h5_direct) then
   value = 1
else
   value = 0
endif
if(allocated(qvalues)) deallocate(qvalues)
allocate(qvalues(out_inc(1), out_inc(2), out_inc(3)))
qvalues = 0.0_PREC_DP
qvalues = d5_data
! END DUMMY SECTION

!
if(allocated(h5_datasets)) deallocate(h5_datasets)
if(allocated(h5_data)) deallocate(h5_data)
if(allocated(d5_data)) deallocate(d5_data)
if(allocated(h5_sigma)) deallocate(h5_sigma)
if(allocated(d5_sigma)) deallocate(d5_sigma)
if(allocated(h5_x)) deallocate(h5_x)
if(allocated(h5_y)) deallocate(h5_y)
if(allocated(h5_z)) deallocate(h5_z)
if(allocated(h5_dx)) deallocate(h5_dx)
if(allocated(h5_dy)) deallocate(h5_dy)
if(allocated(h5_dz)) deallocate(h5_dz)
!
END SUBROUTINE hdf5_read
!
!*******************************************************************************
!
INTEGER FUNCTION op_func(loc_id, name, info, operator_data) bind(C)
     
!   USE allocate_generic
    USE HDF5
    USE ISO_C_BINDING
!use trig_degree_mod
    IMPLICIT NONE
     
    INTEGER, PARAMETER ::MAXSTR = 1024
    INTEGER(HID_T), VALUE :: loc_id
    CHARACTER(LEN=1), DIMENSION(1:MAXSTR) :: name ! must have LEN=1 for bind(C) strings
    TYPE(C_PTR) :: info
    TYPE(C_PTR) :: operator_data
     
    INTEGER   :: status, i, length
 
    TYPE(H5O_info_t), TARGET :: infobuf
!   TYPE(C_PTR) :: ptr
    CHARACTER(LEN=MAXSTR) :: name_string
!
    INTEGER :: all_status ! allocation status
    INTEGER :: ndata      ! upper limit for allocation
 
    !
    ! Get type of the object and display its name and type.
    ! The name of the object is passed to this FUNCTION by
    ! the Library.
    !
 
    DO i = 1, MAXSTR
       name_string(i:i) = name(i)(1:1)
    ENDDO
 
    CALL H5Oget_info_by_name_f(loc_id, name_string, infobuf, status)
 
    ! Include the string up to the C NULL CHARACTER
    length = 0
    DO
       IF(length>=MAXSTR) EXIT
       IF(name_string(length+1:length+1).EQ.C_NULL_CHAR.OR.length.GE.MAXSTR) EXIT
       length = length + 1
    ENDDO
 
    IF(infobuf%type.EQ.H5O_TYPE_GROUP_F)THEN
!      WRITE(*,*) "Group: ", name_string(1:length)
       CONTINUE
    ELSE IF(infobuf%type.EQ.H5O_TYPE_DATASET_F)THEN
!      WRITE(*,*) "Dataset: ", name_string(1:length)
       IF(h5_n_datasets==H5_MAX_DATASETS) THEN ! Need to increase size
          ndata = H5_MAX_DATASETS + 10
          CALL alloc_arr(h5_datasets, 1, ndata, all_status, ' ' )
          H5_MAX_DATASETS = H5_MAX_DATASETS + 10
       ENDIF 

       h5_n_datasets = h5_n_datasets + 1
       h5_datasets(h5_n_datasets) = name_string(1:length)
    ELSE IF(infobuf%type.EQ.H5O_TYPE_NAMED_DATATYPE_F)THEN
!      WRITE(*,*) "Datatype: ", name_string(1:length)
       CONTINUE
    ELSE
!      WRITE(*,*) "Unknown: ", name_string(1:length)
       CONTINUE
    ENDIF
 
    op_func = 0 ! return successful
 
  END FUNCTION op_func
!
!*******************************************************************************
!
subroutine alloc_arr (array, lb, ub, all_status, def_value )
!-
! Simplified subset of DISCUS allocation
!+
character(len=*), dimension(:), allocatable, intent(inout) :: array
integer         , intent(in) :: lb  ! Lower bound
integer         , intent(in) :: ub  ! Lower bound
integer         , intent(inout) :: all_status
character(len=*), intent(in) :: def_value
!
character(len=len(array)), dimension(:), allocatable :: temp
!
allocate(temp(lb:ub))
temp = def_value
!
if(allocated(array)) then
   temp(lb:ubound(array,1)) = array(lb:)
   deallocate(array)
endif
allocate(array(lb:ub))
array = temp
deallocate(temp)
all_status = 0
!
end subroutine alloc_arr
!
!*******************************************************************************
!
SUBROUTINE hdf5_error(infile, file_id, dataname, dset_id, ier_num, ier_typ, &
                      idims,ier_msg,er_nr, er_type, message)
!
IMPLICIT NONE
!
CHARACTER(LEN=*)       , INTENT(IN) :: infile     ! File name
INTEGER(KIND=LIB_HID_T), INTENT(IN) :: file_id    ! File identifier
CHARACTER(LEN=*)       , INTENT(IN) :: dataname   ! dataset name
INTEGER(KIND=LIB_HID_T), INTENT(IN) :: dset_id    ! dataset identifier
INTEGER            , INTENT(OUT) :: ier_num     ! Error number
INTEGER            , INTENT(OUT) :: ier_typ     ! Error type
INTEGER            , INTENT(IN)  :: idims       ! Dimension of ier_msg
CHARACTER(LEN=*), DIMENSION(idims)   , INTENT(OUT) :: ier_msg    ! Error message
INTEGER            , INTENT(IN) :: er_nr      ! Error number
INTEGER            , INTENT(IN) :: er_type    ! Error type
CHARACTER(LEN=*)   , INTENT(IN) :: message    ! Error message
!
INTEGER             :: hdferr      ! Error number
!
ier_num = er_nr             ! Copy error numbers into errlist_mod
ier_typ = er_type
ier_msg(1) = message(1:LEN(ier_msg))
IF(dataname /= ' ') ier_msg(2) = 'Dataset: '//dataname
IF(infile   /= ' ') ier_msg(3) = 'H5 File: '//infile
!
IF(ALLOCATED(h5_datasets)) DEALLOCATE(h5_datasets)
IF(ALLOCATED(h5_data))     DEALLOCATE(h5_data)
!
! CLOSE up, ignore error messages
IF(dset_id /= 0) CALL H5Dclose_f(dset_id , hdferr)                         ! Close dataset
IF(file_id /= 0) CALL h5fclose_f(file_id, hdferr)                          ! Close the input file
CALL H5close_f(hdferr)                                    ! Close HDF interface
!
!
END SUBROUTINE hdf5_error
!
!*******************************************************************************
!
end module lib_hdf5_read_mod
