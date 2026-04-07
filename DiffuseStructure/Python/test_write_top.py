import sys

from pathlib import Path

import numpy as np
import os
from datetime import datetime
#from nexus_structure import NeXusStructureFormat
#from py_write_structure import write_diffuse_structure
from rbn_nexus_structure import write_diffuse_structure

#Define input variables 
file_path = "rbn__nexus_file.hdf5"
file_path = "example_fortran.hdf5"

#Define header variables
#dict_name = "Disorder structure"
#dict_version = "0.0.0"
#creation_date = datetime.today().strftime('%Y-%m-%d')
creation_method = "test_write_straight.py"
author_name  = "R.B.Neder"
#
unit_cell_lengths = np.array([10.0, 10.0, 10.0], dtype=float)
unit_cell_angles  = np.array([90.0, 90.0, 90.0], dtype=float)
#
symmetry_H_M = "P m m 2"
symmetry_origin=1
symmetry_abc = "abc"
#
symmetry_n_mat = 4
symmetry_mat = np.zeros((3,4,symmetry_n_mat), np.float64)
symmetry_mat[0,0,0] = 1.0
symmetry_mat[1,1,0] = 1.0
symmetry_mat[2,2,0] = 1.0
#
symmetry_mat[0,0,1] =-1.0
symmetry_mat[1,1,1] =-1.0
symmetry_mat[2,2,1] = 1.0
symmetry_mat[2,3,1] = 0.0
#
symmetry_mat[0,0,2] = 1.0
symmetry_mat[1,1,2] =-1.0
symmetry_mat[2,2,2] = 1.0
#
symmetry_mat[0,0,3] =-1.0
symmetry_mat[1,1,3] = 1.0
symmetry_mat[2,2,3] = 1.0
symmetry_mat[2,3,3] = 0.0
#
unit_cells = np.array([ 3,1,1], dtype=np.int32)
#
types_number  = 4
types_names   = ['O', 'H' , 'N' , 'H' ]
print("TYPES NAMES ", types_names)
types_ordinal = np.array([8, 1, 7, 1], dtype=np.int32)
types_charge  = np.array([0, 0, 0, 0], dtype=np.int32)
types_isotope = np.array([16, 1, 14, 1], dtype=np.int32)
types_occupancy = np.ones((types_number, ), dtype=np.float64)
#
number_of_atoms = 30
atom_id   = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11,12,13,14,15,16,17,18,19, 10,21,22,23,24,25,26,27,28,29, 30])
atom_type = np.array([1, 2, 2, 2, 2, 3, 4, 4, 4,  4, 1, 2, 2, 2, 2, 3, 4, 4, 4,  4, 1, 2, 2, 2, 2, 3, 4, 4, 4,  4], dtype=np.int32)
atom_pos  = np.array([[-1.00, -0.87, -1.13, -0.87, -1.13, -0.50, -0.37, -0.63, -0.37, -0.63, \
                        0.00,  0.13, -0.13,  0.13, -0.13,  0.50,  0.63,  0.37,  0.63,  0.37, \
                        1.00,  1.13,  0.87,  1.13,  0.87,  1.50,  0.87,  0.87,  0.63,  1.37],\
                      [ 0.00,  0.17, -0.17, -0.17,  0.17,  0.50,  0.67,  0.33,  0.33,  0.67, \
                        0.00,  0.17, -0.17, -0.17,  0.17,  0.50,  0.67,  0.33,  0.33,  0.67, \
                        0.00,  0.17, -0.17, -0.17,  0.17,  0.50,  0.67,  0.33,  0.33,  0.67],\
                      [ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, \
                        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, \
                        0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00] \
                     ], dtype=np.float64)
atom_unit_cell = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],\
                           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],\
                           [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] \
                          ], dtype=np.int32)
atom_site = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], dtype=np.int32)
        

status_flags = {'is_super_structure': 'is_true', 'is_asymmetric_unit': 'is_false', \
                'is_periodic_x': 'is_true', \
                'is_periodic_y': 'is_true', \
                'is_periodic_z': 'is_true', \
                'is_homogeneous': 'is_true' \
}
#
average_number = 10
average_type = np.array([1, 2, 2, 2, 2, 3, 4, 4, 4, 4])
average_pos  = np.array([[ 0.00,  0.13,  0.87,  0.13,  0.87,  0.50,  0.63,  0.37,  0.63,  0.37], \
                         [ 0.00,  0.17,  0.83,  0.83,  0.17,  0.50,  0.67,  0.33,  0.33,  0.67], \
                         [ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00]])
average_occ  = np.array([1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00])
average_adp  = np.zeros((7, len(average_type)))
average_adp[0,0]    = 0.00127 
average_adp[0,1:5]  = 0.00253 
average_adp[0,5]    = 0.00380 
average_adp[0,6:10] = 0.00507 
average_adp[6,0]    = 0.00127 
average_adp[6,1:5]  = 0.00253 
average_adp[6,5]    = 0.00380 
average_adp[6,6:10] = 0.00507 
average_site = np.array([1,2, 3, 4, 5, 6, 7, 8, 9, 10])
average_structure = (average_number,average_type, average_pos, average_occ, average_adp, average_site)
    
# Define molecular identity
#
# Set molecules
molecules_number = 6
molecules_types   = 2
molecules_length  = 5
molecules_int  = np.zeros([3,molecules_number], dtype=np.int32)
molecules_real = np.zeros([3,molecules_types], dtype=np.float64)
molecules_index = np.zeros([molecules_length, molecules_number], dtype=np.int32)
molecules_int[0,0] = 1
molecules_int[1,0] = 0
molecules_int[2,0] = 5
molecules_int[0,1] = 2
molecules_int[1,1] = 0
molecules_int[2,1] = 5
molecules_int[0,2] = 1
molecules_int[1,2] = 0
molecules_int[2,2] = 5
molecules_int[0,3] = 2
molecules_int[1,3] = 0
molecules_int[2,3] = 5
molecules_int[0,4] = 1
molecules_int[1,4] = 0
molecules_int[2,4] = 5
molecules_int[0,5] = 2
molecules_int[1,5] = 0
molecules_int[2,5] = 5
molecules_real[0,0] = 0.00315
molecules_real[1,0] = 0.1
molecules_real[2,0] = 0.2
molecules_real[0,1] = 0.00630
molecules_real[1,1] = 0.0
molecules_real[2,1] = 0.0
molecules_index[0,0] =  1
molecules_index[1,0] =  2
molecules_index[2,0] =  3
molecules_index[3,0] =  4
molecules_index[4,0] =  5
molecules_index[0,1] =  6
molecules_index[1,1] =  7
molecules_index[2,1] =  8
molecules_index[3,1] =  9
molecules_index[4,1] = 10
molecules_index[0,2] = 11
molecules_index[1,2] = 12
molecules_index[2,2] = 13
molecules_index[3,2] = 14
molecules_index[4,2] = 15
molecules_index[0,3] = 16
molecules_index[1,3] = 17
molecules_index[2,3] = 18
molecules_index[3,3] = 19
molecules_index[4,3] = 20
molecules_index[0,4] = 21
molecules_index[1,4] = 22
molecules_index[2,4] = 23
molecules_index[3,4] = 24
molecules_index[4,4] = 25
molecules_index[0,5] = 26
molecules_index[1,5] = 27
molecules_index[2,5] = 28
molecules_index[3,5] = 29
molecules_index[4,5] = 30
molecules = (molecules_number, molecules_types, molecules_int, molecules_real, molecules_index)
# Define magnetic spin
#
magnetic_spins  = np.array([[-1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -0.37, -1.00, \
                             -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, \
                             -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00, -1.00],\
                            [ 1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00, \
                              1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00, \
                              1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00],\
                            [ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, \
                              0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, \
                              0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00] \
                           ], dtype=np.float64)
    
# Define property flags
property_flags = np.ones((number_of_atoms), dtype=np.int32)
property_flags[:] = 3
#
# Define anisotropic ADP
number_of_anis = 4
anis_is_iso = np.zeros((number_of_anis), dtype=np.int32)
anis_adp = np.zeros((7, number_of_anis), dtype=np.float64)
anis_adp[0,0]    = 0.00127
anis_adp[0,1]    = 0.00253
anis_adp[0,2]    = 0.00380
anis_adp[0,3]    = 0.00507
anis_adp[6,0]    = 0.00127
anis_adp[6,1]    = 0.00253
anis_adp[6,2]    = 0.00380
anis_adp[6,3]    = 0.00507
anis_index = np.zeros((number_of_atoms), dtype=np.int32)
anis_index[0]    = 1
anis_index[1:5]  = 2
anis_index[5]    = 3
anis_index[6:10] = 4
anis_index[10]   = 1
anis_index[11:15]= 2
anis_index[15]   = 3
anis_index[16:20]= 4
anis_index[20]   = 1
anis_index[21:25]= 2
anis_index[25]   = 3
anis_index[26:30]= 4
anisotropic_tuple=(number_of_anis, anis_is_iso, anis_adp, anis_index)
#
print("#####################################################################")
print("FilePath          ", file_path)
print("Creation_method   ", creation_method)
print("Author_name       ", author_name)
print("Unit_cell_lengths ", unit_cell_lengths)
print("Unit_cell_angles  ", unit_cell_angles)
print("Symmetry_H_M      ", symmetry_H_M)
print("Symmetry_origin   ", symmetry_origin)
print("Symmetry_abc      ", symmetry_abc)
for i in range(symmetry_n_mat):
   print("Symmetry_mat      ", symmetry_mat[:,0,i])
   print("Symmetry_mat      ", symmetry_mat[:,1,i])
   print("Symmetry_mat      ", symmetry_mat[:,2,i])
   print(" ")
#
print("Unit cells        ", unit_cells[0])
print("Unit cells        ", unit_cells[1])
print("Unit cells        ", unit_cells[2])
#
print("Number of types   ", types_number)
print("Types_names       ", types_names )
print("Types_ordinal     ", types_ordinal )
print("Types_charge      ", types_charge)
print("Types_isotope     ", types_isotope)
print(" ")
print("Number of atoms   ", number_of_atoms)
for i in range(number_of_atoms):
   print(" ID TY POS CL SI: ", atom_id[i], atom_type[i], atom_pos[:,i], atom_unit_cell[:,i], atom_site[i] )
print(" POS CLL  shape ", atom_pos.shape, atom_unit_cell.shape)
print(" ")
print("Status_flags    ", status_flags)
print("Average structure ")
for i in range(len(average_type)):
   print(" ID TY POS CL SI: ", average_type[i], average_pos[:,i], average_occ[i], average_adp[:,i], average_site[i] )
print(" POS ADP  shape ", average_pos.shape, average_adp.shape)
print(" ")
print("Property_Flags  ", property_flags)
for i in range(anis_adp.shape[1]):
   print("ANIS_ADP  :     ",anis_adp[:,i])
print(" ANIS_ADP.shape ", anis_adp.shape)
print(" ")
print("Molecules       ", molecules_number, molecules_types, molecules_length)
for i in range(molecules_number):
   print("TY CHAR LEN CON :   ", molecules_int[:,i], molecules_index[:,i])
print("Molecule types  ")
for i in range(molecules_types):
   print("Real values     :   ", molecules_real[:,i])
print(" I R INDX shape ", molecules_int.shape, molecules_real.shape, molecules_index.shape)
print(" ")
print("OCCUPANCY :     ", types_occupancy)
print("#####################################################################")
    
print(" ")
print(" SAVE FILE ")
print(" ")

ierror = write_diffuse_structure(file_path,
            creation_method,
            author_name,
            unit_cell_lengths,
            unit_cell_angles,
            symmetry_H_M,
            symmetry_origin,
            symmetry_abc,
            symmetry_n_mat,
            symmetry_mat,
            unit_cells,
            types_number, \
            types_names, \
            types_ordinal, \
            types_charge, \
            types_isotope, \
            number_of_atoms,
            atom_type, \
            atom_pos, \
            atom_unit_cell, \
            atom_site, \
            status_flags=status_flags,
            average_structure=average_structure,
            magnetic_spins=magnetic_spins,
            property_flags=property_flags,
            anisotropic_tuple=anisotropic_tuple,
            molecules=molecules,
            occupancy=types_occupancy
)
print(" ")
print(" SAVED FILE ", ierror)
print(" ")

