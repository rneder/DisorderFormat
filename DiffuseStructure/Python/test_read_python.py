from rbn_nexus_structure import  read_diffuse_structure

#file_path = "discus_asym.hdf5"
#file_path = "discus_molecule.hdf5"
#file_path = "discus_stru.hdf5"
#file_path = "mono_p21_c.hdf5"
#file_path = "disordered.hdf5"
file_path = "rbn__nexus_file.hdf5"
file_path = "example_fortran.hdf5"

return_value = read_diffuse_structure(file_path)
unit_cell_lengths   = return_value[ 0]
unit_cell_angles    = return_value[ 1]
symmetry_H_M        = return_value[ 2]
symmetry_origin     = return_value[ 3]
symmetry_abc        = return_value[ 4]
symmetry_n_mat      = return_value[ 5]
symmetry_mat        = return_value[ 6]
unit_cells          = return_value[ 7]
number_of_types     = return_value[ 8]
types_names         = return_value[ 9]
types_ordinal       = return_value[10]
types_charge        = return_value[11]
types_isotope       = return_value[12]
number_of_atoms     = return_value[13]
atom_type           = return_value[14]
atom_position       = return_value[15]
atom_unit_cell      = return_value[16]
atom_site_number    = return_value[17]
status_flags        = return_value[18]
meta_data           = return_value[19]
average             = return_value[20]
magnetic_spins      = return_value[21]
property_flags      = return_value[22]
anisotropic_tuple   = return_value[23]
molecules           = return_value[24]
occupancy           = return_value[25]
print(" ")
print(" READ DIFFUSE STRUCTURE ")
print(" unit_cell_lengths      ", unit_cell_lengths)
print(" unit_cell_angles       ", unit_cell_angles)
print(" symmetry_H_M           ", symmetry_H_M)
print(" symmetry_origin  ", symmetry_origin)
print(" symmetry_abc     ", symmetry_abc   )
print(" symmetry_n_mat   ", symmetry_n_mat )
for i in range(symmetry_n_mat):
   print(" Symmetry_mat      ", symmetry_mat[:,0,i])
   print(" Symmetry_mat      ", symmetry_mat[:,1,i])
   print(" Symmetry_mat      ", symmetry_mat[:,2,i])
   print(" ")
print(" Unit cells       ", unit_cells[0])
print(" Unit cells       ", unit_cells[1])
print(" Unit cells       ", unit_cells[2])
print(" Number of types  ", number_of_types)
for i in range(number_of_types):
   print(" Name Ord Charge  Iso ", types_names[i], types_ordinal[i], types_charge[i], types_isotope[i])
print(" ")
print(" Number of atoms  ", number_of_atoms )
for i in range(number_of_atoms):
   print(" Atom TY POS   ", atom_type[i], atom_position[0:3,i], atom_unit_cell[0:3,i], atom_site_number[i]) 
#print(" Atom data_i      ", atom_data_i.shape    )
#print(" Atom data_f      ", atom_data_f.shape    )
print(" Status Flags     ", status_flags    )
print(" Meta Data        ", meta_data       )
if average is not None:
   print(" Average structure", average[0]         )
   average_type = average[1]
   average_pos  = average[2]
   average_occ  = average[3]
   average_adp  = average[4]
   average_site = average[5]
   for i in range(average[0]):
      print(" Atom        ", average_type[i], average_pos[:,i], average_occ[i], average_adp[:,i], average_site[i])
   print(" Shapes           ", average_type.shape, average_pos.shape, average_occ.shape, average_adp.shape, average_site.shape         )
else:
   print(" Average structure is not present" )
#
if magnetic_spins is not None:
   print(" Magnetic spins ", magnetic_spins)
else:
   print(" Magnetic spin  is not present" )
#
if property_flags is not None:
   print(" property_flags ", property_flags         )
else:
   print(" property_flags is not present" )
#
if anisotropic_tuple is not None:
   #print(" ADP ", anisotropic_adp)
   anisotropic_n_uij  = anisotropic_tuple[0]
   anisotropic_n_type = anisotropic_tuple[1]
   anisotropic_n_index = anisotropic_tuple[2]
   anisotropic_is_iso  = anisotropic_tuple[3]
   anisotropic_adps    = anisotropic_tuple[4]
   anisotropic_index = anisotropic_tuple[5]
   print(" Anisotropic_adp ", anisotropic_n_uij, anisotropic_n_type, anisotropic_n_index)
   for i in range(anisotropic_n_type):
      print("is_iso==0; UIJ Ueqv ", anisotropic_is_iso[i], anisotropic_adps[:,i])
#   
   print("INDEX   ", anisotropic_index[:])
else:
   print(" Anisotropic_tuple is not present" )
#
if molecules is not None:
   molecules_number = molecules[0]
   molecules_types  = molecules[1]
   molecules_int    = molecules[2]
   molecules_real   = molecules[3]
   molecules_index  = molecules[4]
   print(" molecules ", molecules_number, molecules_types )
   for i in range(molecules_number):
      print(" Type Character Length Atoms", molecules_int[0,i], molecules_int[1,i], molecules_int[2,i], molecules_index[0:molecules_int[2,i],i])
   for i in range(molecules_types):
      print(" Uiso, Clin, Cquad", molecules_real[0,i], molecules_real[1,i], molecules_real[2,i])
   print(" Shapes ", molecules_int.shape, molecules_real.shape)
else:
   print(" molecules is not present" )
#
if occupancy is not None:
   print(" occupancy ", occupancy         )
else:
   print(" occupancy is not present" )


