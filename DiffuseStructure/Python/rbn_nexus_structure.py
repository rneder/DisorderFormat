# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 14:33:12 2024

@author: Maksim Eremenko
@author: Reinhard B. Neder
"""
#
#   Write/read  DiffuseDevelopers common crystal structure format into "outfile"
#
#   Arguments to write_diffuse_structure
#
#   Required part ************************************************************************
#   Data type (dims)        Name               Content in argument
#   char                    file_path          Output file name
#
#   char                    creation_method    Program that ran this script 'DISCUS 6.16.02' or similar
#                                              "yell", "rmcprofile", "rmcdiscord", "scatty"
#   char                    author_name        Authors 'R.B. Neder'                         or similar
#   ndarray float64 (3)     unit_cell_lengths  Unit cell parameter a, b, c in Angstroem
#   ndarray float64 (3)     unit_cell_angles   Unit cell parameter al, be, ga  in Angle
#   char                    symmetry_H_M       Hermann_Mauguin symbols like "F m -3 m"
#   int                     symmetry_orgin     Space group origin as in IT: 1 or 2
#   char                    symmetry_abc       Space group permutation for orthorhombic "abc" ...
#   int                     symmetry_n_mat     Number of symmetry elements
#   ndarray float64 (3,4,N) symmetry_mat       Symmetry elements direct space 3x4 notation N = symmetry_n_mat
#                                              To be understood as:
#                                              (W11  W12  W13  w1)    (x)
#                                              (W21  W22  W23  w2) *  (y)
#                                              (W31  W32  W33  w3)    (z)
#                                                                     (1)
#   ndarray int     (3)     unit_cells         Crystal is composed of this many unit cells along a, b, c
#
#   int                     number_of_types    Grand number of atom types in crystal = M
#   char array      (M)     types_names        Name for each atom type
#   int             (M)     types_ordinal      Ordinal number for each atom type H=1 ... U=92
#   int             (M)     types_charge       Ionic charge of this atom type
#   int             (M)     types_isotope      zero or isotpe number for atom type
# 
#   int                     number_of_atoms    Grand number of atoms in crystal = N
#   int             (N)     atom_type          A reference to the atomtypeID identifies the atom type
#   ndarray float64 (3,N)   atom_postion       Actual fractional atom positions in multiples of basic unit cell
#   ndarray int     (3,N)   atom_unit_cell     Atom is in this unit cell of crystal
#   ndarray int     (N)     atom_site_number   Atom is on this site in the unit cell of crystal
#
########################################################################################################
#   Optional part ** Optional arguments are in **kwargs ************************************
#
#   Optional arguments can be "status_flags"  , "average_structure", "magnetic_spin", 
#                             "property_flags", "anisotropic_tuple"  , "molecules"    ,
#                             "occupancy"
#
#   Data type (dims)        Name               Content in argument
#
#===========================================================================================
#   dict                    status_flags       Dictionary with character entries: 'is_true' or 'is_false'
#                                              'is_super_structure'  This is a supercell
#                                              'is_asymmetric_unit'  This is just an asymmetric unit
#                                              'is_periodic_x'       Periodic boundary conditions apply
#                                              'is_periodic_y'
#                                              'is_periodic_z'
#                                              'is_homogeneous'      Crystal is homogeneous enough, an average
#                                                                    structure can be constructed
#
#===========================================================================================
#   tuple                   average_structure  Average crystal structure with entries:
#           int      N      average_number     Atom types in average structure 
#   ndarray int     (N)     average_type       Atom types in average structure 
#   ndarray float64 (3, N)  average_pos        Atom position for each of the N atoms
#   ndarray float64 (   N)  average_occ        Occupancy values for each of the N atoms
#   ndarray float64 (7, N)  average_adp        [U11, U22, U33, U23, U13, U12, Ueqv] for each of the N atoms
#   ndarray int     (N)     average_site       Atom ns on this site in the unit cell in average structure 
#
#===========================================================================================
#   ndarray float64 (3,M)   magnetic_spins     Spin vector for each atom
#
#===========================================================================================
#   ndarray int     (M)     property_flags     A property flag for each atom
#
#===========================================================================================
#   tuple                   anisotropic_tuple  Info on anisotropic ADPs with entries:
#   ndarray int     (   K)  anisotropic_type   0 for isotropic types, 1 for anisotropic types
#   ndarray float64 (7, K)  anisotropic_adp    [U11, U22, U33, U23, U13, U12, Ueqv] for each of the K types
#                                              [Uiso, 0,  0,   0,   0,   0,   Uiso] for iosotropic ADP
#   ndarray int     (M)     atom_index         Index of the ADP for each atom in the structure
#
#===========================================================================================
#   tuple                   molecules          Info on molecules with entries:
#   ndarray int     (3,K)   molecules_int      Integer info on each of the K molecules:
#                                              0: Molecule type
#                                              1: Molecule character (0=atoms, 1=domains)
#                                              2: Molecule length    (Number of atoms in this molecule)
#   ndarray float64 (3,J)   molecules_real     real valued info on each of the J molecule types:
#                                              0: Molecule biso      an Isotropic Uiso for this molecule
#                                              1: Molecule clin      a linear correction term for PDF 
#                                              2: Molecule cquad     a quadratic correction term for PDF 
#   ndarray int     (O,K)   atom_index         This molecule k contains atoms 1 to o; O is the maximum length
#
#===========================================================================================
#   ndarray float64 (N)     occupany           Occupancy for each atom type
#===========================================================================================
#   ndarray float64 (M)     atom_occupany      An individual occupancy for each atom in the structure
#
############################################################################################

import h5py
import numpy as np
from datetime import datetime

class NeXusStructureFormat:
    def __init__(self, unit_cell_lengths, unit_cell_angles,
                       symmetry_H_M,
                       symmetry_origin,
                       symmetry_abc,
                       symmetry_n_mat,
                       symmetry_mat,
                       unit_cells,
                       number_of_types,
                       types_names,
                       types_ordinal,
                       types_charge,
                       types_isotope,
                       number_of_atoms,
                       atom_type,
                       atom_position,
                       atom_unit_cell,
                       atom_site_number,
                       status_flags=None,
                       metadata=None,
                       average_structure=None,
                       magnetic_spins=None,
                       property_flags=None,
                       anisotropic_tuple=None,
                       molecules=None,
                       occupancy=None,
                       atom_occupancy=None):
        
        self.set_unit_cell_lengths(unit_cell_lengths)
        self.set_unit_cell_angles(unit_cell_angles)
        self.set_symmetry(symmetry_H_M)
        self.set_symmetry_origin(symmetry_origin)
        self.set_symmetry_abc(symmetry_abc)
        self.set_symmetry_n_mat(symmetry_n_mat)
        self.set_symmetry_mat(symmetry_mat)
        self.set_unit_cells(unit_cells)
        self.set_number_of_types(number_of_types)
        self.set_types_names(types_names)
        self.set_types_ordinal(types_ordinal)
        self.set_types_charge(types_charge)
        self.set_types_isotope(types_isotope)
        self.set_number_of_atoms(number_of_atoms)
        self.set_atom_type(atom_type)
        self.set_atom_position(atom_position)
        self.set_atom_unit_cell(atom_unit_cell)
        self.set_atom_site_number(atom_site_number)
        self.set_status_flags(status_flags)
        self.set_metadata(metadata)
        self.set_average_structure(average_structure)
        self.set_magnetic_spins(magnetic_spins)
        self.set_property_flags(property_flags)
        self.set_anisotropic(anisotropic_tuple)
        self.set_molecules(molecules)
        self.set_occupancy(occupancy)
        self.set_atom_occupancy(atom_occupancy)

    def set_unit_cell_lengths(self, unit_cell_lengths):
        if not isinstance(unit_cell_lengths, np.ndarray) or unit_cell_lengths.shape != (3,) or unit_cell_lengths.dtype != float:
            raise ValueError("Unit cell lengths must be a vector dimension 3 of type float.")
        self.unit_cell_lengths = unit_cell_lengths

    def set_unit_cell_angles(self, unit_cell_angles):
        if not isinstance(unit_cell_angles, np.ndarray) or unit_cell_angles.shape != (3,) or unit_cell_angles.dtype != float:
            raise ValueError("Unit cell angles must be a vector dimension 3 of type float.")
        self.unit_cell_angles = unit_cell_angles

    def set_symmetry(self, symmetry_H_M):
        if not isinstance(symmetry_H_M, str) or len(symmetry_H_M) > 32:
            raise ValueError("Symmetry must be a string of up to 32 characters.")
        self.symmetry_H_M = symmetry_H_M

    def set_symmetry_origin(self, symmetry_origin):
        if not isinstance(symmetry_origin, int) or symmetry_origin <= 0 or symmetry_origin >=3:
            raise ValueError("Symmetry origin must be a positive integer 1 or 2.")
        self.symmetry_origin = symmetry_origin

    def set_symmetry_abc(self, symmetry_abc):
        if not isinstance(symmetry_abc, str) or len(symmetry_abc) != 3:
            raise ValueError("Symmetry abc permutation must be a string of 3 characters.")
        self.symmetry_abc = symmetry_abc

    def set_symmetry_n_mat(self, symmetry_n_mat):
        if not isinstance(symmetry_n_mat, int) or symmetry_n_mat <= 0:
            raise ValueError("Number of symmetry matrices must be a positive integer.")
        self.symmetry_n_mat = symmetry_n_mat

    def set_symmetry_mat(self, symmetry_mat):
        if not isinstance(symmetry_mat, np.ndarray) or symmetry_mat.shape != (3, 4, self.symmetry_n_mat) or symmetry_mat.dtype != float:
            raise ValueError("Symmetry matrices must be a 3x4xN matrix of type float.")
        self.symmetry_mat = symmetry_mat

    def set_unit_cells(self, unit_cells):
        if not isinstance(unit_cells, np.ndarray) or unit_cells.shape != (3,) or unit_cells.dtype != np.int32:
            raise ValueError("Unit cells must be a vector of dimension 3 of unsigned integers.")
        self.unit_cells = unit_cells

    def set_number_of_types(self, number_of_types):
        if not isinstance(number_of_types, int) or number_of_types <= 0:
            raise ValueError("Number of types must be a positive integer.")
        self.number_of_types = number_of_types

    def set_number_of_atoms(self, number_of_atoms):
        if not isinstance(number_of_atoms, int) or number_of_atoms <= 0:
            raise ValueError("Number of atoms must be a positive integer.")
        self.number_of_atoms = number_of_atoms

    def set_types_names(self, types_names):
        self.types_names = types_names
        #print("Set name    ", self.types_names  )

    def set_types_ordinal(self, types_ordinal):
        if not isinstance(types_ordinal, np.ndarray) or types_ordinal.dtype != np.int32:
            raise ValueError("Types ordinal numbers must be a vector of unsigned integers.")
        self.types_ordinal = types_ordinal
        #print("Set ordinal ", types_ordinal)

    def set_types_charge(self, types_charge):
        if not isinstance(types_charge, np.ndarray) or types_charge.dtype != np.int32:
            raise ValueError("Types charge must be a vector of unsigned integers.")
        self.types_charge = types_charge
        #print("Set charge  ", types_charge )

    def set_types_isotope(self, types_isotope):
        if not isinstance(types_isotope, np.ndarray) or types_isotope.dtype != np.int32:
            raise ValueError("Types isotope numbers must be a vector of unsigned integers.")
        self.types_isotope = types_isotope
        #print("Set isotope ", types_isotope)

    def set_atom_type(self, atom_type):
        if not isinstance(atom_type, np.ndarray) or atom_type.dtype != np.int32:
            raise ValueError("Atom type must be a vector of unsigned integers.")
        self.atom_type = atom_type

    def set_atom_position(self, atom_position):
        if not isinstance(atom_position, np.ndarray) or atom_position.dtype != np.float64:
            raise ValueError("Atom position must be an array of float64.")
        self.atom_position = atom_position

    def set_atom_unit_cell(self, atom_unit_cell):
        if not isinstance(atom_unit_cell, np.ndarray) or atom_unit_cell.dtype != np.int32:
            raise ValueError("Atom unit cell must be an array of unsigned integers.")
        self.atom_unit_cell = atom_unit_cell

    def set_atom_site_number(self, atom_site_number):
        if not isinstance(atom_site_number, np.ndarray) or atom_site_number.dtype != np.int32:
            raise ValueError("Atom unit cell must be an array of unsigned integers.")
        self.atom_site_number = atom_site_number

    def set_status_flags(self, status_flags):
        if status_flags is not None:
           self.status_flags = status_flags
        else:
           self.status_flags= None
           #print("Status flags not set ")

    def set_metadata(self, metadata):
        if metadata is not None:
           self.metadata = metadata
        else:
           self.metadata = None
           #print("Metadata not set ")
 
    def set_average_structure(self, average_structure):
        if average_structure is not None:
           self.average_structure = average_structure
        else:
           self.average_structure = None
           #print("Average structure not set ")
 
    def set_anisotropic(self, anisotropic_tuple):
        if anisotropic_tuple is not None:
            self.anisotropic_tuple = anisotropic_tuple
        else:
            self.anisotropic_tuple = None
            #print(" ANIS_ADP is None")
        #print(" ANIS_ADP is set ")

    def set_molecules(self, molecules):
        if molecules is not None:
           self.molecules = molecules
        else:
           self.molecules = None
           #print("MOLECULES not set ")

    def set_occupancy(self, occupancy):
        #print("IN set_occupancy ")
        if occupancy is not None:
           #print("occupancy ", occupancy)
           self.occupancy = occupancy
        else:
           #print("occupancy not set ")
           self.occupancy = None

    def set_atom_occupancy(self, atom_occupancy):
        #print("IN set_atom_occupancy ")
        if atom_occupancy is not None:
           #print("atom_occupancy ", atom_occupancy)
           self.atom_occupancy = atom_occupancy
        else:
           #print("atom_occupancy not set ")
           self.atom_occupancy = None

    def set_magnetic_spins(self, magnetic_spins):
        if magnetic_spins is not None:
           if not isinstance(magnetic_spins, np.ndarray) or magnetic_spins.dtype != np.float64:
               raise ValueError("Magnetic_spins must be an array of float64.")
           self.magnetic_spins = magnetic_spins
        else:
           #print("magnetic_spins not set ")
           self.magnetic_spins = None

    def set_property_flags(self, property_flags):
        if property_flags is not None:
           self.property_flags = property_flags
        else:
           self.property_flags = None
           #print("property_flags not set ")
 
    def save_to_nexus(self, file_path):
        ##print("In save_to_nexus")
        try:
            #print("In save_to_nexus typ  A")
            with h5py.File(file_path, 'w') as file:
              
                nx_entry = file.create_group('entry')
                nx_entry.attrs['NX_class'] = 'NXentry'
                
                nx_data = nx_entry.create_group('data')
                nx_data.attrs['NX_class'] = 'NXdata'
                
                # Essential data elements
                #
                nx_data.create_dataset('unit_cell_lengths', data=self.unit_cell_lengths)
                nx_data['unit_cell_lengths'].attrs['Units'] = 'angstroem'
                nx_data.create_dataset('unit_cell_angles', data=self.unit_cell_angles)
                nx_data['unit_cell_angles'].attrs['Units'] = 'angle'
                #
                nx_data.create_dataset('symmetry_space_group_name_H-M', data=np.string_(self.symmetry_H_M))
                nx_data.create_dataset('space_group_origin', data=np.array([self.symmetry_origin], dtype=np.uint))
                nx_data.create_dataset('symmetry_space_group_abc', data=np.string_(self.symmetry_abc))
                nx_data.create_dataset('space_group_symop_number', data=np.array([self.symmetry_n_mat], dtype=np.uint))
                nx_data.create_dataset('space_group_symop_operation_mat', data=self.symmetry_mat)
                #
                nx_data.create_dataset('unit_cells', data=self.unit_cells)
                nx_data.create_dataset('number_of_types', data=np.array([self.number_of_types], dtype=np.uint))
                nx_data.create_dataset('number_of_atoms', data=np.array([self.number_of_atoms], dtype=np.uint))
                #
                # Individual data sets for information on Atom types
                types_names_con = ";".join(self.types_names)
                nx_data.create_dataset('types_names', data=np.string_(types_names_con))
                nx_data.create_dataset('types_ordinal', data=self.types_ordinal)
                nx_data.create_dataset('types_charge', data=self.types_charge)
                nx_data.create_dataset('types_isotope', data=self.types_isotope)

                # Save datasets for the information on the individual atoms 
                nx_data.create_dataset('atom_type', data=self.atom_type)
                nx_data.create_dataset('atom_position', data=self.atom_position)
                nx_data.create_dataset('atom_unit_cell', data=self.atom_unit_cell)
                nx_data.create_dataset('atom_site_number', data=self.atom_site_number)

                ###########################################################################################
                # Optional and flags processing
                if self.status_flags is not None:
                    for key, value in self.status_flags.items():
                        nx_data.attrs['status_flag_' + key] = value
                #
                # Write metadata  
                if self.metadata is not None:
                    #metas_names_con = ";".join(self.metadata)
                    nx_data.create_dataset('audit_conform_dict_name', data=np.string_(self.metadata[0]))
                    nx_data.create_dataset('audit_conform_dict_version', data=np.string_(self.metadata[1]))
                    nx_data.create_dataset('audit_creation_date', data=np.string_(self.metadata[2]))
                    nx_data.create_dataset('audit_creation_method', data=np.string_(self.metadata[3]))
                    nx_data.create_dataset('audit_author_name', data=np.string_(self.metadata[4]))
                    ##nx_data.create_dataset('metadata', data=np.string_(metas_names_con))
                    #nx_metadata = nx_data.create_group('metadata')
                    #for key, value in self.metadata.items():
                    #    nx_metadata.attrs[key] = value
                
                # Optional components handling
                # Save average structure
                if self.average_structure is not None:
                    #print(" AVERAGE IS NOT None", self.average_structure)
                    nx_data.create_dataset('average_number', data=self.average_structure[0])
                    nx_data.create_dataset('average_type',   data=self.average_structure[1])
                    nx_data.create_dataset('average_pos',    data=self.average_structure[2])
                    nx_data.create_dataset('average_occ',    data=self.average_structure[3])
                    nx_data.create_dataset('average_adp',    data=self.average_structure[4])
                    nx_data.create_dataset('average_site',   data=self.average_structure[5])
                #else:
                #    #print(" AVERAGE IS     None")
                
                # Save molecules
                if self.molecules is not None:
                     #print(" MOLECULE IS NOT None", self.molecules[0])
                     #print(" MOLECULE IS NOT None", self.molecules[1])
                     #print(" MOLECULE IS NOT None", self.molecules[2], self.molecules[2].shape)
                     #print(" MOLECULE IS NOT None", self.molecules[3], self.molecules[3].shape)
                     #print(" MOLECULE IS NOT None", self.molecules[4], self.molecules[4].shape)
                     nx_data.create_dataset('molecules_number', data=np.array([self.molecules[0]    ], dtype=np.int32))
                     nx_data.create_dataset('molecules_types' , data=np.array([self.molecules[1]    ], dtype=np.int32))
                     nx_data.create_dataset('molecules_int'   , data=          self.molecules[2]                    )
                     nx_data['molecules_int'].attrs['Content'] = 'Types; Characters; Lengths'
                     nx_data.create_dataset('molecules_real'  , data=          self.molecules[3])
                     nx_data['molecules_real'].attrs['Content'] = 'Ueqv; Corr_lin; Corr_quad'
                     nx_data.create_dataset('molecules_index' , data=          self.molecules[4])

                # Save occupancy
                if self.occupancy is not None:
                    nx_data.create_dataset('occupancy', data=self.occupancy)
                
                # Save magnetic spin
                if self.magnetic_spins is not None:
                    nx_data.create_dataset('magnetic_spins', data=         self.magnetic_spins)
                
                # Save property flags
                if self.property_flags is not None:
                    nx_data.create_dataset('property_flags', data=np.array(self.property_flags, dtype=np.int32))
                
                # Save anisotropic ADP
                if self.anisotropic_tuple is not None:
                     nx_data.create_dataset('anisotropic_number',  data=np.array(self.anisotropic_tuple[0], dtype=np.int32))
                     nx_data.create_dataset('anisotropic_is_iso',  data=np.array(self.anisotropic_tuple[1], dtype=np.int32))
                     nx_data.create_dataset('anisotropic_adp',   data=np.array(self.anisotropic_tuple[2], dtype=np.float64))
                     nx_data.create_dataset('anisotropic_index', data=np.array(self.anisotropic_tuple[3], dtype=np.int32))
                     nx_data['anisotropic_is_iso'].attrs['Content'] = '0 for isotropic; 1 for anisotropic'
                     nx_data['anisotropic_adp'].attrs['Content'] = 'U11; U22; U33; U23; U13; U12; Ueqv'

        except IOError:
            raise IOError("Failed to write to the specified file path.")
    @classmethod
    def load_from_nexus(cls, file_path):
        #print("In load_from_nexus", file_path)
        try:
            #print(" READING FILE")
            with h5py.File(file_path, 'r') as file:
                nx_data = file['entry/data']
                
                # Unit cell lengths validation
                unit_cell_lengths = nx_data['unit_cell_lengths'][:]
                if unit_cell_lengths.shape != (3, ) or unit_cell_lengths.dtype != float:
                    raise ValueError("Unit cell lengths must be array of length 3  of floats.")
                
                # Unit cell angles validation
                unit_cell_angles = nx_data['unit_cell_angles'][:]
                if unit_cell_angles.shape != (3, ) or unit_cell_angles.dtype != float:
                    raise ValueError("Unit cell angles must be array of length 3  of floats.")
                
                # Symmetry validation
                symmetry_H_M = nx_data['symmetry_space_group_name_H-M'][()].decode('utf-8')
                if not isinstance(symmetry_H_M, str) or len(symmetry_H_M) > 32:
                    raise ValueError("Symmetry must be a string of up to 32 characters.")

                # Symmetry Origin number
                symmetry_origin_array = nx_data['space_group_origin'][:]  # Load as array
                if symmetry_origin_array.size == 0:
                    raise ValueError("Space group origin array is empty.")
                symmetry_origin = int(symmetry_origin_array[0])  # Convert the first element to int
                
                if symmetry_origin <= 0 or symmetry_origin >=3:
                    raise ValueError("Space group origin must 1 or 2.")
                
                # Symmetry abc validation
                symmetry_abc = nx_data['symmetry_space_group_abc'][()].decode('utf-8')
                if not isinstance(symmetry_abc, str) or len(symmetry_abc) != 3:
                    raise ValueError("Symmetry abc permutation must be a string of 3 characters.")

                # Symmetry Number of matrices
                symmetry_n_mat_array = nx_data['space_group_symop_number'][:]  # Load as array
                if symmetry_n_mat_array.size == 0:
                    raise ValueError("Number of symmetry matrices array is empty.")
                symmetry_n_mat = int(symmetry_n_mat_array[0])  # Convert the first element to int
                
                if symmetry_n_mat <= 0:
                    raise ValueError("Number of symmetry matrices must be a positive integer.")
                
                symmetry_mat = nx_data['space_group_symop_operation_mat'][:]
                if symmetry_mat.shape != (3, 4, symmetry_n_mat) or symmetry_mat.dtype != float:
                    raise ValueError("Symmetry matrices must be a 3x4xN matrix of floats.")

                # Unit cells validation
                unit_cells = nx_data['unit_cells'][:]
                if unit_cells.shape != (3,) or unit_cells.dtype != np.int32:
                    raise ValueError("Unit cells must be a vector of dimension 3 of unsigned integers.")
                
                # Number of types validation
                # Load number of types, ensuring we extract the first element and convert appropriately
                number_of_types_array = nx_data['number_of_types'][:]  # Load as array
                if number_of_types_array.size == 0:
                    raise ValueError("Number of types array is empty.")
                number_of_types = int(number_of_types_array[0])  # Convert the first element to int
                
                if number_of_types <= 0:
                    raise ValueError("Number of types must be a positive integer.")

                # Number of atoms validation
                # Load number of atoms, ensuring we extract the first element and convert appropriately
                number_of_atoms_array = nx_data['number_of_atoms'][:]  # Load as array
                if number_of_atoms_array.size == 0:
                    raise ValueError("Number of atoms array is empty.")
                number_of_atoms = int(number_of_atoms_array[0])  # Convert the first element to int
                
                if number_of_atoms <= 0:
                    raise ValueError("Number of atoms must be a positive integer.")

                # Load type information
                types_names   = nx_data['types_names'][()].decode('utf-8').split(";")
                types_ordinal = nx_data['types_ordinal'][:]
                types_charge  = nx_data['types_charge'][:]
                types_isotope = nx_data['types_isotope'][:]

                ### Load the atom data
                atom_type = nx_data['atom_type'][:]
                atom_position = nx_data['atom_position'][:]
                atom_unit_cell = nx_data['atom_unit_cell'][:]
                atom_site_number = nx_data['atom_site_number'][:]
                
                # Flags and metadata handling
                # Load flags that start with 'status_flag_' and strip the prefix
                status_flags = {key[12:]: nx_data.attrs[key] for key in nx_data.attrs if key.startswith('status_flag_')}
                #
                ## # Load metadata 
                metadata = (nx_data['audit_conform_dict_name'][()].decode('utf-8'),
                            nx_data['audit_conform_dict_version'][()].decode('utf-8'),
                            nx_data['audit_creation_date'][()].decode('utf-8'),
                            nx_data['audit_creation_method'][()].decode('utf-8'),
                            nx_data['audit_author_name'][()].decode('utf-8'))
                
                # Optional parts loading with checks
                # Load average structure                
                if 'average_type' in nx_data and 'average_occ' in nx_data and 'average_pos' in nx_data and \
                   'average_adp' in nx_data and 'average_site' in nx_data:
                    average_type = np.array(nx_data['average_type'], dtype=np.int32)
                    average_pos  = np.array(nx_data['average_pos' ], dtype=np.float64)
                    average_occ  = np.array(nx_data['average_occ' ], dtype=np.float64)
                    average_adp  = np.array(nx_data['average_adp' ], dtype=np.float64)
                    average_site = np.array(nx_data['average_site'], dtype=np.int32)
                    average_n_atoms = average_type.shape[0]
                    average_structure = (average_n_atoms, average_type, average_pos, average_occ, average_adp, average_site)
                    #average_structure = nx_data['average_structure'][:]
                else:
                    average_structure = None
                
                # Load molecules
                if 'molecules_number' in nx_data and 'molecules_types' in nx_data and 'molecules_int' in nx_data \
                    and 'molecules_real' in nx_data and 'molecules_index' in nx_data:
                    molecules_number = int(nx_data['molecules_number'][0])
                    molecules_types = int(nx_data['molecules_types'][0])
                    molecules_int = np.array(nx_data['molecules_int'][:][:], dtype=np.int32)
                    molecules_real = np.array(nx_data['molecules_real'][:][:], dtype=np.float64)
                    molecules_index = np.array(nx_data['molecules_index'][:][:], dtype=np.int32)
                    #molecules_int = np.array(nx_data['molecules_int'][:][:][0], dtype=np.int32)
                    #molecules_real = np.array(nx_data['molecules_real'][:][:][0], dtype=np.float64)
                    #molecules_index = np.array(nx_data['molecules_index'][:][:][0], dtype=np.int32)
                    molecules = (molecules_number, molecules_types, molecules_int, molecules_real, molecules_index)
                    #print("READ MOLECULES ")
                    #print("INT  MOLECULES ", molecules_int)
                    #print("REAL MOLECULES ", molecules_real)
                    #print("INDX MOLECULES ", molecules_index)
                else:
                    molecules_number = 0
                    molecules_types = 0
                    #print("NO   MOLECULES IN FILE")
                    molecules = None

                # Load occupancy
                if 'occupancy' in nx_data:
                    occupancy = np.array(nx_data['occupancy'][:], dtype=np.float64)
                else:
                    occupancy = None

                # Load atom_occupancy flags
                if 'atom_occupancy' in nx_data:
                    atom_occupancy = nx_data['atom_occupancy'][:]
                else:
                    atom_occupancy = None
                    #print("READ atom_occupancy none    "), atom_occupancy


                # Load magnetic spin
                if 'magnetic_spins' in nx_data:
                    magnetic_spins = nx_data['magnetic_spins'][:]
                else:
                    magnetic_spins = None
                
                # Load property flags
                if 'property_flags' in nx_data:
                    property_flags = nx_data['property_flags'][:]
                else:
                    property_flags = None
                    #print("READ PROPERTY none    "), property_flags

                # Load anisotropic ADP
                if 'anisotropic_adp' in nx_data and 'anisotropic_index' in nx_data:
                    anisotropic_is_iso  = np.array(nx_data['anisotropic_is_iso'][:], dtype=np.int32)
                    anisotropic_adp     = np.array(nx_data['anisotropic_adp'][:][:], dtype=np.float64)
                    anisotropic_index   = np.array(nx_data['anisotropic_index'][:], dtype=np.int32)
                    anisotropic_n_uij   = anisotropic_adp.shape[0]
                    anisotropic_n_type  = anisotropic_adp.shape[1]
                    anisotropic_n_index = anisotropic_index.shape[0]
                    anisotropic_tuple = (anisotropic_n_uij, anisotropic_n_type, anisotropic_n_index, anisotropic_is_iso, anisotropic_adp, anisotropic_index)
                else:
                    anisotropic_tuple = None
                   
                return cls(unit_cell_lengths=unit_cell_lengths,
                           unit_cell_angles=unit_cell_angles,
                           symmetry_H_M=symmetry_H_M,
                           symmetry_origin=symmetry_origin,
                           symmetry_abc=symmetry_abc,
                           symmetry_n_mat=symmetry_n_mat,
                           symmetry_mat=symmetry_mat,
                           unit_cells=unit_cells, 
                           number_of_types=number_of_types,
                           types_names=types_names,
                           types_ordinal=types_ordinal,
                           types_charge=types_charge,
                           types_isotope=types_isotope,
                           number_of_atoms=number_of_atoms,
                           atom_type=atom_type,
                           atom_position=atom_position,
                           atom_unit_cell=atom_unit_cell,
                           atom_site_number=atom_site_number,
                           status_flags=status_flags,
                           metadata=metadata,
                           average_structure=average_structure,
                           magnetic_spins=magnetic_spins,
                           property_flags=property_flags,
                           anisotropic_tuple=anisotropic_tuple,
                           molecules=molecules,
                           occupancy=occupancy,
                           atom_occupancy=atom_occupancy)
        
        except KeyError as e:
            raise KeyError(f"Expected dataset not found in the file: {e}")
        except Exception as e:
            raise Exception(f"An error occurred while loading from the NeXus file: {e}")

#
##########################################################################################
#

def write_diffuse_structure(file_path, \
            creation_method, \
            author_name, \
            unit_cell_lengths, \
            unit_cell_angles, \
            symmetry_H_M, \
            symmetry_origin, \
            symmetry_abc, \
            symmetry_n_mat, \
            symmetry_mat, \
            unit_cells, \
            types_number, \
            types_names, \
            types_ordinal, \
            types_charge, \
            types_isotope, \
            number_of_atoms, \
            atom_type, \
            atom_position, \
            atom_unit_cell, \
            atom_site_number, \
            **kwargs):
    #
    # Get optional parameters from **kwargs
    #
    #print("IN WRITE_DIFFUSE_STRUCTURE")
    status_flags = kwargs.get('status_flags', None )
    average_structure = kwargs.get('average_structure', None )
    magnetic_spins = kwargs.get('magnetic_spins', None )
    property_flags = kwargs.get('property_flags', None )
    anisotropic_tuple = kwargs.get('anisotropic_tuple', None )
    molecules = kwargs.get('molecules', None )
    occupancy = kwargs.get('occupancy', None )
    atom_occupancy = kwargs.get('atom_occupancy', None )
    #
    #print("Incoming shapes ")
    #print("Unit cell length, angles ", (unit_cell_lengths.shape), (unit_cell_angles.shape) )
    #print("Symmetry matrices N, sh  ", symmetry_n_mat, (symmetry_mat.shape     ))
    #print("unit_cells               ", (unit_cells.shape       ))
    #print("Types n, o, c, i         ", types_number, (types_ordinal.shape), (types_charge.shape ), (types_isotope.shape  ))
    #print("Atom N type, pos         ", number_of_atoms, (atom_type.shape ), (atom_position.shape)     , (atom_unit_cell.shape ), (atom_site.shape) )
    #if average_structure is not None:
       #print("Average_structure TPOAS  ", average_structure[0].shape, (average_structure[1].shape ), (average_structure[2].shape ), (average_structure[3].shape ), (average_structure[4].shape ))
    #if property_flags    is not None:
       #print("Propery_flags            ", (property_flags.shape))
    #if magnetic_spin     is not None:
       #print("Magnetic_spin            ", (magnetic_spin.shape ))
    #if occupancy         is not None:
       #print("Occupancy                ", (occupancy.shape     ))
    #if anisotropic_tuple   is not None:
       #print("Anisotropic_adp ADP, IND ", (anisotropic_tuple[0].shape), (anisotropic_tuple[1].shape ))
    #if molecules         is not None:
       #print("Molecules N, T, i,f, IND ", molecules[0], molecules[1], (molecules[2].shape     ), (molecules[3].shape       ), (molecules[4].shape) )
    #print("")
    #print(" status_flags      ", status_flags     )
    #print(" average_structure ", average_structure)
    #print(" magnetic_spin     ", magnetic_spin    )
    #print(" property_flags    ", property_flags   )
    #
    # Define metadata structure
    dict_name = "Disorder structure"
    dict_version = "0.0.0"
    creation_date = datetime.today().strftime('%Y-%m-%d')
    metadata1 = ( dict_name,
             dict_version,
             creation_date,
             creation_method,
             author_name
    )
    print("METADATA ::", metadata1)
    metadata = {'audit_conform_dict_name': dict_name,
            'audit_conform_dict_version': dict_version,
            'audit_creation_date': creation_date,
            'audit_creation_method': creation_method,
            'audit_author_name': author_name
    }
    #
    #
    #print(" ABOUT TO INITIALIZE ")
    # Initialize instance with all fields
    nexus = NeXusStructureFormat(
            unit_cell_lengths,
            unit_cell_angles,
            symmetry_H_M,
            symmetry_origin,
            symmetry_abc,
            symmetry_n_mat,
            symmetry_mat,
            unit_cells,
            types_number,
            types_names,
            types_ordinal,
            types_charge,
            types_isotope,
            number_of_atoms,
            atom_type,
            atom_position,
            atom_unit_cell,
            atom_site_number,
            status_flags,
            metadata1,
            average_structure,
            magnetic_spins,
            property_flags,
            anisotropic_tuple,
            molecules,
            occupancy,
            atom_occupancy
    )
    # Save into actual NeXus file
    nexus.save_to_nexus(file_path)
    #
    return 0

#
##########################################################################################
#

def read_diffuse_structure(file_path):
    #
    loaded = NeXusStructureFormat.load_from_nexus(file_path)
    #
    #print(" FINISHED LOADING ")
    #
    #
    #print("Outgoing shapes ")
    #print("Unit cell length, angles ", (loaded.unit_cell_lengths.shape), (loaded.unit_cell_angles.shape) )
    #print("Symmetry matrices N, sh  ", loaded.symmetry_n_mat, (loaded.symmetry_mat.shape     ))
    #print("unit_cells               ", (loaded.unit_cells.shape       ))
    #print("Types N, len             ", loaded.number_of_types, len(loaded.atom_types) )
    #print("Atom N type, pos         ", loaded.number_of_atoms)
    #if loaded.average_structure is not None:
    #   print("Average_structure TPOAS  ", loaded.average_structure[0], loaded.average_structure[1].shape, (loaded.average_structure[2].shape ), (loaded.average_structure[3].shape ), (loaded.average_structure[4].shape ), (loaded.average_structure[5].shape ))
    #if loaded.property_flags    is not None:
    #   print("Propery_flags            ", (loaded.property_flags.shape))
    #if loaded.magnetic_spin     is not None:
    #   print("Magnetic_spins            ", (loaded.magnetic_spins.shape ))
    #if loaded.occupancy         is not None:
    #   print("Occupancy                ", (loaded.occupancy.shape     ))
    #if loaded.anisotropic_tuple   is not None:
    #   print("Anisotropic_tuple ADP, IND ", loaded.anisotropic_tuple[0], loaded.anisotropic_tuple[1], loaded.anisotropic_tuple[2], (loaded.anisotropic_tuple[3].shape), (loaded.anisotropic_tuple[4].shape ))
    #if loaded.molecules         is not None:
    #   print("Molecules N, T, i,f, IND ", loaded.molecules[0], loaded.molecules[1], (loaded.molecules[2].shape     ), (loaded.molecules[3].shape       ), (loaded.molecules[4].shape) )
    #print("")
    #
    return loaded.unit_cell_lengths, \
           loaded.unit_cell_angles, \
           loaded.symmetry_H_M, \
           loaded.symmetry_origin, \
           loaded.symmetry_abc, \
           loaded.symmetry_n_mat, \
           loaded.symmetry_mat, \
           loaded.unit_cells, \
           loaded.number_of_types, \
           loaded.types_names, \
           loaded.types_ordinal, \
           loaded.types_charge, \
           loaded.types_isotope, \
           loaded.number_of_atoms, \
           loaded.atom_type,\
           loaded.atom_position,\
           loaded.atom_unit_cell,\
           loaded.atom_site_number,\
           loaded.status_flags, \
           loaded.metadata, \
           loaded.average_structure, \
           loaded.magnetic_spins, \
           loaded.property_flags, \
           loaded.anisotropic_tuple, \
           loaded.molecules, \
           loaded.occupancy, \
           loaded.atom_occupancy
