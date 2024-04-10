Proposal for a common 3D-Diffuse / 3D PDF data structure
 
The intensity/3D-PDF is essentially
- A 3D voxel of real valued data points
- It is based on a non-cartesian metric
- Usually, it is based on an even spaced grid (different steps sizes along a/b/c are possible)
 
For the Nexus file format:
The only required field is a NXentry “These contain all the data that is required to describe an experimental run or scan” → This would be useful if we have several structures/datasets that are afterwards averaged. The default would only be one entry. 

NXdata is one level below the NXentry. “contain the experimental results in a self-contained way, i.e., it should be possible to generate a sensible plot of the data from the information contained in each NXdata group.”

NXsample and NXinstrument are likely not necessary in our case. We could think about making a NXProgram that dumps important program information.

Necessary entries:

NXentry/title:
Content: title of the data
Type: H5T_STRING

NXentry/NXdata/type: 
Content type of data, options could be, e.g.  experimental data, intensity, real part, imaginary part, structure factor,  …)
Type: H5T_STRING

NXentry/NXdata/space: 
Content: real space or reciprocal space data 
Type: boolean?  

NXentry/NXdata/dimension: 
Content: contains the dimension of the data: 1d, 2d, 3d 
Type: integer

NXentry/NXdata/data: 
Content:The actual data array nd array of  
Type: H5T_NATIVE_DOUBLE

NXentry/NXdata/lower_limits:  
Content: the lower left bottom corner in terms of a*/b*/c*  or a/b/c
Type: 3 H5T_NATIVE_DOUBLE

NXentry/NXdata/step_sizes:  
Content: increment along the n axis in terms of  a*/b*/c*  or a/b/c
Type: n*3 H5T_NATIVE_DOUBLE

NXentry/NXdata/unit_cell:  
Content: Metric tensor, alternatively a,b,c, alpha, beta, gramma
Type: 6 H5T_NATIVE_DOUBLE

NXentry/NXdata/radiation:  
Content: radiation type used: x-ray, electron, neutron
Type: H5T_STRING

NXentry/NXdata/wavelength:  
Content: Wavelength
Type:  H5T_NATIVE_DOUBLE

NXentry/NXprogram/program_version:
Content: program used to generate the data with version number
Type: H5T_STRING

Optional entries:
NXentry/NXdata/symmetry:  
Content: Symmetry that has been applied to the data, specified as Laue symmetry
Type:  H5T_STRING

NXentry/NXdata/adp:  
Content: Specifies if adp have been used in the calculation, possible value: non, isotropic, anisotropic
Type:  H5T_STRING

NXentry/NXdata/dispersion:  
Content: Specifies if dispersion correction was applied
Type:  Boolean

NXentry/NXprogram/formfactors: 
Content: Specifies a link to the form factors that have been used, e.g. which gaussian approximation
Type:  H5T_STRING

NXentry/NXprogram/averaging: 
Content: Specifies the type of averaging used, e.g. Lots or Lanczos filter or similar
Type:  H5T_STRING

NXentry/NXprogram/approximation: 
Content: Specifies approximations used in the calculation, e.g. fft, level of approximation in nufft….
Type:  H5T_STRING
 
NXentry/NXdata/structurefactor: 
Content: complex structure factor
Type: Is there a complex pre-defined h5 type?

NXentry/NXdata/mask: 
Content: mask applied to the data
Type: H5T_NATIVE_DOUBLE

NXentry/NXdata/kinematic: 
Content: Kinematic or dynamical scattering intensity
Type: Boolean

NXentry/NXdata/elastic: 
Content: Elastic or integral scattering intensity
Type: Boolean

