# CyRSoXS Changes History

## Version 1.1.5.2

* Add option to install build dependencies using conda rather than manual
* Changes to CMakeLists to enable Conda packaging of CyRSoXS

## Version 1.1.4.0

* NumMaterial requirement removed from compilation
* Euler Angle sequence changed to ZYZ
* CMake List updated
* Image rotation bug fixed

## Version 1.0.1 - Beta

* Axis label added in all HDF5 file output

## Version 1.0.0 - Beta

* Support for k rotation added
* Morphology dimensions now read from HDF5 file
* HDF5 file format changed to be more rigorous
* Euler angle support added
* Material information now starts from 1 to be consistent with HDF5

## Version 0.9.0

* Energies now accepted as lists
* Size of `BigUINT` dropped to 32 bit. Compilation option added to enable 64 bit.
* FFT made in place.
* kRotation added (still Experimental)
* Scatter Approach made partial
* Pybind function names modified for Electric field rotations
* Config modified. Now rotations accepted as list `[start:increment:end]`
* Morphology can be added as vector morphology or Euler Angles

## Version 0.8.1

* Proper memory destruction for cuda variables.
Some CUDA allocated variables was not properly released at the end of execution.  
* Minimum gcc version required raised to 7.
* The scattering pattern memory is allocated during constructor. The `allocate` function is removed.

## Version 0.8.0

* Pybind support added
* The code has been made consistent with $\vec{k}$ along $\vec{Z}$  
* Fixed a memory error in Interpolation
