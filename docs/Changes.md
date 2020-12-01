Cy-RSoXS Changes History
====================================
Version 0.9.0
=================
* Energies now accepted as lists
* Size of `BigUINT` dropped to 32 bit. Compilation option added to enable 64 bit.
* FFT made in place.
* kRotation added (still Experimental)
* Scatter Approach made partial 
* Pybind function names modified for Electric field rotations
* Config modified. Now rotations accepted as list `[start:increment:end]`
* Morphology can be added as vector morphology or Euler Angles

Version 0.8.1
=================
* Proper memory destruction for cuda variables. 
Some CUDA allocated variables was not properly released at the end of execution.  
* Minimum gcc version required raised to 7.
* The scattering pattern memory is allocated during constructor. The `allocate` function is removed.


Version 0.8.0
=================

* Pybind support added
* The code has been made consistent with $\vec{k}$ along $\vec{Z}$  
* Fixed a memory error in Interpolation
 