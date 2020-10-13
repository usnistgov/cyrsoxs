Cy-RSoXS Changes History
====================================

Version 0.8.1
=================
* Proper memory destruction for cuda variables. 
Some CUDA allocated variables was not properly released at the end of execution.  
* Minimum gcc version required raised to 7.


Version 0.8.0
=================

* Pybind support added
* The code has been made consistent with $\vec{k}$ along $\vec{Z}$  
* Fixed a memory error in Interpolation
 