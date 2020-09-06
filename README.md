GPU enabled RSoXS simulation (Release 1.0 - Beta version)
====================================
This code solves the X-ray scattering on Graphical processing units (GPU).
Currently it supports the operation with 32 bit / 64 bit float datatypes.




Dependencies
=============
#### Required Dependencies
* A C++ compiler with C++11 support is required. We aim for compatibility with gcc and icc (Intel C++ compiler).
* Nvidia - CUDA
* Cuda Toolkit
* HDF5
* OpenMP

#### Additional dependencies for building with Pybind
* Python >= 3.6

#### Additional dependencies for building without Pybind
* libconfig

#### Optional Dependencies
* Doxygen
* Docker
* OpenCV

Running Cy-RSoXS with docker
============================

The docker for Cy-RSoXS is available at [Docker Hub](https://hub.docker.com/repository/docker/maksbh/cy-rsoxs/general).
This comes installed with nvidia toolkit and other dependencies.
See [docs/DOCKER.md](docs/DOCKER.md) for detailed instructions on how to build Cy-RSoXS with the docker image.


Installation Instruction without Docker
======================================
See [docs/INSTALL.md](docs/INSTALL.md) for detailed instructions on how to build Cy-RSoXS and its dependencies.

Running Cy-RSoXS
================
See [docs/RUN.md](docs/RUN.md) for detailed instructions on how to run Cy-RSoXS.

Visualization
==============
See [docs/Visualiztion.md](docs/Visualization.md) for detailed instructions on how to visualize the output of Cy-RSoXS.

Contributors
============
* Kumar Saurabh
* Adarsh Krishnamurthy
* Baskar Ganapathysubramanian
* Eliot Gann
* Dean Delongchamp
* Michael Chabinyc

Acknowledgement
===============
We thank ONR MURI Center for Self-Assembled Organic Electronics for providing the support for this work.

Contact
=======
Questions and comments should be sent to Dr. Baskar Ganapathysubramanian at [baskarg@iastate.edu](mailto:baskarg@iastate.edu) or Dr.  Adarsh Krishnamurthy [adarsh@iastate.edu](mailto:adarsh@iastate.edu).
