Cy-RSoXS: GPU enabled RSoXS simulation
====================================
* This code calculates resonant X-ray scattering in the Born Approximation on Graphical Processing Units (GPU).
Currently it supports the operation with 32 bit / 64 bit float datatypes.

* The code can be executed through the Python interface (enabled through [Pybind11](https://github.com/pybind/pybind11))
or directly through the executable.


Version  : 1.1.4.0
================
#### The developer version of the CyRSoXS code is available at https://bitbucket.org/baskargroup/cy-rsoxs, and all stable future releases will be merged with the current usnistgov repository.

Dependencies
=============
#### Required Dependencies
* A C++ compiler with C++14 support is required.
* gcc >= 7 (CUDA specific versions might have GCC requirements )
* Cuda Toolkit (>=9)
* HDF5
* OpenMP

#### Additional dependencies for building with Pybind
* Python >= 3.6

#### Additional dependencies for building without Pybind
* libconfig

#### Optional Dependencies
* Doxygen
* Docker

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

This software was developed at Iowa State University in collaboration with NIST. The Iowa State team provided expertise in high performance computing, and the NIST team provided scientific expertise on the RSoXS technique.

#### Iowa State University
* Kumar Saurabh
* Adarsh Krishnamurthy
* Baskar Ganapathysubramanian

#### NIST
* Eliot Gann
* Dean M. DeLongchamp
* Peter J. Dudenas
* Tyler B. Martin



Acknowledgement
===============
We thank ONR MURI Center for Self-Assembled Organic Electronics for providing the support for this work.

Frequently Asked Question
============================
See [docs/FAQ.md](docs/FAQ.md) for some of the known issues or Frequently asked questions.

Contact
=======
Questions and comments should be sent to Dr. Baskar Ganapathysubramanian at [baskarg@iastate.edu](mailto:baskarg@iastate.edu), Dr.  Adarsh Krishnamurthy [adarsh@iastate.edu](mailto:adarsh@iastate.edu), or Dr. Dean DeLongchamp at [dean.delongchamp@nist.gov](mailto:dean.delongchamp@nist.gov)
