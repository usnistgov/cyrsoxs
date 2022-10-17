# CyRSoXS: GPU enabled RSoXS simulation

* This code calculates resonant X-ray scattering in the Born Approximation on Graphics Processing Units (GPU).
Currently it supports the operation with 32 bit / 64 bit float datatypes.

* The code can be executed through the Python interface (enabled through [Pybind11](https://github.com/pybind/pybind11))
or directly through the executable.


## Version  : 1.1.5.0

#### The developer version of the CyRSoXS code is available at https://bitbucket.org/baskargroup/CyRSoXS, and all stable future releases will be merged with the current usnistgov repository.

## Dependencies

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

## [Installation Instructions](docs/INSTALL.md)

## [Running CyRSoXS](docs/RUN.md)

## [Visualization](docs/Visualization.md)

## Running CyRSoXS with Docker

The Docker image for CyRSoXS is available at [Docker Hub](https://hub.docker.com/repository/docker/maksbh/CyRSoXS/general).
This comes installed with Nvidia toolkit and other dependencies.
See [docs/DOCKER.md](docs/DOCKER.md) for detailed instructions on how to build CyRSoXS with the Docker image.


## Contributors

This software was developed at Iowa State University in collaboration with NIST. The Iowa State team provided expertise in high performance computing, and the NIST team provided scientific expertise on the RSoXS technique.

#### Iowa State University
* Kumar Saurabh, [Github](https://github.com/KumarSaurabh1992), [ORCiD](https://orcid.org/0000-0003-2503-367X), [Email](maksbh@iastate.edu)
* Adarsh Krishnamurthy, [Webpage](https://web.me.iastate.edu/idealab/p-krishnamurthy.html), [ORCiD](https://orcid.org/0000-0002-5900-1863), [Email](adarsh@iastate.edu)
* Baskar Ganapathysubramanian, [Webpage](https://bitbucket.org/baskargroup/), [ORCiD](https://orcid.org/0000-0002-8931-4852), [Email](baskarg@iastate.edu)

#### NIST
* Eliot Gann, [Github](https://github.com/EliotGann), [ORCiD](https://orcid.org/0000-0001-5570-8880), [Email](eliot.gann@nist.gov)
* Dean M. DeLongchamp, [Webpage](https://www.nist.gov/people/dean-delongchamp), [ORCiD](https://orcid.org/0000-0003-0840-0757), [Email](dean.delongchamp@nist.gov)
* Peter J. Dudenas, [Github](https://github.com/pdudenas), [ORCiD](https://orcid.org/0000-0002-4578-4182), [Email](peter.dudenas@nist.gov)
* Tyler B. Martin, [Github](https://github.com/martintb), [ORCiD](https://orcid.org/0000-0001-7253-6507), [Email](tyler.martin@nist.gov)
* Peter Beaucage, [Github](https://github.com/pbeaucage), [ORCiD](https://orcid.org/0000-0002-2147-0728), [Email](peter.beaucage@nist.gov)



## Acknowledgement

We thank the Office of Naval Research Multidisciplinary University Research Initiative (ONR MURI) Center for [Self-Assembled Organic Electronics](http://www.mri.psu.edu/mri/facilities-and-centers/soe) for providing support for this work.

## Frequently Asked Questions

See [docs/FAQ.md](docs/FAQ.md) for some of the known issues or Frequently asked questions.

## Contact

Questions and comments should be sent to Dr. Baskar Ganapathysubramanian at [baskarg@iastate.edu](mailto:baskarg@iastate.edu), Dr.  Adarsh Krishnamurthy [adarsh@iastate.edu](mailto:adarsh@iastate.edu), or Dr. Dean DeLongchamp at [dean.delongchamp@nist.gov](mailto:dean.delongchamp@nist.gov)
