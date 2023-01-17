# CyRSoXS: GPU enabled RSoXS simulation

* CyRSoXS is a GPU-accelerated codebase that calculates resonant X-ray scattering in the Born Approximation. It takes a voxel-based model and optical constants as input, and returns simulated X-ray scattering patterns. These models can be derived from real space measurements (AFM, TEM, 4DSTEM), or synthetically generated using procedural morphology generators.
Currently, it supports both 32 bit and 64 bit float data types.

* The code can be executed through Python (enabled through [Pybind11](https://github.com/pybind/pybind11))
or directly through the executable.

## Version  : 1.1.5.2

The developer version of the CyRSoXS code is available at <https://bitbucket.org/baskargroup/cy-rsoxs>. All stable future releases will be merged with the current usnistgov repository.

## Dependencies

### Hardware

* CUDA GPU with [Compute Capability](https://docs.nvidia.com/deploy/cuda-compatibility/index.html) >=7
  * CUDA is a parallel computing platform and programming model creaded by NVIDIA

### Required Dependencies for building from source

* A C++ compiler with C++14 support is required.
* gcc >= 7 (CUDA specific versions might have GCC requirements)
* CUDA Toolkit (>=9)
* HDF5
  * Hierarchical Data Format (version 5) is a set of file formats designed to store and organize large amounts of data. CyRSoXS uses it to store models and simulated scattering patterns
* OpenMP
* libconfig==1.7.2

### Additional dependencies for building with Pybind

* Python >= 3.6

### Optional Dependencies

* Doxygen - for building the documentation.
* Docker - if you want to run CyRSoXS in a Docker container.

## [Installation Instructions](docs/INSTALL.md)

## [Simulation Componenets](docs/DATA.md)

## [Running CyRSoXS](docs/RUN.md)

## [Visualization](docs/Visualization.md)

## Documentation

CyRSoXS documentation is available at <https://cyrsoxs.readthedocs.io>, and instructions to build the documentation are in the [Installation Instructions](docs/INSTALL.md).

## Running CyRSoXS with Docker

The Docker image for CyRSoXS is available at [Docker Hub](https://hub.docker.com/r/maksbh/cy-rsoxs).
This comes installed with Nvidia toolkit and other dependencies.
See [docs/DOCKER.md](docs/DOCKER.md) for detailed instructions on how to build CyRSoXS with the Docker image.

## Contributors

This software was developed at Iowa State University in collaboration with NIST. The Iowa State team provided expertise in high performance computing, and the NIST team provided scientific expertise on the RSoXS technique.

### Iowa State University

* Kumar Saurabh, [Github](https://github.com/KumarSaurabh1992), [ORCiD](https://orcid.org/0000-0003-2503-367X), [Email](maksbh@iastate.edu)
* Adarsh Krishnamurthy, [Webpage](https://web.me.iastate.edu/idealab/p-krishnamurthy.html), [ORCiD](https://orcid.org/0000-0002-5900-1863), [Email](adarsh@iastate.edu)
* Baskar Ganapathysubramanian, [Webpage](https://bitbucket.org/baskargroup/), [ORCiD](https://orcid.org/0000-0002-8931-4852), [Email](baskarg@iastate.edu)

### NIST

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

## NIST Disclaimer

Any mention of commercial products is for information only; it does not imply recommendation or endorsement by NIST.
