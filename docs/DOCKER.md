# GPU enabled RSoXS simulation (1.1.6.0)

* **NOTE:** The current Docker image was created from the developer version of CyRSoXS at <https://bitbucket.org/baskargroup/CyRSoXS>, and the following installation instructions reflect this.

## Installing Docker

On Ubuntu you can use the command to install docker:

```bash
sudo apt-get install docker
sudo apt-get install docker.io
```

## Launching Docker

`sudo docker pull maksbh/CyRSoXS:latest`

To get the imageID run:

`sudo docker images`

It will show the output as:

```bash
REPOSITORY          TAG                 IMAGE ID            CREATED             SIZE
maksbh/CyRSoXS     latest              655455b309d5        17 minutes ago      4.77GB
```

## Installing CyRSoXS

First you need to run the docker interactively:

`sudo docker run -it $IMAGE_ID`

Then run the following command:

```bash
cd; ./configure;
```

Note that you might need to change the configuration file. The contents are:

```bash
#clone the git repository
git clone https://bitbucket.org/baskargroup/cy-rsoxs.git

# Go into the directory
cd CyRSoXS;

# Update the PyBind submodule
git submodule update --init --recursive

# Make the build directory
mkdir build; cd build;

#Cmake command
cmake .. -DCMAKE_BUILD_TYPE=Release -DHDF5_DIR=/root/Dependencies/CMake-hdf5-1.10.5/build/_CPack_Packages/Linux/TGZ/HDF5-1.10.5-Linux/HDF_Group/HDF5/1.10.5/share/cmake/hdf5 -DDLEVEL2=Yes -DNUM_MATERIAL=4

#Make the executable
make
```

## Contributors

This software was developed at Iowa State University in collaboration with NIST. The Iowa State team provided expertise in high performance computing, and the NIST team provided scientific expertise on the RSoXS technique.

### Iowa State University

* Kumar Saurabh
* Adarsh Krishnamurthy
* Baskar Ganapathysubramanian

### NIST

* Eliot Gann
* Dean M. DeLongchamp
* Peter J. Dudenas
* Tyler B. Martin
* Peter Beaucage

## Acknowledgement

We thank the Office of Naval Research Multidisciplinary University Research Initiative (ONR MURI) Center for [Self-Assembled Organic Electronics](http://www.mri.psu.edu/mri/facilities-and-centers/soe) for providing support for this work.

## Contact

Questions and comments should be sent to Dr. Baskar Ganapathysubramanian at [baskarg@iastate.edu](mailto:baskarg@iastate.edu) or Dr.  Adarsh Krishnamurthy at [adarsh@iastate.edu](mailto:adarsh@iastate.edu) or Dr. Dean DeLongchamp at [dean.delongchamp@nist.gov](mailto:dean.delongchamp@nist.gov)
