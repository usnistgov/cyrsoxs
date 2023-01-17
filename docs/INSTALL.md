# CyRSoXS Installation Instructions (1.1.5.2)

## Conda

CyRSoXS is available as a pre-built binary on the `conda-forge` channel. To install:

```bash
conda install cyrsoxs -c conda-forge
```

## Building from source

To build CyRSoXS from source, the following dependencies need to be installed:

* A C++ compiler with C++14 support is required.
* gcc >= 7 (CUDA specific versions might have GCC requirements )
* Cuda Toolkit (>=9)
* HDF5
* OpenMP

### Additional dependencies for building with Pybind

* Python >= 3.6

### Additional dependencies for building without Pybind

* libconfig

### Dependencies for building the documentation

* Doxygen
* Graphviz
* latex

## Compiling libconfig

Libconfig's build system is [Autotools](https://www.gnu.org/software/automake/manual/html_node/Autotools-Introduction.html), which means you'll need to run `./configure` and then `make` to build. Hyperrealm is the maintainer of Libconfig

This guide recommends passing ```--prefix=`pwd`/install``` to `./configure`, which will cause `make install` to copy the output files to `[your_libconfig_dir]/install` instead of `/usr`. This way your libconfig install lives completely inside your libconfig folder. This is necessary if you are working on a system where you don't have admin privileges (i.e. an HPC cluster).

```bash
# Download and verify checksum
wget https://hyperrealm.github.io/libconfig/dist/libconfig-1.7.2.tar.gz
sha256sum libconfig-1.7.2.tar.gz
# 7c3c7a9c73ff3302084386e96f903eb62ce06953bb1666235fac74363a16fad9

# Extract
tar xvf libconfig-1.7.2.tar.gz
rm libconfig-1.7.2.tar.gz

# Compile and copy output files to libconfig-1.7.2/install
cd libconfig-1.7.2
./configure --prefix=`pwd`/install
make -j8  # compile with 8 threads
make install

# Permanently set $LIBCONFIG_DIR environment variable, and set LD_LIBRARY_PATH to include the
# libconfig lib directory to prevent dynamic linking errors with libconfig++.so.
echo "export LIBCONFIG_DIR=`pwd`/install" >> ~/.bashrc
echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$LIBCONFIG_DIR/lib" >> ~/.bashrc
source ~/.bashrc
```

**NOTE:** On some HPC clusters (when using the Intel compiler), the `make` step gives a linker error. This is the libconfig example program failing to link with the Intel runtime. This is okay - the libconfig library itself compiles just fine. Just run `make install` and double check that `install/lib` contains some `*.a` files.

## Installing HDF5

CyRSoXS needs [HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format)-based binary input format.

```bash
cd $HDF5_INSTALL_DIRECTORY
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.5/src/CMake-hdf5-1.10.5.tar.gz
sha256sum CMake-hdf5-1.10.5.tar.gz
# 339bbc4594b6d71ed0794b0861af231bdd06bcc71c8a81563763f72455c1c5c2

tar -xzvf CMake-hdf5-1.10.5.tar.gz
rm CMake-hdf5-1.10.5.tar.gz
cd CMake-hdf5-1.10.5
./build-unix.sh
```

This step might take some time. Do not cancel until all the tests have passed.
This step will create a cmake files in location `$HFD5_DIR/build/_CPack_Packages/Linux/TGZ/HDF5-1.10.5-Linux/HDF_Group/HDF5/1.10.5/share/cmake/hdf5`

Export the path for HDF5:

```bash
cd build/_CPack_Packages/Linux/TGZ/HDF5-1.10.5-Linux/HDF_Group/HDF5/1.10.5/share/cmake/hdf5;
echo "export HDF5_DIR=`pwd`" >> ~/.bashrc
source ~/.bashrc
```

## Downloading CyRSoXS

Clone the CyRSoXS Github repo

```bash
git clone https://github.com/usnistgov/cyrsoxs.git
```

### With Pybind

If you want to use the Python support for CyRSoXS, add the submodule by

```bash
git submodule update --init
```

## Building CyRSoXS

One should have a valid C/C++ compiler and CUDA-toolkit in addition to above modules installed, in order to
compile CyRSoXS.

To compile CyRSoXS as an executable that can be called from the command line or from a bash script, run the following commands:

```bash
cd $CyRSoXS_DIR
mkdir build;
cd build;
cmake .. -DCMAKE_BUILD_TYPE=Release
```

This will generate a `CyRSoXS` executable.

### Building with Pybind

CyRSoXS can also be compiled with Pybind so that it can be imported as a Python library. To compile for this functionality, run the following commands:

```bash
cd $CyRSoXS_DIR
mkdir build-pybind
cd build-pybind
cmake .. -DCMAKE_BUILD_TYPE=Release -DPYBIND=Yes -USE_SUBMODULE_PYBIND=Yes
```

This will generate a `CyRSoXS.so` Shared Library file, which can be imported into Python. Note that this does not create an executable.

#### Optional CMake Flags

```console
    
    -DPYBIND=Yes            # Compiling with Pybind
    -DUSE_PYBIND_SUBMODULE  # Choose to compile with the Pybind submodule, or Conda-installed Pybind 
    -DMAX_NUM_MATERIAL=64   # To change the maximum number of materials (default is 32) 
    -DDOUBLE_PRECISION=Yes  # Calculations will performed with double precision numbers
    -DPROFILING=Yes         # Enables profiling of the code
    -DBUILD_DOCS=Yes        # To build documentation
    -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc # Compiling with the Intel compiler (does not work with Pybind)
    -DOUTPUT_BASE_NAME=CyRSoXS # Changes the name of the built output binary or Python module (if using Pybind)
```

## Making CyRSoXS

Once the CMake files has been generated run the following command:

```bash
make
```

If `-DBUILD_DOCS=Yes`, the make command will build the documentation in html and latex located in `$CyRSoXS_DIR/build/html` and `$CyRSoXS_DIR/build/latex`, respectively. A PDF of the documentation is also built as `$CyRSoXS_DIR/build/latex/CyRSoXS_Manual.pdf`.