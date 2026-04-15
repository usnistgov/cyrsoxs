# CyRSoXS Installation Instructions (1.2.0)

## Conda (recommended for binaries)

CyRSoXS is available as a pre-built binary on the `conda-forge` channel. To install:

```bash
conda install cyrsoxs -c conda-forge
```

Pre-built binaries on conda-forge include compatible CUDA, HDF5, and compiler stacks. PyPI may publish source distributions and limited wheels; GPU-capable wheels are constrained by CUDA manylinux policy, so conda-forge remains the primary binary channel for most users.

## pip (experimental)

Linux installs can be build from wheels or from source directly from pip.
```bash
pip install cyrsoxs
# or
uv pip install cyrsoxs
```
For all dependencies needed to run the tutorial notebook, install the `jupyter` extra:
```bash
uv pip install cyrsoxs[jupyter]
```

## Python package build with uv (PEP 517)

The repository is a [scikit-build-core](https://scikit-build-core.readthedocs.io/) project. The native extension module is named `CyRSoXS` (import `import CyRSoXS`).

System requirements for a local build:

- C++17 compiler and NVIDIA CUDA toolkit (`nvcc` on `PATH`)
- OpenMP

**HDF5:** For `uv sync`, `pip install .`, and other **scikit-build-core** (`SKBUILD`) builds, CMake defaults to **`CYRSOXS_FETCH_HDF5=ON`**. If no system HDF5 is found, HDF5 **1.14.6** is downloaded and built once as **static** libraries under the build tree and linked into the extension with **`--whole-archive`** (Linux) or **`-force_load`** (macOS) so every HDF5 object file is pulled into `CyRSoXS*.so` (avoids undefined symbols such as `H5Rcreate_object` from partial `.a` linking). The first build can take several minutes.

If you use **system shared HDF5**, CMake sets **`BUILD_RPATH` / `INSTALL_RPATH`** to the library directory CMake discovered so the loader does not pick an older `libhdf5.so` from another package (for example `h5py`) instead of the version you built against.

For **plain CMake** builds without `SKBUILD`, fetching defaults to **OFF**; install HDF5 C++/HL dev packages, or pass `-DCYRSOXS_FETCH_HDF5=ON`, or point CMake at an existing install:

```bash
export HDF5_ROOT=/path/to/hdf5/prefix
```

Optional: override GPU architectures (semicolon-separated list, or `native` on CMake 3.24+):

```bash
export CYRSOXS_CUDA_ARCHITECTURES="80;86"
```

Install the project into a uv-managed environment (CUDA required; HDF5 is fetched automatically when missing):

```bash
uv sync --all-groups
```

Equivalent extras-based install:

```bash
uv sync --all-extras
```

Build wheels and sdists without installing into the active environment:

```bash
uv build -v
```

If you only need Python dev tools (for example `pytest`) and do not have HDF5 headers installed yet:

```bash
uv sync --all-extras --no-install-project
```

To run the smoke test without building the extension (the import is skipped when the module is absent):

```bash
uv run --no-project pytest tests/
```

## Building from source (CMake only)

To build CyRSoXS from source without uv, the following dependencies need to be installed:

- A C++ compiler with C++17 support
- gcc compatible with your CUDA toolkit
- CUDA Toolkit with nvcc
- HDF5 (C++, HL)
- OpenMP

### Additional dependencies for building the Python extension (Pybind11)

- Python >= 3.10
- pybind11 (supplied automatically when building through uv / pip; for a pure CMake build, install pybind11 and pass `-DUSE_SUBMODULE_PYBIND=OFF`, or place a pybind11 source tree under `external/pybind11` and pass `-DUSE_SUBMODULE_PYBIND=ON`)

### Additional dependencies for building without Pybind (CLI executable)

- libconfig++

### Dependencies for building the documentation

- Doxygen
- Graphviz
- latex

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

On Debian or Ubuntu, prefer the packaged development files:

```bash
sudo apt-get install libhdf5-dev
```

Alternatively, build HDF5 from upstream sources (older workflow):

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
This step will create a cmake files in location `$HDF5_DIR/build/_CPack_Packages/Linux/TGZ/HDF5-1.10.5-Linux/HDF_Group/HDF5/1.10.5/share/cmake/hdf5`

Export the path for HDF5:

```bash
cd build/_CPack_Packages/Linux/TGZ/HDF5-1.10.5-Linux/HDF_Group/HDF5/1.10.5/share/cmake/hdf5;
echo "export HDF5_DIR=`pwd`" >> ~/.bashrc
source ~/.bashrc
```

## Downloading CyRSoXS

Clone the CyRSoXS Github repo:

```bash
git clone https://github.com/usnistgov/cyrsoxs.git
```

## Building CyRSoXS

One should have a valid C/C++ compiler and CUDA-toolkit in addition to above modules installed, in order to
compile CyRSoXS.

To compile CyRSoXS as an executable that can be called from the command line or from a bash script, run the following commands:

```bash
cd $CyRSoXS_DIR
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
```

This will generate a `CyRSoXS` executable.

### Building with Pybind (CMake)

CyRSoXS can also be compiled with Pybind11 so that it can be imported as a Python library. Example:

```bash
cd $CyRSoXS_DIR
mkdir build-pybind
cd build-pybind
cmake .. -DCMAKE_BUILD_TYPE=Release -DPYBIND=ON -DUSE_SUBMODULE_PYBIND=OFF
cmake --build .
```

This produces a shared library that can be imported as `CyRSoXS` when placed on `PYTHONPATH`, or use `uv build` / `pip install .` for a proper install into the active environment.

#### Optional CMake Flags

```console
    -DPYBIND=ON             # Compiling with Pybind11
    -DUSE_SUBMODULE_PYBIND=ON   # Use vendored pybind11 under external/pybind11 (optional)
    -DUSE_SUBMODULE_PYBIND=OFF  # Use pybind11 from CMAKE_PREFIX_PATH (default for uv builds)
    -DMAX_NUM_MATERIAL=64   # To change the maximum number of materials (default is 32)
    -DDOUBLE_PRECISION=Yes  # Calculations will performed with double precision numbers
    -DPROFILING=Yes         # Enables profiling of the code
    -DBUILD_DOCS=Yes        # To build documentation
    -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc # Compiling with the Intel compiler (does not work with Pybind)
    -DOUTPUT_BASE_NAME=CyRSoXS # Changes the name of the built output binary or Python module (if using Pybind)
    -DCYRSOXS_CUDA_ARCHITECTURES=80;86 # Override default GPU architecture list
```

## Making CyRSoXS

Once the CMake files has been generated run the following command:

```bash
cmake --build .
```

If `-DBUILD_DOCS=Yes`, the build also drives Doxygen and can produce documentation under the build tree, including `$CyRSoXS_DIR/build/latex/CyRSoXS_Manual.pdf` when LaTeX is available.
