GPU enabled RSoXS simulation (Release 1.0 - Beta version)
====================================
## Compiling libconfig

Libconfig's build system is [Autotools](https://www.gnu.org/software/automake/manual/html_node/Autotools-Introduction.html), which means you'll need to run `./configure` and then `make` to build.

This guide recommends passing ```--prefix=`pwd`/install``` to `./configure`, which will cause `make install` to copy the output files to `[your_libconfig_dir]/install` instead of `/usr`. This way your libconfig install lives completely inside your libconfig folder. This is necessary if you are working on a system where you don't have admin privileges (i.e. an HPC cluster).

```bash
# Download and extract
wget http://www.hyperrealm.com/packages/libconfig-1.4.10.tar.gz
tar xvf libconfig-1.4.10.tar.gz
rm libconfig-1.4.10.tar.gz

# Compile and copy output files to libconfig-1.4.10/install
cd libconfig-1.4.10
./configure --prefix=`pwd`/install
make -j4  # compile with 4 threads
make install

# Permanently set $LIBCONFIG_DIR environment variable, which is what TalyFEM uses
# to find your libconfig install. Also set LD_LIBRARY_PATH to include the
# libconfig lib directory to prevent dynamic linking errors with libconfig++.so.
echo "export LIBCONFIG_DIR=`pwd`/install" >> ~/.bashrc
echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$LIBCONFIG_DIR/lib" >> ~/.bashrc
source ~/.bashrc
```

**NOTE:** On some HPC clusters (when using the Intel compiler), the `make` step gives a linker error. This is the libconfig example program failing to link with the Intel runtime. This is okay - the libconfig library itself compiles just fine. Just run `make install` and double check that `install/lib` contains some `*.a` files.

## Installing HDF5 

Cy-RSoXS needs [HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format)-based binary input format.

```bash
cd $HDF5_INSTALL_DIRECTORY
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.5/src/CMake-hdf5-1.10.5.tar.gz
tar - xzvf CMake-hdf5-1.10.5.tar.gz
./build-unix.sh
```
This step might take some time. Do not cancel until all the tests have passed.
This step will create a cmake files in location `$HFD5_DIR/build/_CPack_Packages/Linux/TGZ/HDF5-1.10.5-Linux/HDF_Group/HDF5/1.10.5/share/cmake/hdf5`


Building Cy-RSoXS 
==================

**Downloading Cy-RSoXS**

You can download Cy-RSoXS by cloning into the repository

```bash
git clone https://bitbucket.org/baskargroup/cy-rsoxs.git
```

It will ask for your username and password.

**Cmake options**

One should have a valid C/C++ compiler and CUDA-toolkit in addition to above modules installed to 
compile Cy-RSoXS.
```bash
cd $Cy-RSoXS_DIR
mkdir build; 
cd build;
cmake .. -DCMAKE_BUILD_TYPE=Release -DHDF5_DIR=$HFD5_DIR/build/_CPack_Packages/Linux/TGZ/HDF5-1.10.5-Linux/HDF_Group/HDF5/1.10.5/share/cmake/hdf5 -DDLEVEL2=Yes -DNUM_MATERIAL=4 
```

**Compiling with intel compiler**

If you are compiling with intel compiler:
```bash
 -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc
``` 

**Specifying number of material**

Number of material is specified during the compile time by the flag:
```bash
-DNUM_MATERIAL=4 
```  
The above flag will set the number of material to 4.

**2D computation**

If you are running the case with 1 element in the thickness direction

```bash
-DENABLE_2D=Yes
``` 


Optional Cmake Flags can be added:

**For compiling in double precision mode:**
```bash
-DDOUBLE_PRECISION=Yes
```

**For profiling**
```bash
-DPROFILING=Yes
```

**For ASCII writing of output files (This is slow as compared to BINARY)**
```bash
-DVTI_BINARY=No
```
**For dumping all scattering pattern at all energy level and all angle of rotation**
```bash
-DDLEVEL1=Yes
```
**To dump all intermediate files**
```bash
-DUMP_FILES=Yes
```
**Perform Ewald's projection on CPU (Not recommended)**
```bash
-DEOC=Yes
```
**Build Documentation**
```bash
-DBUILD_DOCS=Yes
```

Making Cy-RSoxs
===============
Once the cmake files has been generated run the following command:
```bash
make 
```

In order to generate the documentation, run
```bash
make doc_doxygen
```
