/////////////////////////////////////////////////////////////////////////////////
// MIT License
//
//Copyright (c) 2019 - 2020 Iowa State University
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.
//////////////////////////////////////////////////////////////////////////////////

/*! \mainpage Welcome to Cy-RSoXS.
 *
 * \section Introduction
 * This is the GPU accelerated version of current RSoXS. Currently,
 * the code is tested on NVIDIA V100 GPU. This code can run on multiple GPU on
 * a single node and operates on single - precision floating point (FP32) operation.
 * In the current version of the code, speedup of up to two orders of magnitude is observed
 * as compared to the community accepted Igor based <a href="https://onlinelibrary.wiley.com/iucr/doi/10.1107/S1600577515019074">RSoXS</a> simulator.
 *
 *
 * \section Dependencies
 * The code has following dependencies:
 * <ul>
 *  <li> A GCC / intel compiler with c++ 14 standard
 *  <li> NVIDIA Graphics processing unit
 *  <li> cuda-toolkit 9  and above
 *  <li> HDF5 library
 *  <li> Libconfig
 *  <li> OpenMP
 *  <li> Doxygen (optional)
 *  </ul>
 *
 *
 * \section Limitation
 * The current version of the code assumes that the complete data fits in
 * GPU memory.
 *
 * \section Contributors
 * <ul>
 * <li> Kumar Saurabh
 * <li> Adarsh Krishnamurthy
 * <li> Baskar Ganapathysubramanian
 * <li> Eliot Gann
 * <li> Dean Delongchamp
 * <li> Michael Chabinyc
 * </ul>
 *
 * \section Acknowledgement
 * We thank ONR MURI Center for Self-Assembled Organic Electronics for providing the
 * support for this work.
 *
 * \copyright Copyright 2019-2020 Iowa State University. All rights reserved.
 * This project is released under the MIT license.
 *
 */
//
//#include <cudaMain.h>
#include "H5Cpp.h"
//#include <cstdlib>
#include <Input/InputData.h>
//#include <Output/writeH5.h>
//#include <cstring>
//#include <omp.h>
//#include <iomanip>
//#include <utils.h>
#include "hdf5_hl.h"
//#include "H5File.h"
/**
 * main function
 * @param argc
 * @param argv
 * @return EXIT_SUCCESS on successful completion.
 */
int main(int argc, char **argv) {
  printf("\n*** Checking HDF5 dimension scales.\n");
#define GRP_NAME "simple_scales"
#define DIMSCALE_NAME "dimscale"
#define NAME_ATTRIBUTE "Billy-Bob"
#define VAR1_NAME "var1"
#define VAR2_NAME "var2"
#define VAR3_NAME "var3"
#define DIM1_LEN 3
#define DIM2_LEN 2
#define FIFTIES_SONG "Mamma said they'll be days like this. They'll be days like this, my mamma said."
#define FILE_NAME "try.h5"
#define ERR return 0
  printf("*** Creating simple dimension scales file...");
  {
    hid_t fileid, grpid, dimscaleid;
    hid_t dimscale_spaceid, var1_spaceid, var3_spaceid;
    hid_t var1_datasetid, var2_datasetid, var3_datasetid;
    hsize_t dims[2] = {DIM1_LEN, DIM2_LEN};
    hsize_t dimscale_dims[1] = {DIM1_LEN};

    /* Open file and create group. */
    if ((fileid = H5Fcreate(FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT,
                            H5P_DEFAULT)) < 0)
      ERR;
    if ((grpid = H5Gcreate1(fileid, GRP_NAME, 0)) < 0) ERR;

    /* Create our dimension scale. Use the built-in NAME attribute
     * on the dimscale. */
    if ((dimscale_spaceid = H5Screate_simple(1, dimscale_dims,
                                             dimscale_dims)) < 0)
      ERR;
    if ((dimscaleid = H5Dcreate1(grpid, DIMSCALE_NAME, H5T_NATIVE_INT,
                                 dimscale_spaceid, H5P_DEFAULT)) < 0)
      ERR;
    if (H5DSset_scale(dimscaleid, NAME_ATTRIBUTE) < 0) ERR;

    /* Create a 1D variable which uses the dimscale. Attach a label
     * to this scale. */
    if ((var1_spaceid = H5Screate_simple(1, dims, dims)) < 0) ERR;
    if ((var1_datasetid = H5Dcreate1(grpid, VAR1_NAME, H5T_NATIVE_INT,
                                     var1_spaceid, H5P_DEFAULT)) < 0)
      ERR;
    if (H5DSattach_scale(var1_datasetid, dimscaleid, 0) < 0) ERR;
    if (H5DSset_label(var1_datasetid, 0, FIFTIES_SONG) < 0) ERR;

    /* Create a 1D variabls that doesn't use the dimension scale. */
    if ((var2_datasetid = H5Dcreate1(grpid, VAR2_NAME, H5T_NATIVE_INT,
                                     var1_spaceid, H5P_DEFAULT)) < 0)
      ERR;

    /* Create a 2D dataset which uses the scale for one of its
     * dimensions. */
    if ((var3_spaceid = H5Screate_simple(2, dims, dims)) < 0) ERR;
    if ((var3_datasetid = H5Dcreate1(grpid, VAR3_NAME, H5T_NATIVE_INT,
                                     var3_spaceid, H5P_DEFAULT)) < 0)
      ERR;
    if (H5DSattach_scale(var3_datasetid, dimscaleid, 0) < 0) ERR;

    /* Close up the shop. */
    if (H5Dclose(dimscaleid) < 0 ||
        H5Dclose(var1_datasetid) < 0 ||
        H5Dclose(var2_datasetid) < 0 ||
        H5Dclose(var3_datasetid) < 0 ||
        H5Sclose(var1_spaceid) < 0 ||
        H5Sclose(var3_spaceid) < 0 ||
        H5Sclose(dimscale_spaceid) < 0 ||
        H5Gclose(grpid) < 0 ||
        H5Fclose(fileid) < 0)
      ERR;

    /* HELP! If you are reading this in the future, and time
     * machines have been invented, please come back to July 10,
     * 2005, the Java Java coffee shop in Lafayette, 8:00 am MST +-
     * 20 minutes. Bring back some advanced weapons systems to
     * destroy the sound system here, which is playing 50's rock and
     * roll. Do-op, do-op, la-ma la-ma, ding dong. Save me!!! (Mind
     * you, James Brown is a different story!) */
  }
  printf("*** Getting simple dimension scales file...");
  {
#define NDIMS 1
#define STR_LEN 2



    H5::H5File  file = H5::H5File("foo.h5", H5F_ACC_RDONLY);

    /* Reopen the file and group. */
//    if ((fileid = H5::H5File(FILE_NAME,H5F_ACC_RDONLY))){};
//    H5::Group group = file.openGroup(GRP_NAME);
//    char obj_name[] = ";";
    const char varName[] = "data";
    H5::DataSet dataSet = file.openDataSet(varName);
    /*printf("\nobj_name %s\n", obj_name);*/


    char label[STR_LEN+1];
    H5DSget_label(dataSet.getId(),2,label,STR_LEN);
    std::cout <<  " " << label <<  "\n";
    file.close();
  }
}

