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

#include <cudaMain.h>
#include <Input/readH5.h>
#include <cstdlib>
#include <Input/InputData.h>
#include <Output/writeH5.h>
#include <cstring>
#include <omp.h>
#include <iomanip>
#include <utils.h>

/**
 * main function
 * @param argc
 * @param argv
 * @return EXIT_SUCCESS on successful completion.
 */
int main(int argc, char **argv) {

    if(argc < 2){
        std::cout << "Usage : " << argv[0] << " "<< "HDF5 fileName";
        exit(EXIT_FAILURE);
    }
    std::vector<Material<NUM_MATERIAL> > materialInput;
    InputData inputData(materialInput);
    const UINT voxelSize[3]{inputData.numX, inputData.numY, inputData.numZ};

    std::string fname = argv[1];
    Voxel<NUM_MATERIAL> *voxelData;
    H5::readFile(fname, voxelSize, voxelData);
    Real *projectionGPUAveraged;
    const UINT
      numEnergyLevel = static_cast<UINT>(std::round((inputData.energyEnd - inputData.energyStart) / inputData.incrementEnergy + 1));
    projectionGPUAveraged = new Real[numEnergyLevel * inputData.numX * inputData.numY];
    printCopyrightInfo();
    cudaMain(voxelSize, inputData, materialInput, projectionGPUAveraged, voxelData);
    if(inputData.writeHDF5) {
      writeH5(inputData, voxelSize, projectionGPUAveraged);
    }
    if(inputData.writeVTI) {
      writeVTI(inputData, voxelSize, projectionGPUAveraged);
    }
    printMetaData(inputData);
    delete[] projectionGPUAveraged;
    delete[] voxelData;
    std::cout << "Complete. Exiting \n";

    return EXIT_SUCCESS;

}

