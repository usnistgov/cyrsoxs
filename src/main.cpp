/////////////////////////////////////////////////////////////////////////////////
// MIT License
//
//Copyright (c) 2019 - 2021 Iowa State University
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
//#include <RotationMatrix.h>

/**
 * main function
 * @param argc
 * @param argv
 * @return EXIT_SUCCESS on successful completion.
 */
int main(int argc, char **argv) {


  if (argc < 2) {
    std::cout << "Usage : " << argv[0] << " " << "HDF5FileName" << " HDF5OutputDirname [optional]";
    exit(EXIT_FAILURE);
  }
  std::string fname = argv[1];

  std::vector<Material> materialInput;
  InputData inputData;
  inputData.NUM_MATERIAL = H5::getNumberOfMaterial(fname);
  inputData.readRefractiveIndexData(materialInput);
  inputData.validate();
  const int & NUM_MATERIAL = inputData.NUM_MATERIAL;
  if (argc > 2) {
    inputData.HDF5DirName = argv[2];
  }
  H5::getDimensionAndOrder(fname, (MorphologyType) inputData.morphologyType, inputData.voxelDims,inputData.physSize,
                           inputData.morphologyOrder);
  inputData.check2D();
  inputData.print();
  if(inputData.caseType != CaseTypes::DEFAULT){
    std::cout << BLU << "This is an experimental feature which is not tested. " << NRM <<"\n";
  }
  RotationMatrix matrix(&inputData);
  Voxel *voxelData;
  BigUINT voxelSize = inputData.voxelDims[0] * inputData.voxelDims[1] * inputData.voxelDims[2];

  mallocCPUPinned(voxelData, voxelSize * NUM_MATERIAL);
  H5::readFile(fname, inputData.voxelDims, voxelData, static_cast<MorphologyType>(inputData.morphologyType),
               inputData.morphologyOrder, NUM_MATERIAL, true);
  if(inputData.dumpMorphology){
    H5::writeXDMF(inputData,voxelData);
  }
  if(not(checkMorphology(voxelData,inputData.voxelDims,NUM_MATERIAL))){
    throw std::runtime_error("Nan detected in the morphology");
  }
  Real *projectionGPUAveraged;
  const UINT
    numEnergyLevel = inputData.energies.size();


  projectionGPUAveraged = new Real[numEnergyLevel * (inputData.voxelDims[0] * inputData.voxelDims[1]) *
                                   inputData.kVectors.size()];

  printCopyrightInfo();
  RotationMatrix rotationMatrix(&inputData);
  if (inputData.algorithmType == Algorithm::MemoryMinizing) {
    cudaMainStreams(inputData.voxelDims, inputData, materialInput, projectionGPUAveraged, rotationMatrix, voxelData);
  } else {
    cudaMain(inputData.voxelDims, inputData, materialInput, projectionGPUAveraged, rotationMatrix, voxelData);
  }
  if (inputData.writeHDF5) {
    writeH5(inputData, inputData.voxelDims, projectionGPUAveraged, inputData.HDF5DirName);
  }
  if (inputData.writeVTI) {
    writeVTI(inputData, inputData.voxelDims, projectionGPUAveraged, inputData.VTIDirName);
  }
  printMetaData(inputData, rotationMatrix);
  delete[] projectionGPUAveraged;
  cudaFreeHost(voxelData);
  std::cout << "Complete. Exiting \n";

  return EXIT_SUCCESS;

}

