/////////////////////////////////////////////////////////////////////////////////
// MIT License
//
//Copyright (c) 2019 - 2022 Iowa State University
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

/*! \mainpage Welcome to CyRSoXS
 *
 * \section Introduction
 * CyRSoXS is a GPU-accelerated codebase that calculates resonant X-ray scattering 
 * in the Born Approximation. It takes a voxel-based model and optical constants as 
 * input, and returns simulated X-ray scattering patterns. These models can be derived 
 * from real space measurements (AFM, TEM, 4DSTEM), or synthetically generated using 
 * procedural morphology generators. Currently, the code can run on multiple GPUs on a 
 * single node and it supports both single- and double-precision floating point operations.
 * A speedup of up to two orders of magnitude is observed as compared to the previous 
 * state-of-the-art Igor-based <a href="https://onlinelibrary.wiley.com/iucr/doi/10.1107/S1600577515019074">RSoXS simulator.</a>
 *
 * This core C++/CUDA simulation engine is used in the NIST RSoXS Simulation Suite (NRSS), which 
 * includes additional Python code for creating, simulating, and analyzing RSoXS. For more
 * information on NRSS, please see the <a href="https://nrss.readthedocs.io">documentation.</a>
 * For more information on what RSoXS is and how you can possible apply it in your own 
 * research, check out the <a href="https://www.nist.gov/programs-projects/resonant-soft-x-ray-scattering-rsoxs">NIST RSOXS project page.</a>
 *
 * \section Dependencies
 * The code has following dependencies:
 * <ul>
 *  <li> A GCC / Intel compiler with C++ 14 standard
 *  <li> NVIDIA GPU
 *  <li> cuda-toolkit 9  and above
 *  <li> HDF5 library
 *  <li> Libconfig
 *  <li> OpenMP
 *  <li> Doxygen (optional)
 *  </ul>
 *
 *
 * \section Limitations
 * The current version of the code assumes that the complete data fits in
 * GPU memory.
 *
 * \section Contributors
 * This software was developed at Iowa State University in collaboration with NIST.
 * The Iowa State team provided expertise in high performance computing, and the
 * NIST team provided scientific expertise on the RSoXS technique.
 *
 * <b>Iowa State University</b>
 * <ul>
 * <li> Kumar Saurabh
 * <li> Adarsh Krishnamurthy
 * <li> Baskar Ganapathysubramanian
 * </ul>
 * <b>NIST</b>
 * <ul>
 * <li> Eliot Gann
 * <li> Dean M. Delongchamp
 * <li> Peter J. Dudenas
 * <li> Tyler B. Martin
 * <li> Peter Beaucage
 * </ul>
 *
 * \section Acknowledgement
 * We thank ONR MURI Center for Self-Assembled Organic Electronics for providing 
 * support for this work.
 *
 * \copyright 
 * Copyright 2019-2022 Iowa State University. All rights reserved.
 * This project is released under the MIT and NIST Licenses.
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

  const std::size_t totalArraySize = static_cast<std::size_t>(numEnergyLevel) * static_cast<std::size_t>(inputData.voxelDims[0] * inputData.voxelDims[1]) *
                                     inputData.kVectors.size();
  projectionGPUAveraged = new Real[totalArraySize];

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

