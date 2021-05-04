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
#ifndef CY_RSOXS_UTILS_H
#define CY_RSOXS_UTILS_H
#include <Input/InputData.h>
#include <iomanip>
#include <Output/outputUtils.h>
#include <Output/writeVTI.h>
#include "version.h"

/**
 * @brief writes to HDF5 file
 * @param inputData Input data
 * @param voxelSize Voxel Size
 * @param projectionGPUAveraged The scattering pattern
 */
static void writeH5(const InputData & inputData, const UINT * voxelSize,const Real *projectionGPUAveraged, const std::string dirName = "HDF5"){
    createDirectory(dirName);
    // TODO: Make it work in parallel.
    omp_set_num_threads(1);

    const UINT numEnergyLevel = inputData.energies.size();
    const UINT voxel2DSize = voxelSize[0] * voxelSize[1];
    Real *oneEnergyData = new Real[voxel2DSize];
    UINT chunkSize = static_cast<UINT>(std::ceil(numEnergyLevel * 1.0 / (omp_get_num_threads() * 1.0)));
    UINT threadID = omp_get_thread_num();
    UINT startID = threadID * chunkSize;
    UINT endID = ((threadID + 1) * chunkSize);

    for (UINT csize = startID; csize < std::min(endID, numEnergyLevel); csize++) {

      std::stringstream stream;
      Real energy = inputData.energies[csize];
      stream << std::fixed << std::setprecision(2) << energy;
      std::string s = stream.str();
      const std::string outputFname = dirName + "/Energy_" + s + ".h5";
      H5::H5File file(outputFname.c_str(), H5F_ACC_TRUNC);
      for (UINT kID = 0; kID < inputData.kVectors.size(); kID++) {
        const UINT offset = csize * voxel2DSize * inputData.kVectors.size() + kID*voxel2DSize;
        std::memcpy(oneEnergyData, &projectionGPUAveraged[offset], sizeof(Real) * voxel2DSize);
        const std::string groupname = "K" + std::to_string(kID);
        H5::writeFile2D(file, oneEnergyData, voxelSize,groupname);
      }
      file.close();
    }

    delete[] oneEnergyData;
}
/**
 * @brief writes to VTI file in parallel
 * @param inputData Input data
 * @param voxelSize Voxel Size
 * @param projectionGPUAveraged The scattering pattern
 */
static void writeVTI(const InputData & inputData, const UINT * voxelSize,const Real *projectionGPUAveraged,const std::string dirName = "VTI"){
    std::cout << "Not supported\n";
    return;
    createDirectory(dirName);
    omp_set_num_threads(inputData.num_threads);
    const UINT numEnergyLevel = inputData.energies.size();
#pragma omp parallel
    {
      UINT threadID = omp_get_thread_num();
      UINT chunkSize = static_cast<UINT>(std::ceil(numEnergyLevel * 1.0 / (omp_get_num_threads() * 1.0)));
      const UINT voxel2DSize = voxelSize[0] * voxelSize[1];
      Real *projectionLocal = new Real[voxel2DSize];
      UINT startID = threadID * chunkSize;
      UINT endID = ((threadID + 1) * chunkSize);
      if (startID < numEnergyLevel) {

        for (UINT csize = startID; csize < std::min(endID, numEnergyLevel); csize++) {

          std::memcpy(projectionLocal, &projectionGPUAveraged[csize * voxel2DSize], sizeof(Real) * voxel2DSize);
          std::stringstream stream;
          Real energy = inputData.energies[csize];
          stream << std::fixed << std::setprecision(2) << energy;
          std::string s = stream.str();
          const std::string outputFname = dirName + "/Energy_" + s;
          VTI::writeDataScalar2DFP(projectionLocal, voxelSize, outputFname, "projection");

        }
      }
      delete[] projectionLocal;
    }
}

/**
 * Prints the copyRight Info
 */
static void printCopyrightInfo(){

  std::cout << " __________________________________________________________________________________________________\n";
  std::cout << "|                                 Thanks for using Cy-RSoXS                                        |\n";
  std::cout << "|--------------------------------------------------------------------------------------------------|\n";
  std::cout << "|  Copyright          : Iowa State University                                                      |\n";
  std::cout << "|  License            : MIT                                                                        |\n";
  std::cout << "|  Acknowledgement    : ONR MURI                                                                   |\n";
  std::cout << "|  Developed at Iowa State University in collaboration with NIST                                   |\n";
  std::cout << "|  Please cite the following publication :                                                         |\n";
  std::cout << "|  Comments/Questions :                                                                            |\n";
  std::cout << "|          1. Dr. Baskar GanapathySubramanian (baskarg@iastate.edu)                                |\n";
  std::cout << "|          2. Dr. Adrash Krishnamurthy        (adarsh@iastate.edu)                                 |\n";
  std::cout << " -------------------------------------------------------------------------------------------------- \n";

  std::cout << "Version   : " << VERSION_MAJOR << "."<< VERSION_MINOR << "."<< VERSION_PATCH << "\n";
  std::cout << "Git patch : " << GIT_HASH << "\n";
  std::cout << "\n";
}

/**
 * Prints the copyRight Info to file
 */
static void printCopyrightInfo(std::ofstream & fout){

  fout << " __________________________________________________________________________________________________\n";
  fout << "|                                 Thanks for using Cy-RSoXS                                        |\n";
  fout << "|--------------------------------------------------------------------------------------------------|\n";
  fout << "|  Copyright          : Iowa State University                                                      |\n";
  fout << "|  License            : MIT                                                                        |\n";
  fout << "|  Acknowledgement    : ONR MURI                                                                   |\n";
  fout << "|  Developed at Iowa State University in collaboration with NIST                                   |\n";
  fout << "|  Please cite the following publication :                                                         |\n";
  fout << "|  Comments/Questions :                                                                            |\n";
  fout << "|          1. Dr. Baskar GanapathySubramanian (baskarg@iastate.edu)                                |\n";
  fout << "|          2. Dr. Adrash Krishnamurthy        (adarsh@iastate.edu)                                 |\n";
  fout << " -------------------------------------------------------------------------------------------------- \n";

  fout << "\n";
}
/**
 * @brief Writes meta data
 * @param inputData Input data
 */

static void printMetaData(const InputData & inputData, const RotationMatrix & rotationMatrix){
  std::ofstream file("metadata.txt");
  printCopyrightInfo(file);
  file << "\n\nCyRSoXS: \n";
  file << "=========================================================================================\n";
  file << "Version              : " << VERSION_MAJOR << "."<< VERSION_MINOR << "."<< VERSION_PATCH << "\n";
  file << "Git patch            : " << GIT_HASH << "\n";
  file << "Size of Real         : " << sizeof(Real) << "\n";
  file << "Number of material   : " << NUM_MATERIAL << "\n";

  file << "\nScaling Information:\n";
  file << "=========================================================================================\n";
  file << "Number of pixel      :[" << inputData.numX << "," << inputData.numY << "]\n";
  file << "Q range              :[" << -M_PI / inputData.physSize << "," << M_PI / inputData.physSize << "]\n";
  file << "\n\n";

  file << "\nInputData : \n";
  file << "=========================================================================================\n";
  inputData.printToFile(file);
  rotationMatrix.printToFile(file);

  file.close();
}

#endif //CY_RSOXS_UTILS_H
