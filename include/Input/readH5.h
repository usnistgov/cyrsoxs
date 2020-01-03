/////////////////////////////////////////////////////////////////////////////////
// MIT License
//
//Copyright (c) 2019 Iowa State University
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

#ifndef PRS_READH5_H
#define PRS_READH5_H

#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <string>
#include <sstream>      // std::stringstream

#include <vector>
#include <assert.h>
#include "Input.h"
#include "H5Cpp.h"

namespace H5 {

inline const UINT getNumMaterial(const H5::H5File &file) {
  float numMaterial;
  DataSet dataSet = file.openDataSet("igor_parameters/igormaterialnum");
  dataSet.read(&numMaterial, PredType::NATIVE_FLOAT);
  std::cout << "Number of material = " << numMaterial << "\n";
  assert(numMaterial == NUM_MATERIAL);

  return (static_cast<UINT>(numMaterial));
}

/**
 *
 * @param [in] file HDF5 file pointer
 * @param [in] numMaterial  number of material
 * @param [in] voxelSize voxel size in 3D
 * @param [out] inputData inputData for Mat_alligned
 */

inline void getMatAllignment(const H5::H5File &file,
                             const UINT numMaterial,
                             const UINT *voxelSize,
                             std::vector<std::vector<Real> > &inputData) {
  assert(numMaterial == NUM_MATERIAL);
  std::string groupName = "vector_morphology/";
  BigUINT numVoxel = static_cast<BigUINT>(voxelSize[0] * voxelSize[1] * voxelSize[2]);
  inputData.resize(numMaterial);
  for (UINT i = 0; i < numMaterial; i++) {
    inputData[i].resize(numVoxel * 3);
  }
  int rank = 4;
  hsize_t dimsm[4]{static_cast<hsize_t>(voxelSize[0]),
                   static_cast<hsize_t>(voxelSize[1]),
                   static_cast<hsize_t>(voxelSize[2]),
                   3};
  H5::DataSpace memspace(rank, dimsm);

  for (int i = 1; i <= numMaterial; i++) {
    std::string varname = groupName + "Mat_" + std::to_string(i) + "_alignment";
    H5::DataSet dataSet = file.openDataSet(varname);
    H5::DataSpace dataspace = dataSet.getSpace();
    dataSet.read(inputData[i - 1].data(), H5::PredType::NATIVE_FLOAT, memspace, dataspace);
  }
}
/**
 *
 * @param [in] file HDF5 file pointer
 * @param [in] numMaterial  number of material
 * @param [in] voxelSize voxel size in 3D
 * @param [out] inputData inputData for Mat_unaligned
 */

inline void getMatUnAlligned(const H5::H5File &file,
                             const UINT numMaterial,
                             const UINT *voxelSize,
                             std::vector<std::vector<Real> > &inputData) {

  assert(numMaterial == NUM_MATERIAL);
  std::string groupName = "vector_morphology/";
  BigUINT numVoxel = static_cast<BigUINT>(voxelSize[0] * voxelSize[1] * voxelSize[2]);

  inputData.resize(numMaterial);
  for (UINT i = 0; i < numMaterial; i++) {
    inputData[i].resize(numVoxel);
  }
  int rank = 3;
  hsize_t dimsm[3]{static_cast<hsize_t>(voxelSize[0]),
                   static_cast<hsize_t>(voxelSize[1]),
                   static_cast<hsize_t>(voxelSize[2])};
  H5::DataSpace memspace(rank, dimsm);

  for (int i = 1; i <= numMaterial; i++) {
    std::string varname = groupName + "Mat_" + std::to_string(i) + "_unaligned";
    H5::DataSet dataSet = file.openDataSet(varname);
    H5::DataSpace dataspace = dataSet.getSpace();
    dataSet.read(inputData[i - 1].data(), H5::PredType::NATIVE_FLOAT, memspace, dataspace);
  }

}

/**
 * @brief reads the hdf5 file
 * @param hdf5file hd5 file to read
 * @param voxelSize voxelSize in 3D
 * @param voxelData voxelData
 * @param isAllocated true if the size of voxelData is allocated
 */

void readFile(const std::string hdf5file, const UINT *voxelSize, Voxel<NUM_MATERIAL> *& voxelData, bool isAllocated = false) {
  H5::H5File file(hdf5file, H5F_ACC_RDONLY);
  const UINT numMaterial = getNumMaterial(file);

  BigUINT numVoxel = static_cast<BigUINT>(voxelSize[0] * voxelSize[1] * voxelSize[2]);
  if(not isAllocated){
    voxelData = new Voxel<NUM_MATERIAL>[numVoxel];
  }
  {
    std::vector<std::vector<Real> > alignmentData;
    getMatAllignment(file, numMaterial, voxelSize, alignmentData);
    for(UINT numMat = 0; numMat < numMaterial; numMat++) {
      for (int i = 0; i < numVoxel; i++) {
        voxelData[i].s1[numMat].x = alignmentData[numMat][3*i + 0];
        voxelData[i].s1[numMat].y = alignmentData[numMat][3*i + 1];
        voxelData[i].s1[numMat].z = alignmentData[numMat][3*i + 2];
      }
    }
  }

  {
    std::vector<std::vector<Real> > unalignedData;
    getMatUnAlligned(file, numMaterial, voxelSize, unalignedData);
    for(UINT numMat = 0; numMat < numMaterial; numMat++) {
      for (int i = 0; i < numVoxel; i++) {
        voxelData[i].s1[numMat].w = unalignedData[numMat][i];
      }
    }
  }

}

}
#endif //PRS_READH5_H
