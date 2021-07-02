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
/**
 *
 * @param [in] file HDF5 file pointer
 * @param [in] numMaterial  number of material
 * @param [in] voxelSize voxel size in 3D
 * @param [out] inputData inputData for Mat_alligned
 */

  static inline void getMatAllignment(const H5::H5File &file,
                                      const UINT *voxelSize,
                                      std::vector<std::vector<Real> > &inputData) {

    std::string groupName = "vector_morphology/";
    BigUINT numVoxel = static_cast<BigUINT>((BigUINT) voxelSize[0] * (BigUINT) voxelSize[1] * (BigUINT) voxelSize[2]);
    inputData.resize(NUM_MATERIAL);
    for (UINT i = 0; i < NUM_MATERIAL; i++) {
      inputData[i].resize(numVoxel * 3);
    }
    for (int i = 1; i <= NUM_MATERIAL; i++) {
      std::string varname = groupName + "Mat_" + std::to_string(i) + "_alignment";
      H5::DataSet dataSet;
      try {
        dataSet = file.openDataSet(varname);
      }
      catch (GroupIException not_found_error) {
        std::cout << "Data set = " << varname << "not found";
        throw std::logic_error(" Dataset is not found.");
      }
      H5::DataType dataType = dataSet.getDataType();
      H5::DataSpace space = dataSet.getSpace();
      hsize_t voxelDims[4];
      const int ndims = space.getSimpleExtentDims( voxelDims, NULL);
      // Note HDF5 wrotes dimensions as (Z,Y,X)
      if((ndims != 4) or (voxelDims[0]!=voxelSize[2]) or (voxelDims[1] != voxelSize[1]) or
         (voxelDims[2] != voxelSize[0]) or (voxelDims[3] != 3)) {
        std::cout << "Error in morphology for Material = " << i << " for vector components\n";
        std::cout << "Ndims = " << ndims << "\n";
        std::cout << "Expected dimension from config (X,Y,Z) = " << voxelSize[0] << " " << voxelSize[1] << " " << voxelSize[2] << "\n";
        std::cout << "Dimensions from HDF5 (X,Y,Z)          = "  << voxelDims[2] << " " << voxelDims[1] << " " << voxelDims[0] << "\n";
        throw std::logic_error("Dimension mismatch for morphology");
      }

#ifdef DOUBLE_PRECISION
      if(dataType != PredType::NATIVE_DOUBLE){
         std::cout << "The data format is not supported for double precision \n";
         exit(EXIT_FAILURE);
      }
      dataSet.read(inputData[i - 1].data(), H5::PredType::NATIVE_DOUBLE);

#else
      if (dataType == PredType::NATIVE_DOUBLE) {

        std::vector<double> alignedData(numVoxel * 3);
        dataSet.read(alignedData.data(), H5::PredType::NATIVE_DOUBLE);
        for (int id = 0; id < numVoxel; id++) {
          inputData[i - 1][id * 3 + 0] = static_cast<Real>(alignedData[3 * id + 0]);
          inputData[i - 1][id * 3 + 1] = static_cast<Real>(alignedData[3 * id + 1]);
          inputData[i - 1][id * 3 + 2] = static_cast<Real>(alignedData[3 * id + 2]);
        }
      } else if (dataType == PredType::NATIVE_FLOAT) {
        dataSet.read(inputData[i - 1].data(), H5::PredType::NATIVE_FLOAT);
      }
      dataType.close();
      dataSet.close();
#endif
    }
  }

/**
 *
 * @param [in] file HDF5 file pointer
 * @param [in] numMaterial  number of material
 * @param [in] voxelSize voxel size in 3D
 * @param [out] inputData inputData for Mat_unaligned
 */

  static inline void getScalar(const H5::H5File &file,
                               const std::string & groupName,
                                      const std::string & strName,
                                      const UINT *voxelSize,
                                      std::vector<std::vector<Real> > &inputData) {

    BigUINT numVoxel = static_cast<BigUINT>((BigUINT) voxelSize[0] * (BigUINT) voxelSize[1] * (BigUINT) voxelSize[2]);

    inputData.resize(NUM_MATERIAL);
    for (UINT i = 0; i < NUM_MATERIAL; i++) {
      inputData[i].resize(numVoxel);
    }
    for (int i = 1; i <= NUM_MATERIAL; i++) {
      std::string varname = groupName + "/Mat_" + std::to_string(i) + strName;

      H5::DataSet dataSet;
      try {
        dataSet = file.openDataSet(varname);
      }
      catch (GroupIException not_found_error) {
        std::cout << "Data set = " << varname << "not found";
        throw std::logic_error(" Dataset is not found.");
      }
      H5::DataType dataType = dataSet.getDataType();
      H5::DataSpace space = dataSet.getSpace();
      hsize_t voxelDims[3];
      const int ndims = space.getSimpleExtentDims( voxelDims, NULL);
      // Note HDF5 wrotes dimensions as (Z,Y,X)
      if((ndims != 3) or (voxelDims[0]!=voxelSize[2]) or (voxelDims[1] != voxelSize[1]) or
         (voxelDims[2] != voxelSize[0])) {
        std::cout << "Error in morphology for Material = " << i << "\n";
        std::cout << "Expected dimension from config (X,Y,Z) = " << voxelSize[0] << " " << voxelSize[1] << " " << voxelSize[2] << "\n";
        std::cout << "Dimensions from HDF5 (X,Y,Z)          = "  << voxelDims[2] << " " << voxelDims[1] << " " << voxelDims[0] << "\n";
        throw std::logic_error("Dimension mismatch for morphology");
      }
#ifdef DOUBLE_PRECISION
      if(dataType != PredType::NATIVE_DOUBLE){
         std::cout << "The data format is not supported for double precision \n";
         exit(EXIT_FAILURE);
      }
      dataSet.read(inputData[i - 1].data(), H5::PredType::NATIVE_DOUBLE);
#else
      if (dataType == PredType::NATIVE_DOUBLE) {
        std::vector<double> unalignedData(numVoxel);
        dataSet.read(unalignedData.data(), H5::PredType::NATIVE_DOUBLE);
        for (int id = 0; id < unalignedData.size(); id++) {
          inputData[i - 1][id] = static_cast<Real>(unalignedData[id]);
        }
      } else if (dataType == PredType::NATIVE_FLOAT) {
        dataSet.read(inputData[i - 1].data(), H5::PredType::NATIVE_FLOAT);
      } else {
        std::cout << "This data format is not supported \n";
        exit(EXIT_FAILURE);
      }
      dataType.close();
      dataSet.close();
#endif
    }

  }

/**
 * @brief reads the hdf5 file
 * @param hdf5file hd5 file to read
 * @param voxelSize voxelSize in 3D
 * @param voxelData voxelData
 * @param isAllocated true if the size of voxelData is allocated
 */

  static void readFile(const std::string &hdf5file, const UINT *voxelSize, Voxel<NUM_MATERIAL> *&voxelData,
                      const MorphologyType & morphologyType, bool isAllocated = false) {
    H5::H5File file(hdf5file, H5F_ACC_RDONLY);
    BigUINT numVoxel = static_cast<BigUINT>((BigUINT) voxelSize[0] * (BigUINT) voxelSize[1] * (BigUINT) voxelSize[2]);
    if (not isAllocated) {
      voxelData = new Voxel<NUM_MATERIAL>[numVoxel];
    }
    if (morphologyType == MorphologyType::VECTOR_MORPHOLOGY) {
      {
        std::vector<std::vector<Real> > alignmentData;
        getMatAllignment(file, voxelSize, alignmentData);
        for (UINT numMat = 0; numMat < NUM_MATERIAL; numMat++) {
          for (int i = 0; i < numVoxel; i++) {
            voxelData[i].s1[numMat].x = alignmentData[numMat][3 * i + 0];
            voxelData[i].s1[numMat].y = alignmentData[numMat][3 * i + 1];
            voxelData[i].s1[numMat].z = alignmentData[numMat][3 * i + 2];
          }
        }
      }
      {
        std::vector<std::vector<Real> > unalignedData;
        getScalar(file,"vector_morphology","_unaligned", voxelSize, unalignedData);
        for (UINT numMat = 0; numMat < NUM_MATERIAL; numMat++) {
          for (UINT i = 0; i < numVoxel; i++) {
            voxelData[i].s1[numMat].w = unalignedData[numMat][i];
          }
        }
      }
    }
    else{
      std::vector<std::vector<Real> > s, phi, theta, vfrac;
      getScalar(file,"morphology","_S", voxelSize, s);
      getScalar(file,"morphology","_Phi", voxelSize, phi);
      getScalar(file,"morphology","_Theta", voxelSize, theta);
      getScalar(file,"morphology","_vfrac", voxelSize, vfrac);
      if(morphologyType == MorphologyType::SPHERICAL_COORDINATES){
      for (UINT matID = 0; matID < NUM_MATERIAL; matID++) {
        for (UINT i = 0; i < numVoxel; i++) {
          voxelData[i].s1[matID].x = s[matID][i] * cos(theta[matID][i]);
          voxelData[i].s1[matID].y = s[matID][i] * sin(theta[matID][i]) * sin(phi[matID][i]);
          voxelData[i].s1[matID].z = s[matID][i] * sin(theta[matID][i]) * cos(phi[matID][i]);
          voxelData[i].s1[matID].w = vfrac[matID][i] - s[matID][i];
        }
      }
      }
      else if(morphologyType == MorphologyType::EULER_ANGLES){
        for (UINT matID = 0; matID < NUM_MATERIAL; matID++) {
          for (UINT i = 0; i < numVoxel; i++) {
            voxelData[i].s1[matID].x = s[matID][i] * cos(phi[matID][i]);
            voxelData[i].s1[matID].y = s[matID][i] * sin(phi[matID][i]) * cos(theta[matID][i]);
            voxelData[i].s1[matID].z = s[matID][i] * sin(phi[matID][i]) * sin(theta[matID][i]);
            voxelData[i].s1[matID].w = vfrac[matID][i] - s[matID][i];
          }
        }
      }
    }
  }
}
#endif //PRS_READH5_H
