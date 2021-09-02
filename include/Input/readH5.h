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
#include <hdf5_hl.h>

namespace H5 {


  static constexpr int AXIS_LABEL_LEN = 2;


  static inline void
  getDimensionAndOrder(const std::string &hdf5fileName, const MorphologyType &morphologyType, UINT *voxelSize, Real & physSize,
                       MorphologyOrder &morphologyOrder) {
    H5::H5File file(hdf5fileName, H5F_ACC_RDONLY);
    { // PhysSize
      std::string groupName = "morphology_parameter";
      std::string dataName = "PhysSize";
      bool groupExists = file.nameExists(groupName.c_str());
      if (not groupExists) {
        std::cerr << "Group " << groupName << "not found";
        exit(EXIT_FAILURE);
      }
      Group group = file.openGroup(groupName.c_str());
      bool dataExists = group.nameExists(dataName.c_str());
      if (not(dataExists)) {
        std::cerr << "DataSet " << dataName << "not found";
        exit(EXIT_FAILURE);
      }
      H5::DataSet dataSet = group.openDataSet(dataName.c_str());
      H5::DataType dataType = dataSet.getDataType();
#ifdef DOUBLE_PRECISION
      if(dataType == PredType::NATIVE_DOUBLE) {
        dataSet.read(&physSize,PredType::NATIVE_FLOAT);
      }
      else {
        throw std::runtime_error("Wrong Data type for physSize");
      }
#else
      if(dataType == PredType::NATIVE_FLOAT) {
        dataSet.read(&physSize,PredType::NATIVE_FLOAT);
      }
      else if(dataType == PredType::NATIVE_DOUBLE) {
        double _physSize;
        dataSet.read(&_physSize,PredType::NATIVE_DOUBLE);
        physSize = _physSize;
      }
      else {
        throw std::runtime_error("Wrong Data type for physSize");
      }
#endif
      dataSet.close();
      group.close();
    }
    if (morphologyType == MorphologyType::VECTOR_MORPHOLOGY) {
      std::string groupName = "vector_morphology";
      std::string dataName = "Mat_1_unaligned";
      bool groupExists = file.nameExists(groupName.c_str());
      if (not groupExists) {
        std::cerr << "Group " << groupName << "not found";
        exit(EXIT_FAILURE);
      }

      Group group = file.openGroup(groupName.c_str());
      bool dataExists = group.nameExists(dataName.c_str());

      if (not(dataExists)) {
        std::cerr << "DataSet " << dataName << "not found";
        exit(EXIT_FAILURE);
      }
      H5::DataSet dataSet = group.openDataSet(dataName.c_str());
      H5::DataSpace space = dataSet.getSpace();
      hsize_t voxelDims[3];
      const int ndims = space.getSimpleExtentDims(voxelDims, NULL);
      if (ndims != 3) {
        std::cerr << "Expected 3D array. Found Dim = " << ndims << "for " << dataName << "\n";
        exit(EXIT_FAILURE);
      }
      char label[2][AXIS_LABEL_LEN];
      H5DSget_label(dataSet.getId(), 0, label[0], AXIS_LABEL_LEN);
      H5DSget_label(dataSet.getId(), 2, label[1], AXIS_LABEL_LEN);
      if (((strcmp(label[0], "Z") == 0) and (strcmp(label[1], "X") == 0))) {
        morphologyOrder = MorphologyOrder::ZYX;
        voxelSize[0] = voxelDims[2];
        voxelSize[1] = voxelDims[1];
        voxelSize[2] = voxelDims[0];
      } else if (((strcmp(label[0], "X") == 0) and (strcmp(label[1], "Z") == 0))) {
        morphologyOrder = MorphologyOrder::XYZ;
        voxelSize[0] = voxelDims[0];
        voxelSize[1] = voxelDims[1];
        voxelSize[2] = voxelDims[2];
      } else {
        throw std::runtime_error("Only XYZ/ZYX ordering supported");
      }
      group.close();
    } else {
      throw std::runtime_error("Not supported");
    }

    file.close();
  }

  /**
   * @brief Converts XYZ format data to ZYX format. Cy-RSoXS assumes the data is always in ZYX format
   * @tparam T
   * @param data data for single material
   * @param numComponents number of components (1/3 for scalar / vector)
   * @param voxelSize 3D array of number of voxels in each direction
   */
  template<typename T>
  static void XYZ_to_ZYX(std::vector<T> &data, const int numComponents, const UINT *voxelSize) {
    std::vector<T> _data = data;
    int counter = 0;
    UINT X = voxelSize[0];
    UINT Y = voxelSize[1];
    UINT Z = voxelSize[2];
    for (int k = 0; k < Z; k++) {
      for (int j = 0; j < Y; j++) {
        for (int i = 0; i < X; i++) {
          UINT flattenZYX = k * (X * Y) + j * X + i;
          UINT flattenXYZ = i * (Y * Z) + j * Z + k;
          for (int c = 0; c < numComponents; c++) {
            _data[flattenZYX * numComponents + c] = data[flattenXYZ * numComponents + c];
          }
          counter++;
        }
      }
    }
    std::swap(data, _data);
  }

/**
 *
 * @param [in] file HDF5 file pointer
 * @param [in] numMaterial  number of material
 * @param [in] voxelSize voxel size in 3D
 * @param [out] inputData inputData for Mat_alligned
 * Note : This is read after unaligned data. So, the morphology order and voxelSize is set according to alligned data.
 * Any mismatch will throw an error.
 */

  static inline void getMatAllignment(const H5::H5File &file,
                                      const UINT *voxelSize,
                                      std::vector<Real> &inputData,
                                      const MorphologyOrder &morphologyOrder,
                                      const int materialID) {
    std::string groupName = "vector_morphology";
    BigUINT numVoxel = static_cast<BigUINT>((BigUINT) voxelSize[0] * (BigUINT) voxelSize[1] * (BigUINT) voxelSize[2]);

    int i = materialID;

    std::string dataName = "Mat_" + std::to_string(i) + "_alignment";

    H5::DataSet dataSet;
    bool groupExists = file.nameExists(groupName.c_str());
    if (not groupExists) {
      std::cerr << "Group " << groupName << "not found";
      exit(EXIT_FAILURE);
    }
    Group group = file.openGroup(groupName.c_str());
    bool dataExists = group.nameExists(dataName.c_str());
    if (not(dataExists)) {
      std::cerr << "Dataset = " << dataName << "does not exists";
      exit(EXIT_FAILURE);
    }

    dataSet = group.openDataSet(dataName.c_str());
    H5::DataType dataType = dataSet.getDataType();
    H5::DataSpace space = dataSet.getSpace();
    hsize_t voxelDims[4];
    const int ndims = space.getSimpleExtentDims(voxelDims, NULL);
    if (ndims != 4) {
      std::cerr << "Expected 4D array. Found Dim = " << ndims << "for " << dataName << "\n";
      exit(EXIT_FAILURE);
    }
    char label[2][AXIS_LABEL_LEN];
    H5DSget_label(dataSet.getId(), 0, label[0], AXIS_LABEL_LEN);
    H5DSget_label(dataSet.getId(), 2, label[1], AXIS_LABEL_LEN);

    // We store voxel dimension as (X,Y,Z) irrespective of HDF axis label
    if (morphologyOrder == MorphologyOrder::ZYX) {
      if (not((strcmp(label[0], "Z") == 0) and (strcmp(label[1], "X") == 0))) {
        std::cerr << "Axis label mismatch for morphology for " << dataName << "\n";
        exit(EXIT_FAILURE);
      }
      if ((voxelDims[0] != voxelSize[2]) or (voxelDims[1] != voxelSize[1]) or
          (voxelDims[2] != voxelSize[0]) or (voxelDims[3] != 3)) {
        std::cout << "Error in " << dataName << "\n";
        std::cout << "Error in morphology for Material = " << i << "\n";
        std::cout << "Expected dimension (X,Y,Z) = " << voxelSize[0] << " " << voxelSize[1] << " " << voxelSize[2]
                  << "\n";
        std::cout << "Dimensions from HDF5 (X,Y,Z)          = " << voxelDims[2] << " " << voxelDims[1] << " "
                  << voxelDims[0] << "\n";
        throw std::logic_error("Dimension mismatch for morphology");
      }
    }
    if (morphologyOrder == MorphologyOrder::XYZ) {
      if (not((strcmp(label[0], "X") == 0) and (strcmp(label[1], "Z") == 0))) {
        std::cerr << "Axis label mismatch for morphology for " << dataName << "\n";
        exit(EXIT_FAILURE);
      }
      if ((voxelDims[0] != voxelSize[0]) or (voxelDims[1] != voxelSize[1]) or
          (voxelDims[2] != voxelSize[2]) or (voxelDims[3] != 3)) {
        std::cout << "Error in " << dataName << "\n";
        std::cout << "Error in morphology for Material = " << i << "\n";
        std::cout << "Expected dimension (X,Y,Z)   = " << voxelSize[0] << " " << voxelSize[1] << " " << voxelSize[2]
                  << "\n";
        std::cout << "Dimensions from HDF5 (X,Y,Z) = " << voxelDims[0] << " " << voxelDims[1] << " " << voxelDims[2]
                  << "\n";
        throw std::logic_error("Dimension mismatch for morphology");
      }
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
      if (morphologyOrder == MorphologyOrder::XYZ) {
        XYZ_to_ZYX(alignedData, 3, voxelSize);
      }
      for (BigUINT id = 0; id < numVoxel; id++) {
        inputData[id * 3 + 0] = static_cast<Real>(alignedData[3 * id + 0]);
        inputData[id * 3 + 1] = static_cast<Real>(alignedData[3 * id + 1]);
        inputData[id * 3 + 2] = static_cast<Real>(alignedData[3 * id + 2]);
      }
    } else if (dataType == PredType::NATIVE_FLOAT) {
      std::vector<float> alignedData(numVoxel * 3);
      dataSet.read(alignedData.data(), H5::PredType::NATIVE_FLOAT);
      if (morphologyOrder == MorphologyOrder::XYZ) {
        XYZ_to_ZYX(alignedData, 3, voxelSize);
      }
      for (BigUINT id = 0; id < numVoxel; id++) {
        inputData[id * 3 + 0] = static_cast<Real>(alignedData[3 * id + 0]);
        inputData[id * 3 + 1] = static_cast<Real>(alignedData[3 * id + 1]);
        inputData[id * 3 + 2] = static_cast<Real>(alignedData[3 * id + 2]);
      }
    }
    group.close();
    dataType.close();
    dataSet.close();
#endif
  }

/**
 *
 * @param [in] file HDF5 file pointer
 * @param [in] numMaterial  number of material
 * @param [in] voxelSize voxel size in 3D
 * @param [out] inputData inputData for Mat_unaligned
 */

  static inline bool getScalar(const H5::H5File &file,
                               const std::string &groupName,
                               const std::string &strName,
                               const UINT *voxelSize,
                               const MorphologyOrder &morphologyOrder,
                               std::vector<Real> &morphologyData,
                               const int materialID,
                               const bool isRequired = true) {

    BigUINT numVoxel = static_cast<BigUINT>((BigUINT) voxelSize[0] * (BigUINT) voxelSize[1] * (BigUINT) voxelSize[2]);


    int i = materialID;
    std::string dataName = "Mat_" + std::to_string(i) + strName;

    bool groupExists = file.nameExists(groupName.c_str());
    if (not groupExists) {
      std::cerr << "Group " << groupName << "not found";
      exit(EXIT_FAILURE);
    }

    Group group = file.openGroup(groupName.c_str());
    bool dataExists = group.nameExists(dataName.c_str());

    // Check if dataset exists and required
    if (isRequired and not(dataExists)) {
      std::cerr << "Dataset = " << dataName << "does not exists";
      exit(EXIT_FAILURE);
    }

    // Fill with 0 if not exists
    if (not(dataExists)) {
      std::fill(morphologyData.begin(), morphologyData.end(), 0.0);
      return dataExists;
    }

    H5::DataSet dataSet = group.openDataSet(dataName.c_str());
    H5::DataType dataType = dataSet.getDataType();
    H5::DataSpace space = dataSet.getSpace();
    hsize_t voxelDims[3];
    const int ndims = space.getSimpleExtentDims(voxelDims, NULL);
    if (ndims != 3) {
      std::cerr << "Expected 3D array. Found Dim = " << ndims << "for " << dataName << "\n";
      exit(EXIT_FAILURE);
    }
    char label[2][AXIS_LABEL_LEN];
    H5DSget_label(dataSet.getId(), 0, label[0], AXIS_LABEL_LEN);
    H5DSget_label(dataSet.getId(), 2, label[1], AXIS_LABEL_LEN);


    // We store voxel dimension as (X,Y,Z) irrespective of HDF axis label
    if (morphologyOrder == MorphologyOrder::ZYX) {
      if (not((strcmp(label[0], "Z") == 0) and (strcmp(label[1], "X") == 0))) {
        std::cerr << "Axis label mismatch for morphology for " << dataName << "\n";
        exit(EXIT_FAILURE);
      }
      if ((voxelDims[0] != voxelSize[2]) or (voxelDims[1] != voxelSize[1]) or
          (voxelDims[2] != voxelSize[0])) {
        std::cout << "Error in " << dataName << "\n";
        std::cout << "Error in morphology for Material = " << i << "\n";
        std::cout << "Expected dimension (X,Y,Z) = " << voxelSize[0] << " " << voxelSize[1] << " " << voxelSize[2]
                  << "\n";
        std::cout << "Dimensions from HDF5 (X,Y,Z) = " << voxelDims[2] << " " << voxelDims[1] << " "
                  << voxelDims[0] << "\n";

        throw std::logic_error("Dimension mismatch for morphology");
      }
    }
    if (morphologyOrder == MorphologyOrder::XYZ) {
      if (not((strcmp(label[0], "X") == 0) and (strcmp(label[1], "Z") == 0))) {
        std::cerr << "Axis label mismatch for morphology for " << dataName << "\n";
        exit(EXIT_FAILURE);
      }
      if ((voxelDims[0] != voxelSize[0]) or (voxelDims[1] != voxelSize[1]) or
          (voxelDims[2] != voxelSize[2])) {
        std::cout << "Error in " << dataName << "\n";
        std::cout << "Error in morphology for Material = " << i << "\n";
        std::cout << "Expected dimension (X,Y,Z)   = " << voxelSize[0] << " " << voxelSize[1] << " " << voxelSize[2]
                  << "\n";
        std::cout << "Dimensions from HDF5 (X,Y,Z) = " << voxelDims[0] << " " << voxelDims[1] << " " << voxelDims[2]
                  << "\n";
        throw std::logic_error("Dimension mismatch for morphology");
      }
    }


#ifdef DOUBLE_PRECISION
    if(dataType != PredType::NATIVE_DOUBLE){
       std::cout << "The data format is not supported for double precision \n";
       exit(EXIT_FAILURE);
    }
    dataSet.read(inputData[i - 1].data(), H5::PredType::NATIVE_DOUBLE);
#else
    if (dataType == PredType::NATIVE_DOUBLE) {
      std::vector<double> scalarData(numVoxel);
      dataSet.read(scalarData.data(), H5::PredType::NATIVE_DOUBLE);
      if (morphologyOrder == MorphologyOrder::XYZ) {
        XYZ_to_ZYX(scalarData, 1, voxelSize);
      }
      for (int id = 0; id < scalarData.size(); id++) {
        morphologyData[id] = static_cast<Real>(scalarData[id]);
      }
    } else if (dataType == PredType::NATIVE_FLOAT) {
      dataSet.read(morphologyData.data(), H5::PredType::NATIVE_FLOAT);
      if (morphologyOrder == MorphologyOrder::XYZ) {
        XYZ_to_ZYX(morphologyData, 1, voxelSize);
      }
    } else {
      std::cerr << "This data format is not supported \n";
      exit(EXIT_FAILURE);
    }
    group.close();
    dataType.close();
    dataSet.close();
#endif
    return dataExists;
  }


/**
 * @brief reads the hdf5 file
 * @param hdf5file hd5 file to read
 * @param voxelSize voxelSize in 3D
 * @param voxelData voxelData
 * @param isAllocated true if the size of voxelData is allocated
 */

  static int readFile(const std::string &hdf5file, const UINT *voxelSize, Voxel *&voxelData,
                      const MorphologyType &morphologyType, const MorphologyOrder &morphologyOrder,
                      bool isAllocated = false) {
    H5::H5File file(hdf5file, H5F_ACC_RDONLY);
    BigUINT numVoxel = static_cast<BigUINT>((BigUINT) voxelSize[0] * (BigUINT) voxelSize[1] * (BigUINT) voxelSize[2]);


    if (not isAllocated) {
      voxelData = new Voxel[numVoxel * NUM_MATERIAL];
    }
    if (morphologyType == MorphologyType::VECTOR_MORPHOLOGY) {
      {
        std::vector<Real> unalignedData(numVoxel);
        for (int numMat = 1; numMat < NUM_MATERIAL + 1; numMat++) {
#ifdef PYBIND
          getScalar(file, "vector_morphology", "_unaligned", voxelSize, morphologyOrder, unalignedData,numMat,true, false);
#else
          getScalar(file, "vector_morphology", "_unaligned", voxelSize, morphologyOrder, unalignedData, numMat, true);
#endif
          for (UINT i = 0; i < numVoxel; i++) {
            voxelData[i * NUM_MATERIAL + numMat - 1].s1.w = unalignedData[i];
          }
        }
      }
      {
        std::vector<Real> alignmentData(numVoxel * 3);
        for (UINT numMat = 1; numMat < NUM_MATERIAL + 1; numMat++) {
          getMatAllignment(file, voxelSize, alignmentData, static_cast<const MorphologyOrder>(morphologyOrder), numMat);
          for (BigUINT i = 0; i < numVoxel; i++) {
            voxelData[i * NUM_MATERIAL + numMat - 1].s1.x = alignmentData[3 * i + 0];
            voxelData[i * NUM_MATERIAL + numMat - 1].s1.y = alignmentData[3 * i + 1];
            voxelData[i * NUM_MATERIAL + numMat - 1].s1.z = alignmentData[3 * i + 2];
          }
        }
      }
    } else if (morphologyType == MorphologyType::EULER_ANGLES) {
      /// TODO
    } else {
      throw std::runtime_error("Wrong type of morphology");
    }
    file.close();
    return EXIT_SUCCESS;
  }
}
#endif //PRS_READH5_H
