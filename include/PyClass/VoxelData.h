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

#ifndef CY_RSOXS_VOXELDATA_H
#define CY_RSOXS_VOXELDATA_H

#include <Input/InputData.h>
#include <pybind11/numpy.h>
#include <Output/writeH5.h>
#include <Input/readH5.h>

namespace py = pybind11;
/*
 * Stores the voxel data when passed through the Python.
 *
 * There are two ways of passing the voxel data:
 *
 * 1. Pass directly through HDF5 file.
 * 2. Pass through numpy arrays.
 *
 * Note: The above two options are mutually exclusive. They can not be mixed.
 */
class VoxelData {
private:
    Voxel<NUM_MATERIAL> *voxel = nullptr; /// Voxel data
    const InputData &inputData_;           /// input data
    std::bitset<NUM_MATERIAL> validData_; /// Check that voxel data is correct
    const MorphologyType morphologyType_; /// MorphologyType
public:
    /**
     * @brief constructor
     * @param inputData Input data
     */
    VoxelData(const InputData &inputData)
        : inputData_(inputData),morphologyType_(static_cast<MorphologyType>(inputData.morphologyType)) {
      clear();
      const BigUINT numVoxels = inputData_.numX * inputData_.numY * inputData_.numZ;
      voxel = new Voxel<NUM_MATERIAL>[numVoxels];
      validData_.reset();
    }

    /**
     * @brief add Material Allignment and unalligned data to the input for a given material. \n
     * Once you add the material the bits corresponding to that bit is turned on.
     * @param matAlignementData Alignment data for the material. Must be in the order of Sx,Sy,Sz.
     * @param matUnalignedData The unalignment component for the material
     * @param matID material ID
     */
    void addMaterialDataVectorMorphology(py::array_t<Real, py::array::c_style | py::array::forcecast> &matAlignementData,
                         py::array_t<Real, py::array::c_style | py::array::forcecast> &matUnalignedData,
                         const UINT matID) {

      if(morphologyType_ == MorphologyType::EULER_ANGLES){
        py::print("Error: [Expected]: Vector Morphology [Found:] Euler Angles. Returning\n");
        return;
      }
      if (matID >= NUM_MATERIAL) {
        throw std::logic_error("Number of material does not match with the compiled version");
      }
      if(validData_.test(matID)){
        py::print("The material is already set. Please first reset to add the entries. Returning.");
        return;
      }
      const BigUINT numVoxels = inputData_.numX * inputData_.numY * inputData_.numZ;
      for (BigUINT i = 0; i < numVoxels; i++) {
        voxel[i].s1[matID].x = matAlignementData.data()[i * 3 + 0];
        voxel[i].s1[matID].y = matAlignementData.data()[i * 3 + 1];
        voxel[i].s1[matID].z = matAlignementData.data()[i * 3 + 2];
        voxel[i].s1[matID].w = matUnalignedData.data()[i];
      }
      validData_.set(matID, true);
    }

    void addMaterialDataEulerAngles(py::array_t<Real, py::array::c_style | py::array::forcecast> &matSVector,
                         py::array_t<Real, py::array::c_style | py::array::forcecast> &matThetaVector,
                         py::array_t<Real, py::array::c_style | py::array::forcecast> &matPhiVector,
                         py::array_t<Real, py::array::c_style | py::array::forcecast> &matVfracVector,
                         const UINT matID) {
      pybind11::print("Not supported");
      return;
        if(morphologyType_ != MorphologyType::EULER_ANGLES){
          py::print("Error: [Expected]: Euler Angles / Spherical Coordinates. [Found:] VectorMorphology\n");
          return;
        }
        if (matID >= NUM_MATERIAL) {
            throw std::logic_error("Number of material does not match with the compiled version");
        }
        if(validData_.test(matID)){
            py::print("The material is already set. Please first reset to add the entries. Returning.");
            return;
        }
        const BigUINT numVoxels = inputData_.numX * inputData_.numY * inputData_.numZ;
        if(morphologyType_ == MorphologyType::EULER_ANGLES) {
          for (BigUINT i = 0; i < numVoxels; i++) {
            voxel[i].s1[matID].x = matSVector.data()[i] * cos(matPhiVector.data()[i]);
            voxel[i].s1[matID].y = matSVector.data()[i] * sin(matPhiVector.data()[i]) * cos(matThetaVector.data()[i]);
            voxel[i].s1[matID].z = matSVector.data()[i] * sin(matPhiVector.data()[i]) * sin(matThetaVector.data()[i]);
            voxel[i].s1[matID].w = matVfracVector.data()[i] - matSVector.data()[i];
          }
        }
        validData_.set(matID, true);
    }

    /**
     * @brief read the material information from the filename
     * @param fname HDF5 filename
     */
    void readFromH5(const std::string &fname) {
      if(validData_.any()){
        py::print("Some of the material is already set. Please first reset to continue. Returning.");
        return;
      }
      const UINT voxelDim[3]{inputData_.numX, inputData_.numY, inputData_.numZ};
      H5::readFile(fname, voxelDim, voxel,morphologyType_, true);
      validData_.flip();
    }

    /**
     * @brief resets the bit. Once this is called you need to add all the materials.
     */
    void reset(){
      validData_.reset();
    }

    /**
     * @brief Write the voxel data to HDF5.
     */
    void writeToH5() const {
      if(not(this->validate())){
        return;
      }
      const BigUINT numVoxels = inputData_.numX * inputData_.numY * inputData_.numZ;
      const UINT dim[3]{inputData_.numX, inputData_.numY, inputData_.numZ};
      Real *scalarValues = new Real[numVoxels];
      for (int nMat = 0; nMat < NUM_MATERIAL; nMat++) {
        std::string fname = "Unalligned_" + std::to_string(nMat);
        for (int i = 0; i < numVoxels; i++) {
          scalarValues[i] = voxel[i].s1[nMat].w;
        }
        H5::writeFile3DScalar(fname, scalarValues, dim, "unalligned");
      }
      delete[] scalarValues;

      Real *vectorValues = new Real[numVoxels * 3];
      for (int nMat = 0; nMat < NUM_MATERIAL; nMat++) {
        std::string fname = "Alligned_" + std::to_string(nMat);
        for (int i = 0; i < numVoxels; i++) {
          vectorValues[i * 3 + 0] = voxel[i].s1[nMat].x;
          vectorValues[i * 3 + 1] = voxel[i].s1[nMat].y;
          vectorValues[i * 3 + 2] = voxel[i].s1[nMat].z;
        }
        H5::writeFile3DVector(fname, vectorValues, dim, "S");
      }
      delete[] vectorValues;
    }

    /**
     * Clear the voxel data.
     */
    void clear() {
      if (voxel != nullptr) {
        delete[] voxel;
      }
      voxel = nullptr;
    }

    /**
     * @brief Destructor
     */
    ~VoxelData() {
      clear();
    }

    /**
     * @brief returns the voxel data
     * @return The voxel data
     */
    const Voxel<NUM_MATERIAL> *data() const {
      return voxel;
    }

    /**
     * @brief Checks if the voxel data is correct.
     * @return True if the input is correct. False otherwise
     */
    bool validate() const {
      if (validData_.all()){
        return true;
      }
      else{
        for(UINT i = 0; i < validData_.size(); i++){
          if(not(validData_.test(i))){
            py::print("Voxel Data missing / corrupt for material = ",i);
          }
        }
      }
      return false;
    }

};
#endif //CY_RSOXS_VOXELDATA_H
