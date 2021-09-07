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
    Voxel *voxel = nullptr; /// Voxel data
    const InputData &inputData_;           /// input data
    std::bitset<NUM_MATERIAL> validData_; /// Check that voxel data is correct
public:
    /**
     * @brief constructor
     * @param inputData Input data
     */
    VoxelData(const InputData &inputData)
        : inputData_(inputData) {
      clear();
      const BigUINT numVoxels = inputData_.voxelDims[0] * inputData_.voxelDims[1] * inputData_.voxelDims[2];
      mallocCPUPinned(voxel,numVoxels*NUM_MATERIAL);
      validData_.reset();
    }

    /**
     * @brief add Material Allignment and unalligned data to the input for a given material. \n
     * Once you add the material the bits corresponding to that bit is turned on.
     * @param matAlignementData Alignment data for the material. Must be in the order of Sx,Sy,Sz.
     * @param matUnalignedData The unalignment component for the material
     * @param matID material ID . Start from 1 to NMat
     */
    void addMaterialDataVectorMorphology(py::array_t<Real, py::array::c_style | py::array::forcecast> &matAlignementData,
                         py::array_t<Real, py::array::c_style | py::array::forcecast> &matUnalignedData,
                         const UINT matID) {

      if(inputData_.morphologyType == MorphologyType::EULER_ANGLES){
        py::print("Error: [Expected]: Vector Morphology [Found:] Euler Angles. Returning\n");
        return;
      }
      if (matID > NUM_MATERIAL or (matID == 0)) {
        throw std::logic_error("Number of material does not match with the compiled version. matID must range from 1 to NUM_MATERIAL");
      }
      if(validData_.test(matID-1)){
        py::print("The material is already set. Please first reset to add the entries. Returning.");
        return;
      }
      const BigUINT numVoxels = inputData_.voxelDims[0] * inputData_.voxelDims[1] * inputData_.voxelDims[2];
      std::vector<Real> _matAlignedData(numVoxels*3);
      std::vector<Real> _matUnalignedData(numVoxels);
      for(int i = 0; i < numVoxels; i++){
        _matAlignedData[3*i + 0] = matAlignementData.data()[3*i+0];
        _matAlignedData[3*i + 1] = matAlignementData.data()[3*i+1];
        _matAlignedData[3*i + 2] = matAlignementData.data()[3*i+2];
        _matUnalignedData[i] = matUnalignedData.data()[i];
      }
      if(inputData_.morphologyType == MorphologyOrder::XYZ){
        H5::XYZ_to_ZYX(_matUnalignedData,1,inputData_.voxelDims);
        H5::XYZ_to_ZYX(_matAlignedData,3,inputData_.voxelDims);
      }
      for (BigUINT i = 0; i < numVoxels; i++) {
        voxel[(matID - 1)*numVoxels + i].s1.x = _matAlignedData.data()[i * 3 + 0];
        voxel[(matID - 1)*numVoxels + i].s1.y = _matAlignedData.data()[i * 3 + 1];
        voxel[(matID - 1)*numVoxels + i].s1.z = _matAlignedData.data()[i * 3 + 2];
        voxel[(matID - 1)*numVoxels + i].s1.w = _matUnalignedData.data()[i];
      }
      validData_.set(matID-1, true);
    }

    void addMaterialDataEulerAngles(py::array_t<Real, py::array::c_style | py::array::forcecast> &matSVector,
                         py::array_t<Real, py::array::c_style | py::array::forcecast> &matThetaVector,
                         py::array_t<Real, py::array::c_style | py::array::forcecast> &matPhiVector,
                         py::array_t<Real, py::array::c_style | py::array::forcecast> &matVfracVector,
                         const UINT matID) {

        if(inputData_.morphologyType != MorphologyType::EULER_ANGLES){
          py::print("Error: [Expected]: Euler Angles / Spherical Coordinates. [Found:] VectorMorphology\n");
          return;
        }
        if ((matID > NUM_MATERIAL) or (matID == 0)) {
            throw std::logic_error("Number of material does not match with the compiled version. matID must range from 1 to NUM_MATERIAL");
        }
        if(validData_.test(matID - 1)){
            py::print("The material is already set. Please first reset to add the entries. Returning.");
            return;
        }
        const BigUINT numVoxels = inputData_.voxelDims[0] * inputData_.voxelDims[1] * inputData_.voxelDims[2];
        std::vector<Real> _S(numVoxels);
        std::vector<Real> _Theta(numVoxels);
        std::vector<Real> _Phi(numVoxels);
        std::vector<Real> _Vfrac(numVoxels);
        for(int i = 0; i < numVoxels; i++){
          _S[i]     = matSVector.data()[i];
          _Theta[i] = matThetaVector.data()[i];
          _Phi[i]   = matPhiVector.data()[i];
          _Vfrac[i] = matVfracVector.data()[i];
        }
        if(inputData_.morphologyType == MorphologyOrder::XYZ){
          H5::XYZ_to_ZYX(_S,1,inputData_.voxelDims);
          H5::XYZ_to_ZYX(_Theta,1,inputData_.voxelDims);
          H5::XYZ_to_ZYX(_Phi,1,inputData_.voxelDims);
          H5::XYZ_to_ZYX(_Vfrac,1,inputData_.voxelDims);
        }
        if(inputData_.morphologyType == MorphologyType::EULER_ANGLES) {
          for (BigUINT i = 0; i < numVoxels; i++) {
            voxel[(matID - 1)*numVoxels + i].s1.x = _S[i];
            voxel[(matID - 1)*numVoxels + i].s1.y = _Theta[i];
            voxel[(matID - 1)*numVoxels + i].s1.z = _Phi[i];
            voxel[(matID - 1)*numVoxels + i].s1.w = _Vfrac[i];
          }
        }
        validData_.set(matID - 1, true);
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

      H5::readFile(fname, inputData_.voxelDims, voxel,(MorphologyType)inputData_.morphologyType,inputData_.morphologyOrder, true);
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

      H5::writeXDMF(inputData_,voxel);
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
    const Voxel *data() const {
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
