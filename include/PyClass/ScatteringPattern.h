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
#ifndef CY_RSOXS_PROJECTIONDATA_H
#define CY_RSOXS_PROJECTIONDATA_H
#include <Input/InputData.h>
#include <utils.h>
/**
 * @brief Stores the Scattering pattern data
 */
class ScatteringPattern{

    /// Input data
    const InputData & inputData_;
    /// Pointer to the scattering pattern data
    Real * data_ = nullptr;

public:
    /**
     * @brief Constructor
     * @param InputData inputData
     */
    ScatteringPattern(const InputData & InputData)
    :inputData_(InputData){
        const UINT
            numEnergyLevel = static_cast<UINT>(inputData_.energies.size());
        const std::size_t totalArraySize = static_cast<std::size_t>(numEnergyLevel) * static_cast<std::size_t>(inputData_.voxelDims[0] * inputData_.voxelDims[1]) *
                                         inputData_.kVectors.size();
        data_ = new Real[totalArraySize];
    }

    /**
     * @brief Pointer to the data
     * @return the pointer to the scattering pattern data
     */
    inline Real * data(){
      return data_;
    }

    /**
     * @brief Clears up the memory
     */
    void clear(){
      if(data_ != nullptr) {
        delete[] data_;
        data_ = nullptr;
      }
    }

    /**
     * @brief Destructor
     */
    ~ScatteringPattern(){
      clear();
    }

    /**
     * @brief Writes Scattering pattern data to HDF5 file
     */
    void writeToHDF5(const std::string & dirName = "HDF5") const {
      writeH5(inputData_, inputData_.voxelDims, data_,dirName);
    }

    /**
     * @brief Writes scattering pattern data to VTI paraview format
     */
    void writeToVTI(const std::string & dirName = "VTI") const {
      writeVTI(inputData_, inputData_.voxelDims, data_,dirName);
    }

    /**
     * @brief returns the data in the numpy array. Note that it is the same memory allocation
     * with numpy wrapper.
     * @param energy Energy for which the numpy array is needed
     * @param kID id of k Vector
     * @return numpy numpy array with the scattering pattern data of the energy
     */
    py::array_t<Real> writeToNumpy(const Real energy, const UINT kID = 0) const {
      if(inputData_.caseType == DEFAULT){
        if(kID > 0){
          py::print("kID cannot be greater than 0 for Default");
          return py::array_t<Real>{};
        }
      }
      if(kID >= inputData_.kVectors.size()){
        py::print("[ERROR] kID = ",kID, " must be smaller than kVector Size = ", inputData_.kVectors.size());
        return py::array_t<Real>{};
      }
      const auto & energies = inputData_.energies;
      if((energy < energies[0]) or (energy > energies[energies.size() - 1])){
        py::print("[ERROR]: Wrong EnergyID");
        return py::array_t<Real>{};
      }
      py::print("[INFO] kVector = [",inputData_.kVectors[kID].x,",",inputData_.kVectors[kID].y,",",inputData_.kVectors[kID].z,"]");
      const UINT energyID = std::lower_bound(energies.begin(),energies.end(),energy) - energies.begin();

      if(not(FEQUALS(energies[energyID],energy))){
        py::print("[ERROR]: Wrong EnergyID");
        return py::array_t<Real>{};
      }
      py::capsule free_when_done(this->data_, [](void *f) {
      });
      const std::size_t voxel2DSize = static_cast<std::size_t>(inputData_.voxelDims[0])*static_cast<std::size_t>(inputData_.voxelDims[1]);
      const std::size_t offset = static_cast<std::size_t>(energyID) * voxel2DSize * inputData_.kVectors.size() +
                                 static_cast<std::size_t>(kID*voxel2DSize);
      return (py::array_t<Real>(
      {(int)inputData_.voxelDims[1],(int)inputData_.voxelDims[0]},
      {sizeof(Real)*inputData_.voxelDims[1],sizeof(Real)},
      &this->data_[energyID*(inputData_.voxelDims[0]*inputData_.voxelDims[1])*inputData_.kVectors.size() + kID*(inputData_.voxelDims[0]*inputData_.voxelDims[1])],
          free_when_done));


    }

    // BEGIN code that was contributed by employees of the National Institute of Standards and Technology (NIST), 
    // an agency of the Federal Government and is being made available as a public service.
    // Pursuant to title 17 United States Code Section 105, works of NIST employees are not subject to copyright protection in the United States
    /**
     * @brief returns the all the data at once in the numpy array. Note that it is the same memory allocation
     * with numpy wrapper.
     * @param kID id of k Vector
     * @return numpy array with the scattering pattern data of the energy
     */
    py::array_t<Real> writeAllToNumpy(const UINT kID = 0) const {
      if(inputData_.caseType == DEFAULT){
        if(kID > 0){
          py::print("kID cannot be greater than 0 for Default");
          return py::array_t<Real>{};
        }
      }
      if(kID >= inputData_.kVectors.size()){
        py::print("[ERROR] kID = ",kID, " must be smaller than kVector Size = ", inputData_.kVectors.size());
        return py::array_t<Real>{};
      }
      const auto & energies = inputData_.energies;
      py::print("[INFO] kVector = [",inputData_.kVectors[kID].x,",",inputData_.kVectors[kID].y,",",inputData_.kVectors[kID].z,"]");

      py::capsule free_when_done(this->data_, [](void *f) {
      });

      return (py::array_t<Real>(
      {(int)inputData_.energies.size(), (int)inputData_.voxelDims[1], (int)inputData_.voxelDims[0]},
      {sizeof(Real)*inputData_.voxelDims[1]*inputData_.voxelDims[0], sizeof(Real)*inputData_.voxelDims[1], sizeof(Real)},
      data_, free_when_done));

    }
    // END code that was contributed by employees of NIST


};
#endif //CY_RSOXS_PROJECTIONDATA_H
