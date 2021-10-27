//
// Created by maksbh on 10/26/21.
//

#ifndef CY_RSOXS_POLARIZATION_H
#define CY_RSOXS_POLARIZATION_H

#include <Input/InputData.h>
#include <Output/writeH5.h>
/**
 * @brief Polarization computation for debugging
 */
class Polarization {
  const InputData & inputData_;
  Complex * polarizationX_, * polarizationY_, *polarizationZ_;

public:
  /**
   * @brief Constructor
   * @param inputData input data
   */
  Polarization(const InputData& inputData)
  :inputData_(inputData){
    const BigUINT numVoxels = inputData.voxelDims[0]*inputData.voxelDims[1]*inputData.voxelDims[2];
    polarizationX_ = new Complex[numVoxels];
    polarizationY_ = new Complex[numVoxels];
    polarizationZ_ = new Complex[numVoxels];
  }

  /**
   * @brief clears the memory
   */
  void clear(){
    delete [] polarizationX_;
    delete [] polarizationY_;
    delete [] polarizationZ_;
  }
  /**
   * @brief write to HDF5 in the file Polarization.h5
   */
  void writeToHDF5()const {
    createDirectory("Polarization");
    for(int i = 0; i < inputData_.voxelDims[0]*inputData_.voxelDims[1]*inputData_.voxelDims[2]; i++){
      polarizationX_[i].x = i;
      polarizationX_[i].y = 10*i;
    }
    H5::writePolarization(polarizationX_,inputData_,"Polarization/PolarizationX","pX");
    H5::writePolarization(polarizationY_,inputData_,"Polarization/PolarizationY","pY");
    H5::writePolarization(polarizationZ_,inputData_,"Polarization/PolarizationZ","pZ");
  }

  Complex * getData(int id){
    if(id == 0){
      return polarizationX_;
    }
    else if(id == 1){
      return polarizationY_;
    }
    else if(id == 2){
      return polarizationZ_;
    }
    else{
      throw std::runtime_error("Wrong id");
    }
  }
  /**
     * @brief returns the data in the numpy array. Note that it is the same memory allocation
     * with numpy wrapper.
     * @param id = 0/1/2 for pX/py/pZ
     * @return numpy numpy array with the scattering pattern data of the energy
    */
  py::array_t<std::complex<Real>> writeToNumpy(int id) const {
    for(int i = 0; i < inputData_.voxelDims[0]*inputData_.voxelDims[1]*inputData_.voxelDims[2]; i++){
      polarizationX_[i].x = i;
      polarizationX_[i].y = 10*i;
    }
    std::complex<Real> * data;
    if(id >= 3){
      pybind11::print("Wrong id. Must be less than 2.");
      return py::array_t<std::complex<Real>>{};
    }
    if(id == 0) {
      data = (std::complex<Real> *)polarizationX_;
    }
    else if(id == 1){
      data = (std::complex<Real> *)polarizationY_;
    }
    else if(id == 2){
      data = (std::complex<Real> *)polarizationZ_;
    }

    py::capsule free_when_done(data, [](void *f) {

    });
    return (py::array_t<std::complex<Real>>(
      {(int) inputData_.voxelDims[2], (int) inputData_.voxelDims[1], (int) inputData_.voxelDims[0]}, /// C++ order (X-> fastest)
      {sizeof(std::complex<Real>) * inputData_.voxelDims[0]*inputData_.voxelDims[1],sizeof(std::complex<Real>) * inputData_.voxelDims[1], sizeof(std::complex<Real>)},
      data,
      free_when_done));
  }

  /**
   * @brief Destructor
   */
  ~Polarization(){
    clear();
  }
};
#endif //CY_RSOXS_POLARIZATION_H
