//
// Created by maksbh on 9/5/20.
//

#ifndef CY_RSOXS_PROJECTIONDATA_H
#define CY_RSOXS_PROJECTIONDATA_H
#include <Input/InputData.h>
#include <utils.h>
class ScatteringPattern{

    const InputData & inputData_;
    Real * data_ = nullptr;
    bool isAllocated_ = false;

public:
    ScatteringPattern(const InputData & InputData)
    :inputData_(InputData){

    }

    void allocate(){

      if(isAllocated_){
        py::print("The memory is already allocated. Please deallocate it first by calling clear");
        return;
      }
      isAllocated_ = true;

      const UINT
          numEnergyLevel = static_cast<UINT>(std::round((inputData_.energyEnd - inputData_.energyStart) / inputData_.incrementEnergy + 1));
      data_ = new Real[numEnergyLevel * inputData_.numX * inputData_.numY];
    }

    inline Real * data(){
      if(not(isAllocated_)){
        py::print("[ERROR] Memory not allocated. Please allocate it first");
      }
      return data_;
    }

    void clear(){
      if(data_ != nullptr) {
        delete[] data_;
        data_ = nullptr;
        isAllocated_ = false;
      }
    }

    ~ScatteringPattern(){
      clear();
    }

    void writeToHDF5(){
      const UINT voxelDimensions[3]{inputData_.numX,inputData_.numY,inputData_.numZ};
      writeH5(inputData_, voxelDimensions, data_);
    }

    void writeToVTI(){
      const UINT voxelDimensions[3]{inputData_.numX,inputData_.numY,inputData_.numZ};
      writeVTI(inputData_, voxelDimensions, data_);
    }


    py::array_t<Real> writeToNumpy(const Real energy){
      if((energy < inputData_.energyStart) or (energy > inputData_.energyEnd)){
        py::print("[LOG]: Wrong EnergyID");
        return py::array_t<Real>{};
      }

      const UINT
      energyID = static_cast<UINT>(std::round( (energy - inputData_.energyStart)/ inputData_.incrementEnergy));

      if(not(FEQUALS(inputData_.energyStart + energyID*inputData_.incrementEnergy,energy))){
        py::print("[LOG]: Wrong EnergyID");
        return py::array_t<Real>{};
      }
      py::capsule free_when_done(this->data_, [](void *f) {
      });

      return (py::array_t<Real>(
      {(int)inputData_.numX,(int)inputData_.numY},
      {sizeof(Real)*inputData_.numY,sizeof(Real)},
      &this->data_[energyID*(inputData_.numY*inputData_.numX)],
          free_when_done));


    }

    void initialize(){
      const UINT
          numEnergyLevel = static_cast<UINT>(std::round((inputData_.energyEnd - inputData_.energyStart) / inputData_.incrementEnergy + 1));
      const UINT size = inputData_.numX*inputData_.numY;
      for(int i = 0; i < numEnergyLevel; i++){
        for(int sz = 0; sz < size; sz++){
          data_[i*size+sz] = i;
        }
      }
    }

};
#endif //CY_RSOXS_PROJECTIONDATA_H
