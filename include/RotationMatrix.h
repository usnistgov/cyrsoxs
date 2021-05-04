//
// Created by maksbh on 5/4/21.
//

#ifndef CY_RSOXS_ROTATIONMATRIX_H
#define CY_RSOXS_ROTATIONMATRIX_H

#include <Input/InputData.h>
#include <Rotation.h>

struct BaseConfiguration{
  Matrix matrix;
  Real baseRotAngle;
} ;

class RotationMatrix{
  struct BaseAxis{
    Real3  X;
    Real3  Y;
    Real3  Z;
  };
  const InputData * inputData_;
  std::vector<BaseConfiguration> baseConfig_;
  Matrix detectorRotationMatrix_;
  std::vector<BaseAxis> baseAxis_;
public:
  RotationMatrix(const InputData * inputData);
  void initComputation();

  inline const auto & getBaseConfigurations()const{
    return baseConfig_;
  }

  inline const auto & getDetectorRotationMatrix() const{
    return detectorRotationMatrix_;
  }

  void printToFile(std::ofstream & fout) const;
};

#endif //CY_RSOXS_ROTATIONMATRIX_H
