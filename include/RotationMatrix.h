//
// Created by maksbh on 5/4/21.
//

#ifndef CY_RSOXS_ROTATIONMATRIX_H
#define CY_RSOXS_ROTATIONMATRIX_H

#include <Input/InputData.h>
#include <Rotation.h>

/**
 * @brief Stores the base configuration.
 */
struct BaseConfiguration{
  /// rotation Matrix
  Matrix matrix;
  /// rotation Angle so that E.x $\approx$ 1
  Real baseRotAngle;
} ;

/**
 * @brief Class to compute the rotation matrix for all k Vector
 */
class RotationMatrix{
  struct BaseAxis{
    Real3  X;
    Real3  Y;
    Real3  Z;
  };
  /// input Data
  const InputData * inputData_;
  /// base configuration vector
  std::vector<BaseConfiguration> baseConfig_;
  /// rotation matrix detector
  Matrix detectorRotationMatrix_;
  /// The axis value for X, Y, Z for 0 degree E rotation
  std::vector<BaseAxis> baseAxis_;
public:
  /**
   * @ brief Constructor
   * @param inputData input data
   */
  RotationMatrix(const InputData * inputData);
  /**
   * @brief computes the required rotation matrices
   */
  void initComputation();

  /**
   * @brief Getter
   * @return the base configuration
   */
  inline const auto & getBaseConfigurations()const{
    return baseConfig_;
  }

  /**
   * @brief Getter
   * @return Rotation matrix for detector
   */
  inline const auto & getDetectorRotationMatrix() const{
    return detectorRotationMatrix_;
  }

  /**
   * @brief Prints the information to File
   * @param fout ofstream object
   */
  void printToFile(std::ofstream & fout) const;
};

#endif //CY_RSOXS_ROTATIONMATRIX_H
