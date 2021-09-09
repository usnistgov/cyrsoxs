/////////////////////////////////////////////////////////////////////////////////
// MIT License
//
//Copyright (c) 2019 - 2021 Iowa State University
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
