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

#ifndef CUDA_BASE_INPUT_H
#define CUDA_BASE_INPUT_H

#include <Datatypes.h>
#include <stdio.h>
#include <string.h>
#include <cstring>

#ifndef BIAXIAL
/// Stores the material refractive index data
struct Material {

  /** Parallel component of Refractive index **/
  Complex npara;

  /** Perpendicular component of Refractive index **/
  Complex nperp;

  /**
   * @brief Constructor
   */
  Material(){
      std::memset (&npara, 0, sizeof (Complex));
      std::memset (&nperp, 0, sizeof (Complex));
  }
};

#endif

/**
 * This is the structure for the voxel data
 * @tparam num_component number of component
 */

struct Voxel {

  enum EULER_ANGLE:int{
    S = 0,
    THETA = 1,
    PSI = 2,
    VFRAC = 3
  };
  /// s1 the uniaxial director vector \n for Vector Morphology
  /// s1[i].x = x component of the director vector (\f$s_x\f$) \n
  /// s1[i].y = y component of the director vector (\f$s_y\f$) \n
  /// s1[i].z = z component of the director vector (\f$s_z\f$) \n
  /// s1[i].w = fraction of unaligned component (\f$\phi_{ua}\f$)
  /// s1 stores the Euler Angle information for Euler Angle
  /// s1[i].x = fraction of aligned component \n
  /// s1[i].y = (\f$\theta\f$)rotation angle about X-axis  \n
  /// s1[i].z = (\f$\psi\f$) Second rotation angle about Z-axis \n
  /// s1[i].w = volume fraction of the material
  Real4 s1;

  /**
   * @brief Constructor
   */
  Voxel():
  s1{0.0,0.0,0.0,0.0}{
  }
  /**
   * @brief Getter
   * @param [in] id 0-3 depending to s.x,s.y,s.z or s.w
   * @return values at given id
   */
  Real  __host__ __device__ getValueAt(const int id) const{
    const Real val[4]{s1.x,s1.y,s1.z,s1.w};
    return val[id];
  }
};

#endif //CUDA_BASE_INPUT_H
