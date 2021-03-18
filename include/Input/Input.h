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

#ifndef CUDA_BASE_INPUT_H
#define CUDA_BASE_INPUT_H

#include <Datatypes.h>
#include <stdio.h>
#include <string.h>
#include <cstring>

#ifndef BIAXIAL
template<int num_component>
struct Material {

  /** Parallel component of Refractive index **/
  Complex npara[num_component];

  /** Perpendicular component of Refractive index **/
  Complex nperp[num_component];

  Material(){
      std::memset (npara, 0, sizeof (Complex)*num_component);
      std::memset (nperp, 0, sizeof (Complex)*num_component);
  }
};

#endif

/**
 * This is the structure for the voxel data
 * @tparam num_component number of component
 */
template<int num_component>
struct Voxel {
  /// s1 the uniaxial director vector \n
  /// s1[i].x = x component of the director vector (\f$s_x\f$) \n
  /// s1[i].y = y component of the director vector (\f$s_y\f$) \n
  /// s1[i].z = z component of the director vector (\f$s_z\f$) \n
  /// s1[i].w = fraction of unaligned component (\f$\phi_{ua}\f$)

  Real4 s1[num_component];
  __host__ __device__ inline Real  s(const UINT & id) const{
    return s1[id].x;
  }
  __host__ __device__ inline Real  phi(const UINT & id) const{
    return s1[id].y;
  }
  __host__ __device__ inline Real  theta(const UINT & id) const{
    return s1[id].z;
  }
  __host__ __device__ inline Real  vFrac(const UINT & id) const{
    return s1[id].w;
  }

  Voxel(){
      std::memset (s1, 0, sizeof (Real4)*num_component);
  }
};

struct ElectricField {
  /// Polarized electric field
  Real3 e;
  /// k vector
  Real3 k;
};

#endif //CUDA_BASE_INPUT_H
