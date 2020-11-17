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

#ifndef UNIAXIAL_H
#define UNIAXIAL_H

#include <Datatypes.h>
#include "cudaHeaders.h"
#include "cudaUtils.h"
#include <Input/Input.h>
#include <cmath>
#include <complex>
#include <math.h>

#ifdef EOC
#include <opencv2/opencv.hpp>
#endif

#include <cstring>

#include <omp.h>
#include <Output/writeVTI.h>



/**
 * @brief This function computes the polarization in real space for the uniaxial case.
 * @param [in] material material data for a particular energy level under consideration.
 * @param [in] angle The angle of rotation.
 * @param [in] voxelInput voxel Input data.
 * @param [in] threadID threadID to access the index of the entry.
 * @param [out] polarizationX
 * @param [out] polarizationY
 * @param [out] polarizationZ
 */

__device__ void computePolarizationUniaxial(const Material<NUM_MATERIAL> *material, const Real angle,
                                            const Voxel<NUM_MATERIAL> *voxelInput, const BigUINT threadID,
                                            Complex *polarizationX, Complex *polarizationY, Complex *polarizationZ) {

  Complex pX, pY, pZ, npar[NUM_MATERIAL], nper[NUM_MATERIAL];
  Real4 s1[NUM_MATERIAL];

  pX.x = 0;
  pX.y = 0;
  pY.x = 0;
  pY.y = 0;
  pZ.x = 0;
  pZ.y = 0;
  const Real cos_a = cos(angle);
  const Real sin_a = sin(angle);
  static constexpr Real OneBy4Pi = static_cast<Real> (1.0 / (4.0 * M_PI));
  Real temp;
  for (int i = 0; i < NUM_MATERIAL; i++) {
    s1[i] = voxelInput[threadID].s1[i];
    temp = cos_a * voxelInput[threadID].s1[i].y - sin_a * voxelInput[threadID].s1[i].x;
    s1[i].x = sin_a * s1[i].y + cos_a * s1[i].x;
    s1[i].y = temp;
    npar[i] = material->npara[i];
    nper[i] = material->nperp[i];
  }

  Complex nsum;
  for (int i = 0; i < NUM_MATERIAL; i++) {
    const Real &sx = s1[i].x;
    const Real &sy = s1[i].y;
    const Real &sz = s1[i].z;
    const Real &phi_ui = s1[i].w;

    Real phi = phi_ui + sx * sx + sy * sy + sz * sz;

    nsum.x = npar[i].x + 2 * nper[i].x;
    nsum.y = npar[i].y + 2 * nper[i].y;

    computeComplexSquare(nsum);
    computeComplexSquare(npar[i]);
    computeComplexSquare(nper[i]);

//    pX.x += (npar[i].x - nper[i].x) * s1[i].z * s1[i].x;
//    pX.y += (npar[i].y - nper[i].y) * s1[i].z * s1[i].x;
//
//    pY.x += (npar[i].x - nper[i].x) * s1[i].z * s1[i].y;
//    pY.y += (npar[i].y - nper[i].y) * s1[i].z * s1[i].y;
//
//    pZ.x +=
//        (s1[i].w * nsum.x) / 9.0 + npar[i].x * s1[i].z * s1[i].z + nper[i].x * (s1[i].x * s1[i].x + s1[i].y * s1[i].y)
//            - phi;
//    pZ.y +=
//        (s1[i].w * nsum.y) / 9.0 + npar[i].y * s1[i].z * s1[i].z + nper[i].y * (s1[i].x * s1[i].x + s1[i].y * s1[i].y);

/**
 *  npar^2*sx^2 + nper^2*sy^2 + nper^2*sz^2
                 sx*sy*(npar^2 - nper^2)
                 sx*sz*(npar^2 - nper^2)
 */
    pX.x +=
        (phi_ui * nsum.x) / (Real) 9.0 + npar[i].x * sx * sx + nper[i].x * (sz * sz + sy * sy) - phi;
    pX.y +=
        (phi_ui * nsum.y) / (Real) 9.0 + npar[i].y * sx * sx + nper[i].y * (sz * sz + sy * sy);

    pY.x += (npar[i].x - nper[i].x) * sx * sy;
    pY.y += (npar[i].y - nper[i].y) * sx * sy;

    pZ.x += (npar[i].x - nper[i].x) * sz * sx;
    pZ.y += (npar[i].y - nper[i].y) * sz * sx;

  }

  pX.x *= OneBy4Pi;
  pX.y *= OneBy4Pi;

  pY.x *= OneBy4Pi;
  pY.y *= OneBy4Pi;

  pZ.x *= OneBy4Pi;
  pZ.y *= OneBy4Pi;

  polarizationX[threadID] = pX;
  polarizationY[threadID] = pY;
  polarizationZ[threadID] = pZ;

}

__host__ __device__ inline BigUINT reshape3Dto1D(UINT i, UINT j, UINT k, uint3 voxel) {
  return ((i) + (j) * voxel.x + (k) * voxel.x * voxel.y);
}

__host__ __device__ inline void reshape1Dto3D(BigUINT id, UINT &X, UINT &Y, UINT &Z, uint3 voxel) {
  Z = static_cast<UINT>(id / (voxel.y * voxel.x * 1.0));
  Y = static_cast<UINT>(id - Z * voxel.y * voxel.x) / (voxel.x * 1.0);
  X = static_cast<UINT>(id - Y * voxel.x - Z * voxel.y * voxel.x);
}

template<typename T>
__device__ inline void swap(T &var1, T &var2) {
  T temp;
  temp = var1;
  var1 = var2;
  var2 = temp;
}


template<typename T>
__global__ void FFTIgor(T *polarization, uint3 voxel) {
  UINT threadID = threadIdx.x + blockIdx.x * blockDim.x;

  BigUINT totalSize = (voxel.x * voxel.y * voxel.z);
  if (threadID >= totalSize) {
    return;
  }


  UINT X, Y, Z;
  reshape1Dto3D(threadID, X, Y, Z, voxel);

  uint3 midX;
  midX.x = voxel.x / 2;
  midX.y = voxel.y / 2;
  midX.z = voxel.z / 2;

  uint3 id;
  if (X <= midX.x) {
    id.x = midX.x - X;
  } else {
    id.x = voxel.x + (midX.x - X);
  }

  if (Y <= midX.y) {
    id.y = midX.y - Y;
  } else {
    id.y = voxel.y + (midX.y - Y);
  }

  if (Z <= midX.z) {
    id.z = midX.z - Z;
  } else {
    id.z = voxel.z + (midX.z - Z);
  }


  BigUINT copyID = reshape3Dto1D(id.x, id.y, id.z, voxel);
  if(copyID > threadID){
      swap(polarization[copyID],polarization[threadID]);
  }

}

template<typename T>
__device__ void FFTShift(T *array, uint3 voxel) {
  uint3 mid;
  mid.x = voxel.x / 2;
  mid.y = voxel.y / 2;
  mid.z = voxel.z / 2;
  UINT threadID = threadIdx.x + blockIdx.x * blockDim.x;
  BigUINT totalSize = (voxel.x * voxel.y * voxel.z);
  if (threadID >= totalSize) {
    return;
  }

  UINT X, Y, Z;
  reshape1Dto3D(threadID, X, Y, Z, voxel);
  if (voxel.z == 1) {
    if ((X < mid.x) and (Y < mid.y)) {
      BigUINT copyID = reshape3Dto1D(X + mid.x, Y + mid.y, Z + mid.z, voxel);
      BigUINT currID = threadID;
      swap(array[currID], array[copyID]);

      copyID = reshape3Dto1D(mid.x + X, Y, mid.z + Z, voxel);
      currID = reshape3Dto1D(X, mid.y + Y, Z, voxel);
      swap(array[currID], array[copyID]);
    }
    return;
  }
  if ((X < mid.x) and (Y < mid.y) and (Z < mid.z)) {
    BigUINT copyID = reshape3Dto1D(X + mid.x, Y + mid.y, Z + mid.z, voxel);
    BigUINT currID = threadID;
    swap(array[currID], array[copyID]);

    copyID = reshape3Dto1D(mid.x + X, mid.y + Y, Z, voxel);
    currID = reshape3Dto1D(X, Y, mid.z + Z, voxel);
    swap(array[currID], array[copyID]);

    copyID = reshape3Dto1D(mid.x + X, Y, mid.z + Z, voxel);
    currID = reshape3Dto1D(X, mid.y + Y, Z, voxel);
    swap(array[currID], array[copyID]);

    copyID = reshape3Dto1D(X, mid.y + Y, mid.z + Z, voxel);
    currID = reshape3Dto1D(mid.x + X, Y, Z, voxel);
    swap(array[currID], array[copyID]);

    copyID = reshape3Dto1D(mid.x + X, mid.y + Y, mid.z + Z, voxel);
    currID = reshape3Dto1D(X, Y, mid.z + Z, voxel);
    swap(array[currID], array[copyID]);
  }


}


/**
 * @brief This function computes the Scatter3D computation.
 *
 * @param [in] polarizationX  X component of polarization in Fourier space
 * @param [in] polarizationY  Y component of polarization in Fourier space
 * @param [in] polarizationZ  Z component of polarization in Fourier space
 * @param [out] Scatter3D     Scatter 3D result
 * @param [in] elefield       Electric field
 * @param [in] voxelNum       Number of total voxel
 * @param [in] voxel          Number of voxel in each direciton.
 * @param [in] physSize       Physical Size
 */
__global__ void computeScatter3D(Complex *polarizationX,
                                 Complex *polarizationY,
                                 Complex *polarizationZ,
                                 Real *Scatter3D,
                                 const ElectricField elefield,
                                 const Real eAngle,
                                 const Real kAngle,
                                 const BigUINT voxelNum,
                                 const uint3 voxel,
                                 const Real physSize,
                                 const bool enable2D) {


  UINT threadID = threadIdx.x + blockIdx.x * blockDim.x;
  if (threadID >= voxelNum) {
    return;
  }

  UINT X, Y, Z;
  reshape1Dto3D(threadID, X, Y, Z, voxel);

  Real3 dx;
  dx.x = static_cast<Real>((2 * M_PI / physSize) / ((voxel.x - 1) * 1.0));
  dx.y = static_cast<Real>((2 * M_PI / physSize) / ((voxel.y - 1) * 1.0));
  dx.z = 0;
  if (not(enable2D)) {
    dx.z = static_cast<Real>((2 * M_PI / physSize) / ((voxel.z - 1) * 1.0));
  }
  Real3 q;
  q.x = static_cast<Real>((-M_PI / physSize) + X * dx.x);
  q.y = static_cast<Real>((-M_PI / physSize) + Y * dx.y);
  q.z = 0;
  if (not(enable2D)) {
    q.z = static_cast<Real>((-M_PI / physSize) + Z * dx.z);
  }
 const Real cosPhi   = cos(eAngle);
 const Real cosTheta = cos(kAngle);
 const Real sinPhi   = sin(eAngle);
 const Real sinTheta = sin(kAngle);

 const Real &k = elefield.k.z;
 const Real &qX = q.x;
 const Real &qY = q.y;
 const Real &qZ = q.z;

 Complex p1,p2,p3;

 p1.x = polarizationX[threadID].x*cosPhi + polarizationY[threadID].x*sinPhi;
 p1.y = polarizationX[threadID].y*cosPhi + polarizationY[threadID].y*sinPhi;

 p2.x = -polarizationX[threadID].x*sinPhi*cosTheta + polarizationY[threadID].x*cosPhi*cosTheta -polarizationZ[threadID].x*sinTheta;
 p2.y = -polarizationX[threadID].y*sinPhi*cosTheta + polarizationY[threadID].y*cosPhi*cosTheta -polarizationZ[threadID].y*sinTheta;

 p3.x = -polarizationX[threadID].x*sinPhi*sinTheta + polarizationY[threadID].x*cosPhi*sinTheta +polarizationZ[threadID].x*cosTheta;
 p3.y = -polarizationX[threadID].y*sinPhi*sinTheta + polarizationY[threadID].y*cosPhi*sinTheta +polarizationZ[threadID].y*cosTheta;

 Real q1,q2,q3;

 q1 = qX*cosPhi + qY*sinPhi;
 q2 = -qX*sinPhi*cosTheta + qY*cosPhi*cosTheta -qZ*sinTheta;
 q3 = -qX*sinPhi*sinTheta + qY*cosPhi*sinTheta +qZ*cosTheta;

 const Complex  pVec[3]{p1,p2,p3};
 const Real     qVec[3]{q3,q2,k+q1};

 Scatter3D[threadID] = computeMagVec1TimesVec1TTimesVec2(qVec,pVec,k);
}

/**
 * This function computes the equivalent id of the given position.
 * The id is calculated keeping in mind the effect of polarization FFT
 * shift done in cufft calculation with respect to Igor.
 * @param [in] pos position of which the index is to be found.
 * @param [in] i corresponding X coordinate
 * @param [in] j corresponding Y coordinate
 * @param [in] start start position
 * @param [in] dx the grid podition in each direction
 * @param [in] voxel voxel dimension
 * @return the equivalent id for the given position.
 */
__host__ __device__ BigUINT computeEquivalentID(const Real3 pos,
                                                const UINT i,
                                                const UINT j,
                                                const Real start,
                                                const Real3 dx,
                                                const uint3 voxel,
                                                const bool enable2D) {

  UINT k = 0;
  if (not(enable2D)) {
    k = static_cast<UINT >(round((pos.z - start) / (dx.z)));
  }
  BigUINT id = k * voxel.x * voxel.y + j * voxel.x + i;
  return id;
}


__device__ inline Real computeBilinearInterpolation(const Real *data,
                                                             const Real3 & pos,
                                                             const Real & start,
                                                             const Real3 & dx,
                                                             const UINT & X,
                                                             const UINT & Y,
                                                             const uint3 & voxel
) {

  static constexpr unsigned short num_neigbour = 4;
  Real2 diffX;
  diffX.x = (pos.x - (start + X * dx.x)) / dx.x;
  diffX.y = (pos.y - (start + Y * dx.y)) / dx.y;

  Real _data[num_neigbour];

  _data[0] = data[reshape3Dto1D(X, Y, 0, voxel)];                /// (X,Y)
  _data[1] = data[reshape3Dto1D(X + 1, Y, 0, voxel)];            /// (X + 1,Y)
  _data[2] = data[reshape3Dto1D(X, Y + 1, 0, voxel)];           /// (X,Y + 1)
  _data[3] = data[reshape3Dto1D(X + 1, Y + 1, 0, voxel)];        /// (X + 1,Y + 1)

  Real valx1 = ((1 - diffX.x) * (_data[0]) + (diffX.x) * (_data[1])); /// Interpolate along X1
  Real valx2 = ((1 - diffX.x) * (_data[2]) + (diffX.x) * (_data[3])); /// Interpolate along X2

  Real valy = ((1 - diffX.y) * valx1 + (diffX.y) * valx2);
  return (valy);


}

__device__ inline Real computeTrilinearInterpolation(const Real *data,
                                                              const Real3 & pos,
                                                              const Real & start,
                                                              const Real3 & dx,
                                                              const UINT & X,
                                                              const UINT & Y,
                                                              const UINT & Z,
                                                              const uint3 & voxel
) {
  static constexpr unsigned short num_neigbour = 8;
  Real3 diffX;
  diffX.x = (pos.x - (start + X * dx.x)) / dx.x;
  diffX.y = (pos.y - (start + Y * dx.y)) / dx.y;
  diffX.z = (pos.z - (start + Z * dx.z)) / dx.z;


  Real _data[num_neigbour];

  if ((Z + 1) >= voxel.z) {
    return NAN;
  }

  _data[0] = data[reshape3Dto1D(X, Y, Z, voxel)];                /// (X,Y,Z)
  _data[1] = data[reshape3Dto1D(X + 1, Y, Z, voxel)];            /// (X + 1,Y,Z)
  _data[2] = data[reshape3Dto1D(X, Y + 1, Z, voxel)];           /// (X,Y + 1,Z)
  _data[3] = data[reshape3Dto1D(X + 1, Y + 1, Z, voxel)];        /// (X + 1,Y + 1,Z)

  _data[4] = data[reshape3Dto1D(X, Y, Z + 1, voxel)];             /// (X,Y,Z + 1)
  _data[5] = data[reshape3Dto1D(X + 1, Y, Z + 1, voxel)];           /// (X + 1,Y,Z + 1)
  _data[6] = data[reshape3Dto1D(X, Y + 1, Z + 1, voxel)];          /// (X,Y + 1,Z + 1)
  _data[7] = data[reshape3Dto1D(X + 1, Y + 1, Z + 1, voxel)];       /// (X + 1,Y + 1,Z + 1)

  /*** https://en.wikipedia.org/wiki/Trilinear_interpolation **/
  Real c00 = ((1 - diffX.x) * (_data[0]) + (diffX.x) * (_data[1])); /// Interpolate along X1
  Real c01 = ((1 - diffX.x) * (_data[2]) + (diffX.x) * (_data[3])); /// Interpolate along X2
  Real c10 = ((1 - diffX.x) * (_data[4]) + (diffX.x) * (_data[5])); /// Interpolate along X3 Z is changed
  Real c11 = ((1 - diffX.x) * (_data[6]) + (diffX.x) * (_data[7])); /// Interpolate along X4  Z is chaged

  Real c0 = ((1 - diffX.y) * (c00) + (diffX.y) * (c01)); /// Interpolate along Y1
  Real c1 = ((1 - diffX.y) * (c10) + (diffX.y) * (c11)); /// Interpolate along Y2

  Real c = ((1 - diffX.z) * c0 + (diffX.z) * c1); /// Interpolate along Z

  return (c);

}

/**
 * @brief This function computes the equivalent Projection of 3D Scatter
 * field on Ewald's sphere on GPU. In this case we assume a constant k in
 * the z - direction. In future, we would rotate the k vector in plain
 * to gain more statistics.
 *
 * @param [out] projection The projection result
 * @param [in] scatter3D Scatter3D result.
 * @param [in] voxel Number of voxel in each direction
 * @param [in] k Electric field k.
 * @param [in] physSize Physical Size.
 */
/// TODO: template it
__global__ void computeEwaldProjectionGPU(Real *projection,
                                          const Real *scatter3D,
                                          const uint3 voxel,
                                          const Real k,
                                          const Real eAngle,
                                          const Real kAngle,
                                          const Real physSize,
                                          const Interpolation::EwaldsInterpolation interpolation,
                                          const bool enable2D) {
  UINT threadID = threadIdx.x + blockIdx.x * blockDim.x;
  const UINT totalSize = voxel.x * voxel.y;
  if (threadID >= totalSize) {
    return;
  }
  Real val, start;
  Real3 dx, pos;
  start = -static_cast<Real>(M_PI / physSize);

  const Real & kx = 0;
  const Real & ky = 0;
  const Real & kz = k;

  dx.x = static_cast<Real>((2.0 * M_PI / physSize) / ((voxel.x - 1) * 1.0));
  dx.y = static_cast<Real>((2.0 * M_PI / physSize) / ((voxel.y - 1) * 1.0));
  dx.z = 0;
  if (not(enable2D)) {
    dx.z = static_cast<Real>((2.0 * M_PI / physSize) / ((voxel.z - 1) * 1.0));
  }
  UINT Y = static_cast<UINT>(threadID / (voxel.x * 1.0));
  UINT X = static_cast<UINT>(threadID - Y * voxel.x);
  pos.y = (start + Y * dx.y) ;
  pos.x = (start + X * dx.x) ;

  val = k * k - (pos.x + kx) * (pos.x + kx) - (pos.y + ky) * (pos.y + ky);

  if((val < 0) or (X == (voxel.x - 1)) or (Y == (voxel.y - 1))) {
    projection[threadID] = NAN;
  }
  else
  {
    pos.z = -kz + sqrt(val);
    if (interpolation == Interpolation::EwaldsInterpolation::NEARESTNEIGHBOUR) {
      BigUINT id = computeEquivalentID(pos, X, Y, start, dx, voxel, enable2D);
      projection[threadID] += scatter3D[id];
    }
    else {
      if (enable2D) {
        projection[threadID] += computeBilinearInterpolation(scatter3D, pos, start, dx, X, Y, voxel);
      } else {
        UINT Z = static_cast<UINT >(((pos.z - start) / (dx.z)));
        projection[threadID] += computeTrilinearInterpolation(scatter3D, pos, start, dx, X, Y, Z, voxel);
      }
    }
  }
}
/**
 * @brief This function computes the rotation masks.
 * If during the rotation, some values has been NANs, because they do not belong
 * within the boundary, this function resets it to zero. In addition, it computes a
 * mask counter keeping track of which pixels has been NANs as we do not want to average
 * the final result which has been NANs.
 * @param [in,out] rotProjection the rotated projection
 * @param [out] mask the mask vector
 * @param [in] Number of voxel in each direction
 */
__global__ void computeRotationMask(Real *rotProjection, UINT *mask, const uint3 voxel) {
  UINT threadID = threadIdx.x + blockIdx.x * blockDim.x;
  const UINT totalSize = voxel.x * voxel.y;
  if (threadID >= totalSize) {
    return;
  }
  if (isnan(rotProjection[threadID])) {
    rotProjection[threadID] = 0.0;
  } else {
    mask[threadID]++;
  }

}

/**
 * @brief This functions performs the averaging of the pixels keeping track of how
 * many values are Non NANs for each pixel.
 * @param rotProjection rotation vector
 * @param mask mask vector
 * @param voxel Number of voxel in each direction
 */
__global__ void averageRotation(Real *rotProjection, const UINT *mask, const uint3 voxel) {
  UINT threadID = threadIdx.x + blockIdx.x * blockDim.x;
  const UINT totalSize = voxel.x * voxel.y;
  if (threadID >= totalSize) {
    return;
  }
  if(mask[threadID] == 0){
      rotProjection[threadID] = 0;
  }
  else {
      rotProjection[threadID] = rotProjection[threadID] * static_cast<Real>(1.0 / (mask[threadID] * 1.0));
  }

}

#ifdef EOC

/**
 * @brief This function computes the equivalent Projection of 3D Scatter
 * field on Ewald's sphere on GPU. In this case we assume a constant k in
 * the z - direction. In future, we would rotate the k vector in plain
 * to gain more statistics.
 * @param [out] projection The projection result
 * @param [in] scatter3D Scatter3D result.
 * @param [in] voxel Number of voxel in each direction
 * @param [in] k Electric field k.
 */

__host__ void computeEwaldProjectionCPU(Real *projection, const Real *scatter3D, const uint3 voxel, const Real k) {
  Real val, start;
  Real3 dx, pos;
  start = -static_cast<Real>(M_PI / PHYSSIZE);

  dx.x = static_cast<Real>((2.0 * M_PI / PHYSSIZE) / ((voxel.x - 1) * 1.0));
  dx.y = static_cast<Real>((2.0 * M_PI / PHYSSIZE) / ((voxel.y - 1) * 1.0));
#if (ENABLE_2D)
  dx.z = 0;
#else
  dx.z = static_cast<Real>((2.0 * M_PI / PHYSSIZE) / ((voxel.z - 1) * 1.0));
#endif

//#pragma omp parallel for
  for (UINT j = 0; j < voxel.y; j++) {
    pos.y = start + j * dx.y;
    for (UINT i = 0; i < voxel.x; i++) {
      pos.x = start + i * dx.x;
      val = k * k - pos.x * pos.x - pos.y * pos.y;
      if (val >= 0) {

        pos.z = -k + sqrt(val);
        BigUINT id = computeEquivalentID(pos, i, j, start, dx, voxel);
        projection[j * voxel.x + i] = scatter3D[id];
      } else {
        projection[j * voxel.x + i] = NAN;
      }
    }
  }
}

/**
 * This function rotates the image on CPU. The rotation is performed with the library openCV.
 * @param [in,out] projection The rotated projection image
 * @param [in] voxel voxel data
 * @param [in] rotAngle the rotation angle
 */
__host__ void rotateImage(Real *projection, const uint3 voxel, const Real rotAngle){

  const int size[2]{static_cast<int>(voxel.x), static_cast<int>(voxel.y)};

  cv::Mat image(voxel.x,voxel.y,CV_32FC1);
  memcpy(image.data, projection,size[0]*size[1]*sizeof(Real));
  cv::Point2f pc(image.cols/2., image.rows/2.);
  cv::Mat r = cv::getRotationMatrix2D(pc, rotAngle, 1.0);
  cv::Mat rotatedImage;
  warpAffine(image, rotatedImage, r, image.size());
  if (rotatedImage.isContinuous()) {
    std::memcpy(projection,rotatedImage.data,rotatedImage.rows*rotatedImage.cols);
  }
  else{
    throw "Memory not continuous.";
  }
}

#endif
/**
 * @brief This function performs the integration of the projection data.
 * Currently inactive.
 *
 * @param integrate The integrate value.
 * @param angle (min angle , max angle) in radians
 * @param projection Projection input
 * @param voxel 2D voxel size in x and y
 * @param physSize Physical Size
 */
__global__ void radialIntegrate(Real *integrate,
                                const Real2 angle,
                                const Real *projection,
                                const uint2 voxel,
                                const Real physSize) {

  UINT threadID = threadIdx.x + blockIdx.x * blockDim.x;
  const int totalSize = voxel.x * voxel.y;
  if (threadID >= totalSize) {
    return;
  }

  UINT Y = static_cast<UINT>(threadID / (voxel.x * 1.0));
  UINT X = static_cast<UINT>(threadID - Y * voxel.x);
  Real2 dx, pos;
  Real start, val;
  start = -static_cast<Real>(M_PI / physSize);
  dx.x = static_cast<Real>((2.0 * M_PI / physSize) / ((voxel.x - 1) * 1.0));
  dx.y = static_cast<Real>((2.0 * M_PI / physSize) / ((voxel.y - 1) * 1.0));

  pos.x = start + Y * dx.x;
  pos.y = start + X * dx.y;

  bool check1 = isnan(projection[threadID]); // if projection is not nan

  val = atan(-pos.y / -pos.x);

  bool check2;
  {

    bool check2_1 = val > angle.x;
    bool check2_2 = val < angle.y;

    check2 = check2_1 and check2_2;
  }
  bool check3;
  {
    bool check3_1 = val > (angle.x - M_PI);
    bool check3_2 = val < (angle.y - M_PI);
    check3 = check3_1 and check3_2;
  }

  bool check = not(check1) and (check2 or check3);

  Real radialDist = sqrt(pos.y * pos.y + pos.x * pos.x);

  if (check) {
    integrate[threadID] = radialDist * projection[threadID];
  } else {
    integrate[threadID] = 0;
  }
}

#endif //WRITER_UNIAXIAL_H
