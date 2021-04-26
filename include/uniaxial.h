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
#include <Rotation.h>

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

__device__ void computePolarizationEulerAngles(const Material<NUM_MATERIAL> *material, const Real angle,
                                            const Voxel<NUM_MATERIAL> *voxelInput, const BigUINT threadID,
                                            Complex *polarizationX, Complex *polarizationY, Complex *polarizationZ) {

  Complex pX, pY, pZ;

  pX.x = 0;
  pX.y = 0;
  pY.x = 0;
  pY.y = 0;
  pZ.x = 0;
  pZ.y = 0;
  const Real cosAngle = cos(angle);
  const Real sinAngle = sin(angle);
  static constexpr Real OneBy4Pi = static_cast<Real> (1.0 / (4.0 * M_PI));

  Complex nsum;
  /**
 * [0 1 2]
 * [1 3 4]
 * [2 4 5]
 */
  static constexpr  UINT MATINDEXID[3][3]{{0,1,2}, {1,3,4},{2,4,5}};
  Complex rotatedNr[6]; // Only storing what is required
  memset(rotatedNr,0,sizeof(Complex)*6);
  /// TODO: This loop is redundant and can be pre-computed over all the energy level at the cost of communication.
  for (int i = 0; i < NUM_MATERIAL; i++) {
    const Real & S     = voxelInput[threadID].s(i);
    const Real & Phi   = voxelInput[threadID].phi(i);
    const Real & Theta = voxelInput[threadID].theta(i);
    const Real & Vfrac = voxelInput[threadID].vFrac(i);

    Complex  npar = material->npara[i];
    Complex  nper = material->nperp[i];
    const Real &sx = S*cos(Phi);
    const Real &sy = S*sin(Phi)*cos(Theta);
    const Real &sz = S*sin(Phi)*sin(Theta);
    const Real &phi_ui = Vfrac - S;

    const Real & phi = Vfrac;
    if(threadID == 14669) {
      printf("New: %d MatId = %d, S1 = (%f, %f %f %f &f %f )\n", threadID, i, sx, sy, sz, phi_ui,phi);
    }
    nsum.x = npar.x + 2 * nper.x;
    nsum.y = npar.y + 2 * nper.y;

    computeComplexSquare(nsum);
    computeComplexSquare(npar);
    computeComplexSquare(nper);

    rotatedNr[0].x += npar.x*sx*sx + nper.x*(sy*sy + sz*sz) + ((phi_ui * nsum.x) / (Real) 9.0) - phi;
    rotatedNr[0].y += npar.y*sx*sx + nper.y*(sy*sy + sz*sz) + ((phi_ui * nsum.y) / (Real) 9.0);

    rotatedNr[1].x += (npar.x - nper.x)*sx*sy;
    rotatedNr[1].y += (npar.y - nper.y)*sx*sy;

    rotatedNr[2].x += (npar.x - nper.x)*sx*sz;
    rotatedNr[2].y += (npar.y - nper.y)*sx*sz;

    rotatedNr[3].x += npar.x*sy*sy + nper.x*(sx*sx + sz*sz) + ((phi_ui * nsum.x) / (Real) 9.0) - phi;
    rotatedNr[3].y += npar.y*sy*sy + nper.y*(sx*sx + sz*sz) + ((phi_ui * nsum.y) / (Real) 9.0);

    rotatedNr[4].x += (npar.x - nper.x)*sy*sz;
    rotatedNr[4].y += (npar.y - nper.y)*sy*sz;

    rotatedNr[5].x += npar.x*sz*sz + nper.x*(sx*sx + sy*sy) +  + ((phi_ui * nsum.x) / (Real) 9.0) - phi;
    rotatedNr[5].y += npar.y*sz*sz + nper.y*(sx*sx + sy*sy) +  + ((phi_ui * nsum.y) / (Real) 9.0);
  }
  pX.x = rotatedNr[MATINDEXID[0][0]].x*cosAngle + rotatedNr[MATINDEXID[0][1]].x*sinAngle;
  pX.y = rotatedNr[MATINDEXID[0][0]].y*cosAngle + rotatedNr[MATINDEXID[0][1]].y*sinAngle;

  pY.x = rotatedNr[MATINDEXID[1][0]].x*cosAngle + rotatedNr[MATINDEXID[1][1]].x*sinAngle;
  pY.y = rotatedNr[MATINDEXID[1][0]].y*cosAngle + rotatedNr[MATINDEXID[1][1]].y*sinAngle;

  pZ.x = rotatedNr[MATINDEXID[2][0]].x*cosAngle + rotatedNr[MATINDEXID[2][1]].x*sinAngle;
  pZ.y = rotatedNr[MATINDEXID[2][0]].y*cosAngle + rotatedNr[MATINDEXID[2][1]].y*sinAngle;

  pX.x *= OneBy4Pi;
  pX.y *= OneBy4Pi;

  pY.x *= OneBy4Pi;
  pY.y *= OneBy4Pi;

  pZ.x *= OneBy4Pi;
  pZ.y *= OneBy4Pi;
  RotateZ(pX,pY,pZ,angle);
  polarizationX[threadID] = pX;
  polarizationY[threadID] = pY;
  polarizationZ[threadID] = pZ;

}

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

template<ReferenceFrame referenceFrame>
__device__ void computePolarizationVectorMorphologyOptimized(const Material<NUM_MATERIAL> *material, const Real angle,
                                                    const Voxel<NUM_MATERIAL> *voxelInput, const BigUINT threadID,
                                                    Complex *polarizationX, Complex *polarizationY, Complex *polarizationZ,
                                                    const Matrix rotationMatrix) {

  Complex pX, pY, pZ;
  pX.x = 0;
  pX.y = 0;
  pY.x = 0;
  pY.y = 0;
  pZ.x = 0;
  pZ.y = 0;

  const Real cosAngle = cos(angle);
  const Real sinAngle = sin(angle);
  static constexpr Real OneBy4Pi = static_cast<Real> (1.0 / (4.0 * M_PI));
  /**
 * [0 1 2]
 * [1 3 4]
 * [2 4 5]
 */
  static constexpr  UINT MATINDEXID[3][3]{{0,1,2}, {1,3,4},{2,4,5}};
  Complex rotatedNr[6]; // Only storing what is required
  memset(rotatedNr,0,sizeof(Complex)*6);
  Complex nsum;
  /// TODO: This loop is redundant and can be pre-computed over all the energy level at the cost of communication.
  for (int i = 0; i < NUM_MATERIAL; i++) {
    Complex npar = material->npara[i];
    Complex nper = material->nperp[i];
    const Real &sx = voxelInput[threadID].s1[i].x;
    const Real &sy = voxelInput[threadID].s1[i].y;
    const Real &sz = voxelInput[threadID].s1[i].z;
    const Real &phi_ui = voxelInput[threadID].s1[i].w;

    const Real & phi = phi_ui + sx * sx + sy * sy + sz * sz;

    nsum.x = npar.x + 2 * nper.x;
    nsum.y = npar.y + 2 * nper.y;

    computeComplexSquare(nsum);
    computeComplexSquare(npar);
    computeComplexSquare(nper);

    rotatedNr[0].x += npar.x*sx*sx + nper.x*(sy*sy + sz*sz) + ((phi_ui * nsum.x) / (Real) 9.0) - phi;
    rotatedNr[0].y += npar.y*sx*sx + nper.y*(sy*sy + sz*sz) + ((phi_ui * nsum.y) / (Real) 9.0);

    rotatedNr[1].x += (npar.x - nper.x)*sx*sy;
    rotatedNr[1].y += (npar.y - nper.y)*sx*sy;

    rotatedNr[2].x += (npar.x - nper.x)*sx*sz;
    rotatedNr[2].y += (npar.y - nper.y)*sx*sz;

    rotatedNr[3].x += npar.x*sy*sy + nper.x*(sx*sx + sz*sz) + ((phi_ui * nsum.x) / (Real) 9.0) - phi;
    rotatedNr[3].y += npar.y*sy*sy + nper.y*(sx*sx + sz*sz) + ((phi_ui * nsum.y) / (Real) 9.0);

    rotatedNr[4].x += (npar.x - nper.x)*sy*sz;
    rotatedNr[4].y += (npar.y - nper.y)*sy*sz;

    rotatedNr[5].x += npar.x*sz*sz + nper.x*(sx*sx + sy*sy) +  + ((phi_ui * nsum.x) / (Real) 9.0) - phi;
    rotatedNr[5].y += npar.y*sz*sz + nper.y*(sx*sx + sy*sy) +  + ((phi_ui * nsum.y) / (Real) 9.0);
  }
  static constexpr Real3 eleField{1,0,0};
  Real3 matVec;
  doMatVec<false>(rotationMatrix,eleField,matVec);
  pX.x = rotatedNr[MATINDEXID[0][0]].x*matVec.x + rotatedNr[MATINDEXID[0][1]].x*matVec.y + rotatedNr[MATINDEXID[0][2]].x*matVec.z;
  pX.y = rotatedNr[MATINDEXID[0][0]].y*matVec.x + rotatedNr[MATINDEXID[0][1]].y*matVec.y + rotatedNr[MATINDEXID[0][2]].y*matVec.z;

  pY.x = rotatedNr[MATINDEXID[1][0]].x*matVec.x + rotatedNr[MATINDEXID[1][1]].x*matVec.y + rotatedNr[MATINDEXID[1][2]].x*matVec.z;
  pY.y = rotatedNr[MATINDEXID[1][0]].y*matVec.x + rotatedNr[MATINDEXID[1][1]].y*matVec.y + rotatedNr[MATINDEXID[1][2]].y*matVec.z;

  pZ.x = rotatedNr[MATINDEXID[2][0]].x*matVec.x + rotatedNr[MATINDEXID[2][1]].x*matVec.y + rotatedNr[MATINDEXID[2][2]].y*matVec.z;
  pZ.y = rotatedNr[MATINDEXID[2][0]].y*matVec.x + rotatedNr[MATINDEXID[2][1]].y*matVec.y + rotatedNr[MATINDEXID[2][2]].y*matVec.z;

  pX.x *= OneBy4Pi;
  pX.y *= OneBy4Pi;

  pY.x *= OneBy4Pi;
  pY.y *= OneBy4Pi;

  pZ.x *= OneBy4Pi;
  pZ.y *= OneBy4Pi;
  if(referenceFrame == ReferenceFrame::MATERIAL) {
    rotate<true>(rotationMatrix,pX,pY,pZ);
  }
  polarizationX[threadID] = pX;
  polarizationY[threadID] = pY;
  polarizationZ[threadID] = pZ;
}

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

__device__ void computePolarizationVectorMorphology(const Material<NUM_MATERIAL> *material, const Real angle,
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



inline __device__ Real computeScatter3D(const Complex * polarizationX,
                                        const Complex * polarizationY,
                                        const Complex * polarizationZ,
                                        const Real & k,
                                        const Real & eAngle,
                                        const Real & kAngle,
                                        const Real3 & dX,
                                        const Real & physSize,
                                        const BigUINT & id,
                                        const uint3 & voxel,
                                        const bool enable2D,
                                        const Real3 & kVector
){
    const Real cosPhi   = cos(eAngle);
    const Real cosTheta = cos(kAngle);
    const Real sinPhi   = sin(eAngle);
    const Real sinTheta = sin(kAngle);

    UINT X, Y, Z;
    reshape1Dto3D(id, X, Y, Z, voxel);
    Real3 q;
    q.x = static_cast<Real>((-M_PI / physSize) + X * dX.x);
    q.y = static_cast<Real>((-M_PI / physSize) + Y * dX.y);
    q.z = 0;
    if (not(enable2D)) {
        q.z = static_cast<Real>((-M_PI / physSize) + Z * dX.z);
    }

    Real qVec[3];
    qVec[0] =  k*kVector.x + q.x;
    qVec[1] =  k*kVector.y + q.y;
    qVec[2] =  k*kVector.z + q.z;

    Complex pVec[3]{polarizationX[id],polarizationY[id],polarizationZ[id]};
    return (computeMagVec1TimesVec1TTimesVec2(qVec,pVec,k));
}

/**
* @brief This function computes the Scatter3D computation.
* This is an optimized routine based on the fact that the rotation does not
* change the magnitude.
* @param [in] polarizationX  X component of polarization in Fourier space
* @param [in] polarizationY  Y component of polarization in Fourier space
* @param [in] polarizationZ  Z component of polarization in Fourier space
* @param [out] Scatter3D     Scatter 3D result
* @param [in] elefield       Electric field
* @param [in] voxelNum       Number of total voxel
* @param [in] voxel          Number of voxel in each direciton.
* @param [in] physSize       Physical Size
*/
__global__ void computeScatter3D(const Complex *polarizationX,
                                 const Complex *polarizationY,
                                 const Complex *polarizationZ,
                                 Real *Scatter3D,
                                 const ElectricField elefield,
                                 const Real eAngle,
                                 const Real kAngle,
                                 const BigUINT voxelNum,
                                 const uint3 voxel,
                                 const Real physSize,
                                 const bool enable2D,
                                 const Real3 kVector) {


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

// const Real cosPhi   = cos(eAngle);
// const Real cosTheta = cos(kAngle);
// const Real sinPhi   = sin(eAngle);
// const Real sinTheta = sin(kAngle);
//
// const Real &k = elefield.k.z;
// const Real &qX = q.x;
// const Real &qY = q.y;
// const Real &qZ = q.z;
//
// Real qVec[3];
// qVec[0] = -k*sinPhi*sinTheta + qX;
// qVec[1] =  k*cosPhi*sinTheta + qY;
// qVec[2] =  k*cosTheta + qZ;
//
// Complex pVec[3]{polarizationX[threadID],polarizationY[threadID],polarizationZ[threadID]};

 Scatter3D[threadID] = computeScatter3D(polarizationX,polarizationY,polarizationZ,elefield.k.z,eAngle,kAngle,dx,physSize,threadID,voxel,enable2D,kVector);
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


/**
 * @brief This routine is optimized for the case where pixel location is aligned with X,Y which is the case
 * in these simulations
 * @param data data from which the interpolation happens
 * @param pos position of the desired location
 * @param start start of the domain
 * @param dx spacing
 * @param X X voxel id
 * @param Y Y voxel id
 * @param Z Z voxel id
 * @param voxel Number of voxels
 * @return the interpolated value
 */
__device__ inline Real computeTrilinearInterpolation(const Real *data,
                                                              const Real3 & pos,
                                                              const Real & start,
                                                              const Real3 & dx,
                                                              const UINT & X,
                                                              const UINT & Y,
                                                              const UINT & Z,
                                                              const uint3 & voxel
) {
    if((Z+1) > voxel.z){
        return NAN;
    }
    Real diffZ;
    diffZ = (pos.z - (start + Z * dx.z)) / dx.z;
    Real _data[2];
    _data[0] = data[reshape3Dto1D(X, Y, Z, voxel)];                /// (X,Y,Z)
    _data[1] = data[reshape3Dto1D(X, Y, Z + 1, voxel)];             /// (X,Y,Z + 1)
    Real c = ((1 - diffZ) * _data[0] + (diffZ) * _data[1]); /// Interpolate along Z
    return (c);
}

__device__ inline Real computeTrilinearInterpolation(const Real & data1,
                                                     const Real & data2,
                                                     const Real3 & pos,
                                                     const Real & start,
                                                     const Real3 & dx,
                                                     const UINT & X,
                                                     const UINT & Y,
                                                     const UINT & Z,
                                                     const uint3 & voxel
) {
    Real diffZ;
    diffZ = (pos.z - (start + Z * dx.z)) / dx.z;
    Real c = ((1 - diffZ) * data1 + (diffZ) * data2);
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
                                          const bool enable2D,
                                          const Real3 kVector) {
  UINT threadID = threadIdx.x + blockIdx.x * blockDim.x;
  const UINT totalSize = voxel.x * voxel.y;
  if (threadID >= totalSize) {
    return;
  }
  Real val, start;
  Real3 dx, pos;
  start = -static_cast<Real>(M_PI / physSize);


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

  const Real & kx = k*kVector.x;
  const Real & ky = k*kVector.y;
  const Real & kz = k*kVector.z;
  val = k * k - (kx + pos.x) * (kx + pos.x) - (ky + pos.y ) * (ky + pos.y);

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
        projection[threadID] += scatter3D[reshape3Dto1D(X,Y,0,voxel)];
      } else {
        UINT Z = static_cast<UINT >(((pos.z - start) / (dx.z)));
        projection[threadID] += computeTrilinearInterpolation(scatter3D, pos, start, dx, X, Y, Z, voxel);
      }
    }
  }
}

__global__ void computeEwaldProjectionGPU(Real *projection,
                                          const Complex *polarizationX,
                                          const Complex *polarizationY,
                                          const Complex *polarizationZ,
                                          const uint3 voxel,
                                          const Real k,
                                          const Real eAngle,
                                          const Real kAngle,
                                          const Real physSize,
                                          const Interpolation::EwaldsInterpolation interpolation,
                                          const bool enable2D,
                                          const Real3 kVector) {
    UINT threadID = threadIdx.x + blockIdx.x * blockDim.x;
    const UINT totalSize = voxel.x * voxel.y;
    if (threadID >= totalSize) {
        return;
    }
    Real val, start;
    Real3 dx, pos;
    start = -static_cast<Real>(M_PI / physSize);


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
    const Real & kx = k*kVector.x;
    const Real & ky = k*kVector.y;
    const Real & kz = k*kVector.z;

    val = k * k - (kx + pos.x) * (kx + pos.x) - (ky + pos.y ) * (ky + pos.y);

    if((val < 0) or (X == (voxel.x - 1)) or (Y == (voxel.y - 1))) {
        projection[threadID] = NAN;
    }
    else
    {
        pos.z = -kz + sqrt(val);
        if (interpolation == Interpolation::EwaldsInterpolation::NEARESTNEIGHBOUR) {
            BigUINT id = computeEquivalentID(pos, X, Y, start, dx, voxel, enable2D);
            projection[threadID] += computeScatter3D(polarizationX,polarizationY,polarizationZ,k,eAngle,kAngle,dx,physSize,id,voxel,enable2D,kVector);
        }
        else {
            if (enable2D) {
                BigUINT id = reshape3Dto1D(X,Y,0,voxel);
                projection[threadID] += computeScatter3D(polarizationX,polarizationY,polarizationZ,k,eAngle,kAngle,dx,physSize,id,voxel,enable2D,kVector);

            } else {
                UINT Z = static_cast<UINT >(((pos.z - start) / (dx.z)));
                if(Z > (voxel.z - 1)){
                    projection[threadID] = NAN;
                }
                else {
                    BigUINT id1 = reshape3Dto1D(X, Y, Z, voxel);
                    BigUINT id2 = reshape3Dto1D(X, Y, Z + 1, voxel);

                    Real data1 = computeScatter3D(polarizationX, polarizationY, polarizationZ, k, eAngle, kAngle, dx,
                                                  physSize, id1, voxel, enable2D,kVector);
                    Real data2 = computeScatter3D(polarizationX, polarizationY, polarizationZ, k, eAngle, kAngle, dx,
                                                  physSize, id2, voxel, enable2D,kVector);

                    projection[threadID] += computeTrilinearInterpolation(data1, data2, pos, start, dx, X, Y, Z, voxel);
                }

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
