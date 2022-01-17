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
 * @brief This function computes the polarization in real space for the uniaxial case (Euler angle version).
 * @param [in] material material data for a particular energy level under consideration.
 * @param [in] voxelInput voxel Input data.
 * @param [in] threadID threadID to access the index of the entry.
 * @param [out] polarizationX X polarization
 * @param [out] polarizationY Y polarization
 * @param [out] polarizationZ Z polarization
 * @param [in] numVoxels number of voxels
 * @param [in] rotationMatrix rotation matrix for given k/E
 * @param [in] NUM_MATERIAL number of material
 */
template<ReferenceFrame referenceFrame>
__device__ void computePolarizationEulerAngles(const Material *material,
                                            const Voxel *voxelInput, const BigUINT threadID,
                                            Complex *polarizationX, Complex *polarizationY, Complex *polarizationZ,
                                            const BigUINT & numVoxels, const Matrix & rotationMatrix, int NUM_MATERIAL) {

  Complex pX{0.0,0.0}, pY{0.0,0.0}, pZ{0.0,0.0};

  static constexpr Real OneBy4Pi = static_cast<Real> (1.0 / (4.0 * M_PI));
  /**
 * [0 1 2]
 * [1 3 4]
 * [2 4 5]
 */
  Complex rotatedNr{0.0,0.0}; // Only storing what is required
  static constexpr Real3 eleField{1,0,0};
  Real3 matVec;
  doMatVec<false>(rotationMatrix,eleField,matVec);
  Complex nsum;
  for (int numMaterial = 0; numMaterial < NUM_MATERIAL; numMaterial++) {
    Complex npar = material[numMaterial].npara;
    Complex nper = material[numMaterial].nperp;
    const Voxel matProp = voxelInput[numVoxels * numMaterial + threadID];

    const Real & psiAngle      = matProp.getValueAt(Voxel::EULER_ANGLE::PSI);
    const Real & thetaAngle    = matProp.getValueAt(Voxel::EULER_ANGLE::THETA);
    const Real & Vfrac         = matProp.getValueAt(Voxel::EULER_ANGLE::VFRAC);
    const Real & phi_a         = Vfrac*matProp.getValueAt(Voxel::EULER_ANGLE::S);

    const Real   sx     =  cos(psiAngle)*sin(thetaAngle);
    const Real   sy     =  sin(psiAngle)*sin(thetaAngle);
    const Real   sz     =  cos(thetaAngle);
    const Real   phi_ui = Vfrac - phi_a;

    const Real  & phi =  Vfrac;

    nsum.x = npar.x + 2 * nper.x;
    nsum.y = npar.y + 2 * nper.y;

    computeComplexSquare(nsum);
    computeComplexSquare(npar);
    computeComplexSquare(nper);

    // S = fraction of aligned components
    // (0)
    rotatedNr.x = phi_a*(npar.x*sx*sx + nper.x*(sy*sy + sz*sz)) + ((phi_ui * nsum.x) / (Real) 9.0) - phi;
    rotatedNr.y = phi_a*(npar.y*sx*sx + nper.y*(sy*sy + sz*sz)) + ((phi_ui * nsum.y) / (Real) 9.0);

    pX.x += rotatedNr.x*matVec.x;
    pX.y += rotatedNr.y*matVec.x;

    // (1)
    rotatedNr.x = phi_a*(npar.x - nper.x)*sx*sy;
    rotatedNr.y = phi_a*(npar.y - nper.y)*sx*sy;

    pX.x += rotatedNr.x*matVec.y;
    pX.y += rotatedNr.y*matVec.y;

    pY.x += rotatedNr.x*matVec.x;
    pY.y += rotatedNr.y*matVec.x;

    // (2)
    rotatedNr.x = phi_a*(npar.x - nper.x)*sx*sz;
    rotatedNr.y = phi_a*(npar.y - nper.y)*sx*sz;

    pX.x += rotatedNr.x*matVec.z;
    pX.y += rotatedNr.y*matVec.z;

    pZ.x += rotatedNr.x*matVec.x;
    pZ.y += rotatedNr.y*matVec.x;

    // (3)
    rotatedNr.x = phi_a*(npar.x*sy*sy + nper.x*(sx*sx + sz*sz)) + ((phi_ui * nsum.x) / (Real) 9.0) - phi;
    rotatedNr.y = phi_a*(npar.y*sy*sy + nper.y*(sx*sx + sz*sz)) + ((phi_ui * nsum.y) / (Real) 9.0);

    pY.x += rotatedNr.x*matVec.y;
    pY.y += rotatedNr.y*matVec.y;

    // (4)
    rotatedNr.x = phi_a*(npar.x - nper.x)*sy*sz;
    rotatedNr.y = phi_a*(npar.y - nper.y)*sy*sz;

    pY.x += rotatedNr.x*matVec.z;
    pY.y += rotatedNr.y*matVec.z;

    pZ.x += rotatedNr.x*matVec.y;
    pZ.y += rotatedNr.y*matVec.y;

    // (5)
    rotatedNr.x = phi_a*(npar.x*sz*sz + nper.x*(sx*sx + sy*sy)) +  ((phi_ui * nsum.x) / (Real) 9.0) - phi;
    rotatedNr.y = phi_a*(npar.y*sz*sz + nper.y*(sx*sx + sy*sy)) +  ((phi_ui * nsum.y) / (Real) 9.0);

    pZ.x += rotatedNr.x*matVec.z;
    pZ.y += rotatedNr.y*matVec.z;

  }

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
 * @brief This function computes the polarization in real space for the uniaxial case. (Vector Morphology)
 * @param [in] material material data for a particular energy level under consideration.
 * @param [in] voxelInput voxel Input data.
 * @param [in] threadID threadID to access the index of the entry.
 * @param [out] polarizationX X polarization
 * @param [out] polarizationY Y polarization
 * @param [out] polarizationZ Z polarization
 * @param [in] numVoxels number of voxels
 * @param [in] rotationMatrix rotation matrix for given k/E
 * @param [in] NUM_MATERIAL number of material
 */

template<ReferenceFrame referenceFrame>
__device__ void computePolarizationVectorMorphologyOptimized(const Material *material,
                                                    const Voxel *voxelInput, const BigUINT & threadID,
                                                    Complex *polarizationX, Complex *polarizationY, Complex *polarizationZ,
                                                    const BigUINT & numVoxels, const Matrix & rotationMatrix, int NUM_MATERIAL) {

  Complex pX{0.0,0.0}, pY{0.0,0.0}, pZ{0.0,0.0};

  static constexpr Real OneBy4Pi = static_cast<Real> (1.0 / (4.0 * M_PI));
  /**
 * [0 1 2]
 * [1 3 4]
 * [2 4 5]
 */
  Complex rotatedNr{0.0,0.0}; // Only storing what is required
  static constexpr Real3 eleField{1,0,0};
  Real3 matVec;
  doMatVec<false>(rotationMatrix,eleField,matVec);
  Complex nsum;
  for (int numMaterial = 0; numMaterial < NUM_MATERIAL; numMaterial++) {
    Complex npar = material[numMaterial].npara;
    Complex nper = material[numMaterial].nperp;
    const Voxel matProp = voxelInput[numVoxels * numMaterial + threadID];
    const Real & sx     = matProp.s1.x;
    const Real & sy     = matProp.s1.y;
    const Real & sz     = matProp.s1.z;
    const Real & phi_ui = matProp.s1.w;

    /**
     * According to Eliot Dated June 29,2021:
     * In old morphology generator: Vfrac = |s|^2 + phi_ui
     */
    const Real  phi = phi_ui + sx * sx + sy * sy + sz * sz;

    nsum.x = npar.x + 2 * nper.x;
    nsum.y = npar.y + 2 * nper.y;

    computeComplexSquare(nsum);
    computeComplexSquare(npar);
    computeComplexSquare(nper);

    // (0)
    rotatedNr.x = npar.x*sx*sx + nper.x*(sy*sy + sz*sz) + ((phi_ui * nsum.x) / (Real) 9.0) - phi;
    rotatedNr.y = npar.y*sx*sx + nper.y*(sy*sy + sz*sz) + ((phi_ui * nsum.y) / (Real) 9.0);

    pX.x += rotatedNr.x*matVec.x;
    pX.y += rotatedNr.y*matVec.x;

    // (1)
    rotatedNr.x = (npar.x - nper.x)*sx*sy;
    rotatedNr.y = (npar.y - nper.y)*sx*sy;

    pX.x += rotatedNr.x*matVec.y;
    pX.y += rotatedNr.y*matVec.y;

    pY.x += rotatedNr.x*matVec.x;
    pY.y += rotatedNr.y*matVec.x;

    // (2)
    rotatedNr.x = (npar.x - nper.x)*sx*sz;
    rotatedNr.y = (npar.y - nper.y)*sx*sz;

    pX.x += rotatedNr.x*matVec.z;
    pX.y += rotatedNr.y*matVec.z;

    pZ.x += rotatedNr.x*matVec.x;
    pZ.y += rotatedNr.y*matVec.x;

    // (3)
    rotatedNr.x = npar.x*sy*sy + nper.x*(sx*sx + sz*sz) + ((phi_ui * nsum.x) / (Real) 9.0) - phi;
    rotatedNr.y = npar.y*sy*sy + nper.y*(sx*sx + sz*sz) + ((phi_ui * nsum.y) / (Real) 9.0);

    pY.x += rotatedNr.x*matVec.y;
    pY.y += rotatedNr.y*matVec.y;

    // (4)
    rotatedNr.x = (npar.x - nper.x)*sy*sz;
    rotatedNr.y = (npar.y - nper.y)*sy*sz;

    pY.x += rotatedNr.x*matVec.z;
    pY.y += rotatedNr.y*matVec.z;

    pZ.x += rotatedNr.x*matVec.y;
    pZ.y += rotatedNr.y*matVec.y;

    // (5)
    rotatedNr.x = npar.x*sz*sz + nper.x*(sx*sx + sy*sy) +  ((phi_ui * nsum.x) / (Real) 9.0) - phi;
    rotatedNr.y = npar.y*sz*sz + nper.y*(sx*sx + sy*sy) +  ((phi_ui * nsum.y) / (Real) 9.0);

    pZ.x += rotatedNr.x*matVec.z;
    pZ.y += rotatedNr.y*matVec.z;
  }

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
 * @brief computes Nt for Algorithm 2 for Vector morphology
 * @param [in] material refractive index of the material
 * @param [in] voxelInput voxel data
 * @param [out] Nt  computes Nt = (NR:NR - I)
 * @param [in] offset offset in voxels according to streams
 * @param [in] endID finish ID for this particular stream
 * @param [in] materialID  the material ID for the current loop
 * @param [in] numVoxels number of voxels
 * @param [in] NUM_MATERIAL number of materials
 */

__global__ void computeNtVectorMorphology(const Material * materialConstants,
                          const Voxel * __restrict__ voxelInput,
                          Complex * Nt, const BigUINT  offset, const BigUINT  endID, const UINT materialID,
                          const BigUINT numVoxels, int NUM_MATERIAL) {
  const BigUINT threadID = threadIdx.x + blockIdx.x * blockDim.x;
  if((threadID + offset) >= endID){
    return;
  }
  Complex rotatedNr[6]; // Only storing what is required
  memset(rotatedNr,0,sizeof(Complex)*6);
  Complex nsum;
  Complex npar = materialConstants[materialID].npara;
  Complex nper = materialConstants[materialID].nperp;
  const Real4 matProp = voxelInput[offset + threadID].s1;
  const Real &sx = matProp.x;
  const Real &sy = matProp.y;
  const Real &sz = matProp.z;
  const Real &phi_ui = matProp.w;

  const Real &phi = phi_ui + sx * sx + sy * sy + sz * sz;

  nsum.x = npar.x + 2 * nper.x;
  nsum.y = npar.y + 2 * nper.y;

  computeComplexSquare(nsum);
  computeComplexSquare(npar);
  computeComplexSquare(nper);

  rotatedNr[0].x += npar.x * sx * sx + nper.x * (sy * sy + sz * sz) + ((phi_ui * nsum.x) / (Real) 9.0) - phi;
  rotatedNr[0].y += npar.y * sx * sx + nper.y * (sy * sy + sz * sz) + ((phi_ui * nsum.y) / (Real) 9.0);

  rotatedNr[1].x += (npar.x - nper.x) * sx * sy;
  rotatedNr[1].y += (npar.y - nper.y) * sx * sy;

  rotatedNr[2].x += (npar.x - nper.x) * sx * sz;
  rotatedNr[2].y += (npar.y - nper.y) * sx * sz;

  rotatedNr[3].x += npar.x * sy * sy + nper.x * (sx * sx + sz * sz) + ((phi_ui * nsum.x) / (Real) 9.0) - phi;
  rotatedNr[3].y += npar.y * sy * sy + nper.y * (sx * sx + sz * sz) + ((phi_ui * nsum.y) / (Real) 9.0);

  rotatedNr[4].x += (npar.x - nper.x) * sy * sz;
  rotatedNr[4].y += (npar.y - nper.y) * sy * sz;

  rotatedNr[5].x += npar.x * sz * sz + nper.x * (sx * sx + sy * sy) + ((phi_ui * nsum.x) / (Real) 9.0) - phi;
  rotatedNr[5].y += npar.y * sz * sz + nper.y * (sx * sx + sy * sy) + ((phi_ui * nsum.y) / (Real) 9.0);


  Nt[2*(threadID+offset) + 0+  0*numVoxels].x += rotatedNr[0].x; Nt[2*(threadID+offset) + 0  +  0*numVoxels ].y += rotatedNr[0].y;
  Nt[2*(threadID+offset) + 1 + 0*numVoxels].x += rotatedNr[1].x; Nt[2*(threadID+offset) + 1  +  0*numVoxels ].y += rotatedNr[1].y;
  Nt[2*(threadID+offset) + 0 + 2*numVoxels].x += rotatedNr[2].x; Nt[2*(threadID+offset) + 0  +  2*numVoxels ].y += rotatedNr[2].y;
  Nt[2*(threadID+offset) + 1 + 2*numVoxels].x += rotatedNr[3].x; Nt[2*(threadID+offset) + 1  +  2*numVoxels ].y += rotatedNr[3].y;
  Nt[2*(threadID+offset) + 0 + 4*numVoxels].x += rotatedNr[4].x; Nt[2*(threadID+offset) + 0  +  4*numVoxels ].y += rotatedNr[4].y;
  Nt[2*(threadID+offset) + 1 + 4*numVoxels].x += rotatedNr[5].x; Nt[2*(threadID+offset) + 1  +  4*numVoxels ].y += rotatedNr[5].y;
}

/**
 * @brief computes Nt for Algorithm 2 using Euler angles
 * @param [in] material refractive index of the material
 * @param [in] voxelInput voxel data
 * @param [out] Nt  computes Nt = (NR:NR - I)
 * @param offset offset in voxels according to streams
 * @param endID finish ID for this particular stream
 * @param [in] materialID  the material ID for the current loop
 * @param numVoxels number of voxels
 * @param NUM_MATERIAL number of materials
 */
__global__ void computeNtEulerAngles(const Material  * materialConstants,
                                          const Voxel * __restrict__ voxelInput,
                                          Complex * Nt, const BigUINT  offset, const BigUINT  endID, const UINT materialID,
                                          const BigUINT numVoxels, int NUM_MATERIAL) {

  const BigUINT threadID = threadIdx.x + blockIdx.x * blockDim.x;
  if((threadID + offset) >= endID){
    return;
  }
  Complex rotatedNr[6]; // Only storing what is required
  memset(rotatedNr,0,sizeof(Complex)*6);
  Complex nsum;
  Complex npar = materialConstants[materialID].npara;
  Complex nper = materialConstants[materialID].nperp;
  const Voxel matProp = voxelInput[offset + threadID];
  const Real & psiAngle      = matProp.getValueAt(Voxel::EULER_ANGLE::PSI);
  const Real & thetaAngle    = matProp.getValueAt(Voxel::EULER_ANGLE::THETA);
  const Real & Vfrac         = matProp.getValueAt(Voxel::EULER_ANGLE::VFRAC);
  const Real & phi_a             = Vfrac*matProp.getValueAt(Voxel::EULER_ANGLE::S); //  S = fraction of voxel with aligned component

  const Real   sx     =  cos(psiAngle)*sin(thetaAngle);
  const Real   sy     =  sin(psiAngle)*sin(thetaAngle);
  const Real   sz     =  cos(thetaAngle);
  const Real   phi_ui =  Vfrac - phi_a;

  const Real  & phi =  Vfrac;

  nsum.x = npar.x + 2 * nper.x;
  nsum.y = npar.y + 2 * nper.y;

  computeComplexSquare(nsum);
  computeComplexSquare(npar);
  computeComplexSquare(nper);

  rotatedNr[0].x += phi_a*(npar.x * sx * sx + nper.x * (sy * sy + sz * sz)) + ((phi_ui * nsum.x) / (Real) 9.0) - phi;
  rotatedNr[0].y += phi_a*(npar.y * sx * sx + nper.y * (sy * sy + sz * sz)) + ((phi_ui * nsum.y) / (Real) 9.0);

  rotatedNr[1].x += phi_a*(npar.x - nper.x) * sx * sy;
  rotatedNr[1].y += phi_a*(npar.y - nper.y) * sx * sy;

  rotatedNr[2].x += phi_a*(npar.x - nper.x) * sx * sz;
  rotatedNr[2].y += phi_a*(npar.y - nper.y) * sx * sz;

  rotatedNr[3].x += phi_a*(npar.x * sy * sy + nper.x * (sx * sx + sz * sz)) + ((phi_ui * nsum.x) / (Real) 9.0) - phi;
  rotatedNr[3].y += phi_a*(npar.y * sy * sy + nper.y * (sx * sx + sz * sz)) + ((phi_ui * nsum.y) / (Real) 9.0);

  rotatedNr[4].x += phi_a*(npar.x - nper.x) * sy * sz;
  rotatedNr[4].y += phi_a*(npar.y - nper.y) * sy * sz;

  rotatedNr[5].x += phi_a*(npar.x * sz * sz + nper.x * (sx * sx + sy * sy)) + ((phi_ui * nsum.x) / (Real) 9.0) - phi;
  rotatedNr[5].y += phi_a*(npar.y * sz * sz + nper.y * (sx * sx + sy * sy)) + ((phi_ui * nsum.y) / (Real) 9.0);


  Nt[2*(threadID+offset) + 0+  0*numVoxels].x += rotatedNr[0].x; Nt[2*(threadID+offset) + 0  +  0*numVoxels ].y += rotatedNr[0].y;
  Nt[2*(threadID+offset) + 1 + 0*numVoxels].x += rotatedNr[1].x; Nt[2*(threadID+offset) + 1  +  0*numVoxels ].y += rotatedNr[1].y;
  Nt[2*(threadID+offset) + 0 + 2*numVoxels].x += rotatedNr[2].x; Nt[2*(threadID+offset) + 0  +  2*numVoxels ].y += rotatedNr[2].y;
  Nt[2*(threadID+offset) + 1 + 2*numVoxels].x += rotatedNr[3].x; Nt[2*(threadID+offset) + 1  +  2*numVoxels ].y += rotatedNr[3].y;
  Nt[2*(threadID+offset) + 0 + 4*numVoxels].x += rotatedNr[4].x; Nt[2*(threadID+offset) + 0  +  4*numVoxels ].y += rotatedNr[4].y;
  Nt[2*(threadID+offset) + 1 + 4*numVoxels].x += rotatedNr[5].x; Nt[2*(threadID+offset) + 1  +  4*numVoxels ].y += rotatedNr[5].y;
}
/**
 * @brief computes polarization by ALgorithm 2
 * @tparam referenceFrame reference frame LAB/MATERIAL
 * @param [in] Nt Nt = (NR:NR - I)
 * @param [out] polarizationX pX
 * @param [out] polarizationY pY
 * @param [out] polarizationZ pZ
 * @param [in] rotationMatrix rotation matrix corresponding to E/k
 * @param [in] numVoxels number of voxels
 */
template<ReferenceFrame referenceFrame>
__global__ void computePolarizationVectorMorphologyLowMemory(const Real4 * __restrict__ Nt,Complex *polarizationX,
                                                             Complex *polarizationY, Complex *polarizationZ,
                                                             const Matrix rotationMatrix, const BigUINT numVoxels) {
  const BigUINT threadID = threadIdx.x + blockIdx.x * blockDim.x;
  if(threadID > numVoxels){
    return;
  }
  Complex pX{0,0}, pY{0,0}, pZ{0,0};

  /**
 * [0 1 2]
 * [1 3 4]
 * [2 4 5]
 */
  static constexpr Real OneBy4Pi = static_cast<Real> (1.0 / (4.0 * M_PI));


  static constexpr Real3 eleField{1,0,0};
  Real3 matVec;
  doMatVec<false>(rotationMatrix,eleField,matVec);
  Real4 _rotatedNr;
  _rotatedNr = Nt[threadID + 0*numVoxels]; // (0) and (1)
  Complex  rotatedNr[2];
  rotatedNr[0] = {_rotatedNr.x,_rotatedNr.y};
  rotatedNr[1] = {_rotatedNr.z,_rotatedNr.w};

  // (0)
  pX.x += rotatedNr[0].x*matVec.x;
  pX.y += rotatedNr[0].y*matVec.x;

 // (1)
  pX.x += rotatedNr[1].x*matVec.y;
  pX.y += rotatedNr[1].y*matVec.y;

  pY.x += rotatedNr[1].x*matVec.x;
  pY.y += rotatedNr[1].y*matVec.x;

  _rotatedNr = Nt[threadID + 1*numVoxels]; // (2) and (3)
  rotatedNr[0] = {_rotatedNr.x,_rotatedNr.y};
  rotatedNr[1] = {_rotatedNr.z,_rotatedNr.w};
  // (2)
  pX.x += rotatedNr[0].x*matVec.z;
  pX.y += rotatedNr[0].y*matVec.z;

  pZ.x += rotatedNr[0].x*matVec.x;
  pZ.y += rotatedNr[0].y*matVec.x;

  //(3)
  pY.x += rotatedNr[1].x*matVec.y;
  pY.y += rotatedNr[1].y*matVec.y;

  _rotatedNr = Nt[threadID + 2*numVoxels]; // (4) and (5)
  rotatedNr[0] = {_rotatedNr.x,_rotatedNr.y};
  rotatedNr[1] = {_rotatedNr.z,_rotatedNr.w};
  // (4)
  pY.x += rotatedNr[0].x*matVec.z;
  pY.y += rotatedNr[0].y*matVec.z;

  pZ.x += rotatedNr[0].x*matVec.y;
  pZ.y += rotatedNr[0].y*matVec.y;
  // (5)
  pZ.x += rotatedNr[1].x*matVec.z;
  pZ.y += rotatedNr[1].y*matVec.z;

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
 * @brief flattens 3D array to 1D array
 * @param [in] i X id
 * @param [in] j Y id
 * @param [in] k Z id
 * @param [in] voxel voxel dimensions
 * @return flattened id
 */
__host__ __device__ inline BigUINT reshape3Dto1D(UINT i, UINT j, UINT k, uint3 voxel) {
  return ((i) + (j) * voxel.x + (k) * voxel.x * voxel.y);
}

/**
 * @brief returns 3D index corresponding to 1D array
 * @param [in] id 1D flattened id
 * @param [out] X X id
 * @param [out] Y Y id
 * @param [out] Z Z id
 * @param voxel voxel dimension in each direction
 */
__host__ __device__ inline void reshape1Dto3D(BigUINT id, UINT &X, UINT &Y, UINT &Z, uint3 voxel) {
  Z = static_cast<UINT>(id / (voxel.y * voxel.x * 1.0));
  Y = static_cast<UINT>(id - Z * voxel.y * voxel.x) / (voxel.x * 1.0);
  X = static_cast<UINT>(id - Y * voxel.x - Z * voxel.y * voxel.x);
}

/**
 * @brief swaps two variable
 * @tparam T  template
 * @param [in,out] var1 variable 1
 * @param [in,out] var2 variable 2
 */
template<typename T>
__device__ inline void swap(T &var1, T &var2) {
  T temp;
  temp = var1;
  var1 = var2;
  var2 = temp;
}

/**
 * @brief performs FFT shift. The shift logic is consistent with Igor version
 * @tparam T template
 * @param  [in,out] polarization polarization vector in Fourier space
 * @param  [in] voxel voxel dimensions
 */

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


/**
 * @brief computes X(q) at each voxel
 * @param [in] polarizationX X polarization
 * @param [in] polarizationY Y polarization
 * @param [in] polarizationZ Z polarization
 * @param [in] k magnitude of k vector
 * @param [in] dX spacing in each direction
 * @param [in] physSize physical Size
 * @param [in] id voxel id
 * @param [in] voxel voxel dimensions in each direction
 * @param [in] enable2D whether 2D morphology
 * @param [in] kVector 3D k vector
 * @return X(q) for a given voxel
 */
inline __device__ Real computeScatter3D(const Complex * polarizationX,
                                        const Complex * polarizationY,
                                        const Complex * polarizationZ,
                                        const Real & k,
                                        const Real3 & dX,
                                        const Real & physSize,
                                        const BigUINT & id,
                                        const uint3 & voxel,
                                        const bool enable2D,
                                        const Real3 & kVector
){

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
* @param [in] k              magntude of k vector
* @param [in] voxelNum       Number of total voxel
* @param [in] voxel          Number of voxel in each direciton.
* @param [in] physSize       Physical Size
* @param [in] enable2D       2D morphology or not
* @param [in] kVector        3D k Vector
*/
__global__ void computeScatter3D(const Complex *polarizationX,
                                 const Complex *polarizationY,
                                 const Complex *polarizationZ,
                                 Real *Scatter3D,
                                 const Real k,
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

 Scatter3D[threadID] = computeScatter3D(polarizationX,polarizationY,polarizationZ,k,dx,physSize,threadID,voxel,enable2D,kVector);
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
 * @param [in] enable2D 2D morphology or not
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
 * @param [in] interpolation type of interpolation : Nearest neighbor / Trilinear interpolation
 * @param [in] enable2D 2D morpholgy or not
 * @param [in] kVector 3D k vector
 */
__global__ void computeEwaldProjectionGPU(Real *projection,
                                          const Real *scatter3D,
                                          const uint3 voxel,
                                          const Real k,
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
                                          const Real kMagnitude,
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
    const Real & kx = kMagnitude * kVector.x;
    const Real & ky = kMagnitude * kVector.y;
    const Real & kz = kMagnitude * kVector.z;

    val = kMagnitude * kMagnitude - (kx + pos.x) * (kx + pos.x) - (ky + pos.y ) * (ky + pos.y);

    if((val < 0) or (X == (voxel.x - 1)) or (Y == (voxel.y - 1))) {
        projection[threadID] = NAN;
    }
    else
    {
        pos.z = -kz + sqrt(val);
        if (interpolation == Interpolation::EwaldsInterpolation::NEARESTNEIGHBOUR) {
            BigUINT id = computeEquivalentID(pos, X, Y, start, dx, voxel, enable2D);
            projection[threadID] += computeScatter3D(polarizationX, polarizationY, polarizationZ, kMagnitude, dx, physSize, id, voxel, enable2D, kVector);
        }
        else {
            if (enable2D) {
                BigUINT id = reshape3Dto1D(X,Y,0,voxel);
                projection[threadID] += computeScatter3D(polarizationX, polarizationY, polarizationZ, kMagnitude, dx, physSize, id, voxel, enable2D, kVector);

            } else {
                UINT Z = static_cast<UINT >(((pos.z - start) / (dx.z)));
                if(Z > (voxel.z - 1)){
                    projection[threadID] = NAN;
                }
                else {
                    BigUINT id1 = reshape3Dto1D(X, Y, Z, voxel);
                    BigUINT id2 = reshape3Dto1D(X, Y, Z + 1, voxel);

                    Real data1 = computeScatter3D(polarizationX, polarizationY, polarizationZ, kMagnitude, dx,
                                                  physSize, id1, voxel, enable2D, kVector);
                    Real data2 = computeScatter3D(polarizationX, polarizationY, polarizationZ, kMagnitude, dx,
                                                  physSize, id2, voxel, enable2D, kVector);

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
 * @param [in] voxel Number of voxel in each direction
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
