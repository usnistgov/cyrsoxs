/////////////////////////////////////////////////////////////////////////////////
// NIST License
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

#pragma once

#include <cudaHeaders.h>
#include <Datatypes.h>
#include <Input/Input.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <Input/InputData.h>
#include <Rotation.h>
#include <cufft.h>
#include <RotationMatrix.h>
#ifdef DOUBLE_PRECISION
static constexpr cufftType_t fftType = CUFFT_Z2Z;
#else
static constexpr cufftType_t fftType = CUFFT_C2C;
#endif

/**
 * @brief calls for the cuda kernel are made through this function for Algorithm 1. Do not uses streams
 * @param [in] voxel array of size 3 which states the dimension along each axis
 * @param [in] idata inputData object
 * @param [in] materialInput material Input containing the information of material property
 * @param [out] projectionAverage I(q) projected on Ewalds sphere
 * @param [in] voxelInput  voxel input
 * @param [in] rotationMatrix rotation matrices for k / E vector
 * @return EXIT_SUCCESS on success of execution
 */
int cudaMain(const UINT *voxel, const InputData &idata, const std::vector<Material> &materialInput,
             Real *projectionAverage, RotationMatrix & rotationMatrix, const Voxel *voxelInput);


/**
 * @brief calls for the cuda kernel are made through this function for Algorithm 2. It uses streams
 * @param [in] voxel array of size 3 which states the dimension along each axis
 * @param [in] idata inputData object
 * @param [in] materialInput material Input containing the information of material property
 * @param [out] projectionAverage I(q) projected on Ewalds sphere
 * @param [in] voxelInput  voxel input
 * @param [in] rotationMatrix rotation matrices for k / E vector
 * @return EXIT_SUCCESS on success of execution
 */
int cudaMainStreams(const UINT *voxel, const InputData &idata, const std::vector<Material> &materialInput,
                    Real *projectionAverage, RotationMatrix & rotationMatrix, const Voxel *voxelInput);

/**
 * @brief calls to compute polarization only. Only called with Pybind interface. Used in debugging
 * @param [in] voxel array of size 3 which states the dimension along each axis
 * @param [in] idata inputData object
 * @param [in] materialInput material Input containing the information of material property
 * @param [in] voxelInput  voxel input
 * @param [in] rotationMatrix rotation matrices for k / E vector
 * @param [out] polarizationX pX
 * @param [out] polarizationY pY
 * @param [out] polarizationZ pZ
 * @param [out] EAngle Angle of E rotation
 * @param [in] NUM_MATERIAL Number of material on Device.
 * @return
 */
int computePolarization(const UINT *voxel, const InputData &idata, const std::vector<Material> &materialInput,
                    Complex *polarizationX,Complex *polarizationY,Complex *polarizationZ,
                    RotationMatrix & rotationMatrix, const Voxel *voxelInput, const Real EAngle, const UINT EnergyID,
                    const int NUM_MATERIAL);


/**
 * This is done to warmup the GPU. The first instruction takes usulally
 * a larger amount of time. This is done to remove the effect in computing
 * total time.
 * @return EXIT_SUCCESS on success of warmup or  EXIT_FAILURE if warmp fails
 */
int warmup();


/**
 * @brief GPU kernel called from CPU for computing the polarization
 * @param [in] d_materialConstants : Material optical constants for a given energy.
 * @param [in] voxelInput : Voxel  property.
 * @param [in] voxel : dimension of morphology
 * @param [out] polarizationX: polarization X vector
 * @param [out] polarizationY: polarization Y vector
 * @param [out] polarizationZ: polarization Z vector
 * @param [in] windowing The windowing type for FFT
 * @param [in] enable2D weather the morphology is 2D
 * @param [in] morphologyType type of morphology Euler / Vector
 * @param [in] rotationMatrix rotationMatrix that rotates E field by a given angle
 * @param [in] numVoxels Number of voxel.
 * @param [in] DEVICE_NUM_MATERIAL Number of material on Device.
 */

template<ReferenceFrame referenceFrame>
__global__ void computePolarization(const Material * d_materialConstants,
                                    const Voxel *voxelInput,
                                    const uint3 voxel,
                                    Complex *polarizationX,
                                    Complex *polarizationY,
                                    Complex *polarizationZ,
                                    FFT::FFTWindowing windowing,
                                    const bool enable2D,
                                    const MorphologyType morphologyType,
                                    const Matrix rotationMatrix,
                                    const BigUINT numVoxels, const int DEVICE_NUM_MATERIAL
);

/**
 * @brief CPU function to compute Polarization for Algorithm 1
 * @param [in] d_materialConstants : Material property.
 * @param [in] d_voxelInput : Voxel  property on GPU.
 * @param [in] voxel : dimension of morphology
 * @param [out] d_polarizationX: device polarization X vector
 * @param [out] d_polarizationY: device polarization Y vector
 * @param [out] d_polarizationZ: device polarization Z vector
 * @param [in] windowing The windowing type for FFT
 * @param [in] enable2D weather the morphology is 2D
 * @param [in] morphologyType type of morphology Euler / Vector
 * @param [in] blockSize blocksize for GPU
 * @param [in] referenceFrame reference frame where the P is calculated: LAB/MATERIAL
 * @param [in] rotationMatrix rotationMatrix that rotates E field by a given angle
 * @param [in] numVoxel Number of voxel.
 * @param [in] NUM_MATERIAL Number of material on Device.
 */
__host__ int computePolarization(const Material * d_materialConstants,
                                  const Voxel *d_voxelInput,
                                  const uint3 & voxel,
                                  Complex *d_polarizationX,
                                  Complex *d_polarizationY,
                                  Complex *d_polarizationZ,
                                  const FFT::FFTWindowing & windowing,
                                  const bool & enable2D,
                                  const MorphologyType & morphologyType,
                                  const UINT & blockSize,
                                  const ReferenceFrame & referenceFrame,
                                  const Matrix & rotationMatrix,
                                  const BigUINT & numVoxel,
                                  const int NUM_MATERIAL
);
/**
 * @brief GPU kernel called from CPU for computing the polarization using Algorithm 2
 * @param [in] d_Nt : Nt for a given energy.
 * @param [out] d_pX: polarization X vector
 * @param [out] d_pY: polarization Y vector
 * @param [out] d_pZ: polarization Z vector
 * @param [in] blockSize blocksize for GPU
 * @param [in] referenceFrame reference frame where the P is calculated: LAB/MATERIAL
 * @param [in] rotationMatrix rotationMatrix that rotates E field by a given angle
 * @param [in] numVoxels Number of voxel.
 */
__host__ int computePolarization(const Complex *d_Nt, Complex *d_pX,
                                 Complex *d_pY, Complex *d_pZ,
                                 const UINT &blockSize,
                                 const ReferenceFrame &referenceFrame,
                                 const Matrix &rotationMatrix,const BigUINT &numVoxels);

/**
 * @brief computes FFT on GPU in place
 * @param [in,out] polarization px/py/pz
 * @param [in] plan The plan context
 * @return
 */
__host__ INLINE inline cufftResult  performFFT(Complex *polarization, cufftHandle &plan) {
#ifdef DOUBLE_PRECISION
    return (cufftExecZ2Z(plan, polarization, polarization, CUFFT_FORWARD));
#else
  return (cufftExecC2C(plan, polarization, polarization, CUFFT_FORWARD));
#endif
}

/**
 * @brief performs in place FFT Shift
 * @param [in,out] polarization px/py/pz
 * @param blockSize blocksize of GPU kernel
 * @param vx voxel dims in all direction
 * @param stream cuda stream context
 * @return
 */
__host__ int performFFTShift(Complex *polarization, const UINT & blockSize, const uint3 & vx,  const cudaStream_t stream);

/**
 *
 * @param d_polarizationX device polarization X vector
 * @param d_polarizationY device polarization Y vector
 * @param d_polarizationZ device polarization Z vector
 * @param d_scatter3D  scatter 3D computation
 * @param kMagnitude k magnitude
 * @param voxelSize voxel size
 * @param vx voxel dimensions
 * @param physSize physSize
 * @param enable2D true for 2D morphology
 * @param blockSize blockSize
 * @param kVector kVector
 * @return
 */
__host__ int performScatter3DComputation(const Complex * d_polarizationX, const Complex *d_polarizationY, const Complex * d_polarizationZ,
                                          Real * d_scatter3D,
                                          const Real & kMagnitude,
                                          const BigUINT & voxelSize,
                                          const uint3 & vx,
                                          const Real & physSize,
                                          const bool & enable2D,
                                          const UINT & blockSize,
                                          const Real3 & kVector);

/**
 *
 * @param d_projection
 * @param d_scatter
 * @param kMagnitude
 * @param vx
 * @param physSize
 * @param interpolation
 * @param enable2D
 * @param blockSize
 * @param kVector
 * @return
 */
__host__ int peformEwaldProjectionGPU(Real * d_projection,
                                      const Real * d_scatter,
                                      const Real & kMagnitude,
                                      const uint3 & vx,
                                      const Real & physSize,
                                      const Interpolation::EwaldsInterpolation & interpolation,
                                      const bool & enable2D,
                                      const UINT & blockSize,
                                      const Real3 & kVector);

/**
 *
 * @param projection
 * @param d_polarizationX
 * @param d_polarizationY
 * @param d_polarizationZ
 * @param kMagnitude
 * @param vx
 * @param physSize
 * @param interpolation
 * @param enable2D
 * @param blockSize
 * @param kVector
 * @return
 */
__host__ int peformEwaldProjectionGPU(Real * projection,
                                      const Complex * d_polarizationX,
                                      const Complex * d_polarizationY,
                                      const Complex * d_polarizationZ,
                                      const Real & kMagnitude,
                                      const uint3 & vx,
                                      const Real & physSize,
                                      const Interpolation::EwaldsInterpolation & interpolation,
                                      const bool & enable2D,
                                      const UINT & blockSize,
                                      const Real3 & kVector);


/**
 * @brief CPU function to compute Nt for Algorithm 2
 * @param [in] d_materialConstants : Material property.
 * @param [in] d_voxelInput : Voxel  property on GPU.
 * @param [out] d_Nt :  computes Nt = (NR:NR - I)
 * @param [in] morphologyType type of morphology Euler / Vector
 * @param [in] blockSize  blocksize for GPU
 * @param [in] numVoxels Number of voxel.
 * @param [in] offset offset in voxels according to streams
 * @param [in] endID finish ID for this particular stream
 * @param [in] materialID  the material ID for the current loop
 * @param [in] numStreams number of streams
 * @param [in] stream stream context
 * @param [in] NUM_MATERIAL Number of material on Device.
 * @return EXIT_SUCEESS if completed
 */
__host__ int computeNt(const Material * d_materialConstants,
                       const Voxel *d_voxelInput,
                       Complex * d_Nt,
                       const MorphologyType &morphologyType,
                       const UINT &blockSize,
                       const BigUINT & numVoxels,
                       const BigUINT & offset,
                       const BigUINT & endID,
                       const UINT & materialID,
                       const UINT & numStreams,
                       cudaStream_t stream,
                       const int NUM_MATERIAL);



