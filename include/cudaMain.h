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
 * @brief calls for the cuda kernel are made through this function.
 * @param voxel array of size 3 which states the dimension along each axis
 * @param idata inputData object
 * @param materialInput material Input containing the information of material property
 * @param voxelInput voxelInput
 * @return EXIT_SUCCESS on success of execution
 */
int cudaMain(const UINT *voxel, const InputData &idata, const std::vector<Material<NUM_MATERIAL> > &materialInput,
             Real *projectionAverage, RotationMatrix & rotationMatrix, const Voxel<NUM_MATERIAL> *voxelInput);

int cudaMainStreams(const UINT *voxel, const InputData &idata, const std::vector<Material<NUM_MATERIAL> > &materialInput,
                    Real *projectionAverage, RotationMatrix & rotationMatrix, const Voxel<NUM_MATERIAL> *voxelInput);


/**
 * This is done to warmup the GPU. The first instruction takes usulally
 * a larger amount of time. This is done to remove the effect in computing
 * total time.
 * @return EXIT_SUCCESS on success of warmup or  EXIT_FAILURE if warmp fails
 */
int warmup();


/**
 * @brief GPU kernel called from CPU for computing the polarization
 * @param [in] materialInput : Material property.
 * @param [in] voxelInput : Voxel  property.
 * @param [in] elefield: Electric field.
 * @param [in] voxelNum Number of voxel.
 * @param [in] angle Angle of rotation.
 * @param [out] polarizationX: polarization X vector
 * @param [out] polarizationY: polarization Y vector
 * @param [out] polarizationZ: polarization Z vector
 */
template<ReferenceFrame referenceFrame>
__global__ void computePolarization(Material<NUM_MATERIAL> materialInput,
                                    const Voxel<NUM_MATERIAL> *voxelInput,
                                    const ElectricField elefield,
                                    const Real angle,
                                    const uint3 voxel,
                                    Complex *polarizationX,
                                    Complex *polarizationY,
                                    Complex *polarizationZ,
                                    FFT::FFTWindowing windowing,
                                    const bool enable2D,
                                    const MorphologyType morphologyType
);




__host__ int computePolarization(const Material<NUM_MATERIAL> & materialInput,
                                  const Voxel<NUM_MATERIAL> *voxelInput,
                                  const ElectricField  & elefield,
                                  const Real & angle,
                                  const uint3 & voxel,
                                  Complex *d_polarizationX,
                                  Complex *d_polarizationY,
                                  Complex *d_polarizationZ,
                                  FFT::FFTWindowing windowing,
                                  const bool & enable2D,
                                  const MorphologyType & morphologyType,
                                  const UINT & blockSize,
                                  const ReferenceFrame & referenceFrame,
                                  const Matrix & rotationMatrix
);


__host__ INLINE inline cufftResult  performFFT(Complex *polarization, cufftHandle &plan) {
    return (cufftC2C(plan, polarization, polarization, CUFFT_FORWARD));
}

__host__ int performFFTShift(Complex *polarization, const UINT & blockSize, const uint3 & vx);

__host__ int performScatter3DComputation(const Complex * d_polarizationX, const Complex *d_polarizationY, const Complex * d_polarizationZ,
                                          Real * d_scatter3D,
                                          const ElectricField & eleField,
                                          const Real & eAngle,
                                          const Real & kAngle,
                                          const BigUINT & voxelSize,
                                          const uint3 & vx,
                                          const Real & physSize,
                                          const bool & enable2D,
                                          const UINT & blockSize,
                                          const Real3 & kVector);

__host__ int peformEwaldProjectionGPU(Real * d_projection,
                                      const Real * d_scatter,
                                      const Real & k,
                                      const uint3 & vx,
                                      const Real & eAngle,
                                      const Real & kAngle,
                                      const Real & physSize,
                                      const Interpolation::EwaldsInterpolation & interpolation,
                                      const bool & enable2D,
                                      const UINT & blockSize,
                                      const Real3 & kVector);

__host__ int peformEwaldProjectionGPU(Real * projection,
                                      const Complex * d_polarizationX,
                                      const Complex * d_polarizationY,
                                      const Complex * d_polarizationZ,
                                      const Real & k,
                                      const uint3 & vx,
                                      const Real & eAngle,
                                      const Real & kAngle,
                                      const Real & physSize,
                                      const Interpolation::EwaldsInterpolation & interpolation,
                                      const bool & enable2D,
                                      const UINT & blockSize,
                                      const Real3 & kVector);

__host__ int computeNt(const Material<NUM_MATERIAL> &materialInput,
                       const Voxel<NUM_MATERIAL> *d_voxelInput,
                       Complex * d_Nt,
                       const MorphologyType &morphologyType,
                       const UINT &blockSize);


__host__ int computePolarization(const Complex *d_Nt, Complex *d_pX,
                                 Complex *d_pY, Complex *d_pZ,
                                 const UINT &blockSize,
                                 const ReferenceFrame &referenceFrame,
                                 const Matrix &rotationMatrix);