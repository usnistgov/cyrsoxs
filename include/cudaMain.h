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

#include <cufft.h>

/**
 * @brief calls for the cuda kernel are made through this function.
 * @param voxel array of size 3 which states the dimension along each axis
 * @param idata inputData object
 * @param materialInput material Input containing the information of material property
 * @param voxelInput voxelInput
 * @return EXIT_SUCCESS on success of execution
 */
int cudaMain(const UINT *voxel, const InputData &idata, const std::vector<Material<NUM_MATERIAL> > &materialInput,
             Real *&projectionAverage, const Voxel<NUM_MATERIAL> *voxelInput);

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

__global__ void computePolarization(Material<NUM_MATERIAL> materialInput,
                                    const Voxel<NUM_MATERIAL> *voxelInput,
                                    const ElectricField elefield,
                                    const Real angle,
                                    const uint3 voxel,
                                    Complex *polarizationX,
                                    Complex *polarizationY,
                                    Complex *polarizationZ,
                                    FFT::FFTWindowing windowing
);




__host__ inline cufftResult  performFFT(Complex *polarization, cufftHandle &plan) {
#ifdef DOUBLE_PRECISION
    return (cufftExecZ2Z(plan, polarization, polarization, CUFFT_FORWARD));
#else
    return (cufftExecC2C(plan, polarization, polarization, CUFFT_FORWARD));
#endif
}