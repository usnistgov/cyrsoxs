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

#ifndef CUDA_BASE_CUDAHEADERS_H
#define CUDA_BASE_CUDAHEADERS_H

// #include <cuda_runtime_api.h>
#include <hip/hip_runtime.h>
// #include <driver_types.h>
// #include <device_launch_parameters.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <stdlib.h>

#define INLINE __attribute__((always_inline))

#define CUDA_CHECK_RETURN(value) CheckCudaErrorAux(__FILE__,__LINE__, #value, value)


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(hipError_t code, const char *file, int line, bool abort = true)
{
    if (code != hipSuccess)
    {
        fprintf(stderr, "GPUassert: %s %s %d\n", hipGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

/**
 * Check the return value of the CUDA runtime API call and exit
 * the application if the call has failed.
 *
 */
static void CheckCudaErrorAux (const char *file, unsigned line, const char *statement, hipError_t err)
{
    if (err == hipSuccess)
        return;
    std::cerr << statement<<" returned " << hipGetErrorString(err) << "("<<err<< ") at "<<file<<":"<<line << std::endl;
    exit (EXIT_FAILURE);
}

template <typename T, typename GI>
__host__ INLINE inline void hostDeviceExchange(T * dest, const T * src, const GI & size, const hipMemcpyKind direction){
  CUDA_CHECK_RETURN(hipMemcpy(dest,src,sizeof(T) * size ,direction));
  gpuErrchk(hipPeekAtLastError());

}
template <typename T, typename GI>
__host__ INLINE inline void mallocGPU(T *& d_data, const GI & size){
  CUDA_CHECK_RETURN(hipMalloc((void **) &d_data, sizeof(T) * size));
  gpuErrchk(hipPeekAtLastError());

}

template <typename T, typename GI>
__host__ INLINE inline void mallocCPU(T *& data, const GI & size){
  data = new T[size];
}

template <typename T, typename GI>
__host__ INLINE inline void mallocCPUPinned(T *& data, const GI & size){
  CUDA_CHECK_RETURN(hipMalloc((void**)&data,sizeof(T)*size));
  gpuErrchk(hipPeekAtLastError());
}

template <typename T, typename GI>
__host__ INLINE inline void cudaZeroEntries(T * d_data, const GI & size){
  gpuErrchk(hipMemset(d_data,0,sizeof(T)*size));
}

#define freeCudaMemory(X) CUDA_CHECK_RETURN(hipFree(X)); gpuErrchk(hipPeekAtLastError());
#ifdef DOUBLE_PRECISION
#define cublasScale hiplasDscal
#define cublasAXPY hiplasDaxpy
#define warpAffine rppiWarpAffine_64f_C1R
#define cufftC2C hipfftExecZ2Z
#else
#define cublasScale hipblasSscal
#define cublasAXPY hipblasSaxpy
#define warpAffine rppiWarpAffine_32f_C1R
#define cufftC2C hipfftExecC2C
#endif

#endif //CUDA_BASE_CUDAHEADERS_H
