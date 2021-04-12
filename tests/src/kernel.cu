//
// Created by maksbh on 4/11/21.
//

#include <kernels.h>
#include <cudaHeaders.h>
#include <cudaMain.h>

void computePolarization(Material<NUM_MATERIAL> &refracticeIndexData,
                         const Voxel<NUM_MATERIAL>* voxelInput,
                         const InputData & inputData,
                         Complex *& polarizationX,
                         Complex *& polarizationY,
                         Complex *& polarizationZ
){
  const UINT voxelSize[3]{inputData.numX,inputData.numY,inputData.numZ};
  const BigUINT numVoxels = inputData.numX*inputData.numY*inputData.numZ;
  Voxel<NUM_MATERIAL> *d_voxelInput;
  CUDA_CHECK_RETURN(cudaMalloc((void **) &d_voxelInput, sizeof(Voxel<NUM_MATERIAL>) * numVoxels));
  gpuErrchk(cudaPeekAtLastError());

  CUDA_CHECK_RETURN(cudaMemcpy(d_voxelInput,
                               voxelInput,
                               sizeof(Voxel<NUM_MATERIAL>) * numVoxels,
                               cudaMemcpyHostToDevice));
  gpuErrchk(cudaPeekAtLastError());

  Complex *d_polarizationX,*d_polarizationY,*d_polarizationZ;
  CUDA_CHECK_RETURN(cudaMalloc((void **) &d_polarizationX, sizeof(Complex) * numVoxels));
  gpuErrchk(cudaPeekAtLastError());
  CUDA_CHECK_RETURN(cudaMalloc((void **) &d_polarizationY, sizeof(Complex) * numVoxels));
  gpuErrchk(cudaPeekAtLastError());
  CUDA_CHECK_RETURN(cudaMalloc((void **) &d_polarizationZ, sizeof(Complex) * numVoxels));
  gpuErrchk(cudaPeekAtLastError());

  ElectricField eleField;
  eleField.e.x = 1;
  eleField.e.y = 0;
  eleField.e.z = 0;
  Real wavelength = static_cast<Real>(1239.84197 / inputData.energies[0]);
  eleField.k.x = 0;
  eleField.k.y = 0;
  eleField.k.z = static_cast<Real>(2 * M_PI / wavelength);;
  uint3 vx{voxelSize[0],voxelSize[1],voxelSize[2]};
  UINT BlockSize = static_cast<UINT >(ceil(numVoxels * 1.0 / NUM_THREADS));
  computePolarization <<< BlockSize, NUM_THREADS >>> (refracticeIndexData, d_voxelInput, eleField, 0, vx, d_polarizationX, d_polarizationY, d_polarizationZ,
                                                      static_cast<FFT::FFTWindowing >(inputData.windowingType), inputData.if2DComputation());

  polarizationX = new Complex[numVoxels];
  polarizationY = new Complex[numVoxels];
  polarizationZ = new Complex[numVoxels];

  cudaDeviceSynchronize();
  CUDA_CHECK_RETURN(cudaMemcpy(polarizationX,
                               d_polarizationX,
                               sizeof(Complex) * numVoxels,
                               cudaMemcpyDeviceToHost));
  gpuErrchk(cudaPeekAtLastError());
  CUDA_CHECK_RETURN(cudaMemcpy(polarizationZ,
                               d_polarizationZ,
                               sizeof(Complex) * numVoxels,
                               cudaMemcpyDeviceToHost));
  gpuErrchk(cudaPeekAtLastError());
  CUDA_CHECK_RETURN(cudaMemcpy(polarizationY,
                               d_polarizationY,
                               sizeof(Complex) * numVoxels,
                               cudaMemcpyDeviceToHost));
  gpuErrchk(cudaPeekAtLastError());
  cudaFree(d_polarizationX);
  cudaFree(d_polarizationY);
  cudaFree(d_polarizationZ);
  cudaFree(d_voxelInput);


}
