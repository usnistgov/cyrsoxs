//
// Created by maksbh on 4/11/21.
//

#ifndef CY_RSOXS_UNITTEST_H
#define CY_RSOXS_UNITTEST_H

#include <Input/InputData.h>
#include "testUtils.h"
#include "kernels.h"
//#include "cudaUtils.h"
//template <typename T, typename GI>
//__host__ inline void hostDeviceExcange(T * dest, const T * src, const GI & size, const cudaMemcpyKind direction){
//  CUDA_CHECK_RETURN(cudaMemcpy(dest,src,sizeof(T) * size ,direction));
//  gpuErrchk(cudaPeekAtLastError());
//
//}
//template <typename T, typename GI>
//__host__ inline void mallocGPU(T *& d_data, const GI & size){
//  CUDA_CHECK_RETURN(cudaMalloc((void **) &d_data, sizeof(T) * size));
//  gpuErrchk(cudaPeekAtLastError());
//
//}

#ifdef _WIN32
#include <direct.h>
// MSDN recommends against using getcwd & chdir names
#define cwd _getcwd
#define cd _chdir
#else
#include "unistd.h"
#define cwd getcwd
#define cd chdir
#endif

TEST(CyRSoXS, polarization) {
  const std::string root = CMAKE_ROOT ;
  const std::string fname = root + "/Data/edgespheres32.hd5";
  const std::string configPath = root + "/Data/regressionData/config/";
  if(cd(configPath.c_str()) != 0){
    throw std::runtime_error("Wrong path for config");
  }
  std::vector<Material<NUM_MATERIAL>> refractiveIndexData;
  InputData inputData(refractiveIndexData);
  const UINT voxelSize[3]{32,32,16};
  Voxel<NUM_MATERIAL> * voxelData;


  H5::readFile(fname,voxelSize,voxelData,MorphologyType::VECTOR_MORPHOLOGY,false);
  const BigUINT  numVoxels = voxelSize[0]*voxelSize[1]*voxelSize[2];
  Complex *polarizationX, *polarizationY,*polarizationZ;
  computePolarization(refractiveIndexData[0],voxelData,inputData,polarizationX,polarizationY,polarizationZ);

  Complex *pXOracle = new Complex [numVoxels];
  Complex *pYOracle = new Complex [numVoxels];
  Complex *pZOracle = new Complex [numVoxels];
  std::string pathOfOracle = root+"/Data/regressionData/polarization/";
  readFile(pXOracle,pathOfOracle+"polarizeX.dmp",numVoxels);
  readFile(pYOracle,pathOfOracle+"polarizeY.dmp",numVoxels);
  readFile(pZOracle,pathOfOracle+"polarizeZ.dmp",numVoxels);
  Complex linfError;
  linfError = computeLinfError(pXOracle,polarizationX,numVoxels);
  EXPECT_LE(linfError.x,TOLERANCE_CHECK);
  EXPECT_LE(linfError.y,TOLERANCE_CHECK);
  linfError = computeLinfError(pYOracle,polarizationY,numVoxels);
  EXPECT_LE(linfError.x,TOLERANCE_CHECK);
  EXPECT_LE(linfError.y,TOLERANCE_CHECK);
  linfError = computeLinfError(pZOracle,polarizationZ,numVoxels);
  EXPECT_LE(linfError.x,TOLERANCE_CHECK);
  EXPECT_LE(linfError.y,TOLERANCE_CHECK);

  delete [] pXOracle;
  delete [] pYOracle;
  delete [] pZOracle;

  delete [] voxelData;
  delete [] polarizationX;
  delete [] polarizationY;
  delete [] polarizationZ;
}


TEST(CyRSoXS, FFT) {

  const UINT voxelSize[3]{32,32,16};
  const BigUINT  numVoxels = voxelSize[0]*voxelSize[1]*voxelSize[2];
  const std::string root = CMAKE_ROOT ;
  Complex *pX = new Complex [numVoxels];
  Complex *pY = new Complex [numVoxels];
  Complex *pZ = new Complex [numVoxels];
  const std::string pathOfpolarization = root+"/Data/regressionData/polarization/";
  readFile(pX,pathOfpolarization+"polarizeX.dmp",numVoxels);
  readFile(pY,pathOfpolarization+"polarizeY.dmp",numVoxels);
  readFile(pZ,pathOfpolarization+"polarizeZ.dmp",numVoxels);
  Complex *d_pX, *d_pY,*d_pZ;
  mallocGPU(d_pX,numVoxels);
  mallocGPU(d_pY,numVoxels);
  mallocGPU(d_pZ,numVoxels);
  hostDeviceExcange(d_pX,pX,numVoxels,cudaMemcpyHostToDevice);
  hostDeviceExcange(d_pY,pY,numVoxels,cudaMemcpyHostToDevice);
  hostDeviceExcange(d_pZ,pZ,numVoxels,cudaMemcpyHostToDevice);

  cufftHandle plan;
  cufftPlan3d(&plan, voxelSize[2], voxelSize[1], voxelSize[0], fftType);

  performFFT(d_pX,plan);
  performFFT(d_pY,plan);
  performFFT(d_pZ,plan);
  UINT blockSize = static_cast<UINT >(ceil(numVoxels * 1.0 / NUM_THREADS));
  const uint3 vx{voxelSize[0],voxelSize[1],voxelSize[2]};
  performFFTShift(d_pX,blockSize,vx);
  performFFTShift(d_pY,blockSize,vx);
  performFFTShift(d_pZ,blockSize,vx);
//  CUDA_CHECK_RETURN(cudaMemcpy(data,d_data,sizeof(T) * size ,direction));
//  gpuErrchk(cudaPeekAtLastError());
  hostDeviceExcange(pX,d_pX,numVoxels,cudaMemcpyDeviceToHost);
  hostDeviceExcange(pY,d_pY,numVoxels,cudaMemcpyDeviceToHost);
  hostDeviceExcange(pZ,d_pZ,numVoxels,cudaMemcpyDeviceToHost);

  const std::string pathOfFFT = root+"/Data/regressionData/FFT/";
  Complex *fftPX = new Complex [numVoxels];
  Complex *fftPY = new Complex [numVoxels];
  Complex *fftPZ = new Complex [numVoxels];

  readFile(fftPX,pathOfFFT+"fftpolarizeX.dmp",numVoxels);
  readFile(fftPY,pathOfFFT+"fftpolarizeY.dmp",numVoxels);
  readFile(fftPZ,pathOfFFT+"fftpolarizeZ.dmp",numVoxels);
  Complex linfError;
  linfError = computeLinfError(pX,fftPX,numVoxels);
  EXPECT_LE(linfError.x,TOLERANCE_CHECK);
  EXPECT_LE(linfError.y,TOLERANCE_CHECK);
  linfError = computeLinfError(pY,fftPY,numVoxels);
  EXPECT_LE(linfError.x,TOLERANCE_CHECK);
  EXPECT_LE(linfError.y,TOLERANCE_CHECK);
  linfError = computeLinfError(pZ,fftPZ,numVoxels);
  EXPECT_LE(linfError.x,TOLERANCE_CHECK);
  EXPECT_LE(linfError.y,TOLERANCE_CHECK);

  delete [] pX;
  delete [] pY;
  delete [] pZ;

  delete [] fftPX;
  delete [] fftPY;
  delete [] fftPZ;

}
#endif //CY_RSOXS_UNITTEST_H
