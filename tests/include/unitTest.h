//
// Created by maksbh on 4/11/21.
//

#ifndef CY_RSOXS_UNITTEST_H
#define CY_RSOXS_UNITTEST_H

#include <Input/InputData.h>
#include "testUtils.h"
#include "kernels.h"

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
#endif //CY_RSOXS_UNITTEST_H
