//
// Created by maksbh on 4/11/21.
//

#ifndef CY_RSOXS_UNITTEST_H
#define CY_RSOXS_UNITTEST_H

#include <Input/InputData.h>
#include "testUtils.h"
#include "kernels.h"

TEST(CyRSoXS, polarization) {
  std::vector<Material<NUM_MATERIAL>> refractiveIndexData;
  InputData inputData(refractiveIndexData);
  const UINT voxelSize[3]{32,32,16};
  Voxel<NUM_MATERIAL> * voxelData;
  H5::readFile("edgespheres32.hd5",voxelSize,voxelData,MorphologyType::VECTOR_MORPHOLOGY,false);
  const BigUINT  numVoxels = voxelSize[0]*voxelSize[1]*voxelSize[2];
  Complex *polarizationX, *polarizationY,*polarizationZ;
  computePolarization(refractiveIndexData[0],voxelData,inputData,polarizationX,polarizationY,polarizationZ);
  Complex *pXOracle = new Complex [numVoxels];
  Complex *pYOracle = new Complex [numVoxels];
  Complex *pZOracle = new Complex [numVoxels];
  readFile(pXOracle,"polarizeX.dmp",numVoxels);
  readFile(pYOracle,"polarizeY.dmp",numVoxels);
  readFile(pZOracle,"polarizeZ.dmp",numVoxels);
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
