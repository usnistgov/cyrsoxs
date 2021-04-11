//
// Created by maksbh on 4/10/21.
//

#ifndef CY_RSOXS_BASICFRAMEWORK_H
#define CY_RSOXS_BASICFRAMEWORK_H

#include <Datatypes.h>
#include <cudaMain.h>
#include <Input/readH5.h>

TEST(CyRSoXS, gpuSetup) {
EXPECT_EQ(EXIT_SUCCESS,warmup());
}

TEST(CyRSoXS,H5Reader){
  const std::string root = CMAKE_ROOT ;
  const std::string fname = root + "/Data/edgespheres32.hd5";
  const UINT voxelSize[3]{32,32,16};
  Voxel<NUM_MATERIAL> *voxelData;
  EXPECT_EQ(EXIT_SUCCESS,H5::readFile(fname, voxelSize, voxelData,MorphologyType::VECTOR_MORPHOLOGY));
  delete [] voxelData;
}


#endif //CY_RSOXS_BASICFRAMEWORK_H
