//
// Created by maksbh on 4/14/21.
//

#ifndef CY_RSOXS_FULLFRAMEWORK_H
#define CY_RSOXS_FULLFRAMEWORK_H
TEST(CyRSoXS, fullFramework) {
  const std::string root = CMAKE_ROOT ;
  const std::string fname = root + "/Data/sample.h5";
  const std::string oracleFname = root + "/Data/regressionData/fullData/Energy_280.00.h5";
  const std::string configPath = root + "/Data/config/";
  if(cd(configPath.c_str()) != 0){
    throw std::runtime_error("Wrong path for config");
  }
  std::vector<Material<NUM_MATERIAL>> refractiveIndexData;
  InputData inputData(refractiveIndexData);
  inputData.startAngle = 0;
  inputData.endAngle = 180;
  inputData.incrementAngle = 0.2;
  inputData.ewaldsInterpolation = Interpolation::EwaldsInterpolation::LINEAR;
  const UINT voxelSize[3]{32,32,16};
  Voxel<NUM_MATERIAL> * voxelData;
  H5::readFile(fname,voxelSize,voxelData,MorphologyType::VECTOR_MORPHOLOGY,false);
  Real *projectionGPUAveraged;
  const UINT numEnergyLevel = inputData.energies.size();
  projectionGPUAveraged = new Real[numEnergyLevel * inputData.numX * inputData.numY];
  int suc = cudaMain(voxelSize, inputData, refractiveIndexData, projectionGPUAveraged, voxelData);
  EXPECT_EQ(suc,EXIT_SUCCESS);
  H5::H5File file(oracleFname,H5F_ACC_RDONLY);
  H5::DataSet dataSet = file.openDataSet("projection");
  Real * oracleData = new Real[numEnergyLevel * inputData.numX * inputData.numY];
  dataSet.read(oracleData, H5::PredType::NATIVE_FLOAT);
  dataSet.close();
  file.close();
  BigUINT  numVoxels = voxelSize[0]*voxelSize[1];
  Real linfError = computeLinfError(oracleData,projectionGPUAveraged,numVoxels);
  EXPECT_LE(linfError,TOLERANCE_CHECK);
  delete [] oracleData;
  delete [] projectionGPUAveraged;
}
#endif //CY_RSOXS_FULLFRAMEWORK_H
