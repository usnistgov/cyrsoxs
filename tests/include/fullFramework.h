//
// Created by maksbh on 4/14/21.
//

#ifndef CY_RSOXS_FULLFRAMEWORK_H
#define CY_RSOXS_FULLFRAMEWORK_H
TEST(CyRSoXS, fullFrameworkSingleEnergy) {
  const std::string root = CMAKE_ROOT ;
  const std::string fname = root + "/Data/sample.h5";
  const std::string oracleFname = root + "/Data/regressionData/fullData/Energy_280.00.h5";
  const std::string configPath = root + "/Data/config/";
  if(cd(configPath.c_str()) != 0){
    throw std::runtime_error("Wrong path for config");
  }
  std::vector<Material<NUM_MATERIAL>> refractiveIndexData;
  InputData inputData;
  inputData.startAngle = 0;
  inputData.endAngle = 180;
  inputData.incrementAngle = 0.2;
  inputData.ewaldsInterpolation = Interpolation::EwaldsInterpolation::LINEAR;
  inputData.readRefractiveIndexData(refractiveIndexData);
  for(UINT scatter = 0; scatter < ScatterApproach::MAX_SCATTER_APPROACH; scatter++) {
    inputData.scatterApproach = static_cast<ScatterApproach>(scatter);
    const UINT voxelSize[3]{32, 32, 16};
    Voxel *voxelData;
    H5::readFile(fname, voxelSize, voxelData, MorphologyType::VECTOR_MORPHOLOGY, false);
    Real *projectionGPUAveraged;
    const UINT numEnergyLevel = inputData.energies.size();
    projectionGPUAveraged = new Real[numEnergyLevel * inputData.numX * inputData.numY];
    RotationMatrix rotationMatrix(&inputData);
    int suc = cudaMain(voxelSize, inputData, refractiveIndexData, projectionGPUAveraged, rotationMatrix,voxelData);
    EXPECT_EQ(suc, EXIT_SUCCESS);
    H5::H5File file(oracleFname, H5F_ACC_RDONLY);
    H5::DataSet dataSet = file.openDataSet("projection");
    Real *oracleData = new Real[numEnergyLevel * inputData.numX * inputData.numY];
    dataSet.read(oracleData, H5::PredType::NATIVE_FLOAT);
    dataSet.close();
    file.close();
    BigUINT numVoxels = voxelSize[0] * voxelSize[1];
    Real linfError = computeLinfError(oracleData, projectionGPUAveraged, numVoxels);
    EXPECT_LE(linfError, TOLERANCE_CHECK);
    delete[] oracleData;
    delete[] projectionGPUAveraged;
    delete[] voxelData;
  }
}

TEST(CyRSoXS, fullFrameworkSingleEnergyAlg2) {
  const std::string root = CMAKE_ROOT ;
  const std::string fname = root + "/Data/sample.h5";
  const std::string oracleFname = root + "/Data/regressionData/fullData/Energy_280.00.h5";
  const std::string configPath = root + "/Data/config/";
  if(cd(configPath.c_str()) != 0){
    throw std::runtime_error("Wrong path for config");
  }
  std::vector<Material<NUM_MATERIAL>> refractiveIndexData;
  InputData inputData;
  inputData.startAngle = 0;
  inputData.endAngle = 180;
  inputData.incrementAngle = 0.2;
  inputData.ewaldsInterpolation = Interpolation::EwaldsInterpolation::LINEAR;
  inputData.readRefractiveIndexData(refractiveIndexData);
  for(UINT scatter = 0; scatter < ScatterApproach::MAX_SCATTER_APPROACH; scatter++) {
    inputData.scatterApproach = static_cast<ScatterApproach>(scatter);
    const UINT voxelSize[3]{32, 32, 16};
    Voxel *voxelData;
    mallocCPUPinned(voxelData,inputData.numX*inputData.numY*inputData.numZ*NUM_MATERIAL);
    H5::readFile(fname, voxelSize, voxelData, MorphologyType::VECTOR_MORPHOLOGY, true);
    Real *projectionGPUAveraged;
    const UINT numEnergyLevel = inputData.energies.size();
    projectionGPUAveraged = new Real[numEnergyLevel * inputData.numX * inputData.numY];
    RotationMatrix rotationMatrix(&inputData);
    int suc = cudaMainStreams(voxelSize, inputData, refractiveIndexData, projectionGPUAveraged, rotationMatrix,voxelData);
    EXPECT_EQ(suc, EXIT_SUCCESS);
    H5::H5File file(oracleFname, H5F_ACC_RDONLY);
    H5::DataSet dataSet = file.openDataSet("projection");
    Real *oracleData = new Real[numEnergyLevel * inputData.numX * inputData.numY];
    dataSet.read(oracleData, H5::PredType::NATIVE_FLOAT);
    dataSet.close();
    file.close();
    BigUINT numVoxels = voxelSize[0] * voxelSize[1];
    Real linfError = computeLinfError(oracleData, projectionGPUAveraged, numVoxels);
    EXPECT_LE(linfError, TOLERANCE_CHECK);
    delete[] oracleData;
    delete[] projectionGPUAveraged;
    cudaFreeHost(voxelData);
  }
}

TEST(CyRSoXS, fullFrameworkMultipleEnergy) {
  const std::string root = CMAKE_ROOT ;
  const std::string fname = root + "/Data/sample.h5";

  const std::string configPath = root + "/Data/config/";
  if(cd(configPath.c_str()) != 0){
    throw std::runtime_error("Wrong path for config");
  }



  std::vector<Material<NUM_MATERIAL>> refractiveIndexData;
  InputData inputData;
  inputData.startAngle = 0;
  inputData.endAngle = 180;
  inputData.incrementAngle = 0.2;
  inputData.ewaldsInterpolation = Interpolation::EwaldsInterpolation::LINEAR;
  inputData.energies.clear();
  inputData.energies.resize(4);
  inputData.energies[0] = 280.0;
  inputData.energies[1] = 280.1;
  inputData.energies[2] = 280.2;
  inputData.energies[3] = 280.3;
  inputData.readRefractiveIndexData(refractiveIndexData);
  for(UINT scatter = 0; scatter < ScatterApproach::MAX_SCATTER_APPROACH; scatter++) {
    inputData.scatterApproach = static_cast<ScatterApproach>(scatter);
    const UINT voxelSize[3]{32, 32, 16};
    Voxel *voxelData;
    H5::readFile(fname, voxelSize, voxelData, MorphologyType::VECTOR_MORPHOLOGY, false);
    Real *projectionGPUAveraged;
    const UINT numEnergyLevel = inputData.energies.size();
    projectionGPUAveraged = new Real[numEnergyLevel * inputData.numX * inputData.numY];
    RotationMatrix rotationMatrix(&inputData);
    int suc = cudaMain(voxelSize, inputData, refractiveIndexData, projectionGPUAveraged, rotationMatrix, voxelData);
    EXPECT_EQ(suc, EXIT_SUCCESS);
    static const char *oracleFileName[]{"Energy_280.00.h5", "Energy_280.10.h5", "Energy_280.20.h5", "Energy_280.30.h5"};
    BigUINT numVoxels2D = voxelSize[0] * voxelSize[1];
    Real *oracleData = new Real[numVoxels2D];
    for (int i = 0; i < 4; i++) {
      const std::string oracleName = root + "/Data/regressionData/fullData/" + oracleFileName[i];
      H5::H5File file(oracleName, H5F_ACC_RDONLY);
      H5::DataSet dataSet = file.openDataSet("projection");

      dataSet.read(oracleData, H5::PredType::NATIVE_FLOAT);
      dataSet.close();
      file.close();
      Real linfError = computeLinfError(oracleData, &projectionGPUAveraged[i * numVoxels2D], numVoxels2D);
      EXPECT_LE(linfError, TOLERANCE_CHECK);
    }
    delete [] projectionGPUAveraged;
    delete [] oracleData;
    delete [] voxelData;
  }

}

TEST(CyRSoXS, fullFrameworkMultipleEnergyAlg2) {
  const std::string root = CMAKE_ROOT ;
  const std::string fname = root + "/Data/sample.h5";

  const std::string configPath = root + "/Data/config/";
  if(cd(configPath.c_str()) != 0){
    throw std::runtime_error("Wrong path for config");
  }



  std::vector<Material<NUM_MATERIAL>> refractiveIndexData;
  InputData inputData;
  inputData.startAngle = 0;
  inputData.endAngle = 180;
  inputData.incrementAngle = 0.2;
  inputData.ewaldsInterpolation = Interpolation::EwaldsInterpolation::LINEAR;
  inputData.energies.clear();
  inputData.energies.resize(4);
  inputData.energies[0] = 280.0;
  inputData.energies[1] = 280.1;
  inputData.energies[2] = 280.2;
  inputData.energies[3] = 280.3;
  inputData.readRefractiveIndexData(refractiveIndexData);
  for(UINT scatter = 0; scatter < ScatterApproach::MAX_SCATTER_APPROACH; scatter++) {
    inputData.scatterApproach = static_cast<ScatterApproach>(scatter);
    const UINT voxelSize[3]{32, 32, 16};
    Voxel *voxelData;
    mallocCPUPinned(voxelData,inputData.numX*inputData.numY*inputData.numZ*NUM_MATERIAL);
    H5::readFile(fname, voxelSize, voxelData, MorphologyType::VECTOR_MORPHOLOGY, true);
    Real *projectionGPUAveraged;
    const UINT numEnergyLevel = inputData.energies.size();
    projectionGPUAveraged = new Real[numEnergyLevel * inputData.numX * inputData.numY];
    RotationMatrix rotationMatrix(&inputData);
    int suc = cudaMainStreams(voxelSize, inputData, refractiveIndexData, projectionGPUAveraged, rotationMatrix, voxelData);
    EXPECT_EQ(suc, EXIT_SUCCESS);
    static const char *oracleFileName[]{"Energy_280.00.h5", "Energy_280.10.h5", "Energy_280.20.h5", "Energy_280.30.h5"};
    BigUINT numVoxels2D = voxelSize[0] * voxelSize[1];
    Real *oracleData = new Real[numVoxels2D];
    for (int i = 0; i < 4; i++) {
      const std::string oracleName = root + "/Data/regressionData/fullData/" + oracleFileName[i];
      H5::H5File file(oracleName, H5F_ACC_RDONLY);
      H5::DataSet dataSet = file.openDataSet("projection");

      dataSet.read(oracleData, H5::PredType::NATIVE_FLOAT);
      dataSet.close();
      file.close();
      Real linfError = computeLinfError(oracleData, &projectionGPUAveraged[i * numVoxels2D], numVoxels2D);
      EXPECT_LE(linfError, TOLERANCE_CHECK);
    }
    delete [] projectionGPUAveraged;
    delete [] oracleData;
    cudaFreeHost(voxelData);
  }

}
#endif //CY_RSOXS_FULLFRAMEWORK_H
