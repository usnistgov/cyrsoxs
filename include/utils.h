//
// Created by maksbh on 8/31/20.
//

#ifndef CY_RSOXS_UTILS_H
#define CY_RSOXS_UTILS_H
#include <Input/InputData.h>

void writeH5(){
    createDirectory("HDF5");
    omp_set_num_threads(1);

    const UINT
    numEnergyLevel =
    static_cast<UINT>(std::round(
    (inputData.energyEnd - inputData.energyStart) / inputData.incrementEnergy + 1));
    const UINT voxel2DSize = voxelSize[0] * voxelSize[1];
    Real *oneEnergyData = new Real[voxel2DSize];
    UINT chunkSize = static_cast<UINT>(std::ceil(numEnergyLevel * 1.0 / (omp_get_num_threads() * 1.0)));
    UINT threadID = omp_get_thread_num();
    UINT startID = threadID * chunkSize;
    UINT endID = ((threadID + 1) * chunkSize);

    for (UINT csize = startID; csize < std::min(endID, numEnergyLevel); csize++) {
      std::stringstream stream;
      Real energy = inputData.energyStart + csize * inputData.incrementEnergy;
      stream << std::fixed << std::setprecision(2) << energy;
      std::string s = stream.str();
      std::memcpy(oneEnergyData, &projectionGPUAveraged[csize * voxel2DSize], sizeof(Real) * voxel2DSize);
      const std::string outputFname = "HDF5/Energy_" + s;

      H5::writeFile2D(outputFname, oneEnergyData, voxelSize);
    }
};
#endif //CY_RSOXS_UTILS_H
