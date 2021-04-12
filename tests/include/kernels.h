//
// Created by maksbh on 4/11/21.
//

#ifndef CY_RSOXS_KERNELS_H
#define CY_RSOXS_KERNELS_H
#include <Datatypes.h>
#include <Input/InputData.h>

void computePolarization(Material<NUM_MATERIAL> &refracticeIndexData,
                         const Voxel<NUM_MATERIAL>* voxelInput,
                         const InputData & inputData,
                         Complex *& polarizationX,
                         Complex *& polarizationY,
                         Complex *& polarizationZ
);


#endif //CY_RSOXS_KERNELS_H
