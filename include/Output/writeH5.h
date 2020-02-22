//
// Created by maksbh on 2/22/20.
//

#ifndef CY_RSOXS_WRITEH5_H
#define CY_RSOXS_WRITEH5_H

#include <Datatypes.h>
#include "H5Cpp.h"
namespace H5 {
/**
 * @brief Writes the final scattering pattern data in HDF5 file format
 * @param data the data
 * @param dim the dimension corresponding to the X and Y
 */
void writeFile2D(const std::string fname, const Real *data, const UINT *dim) {
  const std::string filename =fname+".h5";
  H5::H5File *file = new H5::H5File(filename.c_str(), H5F_ACC_TRUNC);
  const int RANK = 2;
  const hsize_t dims[2]{dim[0], dim[1]};

  H5::DataSpace dataspace(RANK, dims);
  H5::DataSet dataset = file->createDataSet("projection", H5::PredType::NATIVE_FLOAT, dataspace);

  dataset.write(data, H5::PredType::NATIVE_FLOAT);
  delete file;

}
}
#endif //CY_RSOXS_WRITEH5_H
