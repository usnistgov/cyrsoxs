/////////////////////////////////////////////////////////////////////////////////
// MIT License
//
//Copyright (c) 2019 - 2020 Iowa State University
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.
//////////////////////////////////////////////////////////////////////////////////

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
  H5::H5File file(filename.c_str(), H5F_ACC_TRUNC);
  try{
  const int RANK = 2;
  const hsize_t dims[2]{dim[0], dim[1]};

  H5::DataSpace dataspace(RANK, dims);
  H5::DataSet dataset = file.createDataSet("projection", H5::PredType::NATIVE_FLOAT, dataspace);
  dataset.write(data, H5::PredType::NATIVE_FLOAT);
  dataset.close();
  }
  catch(H5::FileIException error) {
    error.printErrorStack();
  }
  catch(H5::DataSetIException error)
  {
    error.printErrorStack();

  }
  catch(H5::DataSpaceIException error)
  {
    error.printErrorStack();

  }
}
}
#endif //CY_RSOXS_WRITEH5_H
