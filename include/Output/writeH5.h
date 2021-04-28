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
    void writeFile2D(H5::H5File & file, const Real *data, const UINT *dim, const std::string & groupname) {
        try {
            const int RANK = 2;
            const hsize_t dims[2]{dim[1], dim[0]};
            H5::Group group(file.createGroup(groupname.c_str()));
            H5::DataSpace dataspace(RANK, dims);
#ifdef DOUBLE_PRECISION
            H5::DataSet dataset = file.createDataSet("projection", H5::PredType::NATIVE_DOUBLE, dataspace);
            dataset.write(data, H5::PredType::NATIVE_DOUBLE);
#else
            H5::DataSet dataset = file.createDataSet(groupname+"/projection", H5::PredType::NATIVE_FLOAT, dataspace);
            dataset.write(data, H5::PredType::NATIVE_FLOAT);
#endif
            dataset.close();
        }
        catch (H5::FileIException & error) {
            H5::FileIException::printErrorStack();
        }
        catch (H5::DataSetIException & error) {
            H5::DataSetIException::printErrorStack();

        }
        catch (H5::DataSpaceIException & error) {
            H5::DataSpaceIException::printErrorStack();

        }
    }

    void writeFile3DScalar(const std::string& fname, const Real *data, const UINT *dim, const std::string& varName) {
        const std::string filename = fname + ".h5";
        H5::H5File file(filename.c_str(), H5F_ACC_TRUNC);
        try {
            const int RANK = 3;
            const hsize_t dims[3]{dim[2], dim[1], dim[0]};

            H5::DataSpace dataspace(RANK, dims);
#ifdef DOUBLE_PRECISION
            H5::DataSet dataset = file.createDataSet(varName.c_str(), H5::PredType::NATIVE_DOUBLE, dataspace);
            dataset.write(data, H5::PredType::NATIVE_DOUBLE);
#else
            H5::DataSet dataset = file.createDataSet(varName.c_str(), H5::PredType::NATIVE_FLOAT, dataspace);
            dataset.write(data, H5::PredType::NATIVE_FLOAT);
#endif
            dataset.close();
        }
        catch (H5::FileIException & error) {
            H5::FileIException::printErrorStack();
        }
        catch (H5::DataSetIException & error) {
            H5::DataSetIException::printErrorStack();

        }
        catch (H5::DataSpaceIException & error) {
            H5::DataSpaceIException::printErrorStack();

        }
    }

    void writeFile3DVector(const std::string& fname, const Real *data, const UINT *dim, const std::string& _varName) {

        const std::string filename = fname + ".h5";
        H5::H5File file(filename.c_str(), H5F_ACC_TRUNC);
        const char *dimString[]{"_x", "_y", "_z"};
        const int RANK = 3;
        BigUINT  numVoxel = dim[0]*dim[1]*dim[2];
        Real * temp = new Real [numVoxel];
        for (int d = 0; d < 3; d++) {
            for(BigUINT i = 0; i < numVoxel; i++){
                temp[i] = data[i*3+d];
            }
            std::string varName = _varName + dimString[d];

            try {
                const hsize_t dims[4]{dim[2], dim[1], dim[0]};

                H5::DataSpace dataspace(RANK, dims);
#ifdef DOUBLE_PRECISION
                H5::DataSet dataset = file.createDataSet(varName.c_str(), H5::PredType::NATIVE_DOUBLE, dataspace);
                dataset.write(data, H5::PredType::NATIVE_DOUBLE);
#else
                H5::DataSet dataset = file.createDataSet(varName.c_str(), H5::PredType::NATIVE_FLOAT, dataspace);
                dataset.write(temp, H5::PredType::NATIVE_FLOAT);
#endif
                dataset.close();
            }
            catch (H5::FileIException & error) {
                H5::FileIException::printErrorStack();
            }
            catch (H5::DataSetIException & error) {
                H5::DataSetIException::printErrorStack();

            }
            catch (H5::DataSpaceIException & error) {
                H5::DataSpaceIException::printErrorStack();

            }
        }
        delete [] temp;
    }
}


#endif //CY_RSOXS_WRITEH5_H
