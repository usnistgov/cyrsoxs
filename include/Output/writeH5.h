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
#include <Output/outputUtils.h>

namespace H5 {
/**
 * @brief Writes the final scattering pattern data in HDF5 file format
 * @param data the data
 * @param dim the dimension corresponding to the X and Y
 */
  void writeFile2D(H5::H5File &file, const Real *data, const UINT *dim, const std::string &groupname) {
    try {
      const int RANK = 2;
      const hsize_t dims[2]{dim[1], dim[0]};
      H5::Group group(file.createGroup(groupname.c_str()));
      H5::DataSpace dataspace(RANK, dims);
#ifdef DOUBLE_PRECISION
      H5::DataSet dataset = file.createDataSet("projection", H5::PredType::NATIVE_DOUBLE, dataspace);
      dataset.write(data, H5::PredType::NATIVE_DOUBLE);
#else
      H5::DataSet dataset = file.createDataSet(groupname + "/projection", H5::PredType::NATIVE_FLOAT, dataspace);
      dataset.write(data, H5::PredType::NATIVE_FLOAT);
#endif
      dataset.close();
    }
    catch (H5::FileIException &error) {
      H5::FileIException::printErrorStack();
    }
    catch (H5::DataSetIException &error) {
      H5::DataSetIException::printErrorStack();

    }
    catch (H5::DataSpaceIException &error) {
      H5::DataSpaceIException::printErrorStack();

    }
  }

  void writeMorphologyFile(const std::string &fname, const InputData &inputData, const Voxel *morphologyData) {
    std::array<std::string, 4> morphologyDataSets{};
    const MorphologyType morphologyType = static_cast<MorphologyType>(inputData.morphologyType);
    if (morphologyType == MorphologyType::VECTOR_MORPHOLOGY) {
      morphologyDataSets[0] = "UnalignedFraction";
      morphologyDataSets[1] = "Sx";
      morphologyDataSets[2] = "Sy";
      morphologyDataSets[3] = "Sz";
    } else if (morphologyType == MorphologyType::EULER_ANGLES) {
      morphologyDataSets[0] = "Vfrac";
      morphologyDataSets[1] = "S";
      morphologyDataSets[2] = "phi";
      morphologyDataSets[3] = "theta";
    } else {
      throw std::runtime_error("Unexpected Morphology Type");
    }
    const std::string filename = fname + ".h5";
    H5::H5File file(filename.c_str(), H5F_ACC_TRUNC);
    try {

      const BigUINT numVoxels = inputData.voxelDims[0] * inputData.voxelDims[1] * inputData.voxelDims[2];
      Real *data = new Real[numVoxels];
      const int RANK = 3;
      const hsize_t dims[3]{inputData.voxelDims[2], inputData.voxelDims[1], inputData.voxelDims[0]}; // C++ order

      for (int numMat = 0; numMat < NUM_MATERIAL; numMat++) {
        const auto &group = file.createGroup("Material_" + std::to_string(numMat));
        for (int numComponent = 0; numComponent < 4; numComponent++) {
          for (int numVoxel = 0; numVoxel < numVoxels; numVoxel++) {
            int offset = (numVoxel * NUM_MATERIAL) + numMat;
            data[numVoxel] = morphologyData[offset].getValueAt(numComponent);
          }
          H5::DataSpace dataspace(RANK, dims);
#ifdef DOUBLE_PRECISION
          H5::DataSet dataset = group.createDataSet(morphologyDataSets[numComponent].c_str(), H5::PredType::NATIVE_FLOAT, dataspace);
          dataset.write(data, H5::PredType::NATIVE_DOUBLE);
#else
          H5::DataSet dataSet = group.createDataSet(morphologyDataSets[numComponent].c_str(),
                                                    H5::PredType::NATIVE_FLOAT, dataspace);
          dataSet.write(data, H5::PredType::NATIVE_FLOAT);
#endif
          dataSet.close();
        }
      }
      delete[] data;
    }
    catch (H5::FileIException &error) {
      H5::FileIException::printErrorStack();
    }
    catch (H5::DataSetIException &error) {
      H5::DataSetIException::printErrorStack();

    }
    catch (H5::DataSpaceIException &error) {
      H5::DataSpaceIException::printErrorStack();
    }
  }

  void writeXDMF(const InputData &inputData, const Voxel * voxelData) {
    std::array<std::string, 4> morphologyDataSets{};
    const MorphologyType morphologyType = static_cast<MorphologyType>(inputData.morphologyType);
    if (morphologyType == MorphologyType::VECTOR_MORPHOLOGY) {
      morphologyDataSets[0] = "UnalignedFraction";
      morphologyDataSets[1] = "Sx";
      morphologyDataSets[2] = "Sy";
      morphologyDataSets[3] = "Sz";
    } else if (morphologyType == MorphologyType::EULER_ANGLES) {
      morphologyDataSets[0] = "Vfrac";
      morphologyDataSets[1] = "S";
      morphologyDataSets[2] = "phi";
      morphologyDataSets[3] = "theta";
    } else {
      throw std::runtime_error("Unexpected Morphology Type");
    }
    const std::string dirName = "Morphology";
    const std::string fName = "Morphology";
    const std::string xdmfFileName = dirName+"/"+fName+".xdmf";
    static constexpr int precision = sizeof(Real);
    createDirectory(dirName);
    const UINT voxelSize[3]{inputData.voxelDims[2], inputData.voxelDims[1], inputData.voxelDims[0]}; // C++ format

    FILE *fp = fopen(xdmfFileName.c_str(), "w");
    fprintf(fp, "<?xml version=\"1.0\" ?>\n");
    fprintf(fp, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(fp, "<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.2\">\n");
    // Grid section
    fprintf(fp, "\t<Domain>\n");
    fprintf(fp, "\t<Grid Name=\"morphology\" GridType=\"Uniform\">\n");
    fprintf(fp, "\t\t<Topology TopologyType=\"3DCORECTMesh\" NumberOfElements=\" %d %d %d\" />\n", voxelSize[0],
            voxelSize[1], voxelSize[2]);
    fprintf(fp, "\t\t<Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n");
    fprintf(fp,
            "\t\t<DataItem Name=\"origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\" %d\" Format=\"XML\">\n",
            precision);
    fprintf(fp, "\t\t0.0 0.0 0.0\n");
    fprintf(fp, "\t\t</DataItem>\n");
    fprintf(fp,
            "\t\t<DataItem Name=\"origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\" %d\" Format=\"XML\">\n",
            precision);
    fprintf(fp, "\t\t1.0 1.0 1.0\n");
    fprintf(fp, "\t\t</DataItem>\n");
    fprintf(fp, "\t\\t</Geometry>\n");

    for (int numComponent = 0; numComponent < 4; numComponent++) {
      for (int i = 1; i < NUM_MATERIAL + 1; i++) {
        fprintf(fp, "\t\t<Attribute Name=\"Mat_%d_data\" AttributeType=\"Scalar\" Center=\"Node\">\n", i);
        fprintf(fp, "\t\t<DataItem Dimensions=\" %d %d %d \"NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\"",
                voxelSize[0], voxelSize[1], voxelSize[2], precision);
        fprintf(fp, "\t\t %s:/Material_%d/%s\n", fName.c_str(),i,morphologyDataSets[numComponent].c_str());
        fprintf(fp,"\t\t</DataItem>\n");
        fprintf(fp,"\t\t</Attribute>\n");
      }
    }
    fprintf(fp,"\t</Grid>\n");
    fprintf(fp,"\t</Domain>\n");
    fprintf(fp,"</Xdmf>\n");
    fclose(fp);
    const std::string h5MorphologyName = dirName+"/"+fName;
    writeMorphologyFile(h5MorphologyName,inputData,voxelData);
  }
}


#endif //CY_RSOXS_WRITEH5_H
