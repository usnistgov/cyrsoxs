/////////////////////////////////////////////////////////////////////////////////
// MIT License
//
//Copyright (c) 2019 - 2022 Iowa State University
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
#include <hdf5_hl.h>
namespace H5 {
/**
 * @brief Writes the final scattering pattern data in HDF5 file format
 * @param [in] file HDF5 file name
 * @param [in] data the data
 * @param [in] dim the dimension corresponding to the X and Y
 * @param groupname name of the group
 */
  void writeFile2D(H5::H5File &file, const Real *data, const UINT *dim, const std::string &groupname) {
    try {
      const int RANK = 2;
      const hsize_t dims[2]{dim[1], dim[0]};
      H5::Group group(file.createGroup(groupname.c_str()));
      H5::DataSpace dataspace(RANK, dims);
#ifdef DOUBLE_PRECISION
      H5::DataSet dataset = file.createDataSet(groupname + "/projection", H5::PredType::NATIVE_DOUBLE, dataspace);
      dataset.write(data, H5::PredType::NATIVE_DOUBLE);
#else
      H5::DataSet dataset = file.createDataSet(groupname + "/projection", H5::PredType::NATIVE_FLOAT, dataspace);
      dataset.write(data, H5::PredType::NATIVE_FLOAT);
#endif
      H5DSset_label(dataset.getId(),0,"Qy");
      H5DSset_label(dataset.getId(),1,"Qx");
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

  /**
   *
   * @param value polarization
   * @param inputData input data
   * @param _filename filename
   * @param dataSetName data set name
   */
  void writePolarization(const Complex *value, const InputData &inputData, const std::string &_filename, const std::string &dataSetName) {

    const std::string filename = _filename + ".h5";
    H5::H5File file(filename.c_str(), H5F_ACC_TRUNC);
    try {
      const int RANK = 4;
      const hsize_t dims[4]{inputData.voxelDims[2], inputData.voxelDims[1], inputData.voxelDims[0], 2}; // C++ order
      H5::DataSpace dataspace(RANK, dims);
#ifdef DOUBLE_PRECISION
      H5::DataSet dataSet = file.createDataSet(dataSetName.c_str(), H5::PredType::NATIVE_DOUBLE, dataspace);
      dataSet.write(value, H5::PredType::NATIVE_DOUBLE, dataspace);
#else
      H5::DataSet dataSet = file.createDataSet(dataSetName.c_str(), H5::PredType::NATIVE_FLOAT, dataspace);
      dataSet.write(value,H5::PredType::NATIVE_FLOAT,dataspace);
#endif
      H5DSset_label(dataSet.getId(), 0, "Z");
      H5DSset_label(dataSet.getId(), 1, "Y");
      H5DSset_label(dataSet.getId(), 2, "X");
      dataSet.close();
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
    file.close();
  }
  /**
   * @brief dumps the morphology file after conversion to ZYX order if required
   * @param [in] fname file name
   * @param [in] inputData input data
   * @param [in] morphologyData morphology data
   */
  void writeMorphologyFile(const std::string &fname, const InputData &inputData, const Voxel *morphologyData, const int NUM_MATERIAL) {
    std::array<std::string, 4> morphologyDataSets{};
    const MorphologyType morphologyType = static_cast<MorphologyType>(inputData.morphologyType);
    if (morphologyType == MorphologyType::VECTOR_MORPHOLOGY) {
      morphologyDataSets[0] = "Sx";
      morphologyDataSets[1] = "Sy";
      morphologyDataSets[2] = "Sz";
      morphologyDataSets[3] = "UnalignedFraction";
    } else if (morphologyType == MorphologyType::EULER_ANGLES) {
      morphologyDataSets[0] = "S";
      morphologyDataSets[1] = "Theta";
      morphologyDataSets[2] = "Psi";
      morphologyDataSets[3] = "Vfrac";
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
        const auto &group = file.createGroup("Material_" + std::to_string(numMat + 1));
        for (int numComponent = 0; numComponent < 4; numComponent++) {
          for (int numVoxel = 0; numVoxel < numVoxels; numVoxel++) {
            int offset = numMat*numVoxels + numVoxel ;
            data[numVoxel] = morphologyData[offset].getValueAt(numComponent);
          }
          H5::DataSpace dataspace(RANK, dims);
#ifdef DOUBLE_PRECISION
          H5::DataSet dataSet = group.createDataSet(morphologyDataSets[numComponent].c_str(), H5::PredType::NATIVE_DOUBLE, dataspace);
          dataSet.write(data, H5::PredType::NATIVE_DOUBLE);
#else
          H5::DataSet dataSet = group.createDataSet(morphologyDataSets[numComponent].c_str(),
                                                    H5::PredType::NATIVE_FLOAT, dataspace);
          dataSet.write(data, H5::PredType::NATIVE_FLOAT);
#endif
          H5DSset_label(dataSet.getId(),0,"Z");
          H5DSset_label(dataSet.getId(),1,"Y");
          H5DSset_label(dataSet.getId(),2,"X");
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

  /**
   * @brief Writes the XDMF file for loading HDF into Paraview / Visit
   * @param [in] inputData input data
   * @param [in] voxelData voxel data
   */
  void writeXDMF(const InputData &inputData, const Voxel * voxelData) {
    const int & NUM_MATERIAL = inputData.NUM_MATERIAL;
    std::array<std::string, 4> morphologyDataSets{};
    const MorphologyType morphologyType = static_cast<MorphologyType>(inputData.morphologyType);
    if (morphologyType == MorphologyType::VECTOR_MORPHOLOGY) {
      morphologyDataSets[0] = "Sx";
      morphologyDataSets[1] = "Sy";
      morphologyDataSets[2] = "Sz";
      morphologyDataSets[3] = "UnalignedFraction";
    } else if (morphologyType == MorphologyType::EULER_ANGLES) {
      morphologyDataSets[0] = "S";
      morphologyDataSets[1] = "Theta";
      morphologyDataSets[2] = "Psi";
      morphologyDataSets[3] = "Vfrac";
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
            "\t\t<DataItem Name=\"spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\" %d\" Format=\"XML\">\n",
            precision);
    fprintf(fp, "\t\t1.0 1.0 1.0\n");
    fprintf(fp, "\t\t</DataItem>\n");
    fprintf(fp, "\t\t</Geometry>\n");

    for (int numComponent = 0; numComponent < 4; numComponent++) {
      for (int i = 1; i < NUM_MATERIAL + 1; i++) {
        fprintf(fp, "\t\t<Attribute Name=\"Mat_%d_%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n", i,morphologyDataSets[numComponent].c_str());
        fprintf(fp, "\t\t<DataItem Dimensions=\" %d %d %d \" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",
                voxelSize[0], voxelSize[1], voxelSize[2], precision);
        fprintf(fp, "\t\t %s.h5:/Material_%d/%s\n", fName.c_str(),i,morphologyDataSets[numComponent].c_str());
        fprintf(fp,"\t\t</DataItem>\n");
        fprintf(fp,"\t\t</Attribute>\n");
      }
    }
    fprintf(fp,"\t</Grid>\n");
    fprintf(fp,"\t</Domain>\n");
    fprintf(fp,"</Xdmf>\n");
    fclose(fp);
    const std::string h5MorphologyName = dirName+"/"+fName;
    writeMorphologyFile(h5MorphologyName,inputData,voxelData,NUM_MATERIAL);
  }
}


#endif //CY_RSOXS_WRITEH5_H
