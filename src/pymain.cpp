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


#include <pybind11/pybind11.h>
#include <Input/InputData.h>
#include <Output/writeH5.h>
#include <pybind11/stl.h>
#include <cudaMain.h>
#include <complex>

#include "version.h"
#include <Input/readH5.h>
#include <iomanip>
#include <pybind11/iostream.h>
#include <utils.h>
#include <PyClass/VoxelData.h>
#include <PyClass/RefractiveIndex.h>
#include <PyClass/ScatteringPattern.h>
namespace py = pybind11;


/**
 * Returns the number of materials the code was compiled with
 *
 * @return number of materials
 */
int getNumMaterials() {
    return NUM_MATERIAL;
}




/**
 * @brief Launch the GPU kernel for CyRSoXS.
 * @param inputData InputData
 * @param energyData Energy data
 * @param voxelData Voxel data
 */
void  launch(const InputData &inputData, const RefractiveIndexData &energyData,
            const VoxelData &voxelData,ScatteringPattern & scatteringPattern,
            bool ifWriteMetadata = true) {
  if (not(inputData.validate())) {
    py::print("Issues with Input Data");
    return ;
  }

  if (not(energyData.validate())) {
    py::print("Issues with Energy data");
    return ;
  }

  if (not(voxelData.validate())) {
    py::print("Issues with Voxel Data input");
    return ;
  }


  py::gil_scoped_release release;
  static bool printCopyInfo = true;
  if(printCopyInfo) {
      printCopyrightInfo();
      std::cout << "\n\n[INFO] Additional Cy-RSoXS Details: \n";
      std::cout << "[INFO] Version   = " << VERSION_MAJOR << "."<< VERSION_MINOR << "."<< VERSION_PATCH << "\n";
      std::cout << "[INFO] Git patch = " << GIT_HASH << "\n";
      printCopyInfo = false;
  }



  std::cout << "\n [STAT] Executing: \n\n";
  const UINT voxelDimensions[3]{inputData.numX, inputData.numY, inputData.numZ};
  cudaMain(voxelDimensions, inputData, energyData.getRefractiveIndexData(), scatteringPattern.data(), voxelData.data());
  if(ifWriteMetadata) {
      printMetaData(inputData);
  }


  py::gil_scoped_acquire acquire;

}

void cleanup(RefractiveIndexData &energyData, VoxelData &voxelData, ScatteringPattern & scatteringPattern) {
  voxelData.clear();
  energyData.clear();
  scatteringPattern.clear();

}

PYBIND11_MODULE(CyRSoXS, module) {
  module.doc() = "pybind11  plugin for CyRSoXS";

  py::print("CyRSoXS");
  py::print("============================================================================");
  py::print("Number of materials : ", NUM_MATERIAL);
  py::print("Size of Real        :", sizeof(Real));
  py::add_ostream_redirect(module, "ostream_redirect");
  py::enum_<Interpolation::EwaldsInterpolation>(module, "InterpolationType")
      .value("Linear", Interpolation::EwaldsInterpolation::LINEAR)
      .value("NearestNeighour", Interpolation::EwaldsInterpolation::NEARESTNEIGHBOUR)
      .export_values();

  py::enum_<FFT::FFTWindowing>(module, "FFTWindowing")
      .value("NoPadding", FFT::FFTWindowing::NONE)
      .value("Hanning", FFT::FFTWindowing::HANNING)
      .export_values();

  py::enum_<ScatterApproach>(module, "ScatterApproach")
      .value("Partial", ScatterApproach::PARTIAL)
      .value("Full", ScatterApproach::FULL)
      .export_values();

  py::enum_<MorphologyType>(module, "MorphologyType")
    .value("EulerAngles", MorphologyType::EULER_ANGLES)
    .value("VectorMorphology", MorphologyType::VECTOR_MORPHOLOGY)
    .export_values();

  py::enum_<CaseTypes>(module,"CaseType")
    .value("Default", CaseTypes::DEFAULT)
    .value("BeamDivergence", CaseTypes::BEAM_DIVERGENCE)
    .value("GrazingIncidence", CaseTypes::GRAZING_INCIDENCE)
    .export_values();


  py::class_<InputData>(module, "InputData")
      .def(py::init<>())
      .def("setEnergies", &InputData::setEnergies, "Set the energy data", py::arg("energies"))
      .def("print", &InputData::print, "Print the input data")
      .def("setERotationAngle", &InputData::setEAngles, "Set the rotation for Electric field", py::arg("StartAngle"),
           py::arg("EndAngle"),py::arg("IncrementAngle"))
      .def("setPhysSize", &InputData::setPhysSize, "Set the Physical size (in nm)", py::arg("PhysSize"))
      .def("setDimensions", &InputData::setDimension, "Set the Dimensions", py::arg("X"), py::arg("Y"), py::arg("Z"))
      .def("validate", &InputData::validate, "Validate the input data")
      .def("setCaseType",&InputData::setCaseType,"case Type")
      .def("setMorphologyType",&InputData::setMorphologyType,"morphology Type")
      .def("setKVectors",&InputData::setKVectors,"set K vectors")
      .def_readwrite("interpolationType", &InputData::ewaldsInterpolation, "Ewalds interpolation type")
      .def_readwrite("windowingType", &InputData::windowingType, "Windowing type")
      .def_readwrite("rotMask",&InputData::rotMask,"Rotation Mask")
      .def_readwrite("openMP", &InputData::num_threads, "number of OpenMP threads")
      .def_readwrite("scatterApproach", &InputData::scatterApproach, "sets the scatter approach");


  py::class_<RefractiveIndexData>(module, "RefractiveIndex")
      .def(py::init<const InputData &>(), "Constructor")
      .def("addData", &RefractiveIndexData::addData, "Add Optical constants data ", py::arg("OpticalConstants"),
           py::arg("Energy"))
      .def("validate", &RefractiveIndexData::validate, "Validate the refractive Index data")
      .def("clear",&RefractiveIndexData::clear,"Clears the memory")
      .def("print", &RefractiveIndexData::printEnergyData, "Prints the refractive index data");

  py::class_<VoxelData>(module, "VoxelData")
      .def(py::init<const InputData &>(), "Constructor",py::arg("InputData"))
      .def("addVoxelData", &VoxelData::addMaterialDataVectorMorphology,
           "Adds the Allignment and unaligned component to the given material", py::arg("AlignedData"),
           py::arg("UnalignedData"), py::arg("MaterialID"))
      .def("addVoxelData", &VoxelData::addMaterialDataEulerAngles,
         "Adds the Allignment and unaligned component to the given material", py::arg("S"),
         py::arg("Phi"),py::arg("Theta"), py::arg("vFrac"),py::arg("MaterialID"))
      .def("clear", &VoxelData::clear,"Clears the voxel data")
      .def("validate", &VoxelData::validate,"validate the Voxel data")
      .def("readFromH5", &VoxelData::readFromH5, "Reads from HDF5",py::arg("Filename"))
      .def("reset", &VoxelData::reset,"Resets the voxel data")
      .def("writeToH5", &VoxelData::writeToH5,"Writes voxel data to HDF5 file");

  py::class_<ScatteringPattern>(module,"ScatteringPattern")
      .def(py::init<const InputData &>(), "Constructor",py::arg("InputData"))
      .def("clear",&ScatteringPattern::clear,"Clears the memory")
      .def("writeToHDF5",&ScatteringPattern::writeToHDF5,"Dumps data in  HDF5 file format")
      .def("writeToVTI",&ScatteringPattern::writeToVTI,"Dumps data in  VTI file format")
      .def("dataToNumpy",&ScatteringPattern::writeToNumpy,"Returns data in numpy array",py::arg("Energy"),py::arg("kID"));

  module.def("launch", &launch, "GPU computation", py::arg("InputData"), py::arg("RefractiveIndexData"),
             py::arg("VoxelData"),py::arg("ScatteringPattern"),py::arg("WriteMetaData")=true);
  module.def("cleanup", &cleanup, "Cleanup",  py::arg("RefractiveIndex"), py::arg("VoxelData"),py::arg("ScatteringPattern"));
  module.def("getNumMaterials", &getNumMaterials, "Get number of materials");


}