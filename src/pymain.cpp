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
#include <pybind11/numpy.h>
#include <omp.h>
#include <Input/readH5.h>
#include <iomanip>
#include <pybind11/iostream.h>
#include <utils.h>

namespace py = pybind11;

class EnergyData {
    const Real &energyStart;
    const Real &energyEnd;
    const Real &incrementEnergy;
    std::vector<Material<NUM_MATERIAL> > materialInput;
    std::vector<bool> isValid;
public:
    EnergyData(const InputData &inputData)
        : energyStart(inputData.energyStart), energyEnd(inputData.energyStart),
          incrementEnergy(inputData.incrementEnergy) {
      UINT numEnergyData = std::round(inputData.energyEnd - inputData.energyStart) / (inputData.incrementEnergy) + 1;
      materialInput.resize(numEnergyData);
      isValid.resize(numEnergyData);
      std::fill(isValid.begin(), isValid.end(), false);
    }

    void clear() {
      std::vector<Material<NUM_MATERIAL> > _materialInput;
      std::swap(_materialInput, materialInput);
    }


    void addData(const std::vector<std::vector<Real>> &values, const Real Energy) {
      enum EnergyValues : u_short {
          DeltaPara = 0,
          BetaPara = 1,
          DeltaPerp = 2,
          BetaPerp = 3
      };
      if (not(values.size() == NUM_MATERIAL)) {
        py::print("Wrong input for Energy. Number not matching with number of Materials");
        return;
      }
      for (auto &value:values) {
        if ((value.size() != 4)) {
          py::print("Wrong number of input parameters. Parameters must be in the order of "
                    "(DeltaPara, BetaPara, DeltaPerp, , BetaPerp)");
          return;
        }
      }
      UINT counter = std::round((Energy - energyStart) / incrementEnergy);
      for (UINT i = 0; i < NUM_MATERIAL; i++) {
        materialInput[counter].npara[i].x = 1 - values[i][EnergyValues::DeltaPara];
        materialInput[counter].npara[i].y = values[i][EnergyValues::BetaPara];
        materialInput[counter].nperp[i].x = 1 - values[i][EnergyValues::DeltaPerp];
        materialInput[counter].nperp[i].y = values[i][EnergyValues::BetaPerp];
      }

      isValid[counter] = true;

    }

    void printEnergyData() const {
      UINT count = 0;
      for (auto &values : materialInput) {
        Real currEnegy = energyStart + count * incrementEnergy;
        py::print("Energy = ", currEnegy);
        for (int i = 0; i < NUM_MATERIAL; i++) {
          py::print("Material = ", i, "npara = ", std::complex<Real>(values.npara[i].x, values.npara[i].y),
                    "nperp = ", std::complex<Real>(values.nperp[i].x, values.nperp[i].y));
        }
        count++;

      }
    }

    const std::vector<Material<NUM_MATERIAL>> &getEnergyData() const {
      return materialInput;
    }

    bool validate() const {
      if(std::all_of(isValid.begin(), isValid.end(), [](bool x) { return (x == true); })){
        return true;
      }
      else{
        for(UINT i = 0; i < isValid.size(); i++){
          if(not(isValid[i])){
            Real energy = energyStart + i * incrementEnergy;
            py::print("Optical constants data missing for energy = ",energy,"eV");
          }
        }
      }
      return false;
    }
};


class VoxelData {
private:
    Voxel<NUM_MATERIAL> *voxel = nullptr;
    const InputData &inputData;
    std::bitset<NUM_MATERIAL> validData_;
public:
    VoxelData(const InputData &_inputData)
        : inputData(_inputData) {
      clear();
      const BigUINT numVoxels = inputData.numX * inputData.numY * inputData.numZ;
      voxel = new Voxel<NUM_MATERIAL>[numVoxels];
      validData_.reset();
    }

    void addMaterialData(py::array_t<Real, py::array::c_style | py::array::forcecast> &matAlignementData,
                         py::array_t<Real, py::array::c_style | py::array::forcecast> &matUnalignedData,
                         const UINT matID) {
      if (matID >= NUM_MATERIAL) {
        throw std::logic_error("Number of material does not match with the compiled version");
      }
      if(validData_.test(matID)){
        py::print("The material is already set. Please first reset to add the entries. Returning.");
        return;
      }
      const BigUINT numVoxels = inputData.numX * inputData.numY * inputData.numZ;
      for (BigUINT i = 0; i < numVoxels; i++) {
        voxel[i].s1[matID].x = matAlignementData.data()[i * 3 + 0];
        voxel[i].s1[matID].y = matAlignementData.data()[i * 3 + 1];
        voxel[i].s1[matID].z = matAlignementData.data()[i * 3 + 2];
        std::cout << i <<  " " << numVoxels << "\n";
//        voxel[i].s1[matID].w = matUnalignedData.data()[i];
      }
      validData_.set(matID, true);
    }

    void readFromH5(const std::string &fname) {
      if(validData_.any()){
        py::print("Some of the material is already set. Please first reset to continue. Returning.");
        return;
      }
      const UINT voxelDim[3]{inputData.numX, inputData.numY, inputData.numZ};
      H5::readFile(fname, voxelDim, voxel, true);
      validData_.flip();
    }

    void reset(){
      validData_.reset();
    }
    void writeToH5() const {
      const BigUINT numVoxels = inputData.numX * inputData.numY * inputData.numZ;
      const UINT dim[3]{inputData.numX, inputData.numY, inputData.numZ};
      Real *scalarValues = new Real[numVoxels];
      for (int nMat = 0; nMat < NUM_MATERIAL; nMat++) {
        std::string fname = "Unalligned_" + std::to_string(nMat);
        for (int i = 0; i < numVoxels; i++) {
          scalarValues[i] = voxel[i].s1[nMat].w;
        }
        H5::writeFile3DScalar(fname, scalarValues, dim, "unalligned");
      }
      delete[] scalarValues;

      Real *vectorValues = new Real[numVoxels * 3];
      for (int nMat = 0; nMat < NUM_MATERIAL; nMat++) {
        std::string fname = "Alligned_" + std::to_string(nMat);
        for (int i = 0; i < numVoxels; i++) {
          vectorValues[i * 3 + 0] = voxel[i].s1[nMat].x;
          vectorValues[i * 3 + 1] = voxel[i].s1[nMat].y;
          vectorValues[i * 3 + 2] = voxel[i].s1[nMat].z;
        }
        H5::writeFile3DVector(fname, vectorValues, dim, "S");
      }
      delete[] vectorValues;
    }

    void clear() {
      if (voxel != nullptr) {
        delete[] voxel;
      }
      voxel = nullptr;
    }

    ~VoxelData() {
      if (voxel != nullptr) {
        delete[] voxel;
      }
    }

    const Voxel<NUM_MATERIAL> *data() const {
      return voxel;
    }

    bool validate() const {
      if (validData_.all()){
        return true;
      }
      else{
        for(UINT i = 0; i < validData_.size(); i++){
          if(not(validData_.test(i))){
            py::print("Voxel Data missing / corrupt for material = ",i);
          }
        }
      }
      return false;
    }

};

void launch(const InputData  &inputData, const EnergyData &energyData,
            const  VoxelData &voxelData) {



  if (not(inputData.validate())) {
    py::print("Issues with Input Data");
    return;
  }

  if (not(energyData.validate())) {
    py::print("Issues with Energy data");
    return;
  }

  if (not(voxelData.validate())) {
    py::print("Issues with Voxel Data input");
    return;
  }



  py::gil_scoped_release release;
  printCopyrightInfo();
  std::cout << "\n\n----------- Executing:  -----------------\n\n";
  const UINT voxelDimensions[3]{inputData.numX, inputData.numY, inputData.numZ};
  Real *projectionAveraged;
  cudaMain(voxelDimensions, inputData, energyData.getEnergyData(), projectionAveraged, voxelData.data());

  createDirectory("HDF5");
  omp_set_num_threads(1);
  const UINT numEnergyLevel =
      static_cast<UINT>(std::round(
          (inputData.energyEnd - inputData.energyStart) / inputData.incrementEnergy + 1));
  const UINT voxel2DSize = voxelDimensions[0] * voxelDimensions[1];
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
    std::memcpy(oneEnergyData, &projectionAveraged[csize * voxel2DSize], sizeof(Real) * voxel2DSize);
    const std::string fname = "HDF5/Energy_" + s;

    H5::writeFile2D(fname, oneEnergyData, voxelDimensions);
  }
  delete[] oneEnergyData;
  delete[] projectionAveraged;
  std::ofstream file("metadata.txt");
  file << "---------------- Scaling Information--------------\n";
  file << "Number of pixel = [" << inputData.numX << "," << inputData.numY << "]\n";
  file << "Q range  = [" << -M_PI / inputData.physSize << "," << M_PI / inputData.physSize << "]\n";
  file << "Electric field rotated through = [" << inputData.startAngle << "," << inputData.endAngle << "]\n";
  file << "Increments in electric field rotation " << inputData.incrementEnergy << "\n";
  file << "\n\n";

  file << "-----------------Simulation information -----------------\n";
  file << "Size of Real = " << sizeof(Real) << "\n";
  file << "Number of materials =" << NUM_MATERIAL << "\n";
  file << "Energies simulated from " << inputData.energyStart << " to " << inputData.energyEnd << " with increment of "
       << inputData.incrementEnergy << "\n";
  file.close();
#pragma omp barrier
  std::cout << "Execution finished \n";
  py::gil_scoped_acquire acquire;
}

void cleanup(InputData &inputData, EnergyData &energyData, VoxelData &voxelData) {
  voxelData.clear();
  energyData.clear();

}

PYBIND11_MODULE(CyRSoXS, module) {
  module.doc() = "pybind11  plugin for Cy-RSoXS";
  py::print("----------------Compile time options-------------------");
  py::print("Number of materials : ", NUM_MATERIAL);
  py::print("Size of Real", sizeof(Real));
  py::add_ostream_redirect(module, "ostream_redirect");
#ifdef ENABLE_2D
  py::print("Enable 2D : True"  );
#else
  py::print("Enable 2D : False");
#endif
  py::enum_<Interpolation::EwaldsInterpolation>(module, "InterpolationType")
      .value("Linear", Interpolation::EwaldsInterpolation::LINEAR)
      .value("NearestNeighour", Interpolation::EwaldsInterpolation::NEARESTNEIGHBOUR)
      .export_values();

  py::enum_<FFT::FFTWindowing>(module, "FFTWindowing")
      .value("NoPadding", FFT::FFTWindowing::NONE)
      .value("Hanning", FFT::FFTWindowing::HANNING)
      .export_values();

  py::class_<InputData>(module, "InputData")
      .def(py::init<>())
      .def("setEnergy", &InputData::setEnergy, py::arg("StartEnergy"), py::arg("EndEnergy"), py::arg("IncrementEnergy"))
      .def("print", &InputData::print)
      .def("setRotationAngle", &InputData::setAngles, py::arg("StartAngle"), py::arg("EndAngle"),
           py::arg("IncrementAngle"))
      .def("physSize", &InputData::setPhysSize, py::arg("PhysSize"))
      .def("dimensions", &InputData::setDimension, py::arg("X"), py::arg("Y"), py::arg("Z"))
      .def("validate", &InputData::validate)
      .def_readwrite("interpolationType", &InputData::ewaldsInterpolation)
      .def_readwrite("windowingType", &InputData::windowingType)
      .def_readwrite("writeVTI",&InputData::writeVTI)
      .def_readwrite("openMP",&InputData::num_threads);

  py::class_<EnergyData>(module, "RefractiveIndex")
      .def(py::init<const InputData &>())
      .def("addData", &EnergyData::addData,py::arg("OpticalConstants"),py::arg("Energy"))
      .def("validate",&EnergyData::validate)
      .def("print", &EnergyData::printEnergyData);

  py::class_<VoxelData>(module, "VoxelData")
      .def(py::init<const InputData &>())
      .def("addVoxelData", &VoxelData::addMaterialData,py::arg("AlignedData"),py::arg("UnalignedData"),py::arg("MaterialID"))
      .def("clear", &VoxelData::clear)
      .def("validate",&VoxelData::validate)
      .def("readFromH5", &VoxelData::readFromH5,py::arg("Filename"))
      .def("reset",&VoxelData::reset)
      .def("writeToH5", &VoxelData::writeToH5);


  module.def("launch", &launch, "GPU computation",py::arg("InputData"),py::arg("RefractiveIndexData"),py::arg("VoxelData"));
  module.def("cleanup", &cleanup, "Cleanup",py::arg("InputData"),py::arg("RefractiveIndex"),py::arg("VoxelData"));


}