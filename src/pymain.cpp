//
// Created by maksbh on 8/25/20.
//


#include <pybind11/pybind11.h>
#include <Input/InputData.h>
#include <Output/writeH5.h>
#include <pybind11/stl.h>
#include <cudaMain.h>
#include <complex>
#include <pybind11/numpy.h>
#include <omp.h>
#include <iomanip>


namespace py = pybind11;

struct EnergyData {
    const Real &energyStart;
    const Real &energyEnd;
    const Real &incrementEnergy;
    UINT counter = 0;
    std::vector<Material<NUM_MATERIAL> > materialInput;
public:
    EnergyData(const InputData &inputData)
            : energyStart(inputData.energyStart), energyEnd(inputData.energyStart),
              incrementEnergy(inputData.incrementEnergy) {
        UINT numEnergyData = std::round(inputData.energyEnd - inputData.energyStart) / (inputData.incrementEnergy);
        materialInput.resize(numEnergyData);
    }

    void addData(const std::vector<std::vector<Real>> &values, const Real Energy) {

        enum EnergyValues : u_short {
            DeltaPara = 0,
            DeltaPerp = 1,
            BetaPara = 2,
            BetaPerp = 3
        };
        if (not(values.size() == NUM_MATERIAL)) {
            throw std::runtime_error("Wrong input for Energy. Number not matching with number of Materials");
        }
        if (not(values[0].size() == 4)) {
            throw std::runtime_error("Wrong number of input parameters. Parameters must be in the order of "
                                     "(DeltaPara, DeltaPerp, BetaPara, BetaPerp)");
        }
        Real reqEnergy = energyStart + incrementEnergy * counter;
        if (fabs(reqEnergy - Energy) > 1E-5) {
            throw std::runtime_error("Energy input missing for " + std::to_string(reqEnergy) +
                                     ". Input must be in ascending order of energies");
        }
        for (UINT i = 0; i < NUM_MATERIAL; i++) {
            materialInput[counter].npara[i].x = 1 - values[i][EnergyValues::DeltaPara];
            materialInput[counter].npara[i].y = values[i][EnergyValues::BetaPara];

            materialInput[counter].nperp[i].x = 1 - values[i][EnergyValues::DeltaPerp];
            materialInput[counter].nperp[i].y = values[i][EnergyValues::BetaPerp];
        }
    }

    void printEnergyData() const {
        UINT count = 0;
        for (auto &values : materialInput) {
            Real currEnegy = energyStart + count * incrementEnergy;
            py::print("Energy = ", currEnegy);
            const Material<NUM_MATERIAL> &mat = materialInput[count];
            for (int i = 0; i < NUM_MATERIAL; i++) {
                py::print("Material = ", i, "npara = ", std::complex<Real>(mat.npara[i].x, mat.npara[i].y),
                          "nperp = ", std::complex<Real>(mat.nperp[i].x, mat.nperp[i].y));
            }

        }
    }
};


class VoxelData {
private:
    Voxel<NUM_MATERIAL> *voxel = nullptr;
    const InputData &inputData;
public:
    VoxelData(const InputData &_inputData)
            : inputData(_inputData) {
        clear();
        const BigUINT numVoxels = inputData.numX * inputData.numY * inputData.numZ;
        voxel = new Voxel<NUM_MATERIAL>[numVoxels];
    }

    void addMatAllignment(py::array_t<Real, py::array::c_style | py::array::forcecast> &array, const UINT matID) {
        const BigUINT numVoxels = inputData.numX * inputData.numY * inputData.numZ;
        for (BigUINT i = 0; i < numVoxels; i++) {
            voxel[i].s1[matID].x = array.data()[i * 3 + 0];
            voxel[i].s1[matID].y = array.data()[i * 3 + 1];
            voxel[i].s1[matID].z = array.data()[i * 3 + 2];
        }
    }

    void addMatUnAlligned(py::array_t<Real, py::array::c_style | py::array::forcecast> &array, const UINT matID) {
        const BigUINT numVoxels = inputData.numX * inputData.numY * inputData.numZ;
        for (BigUINT i = 0; i < numVoxels; i++) {
            voxel[i].s1[matID].w = array.data()[i];
        }
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
        delete[] voxel;
        voxel = nullptr;
    }

    ~VoxelData() {
        if (voxel != nullptr) {
            delete[] voxel;
        }
    }

    const Voxel<NUM_MATERIAL> *data() const {
        if (voxel != nullptr) {
            delete[] voxel;
        }
    }

};

void launch(const VoxelData &voxelData, const EnergyData &energyData, Real *&projectionGPUAveraged,
            const InputData &inputData) {
    py::gil_scoped_release release;
    const UINT voxelDimensions[3]{inputData.numX, inputData.numY, inputData.numZ};
    Real *projectionAveraged;
    cudaMain(voxelDimensions, inputData, energyData.materialInput, projectionAveraged, voxelData.data());

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
        std::memcpy(oneEnergyData, &projectionGPUAveraged[csize * voxel2DSize], sizeof(Real) * voxel2DSize);
        const std::string fname = "HDF5/Energy_" + s;

        H5::writeFile2D(fname, oneEnergyData, voxelDimensions);
    }
    delete[] oneEnergyData;
    py::gil_scoped_acquire acquire;
    py::print("Finished Execution\n");
}

PYBIND11_MODULE(CyRSoXS, module) {
    module.doc() = "pybind11  plugin for Cy-RSoXS";
    py::class_<InputData>(module, "InputData")
            .def(py::init<>())
            .def("setEnergy", &InputData::setEnergy)
            .def("print", &InputData::print)
            .def("setRotationAngle", &InputData::setAngles)
            .def("physSize", &InputData::setPhysSize)
            .def("dimensions", &InputData::setDimension)
            .def("Validate", &InputData::validate);

    py::class_<EnergyData>(module, "MaterialProperties")
            .def(py::init<const InputData &>())
            .def("addData", &EnergyData::addData)
            .def("print", &EnergyData::printEnergyData);

    py::class_<VoxelData>(module, "VoxelData")
            .def(py::init<const InputData &>())
            .def("addMatAllignment", &VoxelData::addMatAllignment)
            .def("addMatUnalligned", &VoxelData::addMatUnAlligned)
            .def("clear",&VoxelData::clear)
            .def("writeToH5", &VoxelData::writeToH5);


    module.def("launch", &launch, "GPU computation");

}