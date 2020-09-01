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
#include <Input/readH5.h>
#include <iomanip>
#include <pybind11/iostream.h>

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
        std::fill(isValid.begin(),isValid.end(),false);
    }

    void clear(){
        std::vector<Material<NUM_MATERIAL> > _materialInput;
        std::swap(_materialInput,materialInput);
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
        for(auto & value:values) {
            if ((value.size() != 4)) {
                py::print("Wrong number of input parameters. Parameters must be in the order of "
                          "(DeltaPara, BetaPara, DeltaPerp, , BetaPerp)");
                return;
            }
        }
        UINT counter = std::round((Energy - energyStart)/incrementEnergy);
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
    const std::vector<Material<NUM_MATERIAL>>  & getEnergyData() const{
        return materialInput;
    }
    bool validate () const{
        return std::all_of(isValid.begin(),isValid.end(),[](bool x) {return (x == true); });
    }
};


class VoxelData {
private:
    enum VoxelStructure:u_short {
        UnalignedData = 0,
        AlignedData = 1,
        MAX = 2
    };
    Voxel<NUM_MATERIAL> *voxel = nullptr;
    const InputData &inputData;
    std::bitset<NUM_MATERIAL*(VoxelStructure::MAX)> validData_;
public:
    VoxelData(const InputData &_inputData)
            : inputData(_inputData) {
        clear();
        const BigUINT numVoxels = inputData.numX * inputData.numY * inputData.numZ;
        voxel = new Voxel<NUM_MATERIAL>[numVoxels];
        validData_.reset();
    }

    void addMatAllignment(py::array_t<Real, py::array::c_style | py::array::forcecast> &array, const UINT matID) {
        if(matID >= NUM_MATERIAL) {
            throw std::logic_error("Number of material does not match with the compiled version");
        }
        const BigUINT numVoxels = inputData.numX * inputData.numY * inputData.numZ;
        for (BigUINT i = 0; i < numVoxels; i++) {
            voxel[i].s1[matID].x = array.data()[i * 3 + 0];
            voxel[i].s1[matID].y = array.data()[i * 3 + 1];
            voxel[i].s1[matID].z = array.data()[i * 3 + 2];
        }
        validData_.set(matID*(VoxelStructure::MAX) + VoxelStructure::AlignedData,true);
    }

    void addMatUnAlligned(py::array_t<Real, py::array::c_style | py::array::forcecast> &array, const UINT matID) {
        if(matID >= NUM_MATERIAL) {
            throw std::logic_error("Number of material does not match with the compiled version");
        }
        const BigUINT numVoxels = inputData.numX * inputData.numY * inputData.numZ;
        for (BigUINT i = 0; i < numVoxels; i++) {
            voxel[i].s1[matID].w = array.data()[i];
        }
        validData_.set(matID*(VoxelStructure::MAX) + VoxelStructure::UnalignedData,true);
    }

    void readFromH5(const std::string& fname){
        const UINT voxelDim[3]{inputData.numX,inputData.numY,inputData.numZ};
        H5::readFile(fname,voxelDim,voxel,true);
        validData_.flip();
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

    bool isValid() const{
        return (validData_.all());
    }

};

void launch(const VoxelData &voxelData, const EnergyData &energyData,
            const InputData &inputData) {

    if(not(inputData.validate())){
        py::print("Issues with Input Data");
        return;
    }

    if(not(energyData.validate())){
        py::print("Issues with Energy data");
        return;
    }

    if(not(voxelData.isValid())){
        py::print("Issues with Voxel Data input");
        return;
    }
    py::print("--------------- Input Data ---------------------");
    inputData.print();

    py::print("\n\n--------------- Energy Data ---------------------");
    energyData.printEnergyData();

    py::print("\n\n--------------- Execution Begins ---------------------");
    py::gil_scoped_release release;
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
    delete [] projectionAveraged;
    std::ofstream file("metadata.txt");
    file << "---------------- Scaling Information--------------\n";
    file << "Number of pixel = [" << inputData.numX << "," << inputData.numY << "]\n";
    file << "Q range  = [" << -M_PI/inputData.physSize << "," << M_PI/inputData.physSize << "]\n";
    file << "Electric field rotated through = [" << inputData.startAngle << "," << inputData.endAngle << "]\n";
    file << "Increments in electric field rotation " << inputData.incrementEnergy << "\n";
    file << "\n\n";

    file << "-----------------Simulation information -----------------\n";
    file << "Size of Real = " << sizeof(Real) << "\n";
    file << "Number of materials =" << NUM_MATERIAL << "\n";
    file << "Energies simulated from " << inputData.energyStart << " to " << inputData.energyEnd << " with increment of " << inputData.incrementEnergy << "\n";
    file.close();
#pragma omp barrier
    std::cout << "Execution finished \n";
    py::gil_scoped_acquire acquire;
}

void cleanup(InputData & inputData,  EnergyData &energyData, VoxelData & voxelData){
    voxelData.clear();
    energyData.clear();

}

PYBIND11_MODULE(CyRSoXS, module) {
    module.doc() = "pybind11  plugin for Cy-RSoXS";
    py::print("Credits: ISU");
    py::print("----------------Compile time options-------------------");
    py::print("Number of materials : ", NUM_MATERIAL);
    py::print("Size of Real",sizeof(Real));

#ifdef ENABLE_2D
    py::print("Enable 2D : True"  );
#else
    py::print("Enable 2D : False"  );
#endif

    py::class_<InputData>(module, "InputData")
            .def(py::init<>())
            .def("setEnergy", &InputData::setEnergy)
            .def("print", &InputData::print)
            .def("setRotationAngle", &InputData::setAngles)
            .def("physSize", &InputData::setPhysSize)
            .def("dimensions", &InputData::setDimension)
            .def("validate", &InputData::validate)
            .def_readwrite("interpolationType",&InputData::ewaldsInterpolation)
            .def_readwrite("windowingType",&InputData::windowingType);

    py::class_<EnergyData>(module, "MaterialProperties")
            .def(py::init<const InputData &>())
            .def("addData", &EnergyData::addData)
            .def("print", &EnergyData::printEnergyData);

    py::class_<VoxelData>(module, "VoxelData")
            .def(py::init<const InputData &>())
            .def("addMatAllignment", &VoxelData::addMatAllignment)
            .def("addMatUnalligned", &VoxelData::addMatUnAlligned)
            .def("clear",&VoxelData::clear)
            .def("readFromH5", &VoxelData::readFromH5)
            .def("writeToH5", &VoxelData::writeToH5);


    module.def("launch", &launch, "GPU computation");
    module.def("cleanup", &cleanup, "Cleanup");
    py::add_ostream_redirect(module, "ostream_redirect");

}