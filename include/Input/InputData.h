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

#ifndef PRS_READCONFIG_H
#define PRS_READCONFIG_H
#include <string>
#ifndef PYBIND
#include <libconfig.h++>
#include <array>
#endif
#include <iostream>
#include <Input/Input.h>
#include <vector>
#include <cmath>
#include <fstream>
#include <Rotation.h>
#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <bitset>
#include <algorithm>

namespace ParamChecker{
    enum Parameters: u_short {
        ENERGY = 0,
        DIMENSION = 1,
        PHYSSIZE = 2,
        EANGLE = 3,
        KVECTORS = 4,
        MORPHOLOGY_TYPE = 5,
        CASE_TYPE = 6,
        DETECTOR_COORD = 7,
        NUM_MATERIAL = 8,
        MAX_SIZE = 9
    };
    static const char* paramNames[] = {"ENERGY","DIMENSION","PHYSSIZE","EANGLE","KVectors","MorphologyType","CaseType","DetectorCoord","NUM_MATERIAL"};
    static_assert(sizeof(ParamChecker::paramNames)/sizeof(char*) == ParamChecker::Parameters::MAX_SIZE,
            "sizes dont match");

}
#endif

/// This function reads the input data from the file.
class InputData {
private:
    bool enable2D_;
#ifndef PYBIND
  /**
   * This function reads the value that is compulsory to be read from
   * the input data file.
   * @tparam T
   * @param config
   * @param key
   * @param value
   */
  template <typename T>
  void ReadValueRequired(libconfig::Config & config, const std::string key, T & value ){
    bool res = config.lookupValue(key, value);
    if(res == false){
      std::cerr << "[Input Error] No value corresponding to " << key  << " found. Exiting\n";
      exit(EXIT_FAILURE);
    }
  }
  template<typename T,int size>
  static void ReadArray(const libconfig::Setting &setting,const std::string &key_name, std::array<T,size>  &value) {
    for (int i = 0; i < setting.getLength(); ++i) {
      value[i] = setting[i];
    }
  }
  /**
   *
   * @tparam T
   * @param config
   * @param key
   * @param value
   * @return
   */
  template <typename T>
  bool ReadValue(libconfig::Config & config, const std::string key, T & value ){
    bool res = config.lookupValue(key, value);
    if(res == false){
      std::cout << YLW<<  "[WARNING] : No value corresponding to " << key << " found. Setting to default" << NRM << "\n";
    }
    return res;
  }

  template <const char ** enumVarname,int maxEnums>
  static UINT convertStringToEnums(const char * stringName){
    for(int i = 0; i < maxEnums; i++){
      if(strcmp(stringName,enumVarname[i]) == 0){
        return i;
      }
    }
    throw std::logic_error("String did not match not found");
  }
    /**
   * Reads an array from the input file and stores it in the given vector.
   *
   * The given vector is resized to fit the data. Any existing data in the
   * vector is lost. Resizing will only happen if the input file has an array
   * of proper length associated with the given key.
   *
   * If the key is not found, an error will be printed and the program will
   * exit with failure.
   *
   * @param config libconfig object for file
   * @param key the name of the configuration value to retrieve
   * @param[out] arr vector to store the data (will be allocated here)
   */
    template <typename T>
    void ReadArrayRequired(libconfig::Config & config, const std::string& key,
                           std::vector<T>& arr, const UINT size = 0) {
        // this is here because lookup throws an exception if not found
        if (!config.exists(key)) {
            std::cerr << "[Input Error] No value corresponding to " << key
                      << " found. Exiting\n";
            exit(EXIT_FAILURE);
        }

        libconfig::Setting &setting = config.lookup(key);
        if (!setting.isArray()) {  // confirm this is an array
            std::cerr << "expected array input but found single value "
                      << "(key = " << key << ")\n";
            exit(EXIT_FAILURE);
        }
        if(size > 0){
            if(setting.getLength() != size){
                std::cout << RED << "[Input Error]: The array corresponding to " << key << " must be of length " << size
                << ". But found of length " << setting.getLength() << NRM <<"\n";
                exit(EXIT_FAILURE);
            }
        }
        arr.clear();
        for (int i = 0; i < setting.getLength(); i++) {
            T value;
            std::string item_path = key + ".[" + std::to_string(i) + "]";
            ReadValueRequired(config, item_path, value);
            arr.push_back(value);
        }
    }
      /**
   *
   * @tparam T
   * @param name name of the variable
   * @param val value
   * @param max max allowable for that type
   */
  template<typename T>
  void validate(const std::string & name, const T & val, const UINT & max) const{
    if(val >= max){
        std::cout << RED << "[Error] : " << name << "Wrong value. Max acceptable value " <<  max-1 << ". Value found =  " << val << NRM <<"\n";
        exit(EXIT_FAILURE);
    }
  }

#else
    /// Checks the input data parameter for the required parameters
    std::bitset<ParamChecker::Parameters::MAX_SIZE> paramChecker_;
#endif
  /**
   * @brief sets the flag for 2D computataion.
   */


 public:
  /// start of energy
  std::vector<Real> energies;
  /// startAngle
  Real startAngle;
  /// end Angle
  Real endAngle;
  /// increment in Angle
  Real incrementAngle;
  /// number of threads
  UINT num_threads = 4;
  /// Number of voxels in X direction.
  UINT voxelDims[3];
  /// Physical Size
  Real physSize;
  /// Write HDF5 file
  bool writeHDF5 = true;
  /// Write VTI file
  bool writeVTI = false;
  /// whether to do masking for rotation or not
  bool rotMask = false;
  /// Type of Ewalds interpolation
  int ewaldsInterpolation = Interpolation::EwaldsInterpolation::LINEAR;
  /// Windowing Type
  UINT windowingType = FFT::FFTWindowing::NONE;
  /// k
  std::vector<Real3> kVectors;

  /// scatter Approach
  UINT scatterApproach = ScatterApproach::PARTIAL;

  std::string VTIDirName = "VTI";
  std::string HDF5DirName = "HDF5";

  UINT caseType;
  UINT morphologyType;
  bool referenceFrame = ReferenceFrame::LAB;

  Real3 detectorCoordinates{0,0,1};
  int algorithmType = Algorithm::CommunicationMinimizing;
  int numMaxStreams = 1;


  MorphologyOrder morphologyOrder = MorphologyOrder::INVALID;

  bool dumpMorphology = false;

  int NUM_MATERIAL;

  /**
   *
   * @return gets the 2D computation flags
   */
  inline bool if2DComputation() const {
     return enable2D_;
  }
  inline void check2D() {
    if (voxelDims[2] != 1) {
      enable2D_ = false;
    } else {
      std::cout << "[INFO] 2D computation enabled\n";
      enable2D_ = true;
    }
  }


#ifndef PYBIND

  /// Morphology type
  /**
   * Constructor to read input data
   * @param filename filename, default is config.txt
   */
  InputData(std::string filename = "config.txt") {
    libconfig::Config cfg;
    try {
      cfg.readFile(filename.c_str());
    }
    catch (libconfig::FileIOException & e) {
      std::cerr << " Cannot read " << filename << "\n";
      exit(EXIT_FAILURE);
    }catch (libconfig::ParseException &e) {
      std::cerr << "Parse error at " << e.getFile() << ":" << e.getLine()
                << " - " << e.getError() <<"\n";
      exit(EXIT_FAILURE);
    }
    ReadValueRequired(cfg, "CaseType",caseType);
    ReadArrayRequired(cfg, "Energies", energies);
    std::vector<Real> _temp;
    ReadArrayRequired(cfg, "EAngleRotation", _temp,3);

    startAngle = _temp[0]; incrementAngle = _temp[1]; endAngle = _temp[2];
    if(startAngle > endAngle){
      throw std::runtime_error("Start EAngle can not be greater than end EAngle");
    }
    if(incrementAngle == 0){
      if(not(FEQUALS(startAngle,endAngle))){
        throw std::runtime_error("Increment angle = 0");
      }
      else {
        incrementAngle = 1.0;
      }
    }
    ReadValueRequired(cfg, "MorphologyType", morphologyType);

    if(ReadValue(cfg, "NumThreads", num_threads)){}
    if(ReadValue(cfg, "RotMask",rotMask)){}
    if(ReadValue(cfg, "EwaldsInterpolation",ewaldsInterpolation)){}
//    if(ReadValue(cfg, "WriteVTI",writeVTI)){}
//    if(ReadValue(cfg, "VTIDirName",VTIDirName)){}
    if(ReadValue(cfg, "HDF5DirName",HDF5DirName)){}
    if(ReadValue(cfg, "WindowingType",windowingType)){}
    if(ReadValue(cfg, "Algorithm",algorithmType)){}
    if(ReadValue(cfg,"ScatterApproach",scatterApproach)){}
    if(ReadValue(cfg,"DumpMorphology",dumpMorphology)){}
    if(ReadValue(cfg,"MaxStreams",numMaxStreams)){}
    int _temp1;
    if(ReadValue(cfg,"ReferenceFrame",_temp1)){
       referenceFrame = static_cast<bool>(_temp1);
    }

    if(caseType == CaseTypes::DEFAULT) {
      kVectors.resize(1,{0,0,1});
    }
    else {
      const std::string key = "listKVectors";
      if (!cfg.exists(key)) {
        std::cerr << "[Input Error] No value corresponding to " << key
                  << " found. Exiting\n";
        exit(EXIT_FAILURE);
      }
      const libconfig::Setting & listOfKVectors = cfg.getRoot()[key];
      kVectors.resize(listOfKVectors.getLength());
      for(int i = 0; i < kVectors.size(); i++) {
        const libconfig::Setting & kRoot = listOfKVectors[i].lookup("k");
        std::array<Real,3> k;
        ReadArray<Real,3>(kRoot,"k",k);
        kVectors[i].x = k[0];
        kVectors[i].y = k[1];
        kVectors[i].z = k[2];
        normalizeVec(kVectors[i]);
      }
    }
    if(caseType == CaseTypes::GRAZING_INCIDENCE){
      ReadArrayRequired(cfg, "DetectorCoordinates", _temp,3);
      detectorCoordinates.x = _temp[0];
      detectorCoordinates.y = _temp[1];
      detectorCoordinates.z = _temp[2];
      normalizeVec(detectorCoordinates);
    }
  }

  void readRefractiveIndexData(std::vector<Material> &refractiveIndex) const {
    libconfig::Config cfg;
    refractiveIndex.resize(energies.size()*NUM_MATERIAL);
    const UINT & numEnergy = energies.size();
    for (int numMaterial = 0; numMaterial < NUM_MATERIAL; numMaterial++) {
      std::string fname = "Material" + std::to_string(numMaterial+1) + ".txt";
      try {
        cfg.readFile(fname.c_str());
      }
      catch (libconfig::FileIOException & e) {
        std::cerr << " Cannot read " << fname << "\n";
        exit(EXIT_FAILURE);
      }catch (libconfig::ParseException &e) {
        std::cerr << "Parse error at " << e.getFile() << ":" << e.getLine()
                  << " - " << e.getError() <<"\n";
        exit(EXIT_FAILURE);
      }
      for (int i = 0; i < numEnergy; i++) {
        const auto &global = cfg.getRoot()["EnergyData" + std::to_string(i)];
        Real energy = global["Energy"];
        Real currEnergy =energies[i];
        Real diff = fabs(energy - currEnergy);

        if (diff > 1E-3) {
          std::cout << "[Input Error] No energy found for " << currEnergy << "\n";
          exit(EXIT_FAILURE);
        }

        Real deltaPara = global["DeltaPara"];
        Real betaPara = global["BetaPara"];
        Real deltaPerp = global["DeltaPerp"];
        Real betaPerp = global["BetaPerp"];

        refractiveIndex[i*NUM_MATERIAL + numMaterial].npara.x = 1 - deltaPara;
        refractiveIndex[i*NUM_MATERIAL + numMaterial].npara.y = betaPara;

        refractiveIndex[i*NUM_MATERIAL + numMaterial].nperp.x = 1 - deltaPerp;
        refractiveIndex[i*NUM_MATERIAL + numMaterial].nperp.y = betaPerp;
      }

    }
  }
  void validate() const{
      validate("FFT Windowing",windowingType,FFT::FFTWindowing::MAX_SIZE);
      validate("Scatter Approach",scatterApproach,ScatterApproach::MAX_SCATTER_APPROACH);
      validate("Ewalds Interpolation",ewaldsInterpolation,Interpolation::EwaldsInterpolation::MAX_SIZE);
      validate("Morphology Type",ewaldsInterpolation,MorphologyType::MAX_MORPHOLOGY_TYPE);
      validate("Case Type",caseType,CaseTypes::MAX_CASE_TYPE);
      std::cout << GRN << "Input Data : [OK] " << NRM << "\n";
  }
    /**
     * @brief prints the input data
     */
    void print() const{
        std::cout << "NumMaterial          : " << NUM_MATERIAL << "\n";
        if(morphologyOrder == MorphologyOrder::XYZ){
        std::cout << "Dimensions [X Y Z]   : ["<< voxelDims[0] << " " <<  voxelDims[1] << " " << voxelDims[2] << "]\n";
        }
        else{
        std::cout << "Dimensions [Z Y X]   : ["<< voxelDims[2] << " " <<  voxelDims[1] << " " << voxelDims[0] << "]\n";
        }

        std::cout << "PhysSize             : " << physSize << " nm \n";
        std::cout << "E Rotation Angle     : " << startAngle << " : " << incrementAngle << " : " <<endAngle << "\n";
        std::cout << "Morphology Type      : " << morphologyTypeName[morphologyType] << "\n";
        std::cout << "Morphology Order     : " << morphologyOrderName[morphologyOrder] << "\n";
        std::cout << "Energies simulated   : [";
        for (const auto & energy: energies) {
            std::cout << energy << " " ;
        }
        std::cout << "]\n";
        std::cout << "Windowing Type       : " << FFT::windowingName[windowingType] << "\n";
        std::cout << "Rotation Mask        : " << rotMask << "\n";
        std::cout << "Interpolation Type   : " << Interpolation::interpolationName[ewaldsInterpolation] << "\n";
        std::cout << "HDF Output Directory : " << HDF5DirName << "\n";
        std::cout << "Scatter Approach     : " << scatterApproachName[scatterApproach] << "\n";
        std::cout << "Algorithm            : " << algorithmName[algorithmType] << "\n";
	      std::cout << "Reference Frame      : " << referenceFrameName[(UINT)referenceFrame] << "(" << referenceFrame << ")\n";
         if(algorithmType==Algorithm::MemoryMinizing) {
          std::cout  << "MaxStreams           : " << numMaxStreams << "\n";
        }


    }
#else
    /**
    * @brief Constructor
     * @param [in] numberOfMaterial number of materials
    */
    InputData(int numberOfMaterial){
      writeHDF5 = false;
      paramChecker_.reset();
      NUM_MATERIAL = numberOfMaterial;
      paramChecker_[ParamChecker::Parameters::NUM_MATERIAL] = true;
    }
    /**
     * @brief Adds the energy data
     * @param _energies list of energies (in eV)
     */
    void setEnergies(const std::vector<Real> & _energies){
        energies.clear();
        energies = _energies;
        std::sort(energies.begin(),energies.end());
        const bool hasDuplicates = std::adjacent_find(energies.begin(), energies.end()) != energies.end();
        if(hasDuplicates) {
            pybind11::print("energies has duplicate entries . Can not add");
            return;
        }
        paramChecker_.set(ParamChecker::Parameters::ENERGY,true);
    }
    /**
     * @brief  Set the dimensions.
     * @param shape of numpy / related Array
     * @param _morphologyOrder order of morphology XYZ/ZYX
     */
    void setDimension(const std::array<UINT,3> & shape,const MorphologyOrder & _morphologyOrder) {
      morphologyOrder = _morphologyOrder;
      if(morphologyOrder == MorphologyOrder::XYZ) {
        voxelDims[0] = shape[0];
        voxelDims[1] = shape[1];
        voxelDims[2] = shape[2];
      }
      else if(morphologyOrder == MorphologyOrder::ZYX) {
        voxelDims[2] = shape[0];
        voxelDims[1] = shape[1];
        voxelDims[0] = shape[2];
      }
      else {
        pybind11::print("Invalid morphology order.");
        return;
      }
      check2D();
      paramChecker_.set(ParamChecker::Parameters::DIMENSION,true);
    }
    /**
     * @brief set Physical size
     * @param _physSize PhysSize (in nm)
     */
    void setPhysSize(const Real & _physSize) {
      physSize = _physSize;
      paramChecker_.set(ParamChecker::Parameters::PHYSSIZE,true);
    }

    /**
     * @brief Set the angles for rotation for Electric field
     * @param _startAngle start Angle (in degrees)
     * @param _endAngle   end Angle (in degrees)
     * @param _incrementAngle increment in Angle (in degrees)
     */
    void setEAngles(const Real & _startAngle, const Real & _endAngle, const Real & _incrementAngle) {
      startAngle = _startAngle;
      endAngle = _endAngle;
      incrementAngle = _incrementAngle;
      if(startAngle > endAngle) {
        pybind11::print("[ERROR]: Start EAngle cannot be greater than end EAngle");
        return;
      }
      if(FEQUALS(startAngle,endAngle)) {
          incrementAngle = 1.0;
      }
      else {
          if(FEQUALS(incrementAngle,0.0)) {
              pybind11::print("[ERROR] :  Increment angle for Electric field rotation cannot be 0");
              return;
          }
      }
      paramChecker_.set(ParamChecker::Parameters::EANGLE,true);
    }

    void setMorphologyType(const MorphologyType _morphologyType) {
      morphologyType = _morphologyType;
      paramChecker_.set(ParamChecker::Parameters::MORPHOLOGY_TYPE,true);
    }

    void setCaseType(const CaseTypes _caseType) {
      caseType = _caseType;
      paramChecker_.set(ParamChecker::Parameters::CASE_TYPE,true);
      if(caseType == CaseTypes::DEFAULT) {
        Real3 kVec({0,0,1});
        kVectors.resize(1,kVec);
        paramChecker_.set(ParamChecker::Parameters::KVECTORS,true);
        paramChecker_.set(ParamChecker::Parameters::DETECTOR_COORD,true);
      }
    }

    void setKVectors(const std::array<Real,3> & _kVector) {
      if(caseType == CaseTypes::DEFAULT) {
        pybind11::print("Cannot add kVectors for default case type");
        return;
      }
      Real3 kVec({_kVector[0],_kVector[1],_kVector[2]});
      normalizeVec(kVec);
      kVectors.push_back(kVec);
      paramChecker_.set(ParamChecker::Parameters::KVECTORS,true);
      if(caseType == CaseTypes::BEAM_DIVERGENCE) {
        paramChecker_.set(ParamChecker::Parameters::DETECTOR_COORD,true);
      }

    }

    void setAlgorithm(const int & algID, int _numMaxStream = 1) {
      if(algID >= Algorithm::MAXAlgorithmType) {
        pybind11::print("Incorrect AlgID.");
        return;
      }
      algorithmType = algID;
      numMaxStreams = _numMaxStream;

    }

    void setDetectorCoordinates(const std::array<Real,3> & _detectorCoordinates) {
      if(caseType != CaseTypes::GRAZING_INCIDENCE) {
        pybind11::print("Cannot add Detector for default or Beam Divergence type");
        return;
      }
      detectorCoordinates = {_detectorCoordinates[0],_detectorCoordinates[1],_detectorCoordinates[2]};
      paramChecker_.set(ParamChecker::Parameters::DETECTOR_COORD,true);
    }

    /**
     * @brief prints the input data
     */
    void print() const{
        pybind11::print("Required options:");
        pybind11::print("==================================================");
        pybind11::print("NumMaterial          : ",NUM_MATERIAL);
        pybind11::print("CaseType             : ",caseTypenames[caseType]);
        pybind11::print("MorphologyType       : ",morphologyTypeName[morphologyType]);
        if(morphologyOrder == MorphologyOrder::XYZ) {
          pybind11::print("Dimensions [X Y Z]   :  [", voxelDims[0], ",", voxelDims[1], ",", voxelDims[2], "]");
        }
        else {
          pybind11::print("Dimensions [Z Y X]   :  [", voxelDims[2], ",", voxelDims[1], ",", voxelDims[0], "]");
        }
        pybind11::print("PhysSize             : ", physSize , "nm");
        pybind11::print("Energy               : ",energies);
        pybind11::print("ERotation Angle      : ",startAngle , " : ", incrementAngle, " : ",endAngle);
        pybind11::print("Morphology Order     : ",morphologyOrderName[morphologyOrder]);
        pybind11::print("KVectorList           ");
        for(const auto & kVec:kVectors) {
          pybind11::print("                     :  [",kVec.x,",",kVec.y,",",kVec.z,"]");
        }

        pybind11::print("\n");
        pybind11::print("Optional options:");
        pybind11::print("==================================================");
        pybind11::print("Number of openMP threads : ",num_threads);
        pybind11::print("Interpolation Type       : ",Interpolation::interpolationName[ewaldsInterpolation]);
        pybind11::print("Windowing Type           : ",FFT::windowingName[windowingType]);
        pybind11::print("Rotation Mask            : ",rotMask);
        pybind11::print("Reference Frame          : ",referenceFrameName[(UINT)referenceFrame]);
        pybind11::print("Algorithm                : ",algorithmName[algorithmType]);
        if(algorithmType==Algorithm::MemoryMinizing) {
        pybind11::print("NumMaxStreams            : ",numMaxStreams);
        }

    }

    /**
     * @brief validate the input data
     * @return  True if the input data is correct. False otherwise.
     */
    bool validate() const {
      if(voxelDims[2] == 1) {assert(enable2D_);}
      if(voxelDims[0] != 1) {assert(not(enable2D_));}

        if(not(paramChecker_.all())) {
          for(int i = 0; i < paramChecker_.size(); i++) {
              if (not(paramChecker_.test(i))) {
                  pybind11::print("The value for ",ParamChecker::paramNames[i], "is not set");
              }
          }
      }
      return(paramChecker_.all());
    }
#endif
    void printToFile(std::ofstream & fout) const{
        fout << "CaseType             : " << caseTypenames[caseType] << "\n";
        fout << "Morphology Type      : " << morphologyTypeName[morphologyType] << "\n";
        fout << "Morphology Order     : " << morphologyOrderName[morphologyOrder] << "\n";
        fout << "Reference Frame      : " << referenceFrameName[(UINT)referenceFrame] << "(" << referenceFrame << ")\n";
        if(morphologyOrder == MorphologyOrder::XYZ){
          fout << "Dimensions [X Y Z]   : ["<< voxelDims[0] << " " <<  voxelDims[1] << " " << voxelDims[2] << "]\n";
        }
        else{
          fout << "Dimensions [Z Y X]   : ["<< voxelDims[2] << " " <<  voxelDims[1] << " " << voxelDims[0] << "]\n";
        }
        fout << "PhysSize             : " << physSize << "nm \n";
        fout << "E Rotation Angle     : " << startAngle << " : " << incrementAngle << " : " <<endAngle << "\n";
        fout << "Energies simulated   : [";
        for (const auto & energy: energies) {
            fout << energy << " " ;
        }
        fout << "]\n";
        fout << "Windowing Type       : " << FFT::windowingName[windowingType] << "\n";
        fout << "Rotation Mask        : " << rotMask << "\n";
        fout << "Interpolation Type   : " << Interpolation::interpolationName[ewaldsInterpolation] << "\n";
        fout << "kVectors             : (";
        for (const auto & kVec: kVectors) {
          fout << "[" << kVec.x << "," << kVec.y << "," << kVec.z << "]" ;
        }
        fout << ")\n";
#ifndef PYBIND
        fout << "HDF Output Directory : " << HDF5DirName << "\n";
#endif
        fout << "Scatter Approach     : " << scatterApproachName[scatterApproach] << "\n";
        fout << "Algorithm            : " << algorithmName[algorithmType] << "\n";
        if(algorithmType==Algorithm::MemoryMinizing) {
          fout << "MaxStreams           : " << numMaxStreams << "\n";
        }


    }

};

#endif //PRS_READCONFIG_H
