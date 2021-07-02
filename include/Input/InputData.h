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

#ifndef PRS_READCONFIG_H
#define PRS_READCONFIG_H
#include <string>
#ifndef PYBIND
#include <libconfig.h++>
#endif
#include <iostream>
#include <Input/Input.h>
#include <vector>
#include <cmath>
#include <fstream>

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
        KANGLE = 4,
        MAX_SIZE = 5
    };
    static const char* paramNames[] = {"ENERGY","DIMENSION","PHYSSIZE","EANGLE","KANGLE"};
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
   * @param key the name of the configuation value to retrieve
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
  inline void check2D() {
    if (numZ != 1) {
      enable2D_ = false;
    } else {
      std::cout << "[INFO] 2D computation enabled\n";
      enable2D_ = true;
    }
  }

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
  UINT numX;
  /// Number of voxels in Y direction.
  UINT numY;
  /// Number of voxels in Z direction.
  UINT numZ;
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
  /// Start of k rotation
  Real kStart = 0;
  /// End of k rotation
  Real kEnd = 0;
  /// k Intcrement
  Real kIncrement = 1.0;
  /// kRotationType
  UINT kRotationType = KRotationType::NOROTATION;
  /// scatter Approach
  UINT scatterApproach = ScatterApproach::PARTIAL;

  std::string VTIDirName = "VTI";
  std::string HDF5DirName = "HDF5";

  /**
   *
   * @return gets the 2D computation flags
   */
  inline bool if2DComputation() const {
     return enable2D_;
  }


#ifndef PYBIND

  /// Morphology type
  UINT morphologyType;
  /**
   * Constructor to read input data
   * @param refractiveIndex material input
   * @param filename filename, default is config.txt
   */
  InputData(std::vector<Material<NUM_MATERIAL> > &refractiveIndex, std::string filename = "config.txt") {
    libconfig::Config cfg;
    cfg.readFile(filename.c_str());
    ReadArrayRequired(cfg, "Energies", energies);
    std::vector<Real> _temp;
    ReadArrayRequired(cfg, "EAngleRotation", _temp,3);
    startAngle = _temp[0]; incrementAngle = _temp[1]; endAngle = _temp[2];

    ReadValueRequired(cfg, "NumThreads", num_threads);
    ReadValueRequired(cfg, "NumX", numX);
    ReadValueRequired(cfg, "NumY", numY);
    ReadValueRequired(cfg, "NumZ", numZ);
    ReadValueRequired(cfg, "PhysSize", physSize);
    ReadValueRequired(cfg, "MorphologyType", morphologyType);
    if(ReadValue(cfg, "RotMask",rotMask)){}
    if(ReadValue(cfg, "EwaldsInterpolation",ewaldsInterpolation)){}
    if(ReadValue(cfg, "WriteVTI",writeVTI)){}
    if(ReadValue(cfg, "VTIDirName",VTIDirName)){}
    if(ReadValue(cfg, "HDF5DirName",HDF5DirName)){}
    if(ReadValue(cfg, "WindowingType",windowingType)){}
    if(ReadValue(cfg,"KRotationType",kRotationType)){}
    if(ReadValue(cfg,"ScatterApproach",scatterApproach)){}
    else{
        if(numZ < 4){
            scatterApproach = ScatterApproach::FULL;
        }
    }
    if(kRotationType == KRotationType::ROTATION){
        ReadArrayRequired(cfg, "KAngleRotation", _temp,3);
        kStart = _temp[0];kIncrement = _temp[1];kEnd = _temp[2];
        if(FEQUALS(kStart,kEnd)){
                kIncrement = 1.0;
        }
        else{
            if(FEQUALS(kIncrement,0.0)){
                throw std::logic_error("kIncrement can not be 0\n");
            }
        }
        std::cout << MAG << "[WARNING] : This is an experimental routine that you are trying to use " << NRM << "\n";
    } else{
        kStart = 0.0;
        kEnd = 0.0;
        kIncrement = 1.0;
    }

    check2D();

    refractiveIndex.resize(energies.size());
    const UINT & numEnergy = energies.size();
    for (int numMaterial = 0; numMaterial < NUM_MATERIAL; numMaterial++) {
      std::string fname = "Material" + std::to_string(numMaterial) + ".txt";
      try {
        cfg.readFile(fname.c_str());
      }
      catch (libconfig::FileIOException & e) {
        std::cout << "Refractive index config not found for Material " << numMaterial << "\n" << e.what() << "\n";
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

          /** Diagonal enteries **/
          Real deltaPara = global["DeltaPara"];
          Real betaPara = global["BetaPara"];
          Real deltaPerp = global["DeltaPerp"];
          Real betaPerp = global["BetaPerp"];
          /** Diagonal enteries **/
          refractiveIndex[i].npara[numMaterial].x = 1 - deltaPara;
          refractiveIndex[i].npara[numMaterial].y = betaPara;

          refractiveIndex[i].nperp[numMaterial].x = 1 - deltaPerp;
          refractiveIndex[i].nperp[numMaterial].y = betaPerp;
      }

    }
  }

  void validate() const{
      validate("FFT Windowing",windowingType,FFT::FFTWindowing::MAX_SIZE);
      validate("K Rotation",kRotationType,KRotationType::MAX_ROTATION_TYPE);
      validate("Scatter Approach",scatterApproach,ScatterApproach::MAX_SCATTER_APPROACH);
      validate("Ewalds Interpolation",ewaldsInterpolation,Interpolation::EwaldsInterpolation::MAX_SIZE);
      validate("Morphology Type",morphologyType,MorphologyType::MAX_MORPHOLOGY_TYPE);
      if(morphologyType != MorphologyType::VECTOR_MORPHOLOGY) {
        throw std::logic_error("Only vector morphology type is supported");
      }
      std::cout << GRN << "Input Data : [OK] " << NRM << "\n";
  }
    /**
     * @brief prints the input data
     */
    void print() const{

        std::cout << "Dimensions           : ["<< numX << " " <<  numY << " " << numZ << "]\n";
        std::cout << "PhysSize             : " << physSize << " nm \n";
        std::cout << "E Rotation Angle     : " << startAngle << " : " << incrementAngle << " : " <<endAngle << "\n";
        std::cout << "K RotationType       : " << kRotationTypeName[kRotationType] << "\n";
        if(kRotationType == KRotationType::ROTATION) {
            std::cout << "kRotationAngle       : " << kStart << " : " << kIncrement << " : " << kEnd << "\n";
        }
        std::cout << "MorphologyType       : " << morphologyTypeName[morphologyType] << "\n";
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

    }
#else
    /**
    * @brief Constructor
    */
    InputData(){
      writeHDF5 = false;
      paramChecker_.reset();
      paramChecker_.set(ParamChecker::Parameters::KANGLE,true);
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
     * @brief Set the dimensions. Note that HDF5 file dimensions are written in (Z,Y,X)
     * @param _numX number of voxels in X dimensions
     * @param _numY number of voxels in Y dimensions
     * @param _numZ number of voxels in Z dimensions
     */
    void setDimension(const UINT & _numX, const  UINT & _numY,const  UINT & _numZ) {
      numX = _numX;
      numY = _numY;
      numZ = _numZ;
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
    /**
     * @brief Set the angles for rotation for Electric field
     * @param _startAngle start Angle (in degrees)
     * @param _endAngle   end Angle (in degrees)
     * @param _incrementAngle increment in Angle (in degrees)
     */
    void setKAngles(const Real & _startAngle, const Real & _endAngle, const Real & _incrementAngle) {
        if(kRotationType == KRotationType::NOROTATION) {
            pybind11::print("[ERROR] : Trying to set angles with K rotation type set to NONE. First change the KRotationType. Returning." );
            return;
        }
        kStart = _startAngle;
        kEnd = _endAngle;
        kIncrement = _incrementAngle;
        if(FEQUALS(kStart,kEnd)) {
            kIncrement = 1.0;
        }
        else {
            if(FEQUALS(kIncrement,0.0)) {
                pybind11::print("[ERROR] :  Increment angle for K rotation cannot be 0");
                return;
            }
        }
        paramChecker_.set(ParamChecker::Parameters::KANGLE,true);
    }

    void setKRotationType(const KRotationType & _rotationType) {
        kRotationType = _rotationType;
        if(kRotationType == KRotationType::ROTATION) {
            paramChecker_.set(ParamChecker::Parameters::KANGLE,false);
        }

    }
    /**
     * @brief prints the input data
     */
    void print() const{
        pybind11::print("Required options:");
        pybind11::print("==================================================");
        pybind11::print("Dimensions           :  [",numX,",",numY,",",numZ,"]");
        pybind11::print("PhysSize             : ", physSize , "nm");
        pybind11::print("Energy               : ",energies);
        pybind11::print("Rotation Angle       : ",startAngle , " : ", incrementAngle, " : ",endAngle);
        pybind11::print("K Rotation Type      : ",kRotationTypeName[kRotationType]);
        if(kRotationType == KRotationType::ROTATION) {
        pybind11::print("Rotation Angle       : ",kStart , " : ", kIncrement, " : ",kEnd);
        }
        pybind11::print("\n");
        pybind11::print("Optional options:");
        pybind11::print("==================================================");
        pybind11::print("Number of openMP threads : ",num_threads);
        pybind11::print("Interpolation Type       : ",Interpolation::interpolationName[ewaldsInterpolation]);
        pybind11::print("Windowing Type           : ",FFT::windowingName[windowingType]);
        pybind11::print("Rotation Mask            : ",rotMask);

    }

    /**
     * @brief validate the input data
     * @return  True if the input data is correct. False otherwise.
     */
    bool validate() const {
      if(numZ == 1) {assert(enable2D_);}
      if(numZ != 1) {assert(not(enable2D_));}

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
        fout << "Dimensions           : ["<< numX << " " <<  numY << " " << numZ << "]\n";
        fout << "PhysSize             : " << physSize << "nm \n";
        fout << "E Rotation Angle     : " << startAngle << " : " << incrementAngle << " : " <<endAngle << "\n";
        fout << "K RotationType       : " << kRotationTypeName[kRotationType] << "\n";
        if(kRotationType == KRotationType::ROTATION) {
            std::cout << "kRotationAngle       : " << kStart << " : " << kIncrement << " : " << kEnd << "\n";
        }
        fout << "Energies simulated   : [";
        for (const auto & energy: energies) {
            fout << energy << " " ;
        }
        fout << "]\n";
        fout << "Windowing Type       : " << FFT::windowingName[windowingType] << "\n";
        fout << "Rotation Mask        : " << rotMask << "\n";
        fout << "Interpolation Type   : " << Interpolation::interpolationName[ewaldsInterpolation] << "\n";
#ifndef PYBIND
        fout << "HDF Output Directory : " << HDF5DirName << "\n";
#endif
        fout << "Scatter Approach     : " << scatterApproachName[scatterApproach] << "\n";


    }

};

#endif //PRS_READCONFIG_H
