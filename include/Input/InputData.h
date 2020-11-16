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
#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <bitset>
#include <algorithm>

namespace ParamChecker{
    enum Parameters: u_short {
        ENERGY = 0,
        DIMENSION = 1,
        PHYSSIZE = 2,
        ANGLE = 3,
        MAX_SIZE = 4
    };
    static const char* paramNames[] = {"ENERGY","DIMENSION","PHYSSIZE","ANGLE"};
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
      std::cout << "[WARNING] : No value corresponding to " << key << " found. Setting to default\n";
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
                           std::vector<T>& arr) {
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

        arr.clear();
        for (int i = 0; i < setting.getLength(); i++) {
            T value;
            std::string item_path = key + ".[" + std::to_string(i) + "]";
            ReadValueRequired(config, item_path, value);
            arr.push_back(value);
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

  /**
   * Constructor to read input data
   * @param refractiveIndex material input
   * @param filename filename, default is config.txt
   */
  InputData(std::vector<Material<NUM_MATERIAL> > &refractiveIndex, std::string filename = "config.txt") {
    libconfig::Config cfg;
    cfg.readFile(filename.c_str());
    ReadArrayRequired(cfg, "Energies", energies);
    ReadValueRequired(cfg, "StartAngle", startAngle);
    ReadValueRequired(cfg, "EndAngle", endAngle);
    ReadValueRequired(cfg, "IncrementAngle", incrementAngle);
    ReadValueRequired(cfg, "NumThreads", num_threads);
    ReadValueRequired(cfg, "NumX", numX);
    ReadValueRequired(cfg, "NumY", numY);
    ReadValueRequired(cfg, "NumZ", numZ);
    ReadValueRequired(cfg, "PhysSize", physSize);
    if(ReadValue(cfg, "RotMask",rotMask)){}
    if(ReadValue(cfg, "EwaldsInterpolation",ewaldsInterpolation)){}
    if(ReadValue(cfg, "WriteVTI",writeVTI)){}
    if(ReadValue(cfg, "VTIDirName",VTIDirName)){}
    if(ReadValue(cfg, "HDF5DirName",HDF5DirName)){}
    if(ReadValue(cfg, "WindowingType",windowingType)){}
    if(ReadValue(cfg,"kRotationType",kRotationType)){}
    if(kRotationType == KRotationType::ROTATION){
        ReadValueRequired(cfg,"kStart",kStart);
        ReadValueRequired(cfg,"kEnd",kEnd);
        ReadValueRequired(cfg,"kIncrement",kIncrement);
        if(ReadValue(cfg, "kIncrement",kIncrement)){
            if(FEQUALS(kStart,kEnd)){
                kIncrement = 1.0;
            }
            else{
                if(FEQUALS(kIncrement,0.0)){
                    throw std::logic_error("kIncrement can not be 0\n");
                }
            }
        }
        std::cout << "[WARNING:] This is an experimental routine that you are trying to use\n";
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
      cfg.readFile(fname.c_str());
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

#else
    /**
    * @brief Constructor
    */
    InputData() {
      writeHDF5 = false;
      paramChecker_.reset();
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
    void setAngles(const Real & _startAngle, const Real & _endAngle, const Real & _incrementAngle) {
      startAngle = _startAngle;
      endAngle = _endAngle;
      incrementAngle = _incrementAngle;
      paramChecker_.set(ParamChecker::Parameters::ANGLE,true);
    }

    /**
     * @brief prints the input data
     */
    void print() const{
        pybind11::print("--------Required options------------------");
        pybind11::print("Dimensions           : [",numX,",",numY,",",numZ,"]");
        pybind11::print("PhysSize             : ", physSize);
        pybind11::print("Energy from          : ",energies);
        pybind11::print("Rotation Angle  from : ",startAngle , "to",endAngle,"with increment of",incrementAngle);

        pybind11::print("--------Optional options------------------");
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


};

#endif //PRS_READCONFIG_H
