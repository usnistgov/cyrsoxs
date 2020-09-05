/////////////////////////////////////////////////////////////////////////////////
// MIT License
//
//Copyright (c) 2019 - 2020Iowa State University
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
#else
    /// Checks the input data parameter for the required parameters
    std::bitset<ParamChecker::Parameters::MAX_SIZE> paramChecker_;
#endif

 public:
  /// start of energy
  Real energyStart;
  /// end of energy
  Real energyEnd;
  /// increment in Energy
  Real incrementEnergy;
  /// startAngle
  Real startAngle;
  /// end Angle
  Real endAngle;
  /// increment in Angle
  Real incrementAngle;
  /// number of threads
  UINT num_threads=4;
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

#ifndef PYBIND

  /**
   * Constructor to read input data
   * @param refractiveIndex material input
   * @param filename filename, default is config.txt
   */
  InputData(std::vector<Material<NUM_MATERIAL> > &refractiveIndex, std::string filename = "config.txt") {
    libconfig::Config cfg;
    cfg.readFile(filename.c_str());
    ReadValueRequired(cfg, "StartEnergy", energyStart);
    ReadValueRequired(cfg, "EndEnergy", energyEnd);
    ReadValueRequired(cfg, "IncrementEnergy", incrementEnergy);
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
    if(ReadValue(cfg, "WindowingType",windowingType)){}
#if ENABLE_2D

    if(numZ != 1){
      std::cout << "Number of elements in Z - direction is more than 1" << "\n";
      std::cout << "Please re-compile the code by switching DENBLE_2D to No" << "\n";
      exit(EXIT_FAILURE);
    }
#else
    if(numZ == 1){
      std::cout << "Number of elements in Z - direction is  1" << "\n";
      std::cout << "Please re-compile the code by switching DENBLE_2D to Yes" << "\n";
      exit(EXIT_FAILURE);
    }
#endif
    UINT numEnergy = std::round((energyEnd - energyStart) / incrementEnergy + 1);

    refractiveIndex.resize(numEnergy);

    for (int numMaterial = 0; numMaterial < NUM_MATERIAL; numMaterial++) {
      std::string fname = "Material" + std::to_string(numMaterial) + ".txt";
      cfg.readFile(fname.c_str());
      for (int i = 0; i < numEnergy; i++) {
        const auto &global = cfg.getRoot()["EnergyData" + std::to_string(i)];
        Real energy = global["Energy"];
        Real currEnergy = (energyStart + i * incrementEnergy);
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
     * @param _energyStart Start energy (in eV)
     * @param _energyEnd   End energy (in eV)
     * @param _energyIncrement Increment in energy from energy start and energy end (in eV)
     */
    void setEnergy(const Real & _energyStart, const  Real & _energyEnd,const  Real & _energyIncrement){
        if(FEQUALS(_energyIncrement,0.0)) {
          if(FEQUALS(_energyStart,_energyEnd)) {
            pybind11::print("INFO : Adjusting increment Energy to 0.1 for preventing 0/0 error.");
            incrementEnergy = 0.1;
          }
          else {
            pybind11::print("ERROR: Cannot add 0 increment Energy");
            return;
          }
        }
        else {
          incrementEnergy = _energyIncrement;
        }
        energyStart = _energyStart;
        energyEnd = _energyEnd;

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
        pybind11::print("Dimensions           :  [",numX,",",numY,",",numZ,"]");
        pybind11::print("PhysSize             : ", physSize);
        pybind11::print("Energy from          : ",energyStart , "to",energyEnd,"with increment of",incrementEnergy);
        pybind11::print("Rotation Angle  from : ",startAngle , "to",endAngle,"with increment of",incrementAngle);

        pybind11::print("--------Optional options------------------");
        pybind11::print("Number of openMP threads : ",num_threads);
        pybind11::print("Write HDF5               : ",writeHDF5);
        pybind11::print("Interpolation Type       : ",Interpolation::interpolationName[ewaldsInterpolation]);
        pybind11::print("Windowing Type           : ",FFT::windowingName[windowingType]);
    }

    /**
     * @brief validate the input data
     * @return  True if the input data is correct. False otherwise.
     */
    bool validate() const {
#ifdef ENABLE_2D
        if(numZ > 1) {
            throw std::logic_error("Wrong Compilation. Compile with 2D = No");
        }
#else
        if(numZ == 1) {
            throw std::logic_error("Wrong Compilation. Compile with 2D = Yes");
        }

#endif

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
