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
#endif
/// This function reads the input data from the file.
class InputData {
#ifndef PYBIND

 private:
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
  UINT num_threads=1;
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
  UINT ewaldsInterpolation = EwaldsInterpolation::LINEAR;
  /// Windowing Type
  UINT windowingType = FFTWindowing::NONE;

  void setEnergy(double energy){
      energyStart = energy;
  }
  void print() const{

  }

#ifndef PYBIND

  /**
   * Constructor to read input data
   * @param materialInput material input
   * @param filename filename, default is config.txt
   */
  InputData(std::vector<Material<NUM_MATERIAL> > &materialInput, std::string filename = "config.txt") {
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

    materialInput.resize(numEnergy);

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
          materialInput[i].npara[numMaterial].x = 1 - deltaPara;
          materialInput[i].npara[numMaterial].y = betaPara;

          materialInput[i].nperp[numMaterial].x = 1 - deltaPerp;
          materialInput[i].nperp[numMaterial].y = betaPerp;
      }

    }
  }
#else

#endif


};

#endif //PRS_READCONFIG_H
