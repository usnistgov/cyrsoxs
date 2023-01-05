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

#ifndef PRS_DATATYPES_H
#define PRS_DATATYPES_H

#include <cstdio>
#include <cinttypes>
#include <vector_types.h>
#include <cassert>

#ifdef DOUBLE_PRECISION
typedef double Real;
typedef double2 Real2;
typedef double3 Real3;
typedef double4 Real4;
typedef double2 Complex;
#else
typedef float Real;
typedef float2 Real2;
typedef float3 Real3;
typedef float4 Real4;
typedef float2 Complex;

#endif
#ifdef USE_64_BIT_INDICES
typedef uint64_t BigUINT;
#else
typedef uint32_t BigUINT;
#endif
typedef uint32_t UINT;

#define NUM_THREADS 128



/// Interpolation for I(q)
namespace Interpolation{
    /// Ewalds interpolation type
    enum EwaldsInterpolation: UINT{
         /// Nearest neighbor interpolation
        NEARESTNEIGHBOUR = 0,
        /// Linear interpolation
        LINEAR = 1,
        /// Maximum Type of interpolation
        MAX_SIZE=2
    };
    static const char *interpolationName[]{"Nearest Neighbour", "Trilinear interpolation"};
    static_assert(sizeof(interpolationName) / sizeof(char*) == EwaldsInterpolation::MAX_SIZE,
                  "sizes dont match");
}

/// Reference Frame for computing p vector
enum ReferenceFrame:bool{
  /// Material : Rotate p along with E
  MATERIAL = false,
  /// LAB: Compute in Lab axis frame
  LAB = true
};
static const char *referenceFrameName[]{"MATERIAL", "LAB"};

/// Case Types based on the formation
enum CaseTypes : UINT {
  /// Default : k = (0,0,1) Detector = (0,0,1)
  DEFAULT = 0,
  /// BEAM_DIVERGENCE : k = Arbitrary Detector = (0,0,1)
  BEAM_DIVERGENCE = 1,
  /// GRAZING_INCIDENCE : k = Arbitrary Detector = Arbitrary
  GRAZING_INCIDENCE = 2,
  /// MAX_CASE_TYPE : Maximum case type
  MAX_CASE_TYPE = 3
};

static const char *caseTypenames[]{"Default", "BeamDivergence","GrazingIncidence"};
static_assert(sizeof(caseTypenames) / sizeof(char*) == CaseTypes::MAX_CASE_TYPE,
              "sizes dont match");



/// Fast Fourier Transform
namespace FFT {

    /// Windowing type
    enum FFTWindowing : UINT {
        /// Hanning window
        HANNING = 1,
        /// No window
        NONE = 0,
        /// Maximum size
        MAX_SIZE = 2
    };
    static const char *windowingName[]{"NONE","HANNING"};
    static_assert(sizeof(windowingName)/sizeof(char*) == FFTWindowing::MAX_SIZE,
                  "sizes dont match");
}

/// X(q) computation
enum ScatterApproach:UINT{
    /// Partially compute X(q)
    PARTIAL = 0,
    /// Compute Full X(q)
    FULL = 1,
    /// Maximum type
    MAX_SCATTER_APPROACH = 2
};

/// Algorithm for CyRSoXS
enum Algorithm:UINT{
  /// Minimizes communication at cost of memory
  CommunicationMinimizing = 0,
  /// Minimizes memory at cost of communication
  MemoryMinizing = 1,
  /// Maximum type of algorithm
  MAXAlgorithmType = 2
};
static const char *algorithmName[]{"CommunicationMinimizing","MemoryMinimizing"};
static_assert(sizeof(algorithmName)/sizeof(char*) == Algorithm::MAXAlgorithmType,
              "sizes dont match");

static const char *scatterApproachName[]{"Partial","Full"};
static_assert(sizeof(scatterApproachName)/sizeof(char*) == ScatterApproach::MAX_SCATTER_APPROACH,
              "sizes dont match");
/// Comparison for floating operation
#define FEQUALS(x, y) fabs((x) - (y)) < 1E-6 ? true : false

#define RED "\e[1;31m"
#define BLU "\e[2;34m"
#define GRN "\e[0;32m"
#define YLW "\e[0;33m"
#define MAG "\e[0;35m"
#define CYN "\e[0;36m"
#define NRM "\e[0m"

/// Type of morphology
enum MorphologyType:UINT{
    /// Euler angles
    EULER_ANGLES = 0,
    /// Vector morphology
    VECTOR_MORPHOLOGY = 1,
    /// Maximum type of morphology
    MAX_MORPHOLOGY_TYPE = 2
};
static const char *morphologyTypeName[]{"EulerAngles","VectorMorphology"};
static_assert(sizeof(morphologyTypeName)/sizeof(char*) == MorphologyType::MAX_MORPHOLOGY_TYPE,
              "sizes dont match");


/// Morphology order
enum MorphologyOrder:int{
  /// First axis corresponds to Z and second to Y and third to
  ZYX = 0,
  /// First axis corresponds to X and second to Y and third to Z
  XYZ = 1,
  /// Invalid
  INVALID = -1,
  /// Maximum type of order
  MAX_ORDER = 2
};
static const char *morphologyOrderName[]{"ZYX","XYZ"};
static_assert(sizeof(morphologyOrderName)/sizeof(char*) == MorphologyOrder::MAX_ORDER,
              "sizes dont match");
#endif

