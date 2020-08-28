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

#ifndef PRS_DATATYPES_H
#define PRS_DATATYPES_H

#include <cstdio>
#include <cinttypes>
#include <vector_types.h>

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
typedef uint64_t BigUINT;
typedef uint32_t UINT;

#define NUM_THREADS 128

enum EwaldsInterpolation: UINT{
  NEARESTNEIGHBOUR = 0,
  LINEAR = 1
};

enum FFTWindowing:UINT{
  HANNING =1,
  NONE = 0
};
#define FEQUALS(x, y) fabs((x) - (y)) < 1E-10 ? true : false

#endif
