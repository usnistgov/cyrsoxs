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
#ifndef CUDAUTILS_H
#define CUDAUTILS_H

#include "cudaHeaders.h"
#include "Datatypes.h"
/**
 *
 * @brief Computes the A:A, where A is a complex matrix
 * @param [in] Matrix The input matrix
 * @param [out] MatrixTimesMatrix The resultant matrix
 */
__device__ void compute3X3ComplexMultiplication(const Complex Matrix[][3], Complex MatrixTimesMatrix[][3]) {

/** Completely unrolled matrix X matrix with complex enteries of the form:
     [ a00, a01, a02]         [ b00, b01, b02]
     [ a10, a11, a12]   + i   [ b10, b11, b12]
     [ a20, a21, a22]         [ b20, b21, b22]
  **/
  /**************************** Computation of 1st row ******************************************************/
  // (0,0) = a00^2 - b00^2 + a01*a10 + a02*a20 - b01*b10 - b02*b20 + i (2*a00*b00 + a01*b10 + a10*b01 + a02*b20 + a20*b02)

  MatrixTimesMatrix[0][0].x
      = Matrix[0][0].x * Matrix[0][0].x
      - Matrix[0][0].y * Matrix[0][0].y
      + Matrix[0][1].x * Matrix[1][0].x
      + Matrix[0][2].x * Matrix[2][0].x
      - Matrix[0][1].y * Matrix[1][0].y
      - Matrix[0][2].y * Matrix[2][0].y;

  MatrixTimesMatrix[0][0].y
      = 2 * Matrix[0][0].x * Matrix[0][0].y
      + Matrix[0][1].x * Matrix[1][0].y
      + Matrix[1][0].x * Matrix[0][1].y
      + Matrix[0][2].x * Matrix[2][0].y
      + Matrix[2][0].x * Matrix[0][2].y;

  // (0,1) = a00*a01 + a01*a11 + a02*a21 - b00*b01 - b01*b11 - b02*b21 + i (a00*b01 + a01*b00 + a01*b11 + a11*b01 + a02*b21 + a21*b02 )
  MatrixTimesMatrix[0][1].x
      = Matrix[0][0].x * Matrix[0][1].x
      + Matrix[0][1].x * Matrix[1][1].x
      + Matrix[0][2].x * Matrix[2][1].x
      - Matrix[0][0].y * Matrix[0][1].y
      - Matrix[0][1].y * Matrix[1][1].y
      - Matrix[0][2].y * Matrix[2][1].y;

  MatrixTimesMatrix[0][1].y
      = 2 * Matrix[0][0].x * Matrix[0][1].y
      + Matrix[0][1].x * Matrix[0][0].y
      + Matrix[0][1].x * Matrix[1][1].y
      + Matrix[1][1].x * Matrix[0][1].y
      + Matrix[0][2].x * Matrix[2][1].y
      + Matrix[2][1].x * Matrix[0][2].y;

  // (0,2) = a00*a02 + a01*a12 + a02*a22 - b00*b02 - b01*b12 - b02*b22 + i(a00*b02 + a02*b00 + a01*b12 + a12*b01 + a02*b22 + a22*b02)
  MatrixTimesMatrix[0][2].x
      = Matrix[0][0].x * Matrix[0][2].x
      + Matrix[0][1].x * Matrix[1][2].x
      + Matrix[0][2].x * Matrix[2][2].x
      - Matrix[0][0].y * Matrix[0][2].y
      - Matrix[0][1].y * Matrix[1][2].y
      - Matrix[0][2].y * Matrix[2][2].y;

  MatrixTimesMatrix[0][2].y
      = 2 * Matrix[0][0].x * Matrix[0][2].y
      + Matrix[0][2].x * Matrix[0][0].y
      + Matrix[0][1].x * Matrix[1][2].y
      + Matrix[1][2].x * Matrix[0][1].y
      + Matrix[0][2].x * Matrix[2][2].y
      + Matrix[2][2].x * Matrix[0][2].y;

  /**************************** Computation of 2nd row ******************************************************/
  // (1,0) = a00*a10 + a10*a11 + a12*a20 - b00*b10 - b10*b11 - b12*b20 + i (a00*b10 + a10*b00 + a10*b11 + a11*b10 + a12*b20 + a20*b12)
  MatrixTimesMatrix[1][0].x
      = Matrix[1][0].x * Matrix[0][0].x
      + Matrix[1][1].x * Matrix[1][0].x
      + Matrix[1][2].x * Matrix[2][0].x
      - Matrix[1][0].y * Matrix[0][0].y
      - Matrix[1][1].y * Matrix[1][0].y
      - Matrix[1][2].y * Matrix[2][0].y;

  MatrixTimesMatrix[1][0].y
      = Matrix[0][0].x * Matrix[1][0].y
      + Matrix[1][0].x * Matrix[0][0].y
      + Matrix[1][0].x * Matrix[1][1].y
      + Matrix[1][1].x * Matrix[1][0].y
      + Matrix[1][2].x * Matrix[2][0].y
      + Matrix[2][0].x * Matrix[1][2].y;

  // (1,1) = a11^2 - b11^2 + a01*a10 + a12*a21 - b01*b10 - b12*b21 + i (a01*b10 + a10*b01 + 2*a11*b11 + a12*b21 + a21*b12)
  MatrixTimesMatrix[1][1].x
      = Matrix[1][1].x * Matrix[1][1].x
      + Matrix[0][1].x * Matrix[1][0].x
      + Matrix[1][2].x * Matrix[2][1].x
      - Matrix[1][1].y * Matrix[1][1].y
      - Matrix[0][1].y * Matrix[1][0].y
      - Matrix[1][2].y * Matrix[2][1].y;

  MatrixTimesMatrix[1][1].y
      = Matrix[0][1].x * Matrix[1][0].y
      + Matrix[1][0].x * Matrix[0][1].y
      + 2 * Matrix[1][1].x * Matrix[1][1].y
      + Matrix[1][2].x * Matrix[2][1].y
      + Matrix[2][1].x * Matrix[1][2].y;

  // (1,2) = a02*a10 + a11*a12 + a12*a22 - b02*b10 - b11*b12 - b12*b22 + i (a02*b10 + a10*b02 + a11*b12 + a12*b11 + a12*b22 + a22*b12)
  MatrixTimesMatrix[1][2].x
      = Matrix[0][2].x * Matrix[1][0].x
      + Matrix[1][1].x * Matrix[1][2].x
      + Matrix[1][2].x * Matrix[2][2].x
      - Matrix[0][2].y * Matrix[1][0].y
      - Matrix[1][1].y * Matrix[1][2].y
      - Matrix[1][2].y * Matrix[2][2].y;

  MatrixTimesMatrix[1][2].y
      = Matrix[0][2].x * Matrix[1][0].y
      + Matrix[1][0].x * Matrix[0][2].y
      + Matrix[1][1].x * Matrix[1][2].y
      + Matrix[1][2].x * Matrix[1][1].y
      + Matrix[1][2].x * Matrix[2][2].y
      + Matrix[2][2].x * Matrix[1][2].y;


  /**************************** Computation of 3rd row ******************************************************/
  //(2,0):  a00*a20 + a10*a21 + a20*a22 - b00*b20 - b10*b21 - b20*b22 + i (a00*b20 + a20*b00 + a10*b21 + a21*b10 + a20*b22 + a22*b20)
  MatrixTimesMatrix[2][0].x
      = Matrix[0][0].x * Matrix[2][0].x
      + Matrix[1][0].x * Matrix[2][1].x
      + Matrix[2][0].x * Matrix[2][2].x
      - Matrix[0][0].y * Matrix[2][0].y
      - Matrix[1][0].y * Matrix[2][1].y
      - Matrix[2][0].y * Matrix[2][2].y;

  MatrixTimesMatrix[2][0].y
      = Matrix[0][0].x * Matrix[2][0].y
      + Matrix[2][0].x * Matrix[0][0].y
      + Matrix[1][0].x * Matrix[2][1].y
      + Matrix[2][1].x * Matrix[1][0].y
      + Matrix[2][0].x * Matrix[2][2].y
      + Matrix[2][2].x * Matrix[2][0].y;

  //(2,1) a01*a20 + a11*a21 + a21*a22 - b01*b20 - b11*b21 - b21*b22 + i (a01*b20 + a20*b01 + a11*b21 + a21*b11 + a21*b22 + a22*b21)
  MatrixTimesMatrix[2][1].x
      = Matrix[0][1].x * Matrix[2][0].x
      + Matrix[1][1].x * Matrix[2][1].x
      + Matrix[2][1].x * Matrix[2][2].x
      - Matrix[0][1].y * Matrix[2][0].y
      - Matrix[1][1].y * Matrix[2][1].y
      - Matrix[2][1].y * Matrix[2][2].y;

  MatrixTimesMatrix[2][1].y
      = Matrix[0][1].x * Matrix[2][0].y
      + Matrix[2][0].x * Matrix[0][1].y
      + Matrix[1][1].x * Matrix[2][1].y
      + Matrix[2][1].x * Matrix[1][1].y
      + Matrix[2][1].x * Matrix[2][2].y
      + Matrix[2][2].x * Matrix[2][1].y;

  // (2,2) = a22^2 - b22^2 + a02*a20 + a12*a21 - b02*b20 - b12*b21 + i (a02*b20 + a20*b02 + a12*b21 + a21*b12 + 2*a22*b22)
  MatrixTimesMatrix[2][2].x
      = Matrix[2][2].x * Matrix[2][2].x
      + Matrix[0][2].x * Matrix[2][0].x
      + Matrix[1][2].x * Matrix[2][1].x
      - Matrix[2][2].y * Matrix[2][2].y
      - Matrix[0][2].y * Matrix[2][0].y
      - Matrix[1][2].y * Matrix[2][1].y;

  MatrixTimesMatrix[2][2].y
      = Matrix[0][2].x * Matrix[2][0].y
      + Matrix[2][0].x * Matrix[0][2].y
      + Matrix[1][2].x * Matrix[2][1].y
      + Matrix[2][1].x * Matrix[1][2].y
      + 2 * Matrix[2][2].x * Matrix[2][2].y;

}

/**
 * This function computes the inplace square of the complex variable
 * @param [in,out] val Complex variable
 */
__host__ __device__ inline void computeComplexSquare(Complex & val){
  Real realPart = val.x;
  Real imgPart = val.y;
  val.x = realPart*realPart - imgPart*imgPart;
  val.y = 2*realPart*imgPart;

}

/**
 * @brief Computes the Matrix X Vector
 * @param Matrix [in] Complex matrix
 * @param Vec [in] Vector
 * @param matVec  [out] Result
 */
__device__ void computeMatrixTimesVector(const Complex Matrix[][3],const Real * Vec, Complex * matVec) {
  matVec[0].x = Matrix[0][0].x*Vec[0] + Matrix[0][1].x*Vec[1] + Matrix[0][2].x*Vec[2];
  matVec[0].y = Matrix[0][0].y*Vec[0] + Matrix[0][1].y*Vec[1] + Matrix[0][2].y*Vec[2];

  matVec[1].x = Matrix[1][0].x*Vec[0] + Matrix[1][1].x*Vec[1] + Matrix[1][2].x*Vec[2];
  matVec[1].y = Matrix[1][0].y*Vec[0] + Matrix[1][1].y*Vec[1] + Matrix[1][2].y*Vec[2];

  matVec[2].x = Matrix[2][0].x*Vec[0] + Matrix[2][1].x*Vec[1] + Matrix[2][2].x*Vec[2];
  matVec[2].y = Matrix[2][0].y*Vec[0] + Matrix[2][1].y*Vec[1] + Matrix[2][2].y*Vec[2];
}

/**
 * @brief This function computes the matrix vector product.
 * @param Matrix The matrix of complex type.
 * @param Vec  The vector.
 * @param matVec The matrix vector product
 */
__device__ void computeMatrixTimesVector(const Complex Matrix[][3],const Real3 Vec, Complex * matVec) {
  matVec[0].x = Matrix[0][0].x*Vec.x + Matrix[0][1].x*Vec.y + Matrix[0][2].x*Vec.z;
  matVec[0].y = Matrix[0][0].y*Vec.x + Matrix[0][1].y*Vec.y + Matrix[0][2].y*Vec.z;

  matVec[1].x = Matrix[1][0].x*Vec.x + Matrix[1][1].x*Vec.y + Matrix[1][2].x*Vec.z;
  matVec[1].y = Matrix[1][0].y*Vec.x + Matrix[1][1].y*Vec.y + Matrix[1][2].y*Vec.z;

  matVec[2].x = Matrix[2][0].x*Vec.x + Matrix[2][1].x*Vec.y + Matrix[2][2].x*Vec.z;
  matVec[2].y = Matrix[2][0].y*Vec.x + Matrix[2][1].y*Vec.y + Matrix[2][2].y*Vec.z;
}

__global__ void  getMaxandMinimum(const Real *val, const int2 idx, const uint2 voxel, Real *values) {
  const int blockID = blockIdx.x + blockIdx.y * gridDim.x;
  const int threadID = blockID * (blockDim.x * blockDim.y) + threadIdx.y * blockDim.x + threadIdx.x;

  // Minimum
  if (threadID == idx.x) {
    values[0] = val[threadID];
  }

  // Minimum
  if (threadID == idx.y) {
    values[1] = val[threadID];
  }

}

/**
 * @brief This function computes the matrix vector product.
 * @param Matrix The matrix of complex type.
 * @param Vec  The vector.
 * @param matVec The matrix vector product
 */
__device__ Real computeMagVec1TimesVec1TTimesVec2(const Real *vec1 ,const Complex *vec2, const Real & k) {
    const Real d = k*k;

    const Real & a = vec1[0];
    const Real & b = vec1[1];
    const Real & c = vec1[2];

    const Complex & p1 = vec2[0];
    const Complex & p2 = vec2[1];
    const Complex & p3 = vec2[2];

    Real res = 0;
    Complex _temp;
    _temp.x = -a*a*p1.x + d*p1.x - a*(b*p2.x + c*p3.x);
    _temp.y = -a*a*p1.y + d*p1.y - a*(b*p2.y + c*p3.y);
    res = _temp.x*_temp.x + _temp.y*_temp.y;

    _temp.x = -a*b*p1.x - b*b*p2.x + d*p2.x - b*c*p3.x;
    _temp.y = -a*b*p1.y - b*b*p2.y + d*p2.y - b*c*p3.y;
    res += _temp.x*_temp.x + _temp.y*_temp.y;

    _temp.x = -a*c*p1.x - b*c*p2.x + (-c*c + d)*p3.x;
    _temp.y = -a*c*p1.y - b*c*p2.y + (-c*c + d)*p3.y;
    res += _temp.x*_temp.x + _temp.y*_temp.y;

    return res;


}


#endif //WRITER_CUDAUTILS_H
