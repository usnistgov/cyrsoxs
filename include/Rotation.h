/////////////////////////////////////////////////////////////////////////////////
// MIT License
//
//Copyright (c) 2019 - 2021 Iowa State University
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

#ifndef CY_RSOXS_ROTATION_H
#define CY_RSOXS_ROTATION_H

#include <Datatypes.h>
#include <limits>
#include <cudaHeaders.h>

/// Class to store 3D matrix
struct Matrix{
private:
  /// matrix object
  Real matrix[9]{};
public:
  /**
   * @brief Default Constructor
   */
  __host__ __device__ Matrix(){
    std::memset(matrix, 0, sizeof(Real) * 9);
  }
  /**
   * @brief Copies the data from array
   * @param [in] A 3X3 matrix
   */
  __host__  Matrix(const Real * A){
    std::memcpy(matrix, A, sizeof(Real) * 9);
  }
  /**
   * @brief returns the value at mat[id1,id2]
   * @tparam id1  first id
   * @tparam id2 second id
   * @return the value at matrix location
   */
  template<UINT id1, UINT id2>
  __host__ __device__ INLINE inline Real  getValue() const{
    static constexpr UINT matID = id1*3 +id2;
    static_assert((matID < 9),"Wrong matID");
    return matrix[matID];
  }

  /**
   * @brief Getter
   * @param [out] mat returns the value of matrix
   */
  __host__ INLINE inline void  getValue(Real * mat) const{
    std::memcpy(mat,matrix, sizeof(Real)*9);
  }

  /**
   * @brief sets the value at mat[id1,id2]
   * @tparam [in] id1  first id
   * @tparam [in] id2 second id
   * @param [in] value value
   */
  template<UINT id1, UINT id2>
  __host__ __device__ INLINE inline void  setValue(const Real & value) {
    static constexpr UINT matID = id1*3 +id2;
    static_assert((matID < 9),"Wrong matID");
    matrix[matID] = value;
  }

  /**
   * @brief sets matrix to identity matrix
   */
  __host__ __device__ INLINE inline void  setIdentity(){
    std::memset(matrix, 0, sizeof(Real) * 9);
    this->setValue<0,0>(1);
    this->setValue<1,1>(1);
    this->setValue<2,2>(1);
  }
  /**
   * @brief resets all entries to 0
   */
  __host__ __device__ INLINE inline void  reset(){
    std::memset(matrix, 0, sizeof(Real) * 9);
  }

  /**
   * @brief performs matrix multiplication
   * @tparam transpose1 whether to compute the transpose of matrix A
   * @tparam transpose2 whether to compute the transpose of matrix B
   * @param A matrix A
   * @param B matrix B
   */
  template<bool transpose1, bool transpose2>
  __host__ __device__ INLINE inline void performMatrixMultiplication(const Matrix &A, const Matrix &B) {
    this->reset();
    if (not(transpose1) and not(transpose2)) {
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          for (int k = 0; k < 3; ++k) {
            this->matrix[i * 3 + j] += A.matrix[i * 3 + k] * B.matrix[k * 3 + j];
          }
        }
      }
    }

    if(not(transpose1) and transpose2){
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          for (int k = 0; k < 3; ++k) {
            this->matrix[i * 3 + j] += A.matrix[i * 3 + k] * B.matrix[j * 3 + k];
          }
        }
      }
    }

    if((transpose1) and not(transpose2)){
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          for (int k = 0; k < 3; ++k) {
            this->matrix[i * 3 + j] += A.matrix[k * 3 + i] * B.matrix[k * 3 + j];
          }
        }
      }
    }
    if((transpose1) and (transpose2)){
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          for (int k = 0; k < 3; ++k) {
            this->matrix[i * 3 + j] += A.matrix[k * 3 + i] * B.matrix[j * 3 + k];
          }
        }
      }
    }
  }
  __host__ void print()const{
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        std::cout << this->matrix[i * 3 + j] << " ";
      }
      std::cout << "\n";
    }
  }
  __host__ void printToFile(std::ofstream & fout)const{
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        fout << this->matrix[i * 3 + j] << " ";
      }
      fout << "\n";
    }
  }
};
/**
 * @brief computes dot product of two vectors
 * @param [in] vecA 3D vector A
 * @param [in] vecB 3D vector B
 * @return resultant dot product
 */
__host__ static inline Real computeDotProduct(const Real3 & vecA, const Real3 & vecB){
  return (vecA.x*vecB.x + vecA.y*vecB.y + vecA.z*vecB.z);
}

/**
 * @brief computes cross product of two vectors
 * @param [in] vecA 3D vector A
 * @param [in] vecB 3D vector B
 * @return resultant cross product
 */
__host__ static inline Real3 computeCrossProduct(const Real3 & vecA, const Real3 & vecB){
  Real3 crossProduct;
  crossProduct.x =   vecA.y*vecB.z - vecA.z*vecB.y;
  crossProduct.y = -(vecA.x*vecB.z - vecA.z*vecB.x);
  crossProduct.z =   vecA.x*vecB.y - vecA.y*vecB.x;
  return crossProduct;
}

/**
 * @brief computes norm of a vector
 * @param vec vector
 * @return resultant norm
 */
__host__ static inline Real computeVecNorm(const Real3 & vec){
  return (sqrt(computeDotProduct(vec,vec)));
}
/**
 * @brief computes || A X B ||
 * @param [in] vecA 3D vector A
 * @param [in] vecB 3D vector B
 * @return norm of cross product
 */
__host__ static inline Real computeNormCrossproduct(const Real3 & vecA, const Real3 & vecB){
  const Real3 & crossProduct = computeCrossProduct(vecA,vecB);
  return (sqrt(computeDotProduct(crossProduct,crossProduct)));
}

/**
 *
 * @param [in] vecIn input 3D vector
 * @param [out] scaledVec scaled vector
 * @param [in] scaleFactor scale Factor
 */
__host__ static inline void scaleVec(const Real3 & vecIn, Real3 & scaledVec, const Real & scaleFactor){
  scaledVec.x = scaleFactor*vecIn.x;
  scaledVec.y = scaleFactor*vecIn.y;
  scaledVec.z = scaleFactor*vecIn.z;
}

/**
 * @brief returns the inverse of the matrix
 * @param [in] matrixA matrix
 * @param [out] inverseMatrix inverse of the matrix A
 */
__host__ static inline void computeInverseMatrix(const Matrix &  matrixA, Matrix & inverseMatrix){
  double det = matrixA.getValue<0,0>() * (matrixA.getValue<1,1>() * matrixA.getValue<2,2>() - matrixA.getValue<2,1>() * matrixA.getValue<1,2>()) -
               matrixA.getValue<0,1>() * (matrixA.getValue<1,0>() * matrixA.getValue<2,2>() - matrixA.getValue<1,2>() * matrixA.getValue<2,0>()) +
               matrixA.getValue<0,2>() * (matrixA.getValue<1,0>() * matrixA.getValue<2,1>() - matrixA.getValue<1,1>() * matrixA.getValue<2,0>());
  bool val = fabs(det) < 1E-12;
  assert(not(val));

  Real invdet = 1. / det;

  inverseMatrix.setValue<0,0>((matrixA.getValue<1,1>() * matrixA.getValue<2,2>() - matrixA.getValue<2,1>() * matrixA.getValue<1,2>()) * invdet) ;
  inverseMatrix.setValue<0,1>((matrixA.getValue<0,2>() * matrixA.getValue<2,1>() - matrixA.getValue<0,1>() * matrixA.getValue<2,2>()) * invdet) ;
  inverseMatrix.setValue<0,2>((matrixA.getValue<0,1>() * matrixA.getValue<1,2>() - matrixA.getValue<0,2>() * matrixA.getValue<1,1>()) * invdet) ;
  inverseMatrix.setValue<1,0>((matrixA.getValue<1,2>() * matrixA.getValue<2,0>() - matrixA.getValue<1,0>() * matrixA.getValue<2,2>()) * invdet) ;
  inverseMatrix.setValue<1,1>((matrixA.getValue<0,0>() * matrixA.getValue<2,2>() - matrixA.getValue<0,2>() * matrixA.getValue<2,0>()) * invdet) ;
  inverseMatrix.setValue<1,2>((matrixA.getValue<1,0>() * matrixA.getValue<0,2>() - matrixA.getValue<0,0>() * matrixA.getValue<1,2>()) * invdet) ;
  inverseMatrix.setValue<2,0>((matrixA.getValue<1,0>() * matrixA.getValue<2,1>() - matrixA.getValue<2,0>() * matrixA.getValue<1,1>()) * invdet) ;
  inverseMatrix.setValue<2,1>((matrixA.getValue<2,0>() * matrixA.getValue<0,1>() - matrixA.getValue<0,0>() * matrixA.getValue<2,1>()) * invdet) ;
  inverseMatrix.setValue<2,2>((matrixA.getValue<0,0>() * matrixA.getValue<1,1>() - matrixA.getValue<1,0>() * matrixA.getValue<0,1>()) * invdet) ;

}

/**
  * @brief performs matrix multiplication
  * @tparam transpose1 whether to compute the transpose of matrix A
  * @tparam transpose2 whether to compute the transpose of matrix B
  * @param [in] matA matrix A
  * @param [in] matB matrix B
  * @param [out] mat matrix product
  */
template<bool transpose1, bool transpose2>
__host__ static void inline performMatrixMultiplication(const Matrix &  matA, const Matrix &  matB, Matrix &  mat) {
  mat.reset();
  mat.performMatrixMultiplication<transpose1,transpose2>(matA,matB);

}



/**
 * @brief computes Matrix that will transform originalVec into transformed Vec
 *        [RotationMatrix X orignalVec = transformedVec]
 * @param [in/out] originalVec : the orginal 3D vector
 * @param [in] transformedVec  : 3D vector that the original Vec is transformed to.
 * @param [out] RotationMatrix : The resultant 3 X 3 rotation matrix
 */
__host__ static void computeRotationMatrix(const Real3 & originalVec,const Real3 & transformedVec ,Matrix & RotationMatrix){
  /**
   * Matlab code: (https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d)
   * G =  [ dot(A,B) -norm(cross(A,B)) 0;...
           norm(cross(A,B)) dot(A,B)  0;...
           0              0           1];

      F = [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];

      UU =  F*G*inv(F);
   */
  RotationMatrix.reset();
  if((FEQUALS(originalVec.x,transformedVec.x)) and (FEQUALS(originalVec.y,transformedVec.y)) and (FEQUALS(originalVec.z,transformedVec.z))){
    RotationMatrix.setIdentity();
    return;
  }
  const Real dotProduct = computeDotProduct(originalVec,transformedVec);
  const Real3 & crossVec = computeCrossProduct(originalVec,transformedVec);
  const Real normCrossProduct = computeVecNorm(crossVec);

  Matrix G;

  G.setValue<0,0>(dotProduct);       G.setValue<0,1>(-normCrossProduct);  G.setValue<0,2>(0);
  G.setValue<1,0>(normCrossProduct); G.setValue<1,1>(dotProduct)       ;  G.setValue<1,2>(0);
  G.setValue<2,0>(0);                G.setValue<2,1>(0)                ;  G.setValue<2,2>(1);



  Real3 scaledVec;
  scaleVec(originalVec,scaledVec,dotProduct);
  const Real3 subVec{transformedVec.x - scaledVec.x,transformedVec.y - scaledVec.y,transformedVec.z - scaledVec.z};
  const Real normSubVec = computeVecNorm(subVec);

  Matrix F, invF,temp1;

  F.setValue<0,0>( originalVec.x);  F.setValue<0,1>( subVec.x/normSubVec); F.setValue<0,2>(-crossVec.x);
  F.setValue<1,0>( originalVec.y);  F.setValue<1,1>( subVec.y/normSubVec); F.setValue<1,2>(-crossVec.y);
  F.setValue<2,0>( originalVec.z);  F.setValue<2,1>( subVec.z/normSubVec); F.setValue<2,2>(-crossVec.z);

  computeInverseMatrix(F,invF);
  performMatrixMultiplication<false,false>(G,invF,temp1);
  performMatrixMultiplication<false,false>(F,temp1,RotationMatrix);

}

/**
 * @brief compute Rodrigues rotation: Rotates a vector inVec about an axis with an angle
 * (https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula)
 * @param [out] rotatedVec rotated Vector
 * @param [in] inVec input vector
 * @param [in] axis
 * @param [in] angle
 *
 */
__host__ static void inline performRodriguesRotation(Real3  & rotatedVec, const Real3 & inVec,const Real3 & axis, const Real & angle){
  const Real cosAlpha = cos(angle);
  const Real sinAlpha = sin(angle);
  const Real3 & crossProduct = computeCrossProduct(axis,inVec);
  const Real & dotProduct = computeDotProduct(axis,inVec);
  Real3 term1{cosAlpha*inVec.x,cosAlpha*inVec.y,cosAlpha*inVec.z};
  Real3 term2{sinAlpha*crossProduct.x,sinAlpha*crossProduct.y,sinAlpha*crossProduct.z};
  Real3 term3{axis.x,axis.y,axis.z};
  scaleVec(term3,term3,(1-cosAlpha)*dotProduct);
  rotatedVec.x = term1.x + term2.x + term3.x;
  rotatedVec.y = term1.y + term2.y + term3.y;
  rotatedVec.z = term1.z + term2.z + term3.z;
}

/**
 *
 * @tparam transpose
 * @param matrix
 * @param vec
 * @param matVec
 */
template <bool transpose>
__host__ __device__ static void doMatVec(const Matrix & matrix, const Real3 & vec, Real3 & matVec){
  if(transpose) {
    matVec.x = matrix.template getValue<0,0>() * vec.x + matrix.template getValue<1,0>() * vec.y + matrix.template getValue<2,0>() * vec.z;
    matVec.y = matrix.template getValue<0,1>() * vec.x + matrix.template getValue<1,1>() * vec.y + matrix.template getValue<2,1>() * vec.z;
    matVec.z = matrix.template getValue<0,2>() * vec.x + matrix.template getValue<1,2>() * vec.y + matrix.template getValue<2,2>() * vec.z;
  }
  else{
    matVec.x = matrix.template getValue<0,0>() * vec.x + matrix.template getValue<0,1>() * vec.y + matrix.template getValue<0,2>() * vec.z;
    matVec.y = matrix.template getValue<1,0>() * vec.x + matrix.template getValue<1,1>() * vec.y + matrix.template getValue<1,2>() * vec.z;
    matVec.z = matrix.template getValue<2,0>() * vec.x + matrix.template getValue<2,1>() * vec.y + matrix.template getValue<2,2>() * vec.z;
  }

}

/**
 * @brief rotates the  complex vector according to rotation matrix
 * @tparam transpose weather to use transpose of the rotation matrix
 * @param [in] rotationMatrix rotation matrix
 * @param [in,out] vecX rotated vector
 * @param [in,out] vecY rotated vector
 * @param [in,out] vecZ rotated vector
 */
template <bool transpose>
__host__ __device__ static void rotate(const Matrix & rotationMatrix,  Complex & vecX, Complex & vecY, Complex & vecZ){
  Complex tempX, tempY, tempZ;
  if(transpose) {
    tempX.x = rotationMatrix.template getValue<0,0>() * vecX.x + rotationMatrix.template getValue<1,0>() * vecY.x + rotationMatrix.template getValue<2,0>() * vecZ.x;
    tempY.x = rotationMatrix.template getValue<0,1>() * vecX.x + rotationMatrix.template getValue<1,1>() * vecY.x + rotationMatrix.template getValue<2,1>() * vecZ.x;
    tempZ.x = rotationMatrix.template getValue<0,2>() * vecX.x + rotationMatrix.template getValue<1,2>() * vecY.x + rotationMatrix.template getValue<2,2>() * vecZ.x;

    tempX.y = rotationMatrix.template getValue<0,0>() * vecX.y + rotationMatrix.template getValue<1,0>() * vecY.y + rotationMatrix.template getValue<2,0>() * vecZ.y;
    tempY.y = rotationMatrix.template getValue<0,1>() * vecX.y + rotationMatrix.template getValue<1,1>() * vecY.y + rotationMatrix.template getValue<2,1>() * vecZ.y;
    tempZ.y = rotationMatrix.template getValue<0,2>() * vecX.y + rotationMatrix.template getValue<1,2>() * vecY.y + rotationMatrix.template getValue<2,2>() * vecZ.y;
  }
  else {
    tempX.x = rotationMatrix.template getValue<0,0>() * vecX.x + rotationMatrix.template getValue<0,1>() * vecY.x + rotationMatrix.template getValue<0,2>() * vecZ.x;
    tempY.x = rotationMatrix.template getValue<1,0>() * vecX.x + rotationMatrix.template getValue<1,1>() * vecY.x + rotationMatrix.template getValue<1,2>() * vecZ.x;
    tempZ.x = rotationMatrix.template getValue<2,0>() * vecX.x + rotationMatrix.template getValue<2,1>() * vecY.x + rotationMatrix.template getValue<2,2>() * vecZ.x;

    tempX.y = rotationMatrix.template getValue<0,0>() * vecX.y + rotationMatrix.template getValue<0,1>() * vecY.y + rotationMatrix.template getValue<0,2>() * vecZ.y;
    tempY.y = rotationMatrix.template getValue<1,0>() * vecX.y + rotationMatrix.template getValue<1,1>() * vecY.y + rotationMatrix.template getValue<1,2>() * vecZ.y;
    tempZ.y = rotationMatrix.template getValue<2,0>() * vecX.y + rotationMatrix.template getValue<2,1>() * vecY.y + rotationMatrix.template getValue<2,2>() * vecZ.y;
  }
  vecX = tempX;
  vecY = tempY;
  vecZ = tempZ;
}

/**
 * @brief Computes rotation matrix by Rodrigues rotation that maps  (0,0,1) to given k
 * @param [in] k k vector
 * @param [out] rotationMatrixK rotation matrix
 * @return true if rotation is successful
 */
__host__ bool static computeRotationMatrixK(const Real3 & k, Matrix & rotationMatrixK){
  static constexpr Real3 origK{0,0,1};
  computeRotationMatrix(origK,k,rotationMatrixK);
#if DEBUG
  {
    Real3 shiftedK;
    doMatVec<false>(rotationMatrixK, origK, shiftedK);
    assert((FEQUALS(shiftedK.x, k.x)) and (FEQUALS(shiftedK.y, k.y)) and (FEQUALS(shiftedK.z, k.z)));
  }
#endif
  return true;
}
/**
 * @brief finds the base configuration that corresponds to 0 degree of E rotataion
 * @param k k vector
 * @param [in] rotationMatrixK rotation matrix
 * @param [out] rotationMatrix rotation matrix  for E
 * @param [out] rotAngle rotation angle
 * @return true if the rotation is successful
 */
__host__ bool static computeRotationMatrixBaseConfiguration(const Real3 & k, const Matrix & rotationMatrixK, Matrix & rotationMatrix, Real & rotAngle){
  static constexpr Real3 X{1,0,0};
  static constexpr UINT numInterval = 100000;
  static constexpr Real dTheta = M_PI/(numInterval*1.0);
  Real3 shiftedX;
  doMatVec<false>(rotationMatrixK,X,shiftedX);
  Real3 rotatedX;
  Real maxDiff = std::numeric_limits<Real>::infinity();
  rotAngle = 0;
  for(UINT i = 0; i < numInterval; i++){
    performRodriguesRotation(rotatedX,shiftedX,k,i*dTheta);
    Real diff = fabs(rotatedX.y);
    if(maxDiff > diff){
      maxDiff = diff;
      rotAngle = i*dTheta;
    }
  }
  assert(fabs(maxDiff) < 1E-4);
  performRodriguesRotation(rotatedX,shiftedX,k,rotAngle);
  rotAngle = rotatedX.x > 0 ? (rotAngle): (rotAngle)+M_PI;
  performRodriguesRotation(rotatedX,shiftedX,k,rotAngle);
#if DEBUG
  {
    static constexpr Real3 Y{0,1,0};
    static constexpr Real3 Z{0,0,1};
    Real3 shiftedY,shiftedZ;
    doMatVec<false>(rotationMatrixK, Y, shiftedY);
    doMatVec<false>(rotationMatrixK, Z, shiftedZ);
    Real3 rotatedY, rotatedZ;
    performRodriguesRotation(rotatedY,shiftedY,k,rotAngle);
    performRodriguesRotation(rotatedZ,shiftedZ,k,rotAngle);
    assert(FEQUALS(computeDotProduct(shiftedX,shiftedY),0));
    assert(FEQUALS(computeDotProduct(shiftedX,shiftedZ),0));
    assert(FEQUALS(computeDotProduct(shiftedY,shiftedZ),0));
  };
#endif
  Matrix rotationMatrixX;
  computeRotationMatrix(shiftedX,rotatedX,rotationMatrixX);

  performMatrixMultiplication<false,false>(rotationMatrixX,rotationMatrixK,rotationMatrix);


  return true;
}
/**
 *
 * @param [in] k k vector
 * @param [in] rotationMatrixK K rotation matrix that maps (0,0,1) to given K
 * @param [out] rotationMatrix rotation matrix  for E
 * @param [out] rotAngle rotation angle
 * @return  true if the rotation is successful
 */
__host__ bool static computeRotationMatrix(const Real3 & k, const Matrix & rotationMatrixK, Matrix & rotationMatrix, const Real & rotAngle){
  static constexpr Real3 X{1,0,0};
  Real3 shiftedX;
  doMatVec<false>(rotationMatrixK,X,shiftedX);
  Real3 rotatedX;
  Matrix rotationMatrixX;
  performRodriguesRotation(rotatedX,shiftedX,k,rotAngle);
  computeRotationMatrix(shiftedX,rotatedX,rotationMatrixX);
  performMatrixMultiplication<false,false>(rotationMatrixX,rotationMatrixK,rotationMatrix);
  return true;
}
/**
 * @brief normalize the vector
 * @param [in,out] vec the vector
 */
__host__ inline static void normalizeVec(Real3 & vec){
  Real normVec =  computeVecNorm(vec);
  scaleVec(vec,vec,1./normVec);
}
/**
 * @brief Computes matrix for warp affine
 * @param [in] srcPoint src point
 * @param [in] dstPoint destination point
 * @param [out]  warpAffineMatrix the rotation matrix
 */
__host__ static void INLINE inline computeWarpAffineMatrix(const double  srcPoint[][2], const double  dstPoint[][2], double  warpAffineMatrix[][3]){
  const double & x0 = srcPoint[0][0]; const double & y0 = srcPoint[0][1];
  const double & x1 = srcPoint[1][0]; const double & y1 = srcPoint[1][1];
  const double & x2 = srcPoint[2][0]; const double & y2 = srcPoint[2][1];

  const double & u0 = dstPoint[0][0]; const double & v0 = dstPoint[0][1];
  const double & u1 = dstPoint[1][0]; const double & v1 = dstPoint[1][1];
  const double & u2 = dstPoint[2][0]; const double & v2 = dstPoint[2][1];

  double denom = x0*y1-x1*y0-x0*y2+x2*y0+x1*y2-x2*y1;

  warpAffineMatrix[0][0] = (u0*y1-u1*y0-u0*y2+u2*y0+u1*y2-u2*y1)/(denom);
  warpAffineMatrix[0][1] = -(u0*x1-u1*x0-u0*x2+u2*x0+u1*x2-u2*x1)/(denom);
  warpAffineMatrix[0][2] = (u0*x1*y2-u0*x2*y1-u1*x0*y2+u1*x2*y0+u2*x0*y1-u2*x1*y0)/(denom);
  warpAffineMatrix[1][0] = (v0*y1-v1*y0-v0*y2+v2*y0+v1*y2-v2*y1)/(denom);
  warpAffineMatrix[1][1] = -(v0*x1-v1*x0-v0*x2+v2*x0+v1*x2-v2*x1)/(denom);
  warpAffineMatrix[1][2] = (v0*x1*y2-v0*x2*y1-v1*x0*y2+v1*x2*y0+v2*x0*y1-v2*x1*y0)/(denom);
}

#endif //CY_RSOXS_ROTATION_H
