//
// Created by maksbh on 3/30/21.
//

#ifndef CY_RSOXS_ROTATION_H
#define CY_RSOXS_ROTATION_H

#include <Datatypes.h>
#include <limits>

struct Matrix{
private:
  Real matrix[9]{};
public:
  __host__ __device__ Matrix(){
    std::memset(matrix, 0, sizeof(Real) * 9);
  }

  template<UINT id1, UINT id2>
  __host__ __device__ INLINE inline Real  getValue() const{
    static constexpr UINT matID = id1*3 +id2;
    static_assert((matID < 9),"Wrong matID");
    return matrix[matID];
  }
  template<UINT id1, UINT id2>
  __host__ __device__ INLINE inline void  setValue(const Real & value) {
    static constexpr UINT matID = id1*3 +id2;
    static_assert((matID < 9),"Wrong matID");
    matrix[matID] = value;
  }
  __host__ __device__ INLINE inline void  setIdentity(){
    std::memset(matrix, 0, sizeof(Real) * 9);
    this->setValue<0,0>(1);
    this->setValue<1,1>(1);
    this->setValue<2,2>(1);
  }
  __host__ __device__ INLINE inline void  reset(){
    std::memset(matrix, 0, sizeof(Real) * 9);
  }

  __host__ __device__ INLINE inline void performMatrixMultiplication(const Matrix &A, const Matrix &B) {
    this->reset();
#pragma unroll 3
    for (int i = 0; i < 3; ++i) {
#pragma unroll 3
      for (int j = 0; j < 3; ++j) {
#pragma unroll 3
        for (int k = 0; k < 3; ++k) {
          this->matrix[i * 3 + j] += A.matrix[i * 3 + k] * B.matrix[k * 3 + j];
        }
      }
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

[[deprecated]]
__host__ static inline void computeInverseMatrix(const Real  matrixA [][3], Real  inverseMatrix [][3]){
  double det = matrixA[0][0] * (matrixA[1][1] * matrixA[2][2] - matrixA[2][1] * matrixA[1][2]) -
               matrixA[0][1] * (matrixA[1][0] * matrixA[2][2] - matrixA[1][2] * matrixA[2][0]) +
               matrixA[0][2] * (matrixA[1][0] * matrixA[2][1] - matrixA[1][1] * matrixA[2][0]);

  assert(not(FEQUALS(det,0.0)));
  Real invdet = 1. / det;

  inverseMatrix[0][0] = (matrixA[1][1] * matrixA[2][2] - matrixA[2][1] * matrixA[1][2]) * invdet;
  inverseMatrix[0][1] = (matrixA[0][2] * matrixA[2][1] - matrixA[0][1] * matrixA[2][2]) * invdet;
  inverseMatrix[0][2] = (matrixA[0][1] * matrixA[1][2] - matrixA[0][2] * matrixA[1][1]) * invdet;
  inverseMatrix[1][0] = (matrixA[1][2] * matrixA[2][0] - matrixA[1][0] * matrixA[2][2]) * invdet;
  inverseMatrix[1][1] = (matrixA[0][0] * matrixA[2][2] - matrixA[0][2] * matrixA[2][0]) * invdet;
  inverseMatrix[1][2] = (matrixA[1][0] * matrixA[0][2] - matrixA[0][0] * matrixA[1][2]) * invdet;
  inverseMatrix[2][0] = (matrixA[1][0] * matrixA[2][1] - matrixA[2][0] * matrixA[1][1]) * invdet;
  inverseMatrix[2][1] = (matrixA[2][0] * matrixA[0][1] - matrixA[0][0] * matrixA[2][1]) * invdet;
  inverseMatrix[2][2] = (matrixA[0][0] * matrixA[1][1] - matrixA[1][0] * matrixA[0][1]) * invdet;

}

__host__ static inline void computeInverseMatrix(const Matrix &  matrixA, Matrix & inverseMatrix){
  double det = matrixA.getValue<0,0>() * (matrixA.getValue<1,1>() * matrixA.getValue<2,2>() - matrixA.getValue<2,1>() * matrixA.getValue<1,2>()) -
               matrixA.getValue<0,1>() * (matrixA.getValue<1,0>() * matrixA.getValue<2,2>() - matrixA.getValue<1,2>() * matrixA.getValue<2,0>()) +
               matrixA.getValue<0,2>() * (matrixA.getValue<1,0>() * matrixA.getValue<2,1>() - matrixA.getValue<1,1>() * matrixA.getValue<2,0>());

  assert(not(FEQUALS(det,0.0)));
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

[[deprecated]]
__host__ static void inline performMatrixMultiplication(const Real  matA[][3], const Real  matB[][3], Real  mat[][3]){
  std::memset(mat,0, sizeof(Real)*9);
#pragma unroll 3
  for(int i = 0; i < 3; ++i)
#pragma unroll 3
    for(int j = 0; j < 3; ++j)
#pragma unroll 3
      for(int k = 0; k < 3; ++k)
      {
        mat[i][j] += matA[i][k] * matB[k][j];
      }
}

__host__ static void inline performMatrixMultiplication(const Matrix &  matA, const Matrix &  matB, Matrix &  mat) {
  mat.reset();
  mat.performMatrixMultiplication(matA,matB);

}



/**
 * @brief computes Matrix that will transform originalVec into transformed Vec
 *        [RotationMatrix X orignalVec = transformedVec]
 * @param [in/out] originalVec : the orginal 3D vector
 * @param [in] transformedVec  : 3D vector that the original Vec is transformed to.
 * @param [out] RotationMatrix : The resultant 3 X 3 rotation matrix
 */
[[deprecated]]
__host__ static void computeRotationMatrix(const Real3 & originalVec,const Real3 & transformedVec , Real  RotationMatrix [][3]){
  /**
   * Matlab code: (https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d)
   * G =  [ dot(A,B) -norm(cross(A,B)) 0;...
           norm(cross(A,B)) dot(A,B)  0;...
           0              0           1];

      F = [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];

      UU =  F*G*inv(F);
   */
  std::memset(RotationMatrix,0, sizeof(Real)*9);
  if((FEQUALS(originalVec.x,transformedVec.x)) and (FEQUALS(originalVec.y,transformedVec.y)) and (FEQUALS(originalVec.z,transformedVec.z))){
    RotationMatrix[0][0] = 1;
    RotationMatrix[1][1] = 1;
    RotationMatrix[2][2] = 1;
    return;
  }
  const Real dotProduct = computeDotProduct(originalVec,transformedVec);
  const Real3 & crossVec = computeCrossProduct(originalVec,transformedVec);
  const Real normCrossProduct = computeVecNorm(crossVec);

  Real G[3][3];
  G[0][0] = dotProduct;       G[0][1] = -normCrossProduct; G[0][2] = 0;
  G[1][0] = normCrossProduct; G[1][1] =  dotProduct;       G[1][2] = 0;
  G[2][0] = 0;                G[2][1] = 0;                 G[2][2] = 1;



  Real3 scaledVec;
  scaleVec(originalVec,scaledVec,dotProduct);
  const Real3 subVec{transformedVec.x - scaledVec.x,transformedVec.y - scaledVec.y,transformedVec.z - scaledVec.z};
  const Real normSubVec = computeVecNorm(subVec);

  Real F[3][3];
  F[0][0] = originalVec.x;  F[0][1] = subVec.x/normSubVec; F[0][2] = -crossVec.x;
  F[1][0] = originalVec.y;  F[1][1] = subVec.y/normSubVec; F[1][2] = -crossVec.y;
  F[2][0] = originalVec.z;  F[2][1] = subVec.z/normSubVec; F[2][2] = -crossVec.z;


  Real invF[3][3], temp1[3][3];
  computeInverseMatrix(F,invF);
  performMatrixMultiplication(G,invF,temp1);
  performMatrixMultiplication(F,temp1,RotationMatrix);

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
  performMatrixMultiplication(G,invF,temp1);
  performMatrixMultiplication(F,temp1,RotationMatrix);

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
 * @brief computes matrix times vector
 * @tparam transpose weather to use transpose of matrix
 * @param [in] matrix
 * @param [in] vec
 * @param [out] matVec
 */

template <bool transpose>
[[deprecated]]
__host__ static void doMatVec(const Real matrix[][3], const Real3 & vec, Real3 & matVec){
  if(transpose) {
    matVec.x = matrix[0][0] * vec.x + matrix[1][0] * vec.y + matrix[2][0] * vec.z;
    matVec.y = matrix[0][1] * vec.x + matrix[1][1] * vec.y + matrix[2][1] * vec.z;
    matVec.z = matrix[0][2] * vec.x + matrix[1][2] * vec.y + matrix[2][2] * vec.z;
  }
  else{
    matVec.x = matrix[0][0] * vec.x + matrix[0][1] * vec.y + matrix[0][2] * vec.z;
    matVec.y = matrix[1][0] * vec.x + matrix[1][1] * vec.y + matrix[1][2] * vec.z;
    matVec.z = matrix[2][0] * vec.x + matrix[2][1] * vec.y + matrix[2][2] * vec.z;
  }

}
template <bool transpose>
__host__ static void doMatVec(const Matrix & matrix, const Real3 & vec, Real3 & matVec){
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
}
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

  performMatrixMultiplication(rotationMatrixX,rotationMatrixK,rotationMatrix);


  return true;
}
__host__ bool static computeRotationMatrix(const Real3 & k, const Matrix & rotationMatrixK, Matrix & rotationMatrix, const Real & rotAngle){
  static constexpr Real3 X{1,0,0};
  Real3 shiftedX;
  doMatVec<false>(rotationMatrixK,X,shiftedX);
  Real3 rotatedX;
  Matrix rotationMatrixX;
  performRodriguesRotation(rotatedX,shiftedX,k,rotAngle);
  computeRotationMatrix(shiftedX,rotatedX,rotationMatrixX);
  performMatrixMultiplication(rotationMatrixX,rotationMatrixK,rotationMatrix);
}

__host__ inline static void normalizeVec(Real3 & vec){
  Real normVec =  computeVecNorm(vec);
  scaleVec(vec,vec,1./normVec);
}
[[deprecated]]
__host__ bool static computeRotationMatrixBaseConfiguration(const Real3 & k, Real rotationMatrix[][3]){

  static constexpr Real3 origK{0,0,1};
  Real rotationMatrixK[3][3];
  computeRotationMatrix(origK,k,rotationMatrixK);

#if DEBUG
  {
    Real3 shiftedK;
    doMatVec<false>(rotationMatrixK, origK, shiftedK);
    assert((FEQUALS(shiftedK.x, k.x)) and (FEQUALS(shiftedK.y, k.y)) and (FEQUALS(shiftedK.z, k.z)));
  }
#endif
  static constexpr Real3 X{1,0,0};
  static constexpr UINT numInterval = 100000;
  static constexpr Real dTheta = M_PI/(numInterval*1.0);
  Real3 shiftedX;
  doMatVec<false>(rotationMatrixK,X,shiftedX);
  Real3 rotatedX;
  Real maxDiff = std::numeric_limits<Real>::infinity();
  Real rotAngle = 0;
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
  Real rotationMatrixX[3][3];
  computeRotationMatrix(shiftedX,rotatedX,rotationMatrixX);

  performMatrixMultiplication(rotationMatrixX,rotationMatrixK,rotationMatrix);


  return true;

}

#endif //CY_RSOXS_ROTATION_H
