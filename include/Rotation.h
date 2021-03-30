//
// Created by maksbh on 3/30/21.
//

#ifndef CY_RSOXS_ROTATION_H
#define CY_RSOXS_ROTATION_H

#include <Datatypes.h>

/**
 * @brief computes dot product of two vectors
 * @param [in] vecA 3D vector A
 * @param [in] vecB 3D vector B
 * @return resultant dot product
 */
__host__ inline Real computeDotProduct(const Real3 & vecA, const Real3 & vecB){
  return (vecA.x*vecB.x + vecA.y*vecB.y + vecA.z*vecB.z);
}

/**
 * @brief computes cross product of two vectors
 * @param [in] vecA 3D vector A
 * @param [in] vecB 3D vector B
 * @return resultant cross product
 */
__host__ inline Real3 computeCrossProduct(const Real3 & vecA, const Real3 & vecB){
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
__host__ inline Real computeVecNorm(const Real3 & vec){
  return (sqrt(computeDotProduct(vec,vec)));
}
/**
 * @brief computes || A X B ||
 * @param [in] vecA 3D vector A
 * @param [in] vecB 3D vector B
 * @return norm of cross product
 */
__host__ inline Real computeNormCrossproduct(const Real3 & vecA, const Real3 & vecB){
  const Real3 & crossProduct = computeCrossProduct(vecA,vecB);
  return (sqrt(computeDotProduct(crossProduct,crossProduct)));
}

/**
 *
 * @param [in] vecIn input 3D vector
 * @param [out] scaledVec scaled vector
 * @param [in] scaleFactor scale Factor
 */
__host__ inline void scaleVec(const Real3 & vecIn, Real3 & scaledVec, const Real & scaleFactor){
  scaledVec.x = scaleFactor*vecIn.x;
  scaledVec.y = scaleFactor*vecIn.y;
  scaledVec.z = scaleFactor*vecIn.z;
}

__host__ inline void computeInverseMatrix(const Real  matrixA [][3], Real  inverseMatrix [][3]){
  double det = matrixA[0][0] * (matrixA[1][1] * matrixA[2][2] - matrixA[2][1] * matrixA[1][2]) -
               matrixA[0][1] * (matrixA[1][0] * matrixA[2][2] - matrixA[1][2] * matrixA[2][0]) +
               matrixA[0][2] * (matrixA[1][0] * matrixA[2][1] - matrixA[1][1] * matrixA[2][0]);

  assert(not(FEQUALS(det,0.0)));
  Real invdet = 1 / det;

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

__host__ void inline performMatrixMultiplication(const Real  matA[][3], const Real  matB[][3], Real  mat[][3]){
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



/**
 * @brief computes RotationMatrix that will transform originalVec into transformed Vec
 *        [RotationMatrix X orignalVec = transformedVec]
 * @param [in/out] originalVec : the orginal 3D vector
 * @param [in] transformedVec  : 3D vector that the original Vec is transformed to.
 * @param [out] RotationMatrix : The resultant 3 X 3 rotation matrix
 */
__host__ void computeRotationMatrix(const Real3 & originalVec,const Real3 & transformedVec , Real  RotationMatrix [][3]){
  /**
   * Matlab code: (https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d)
   * G =  [ dot(A,B) -norm(cross(A,B)) 0;...
           norm(cross(A,B)) dot(A,B)  0;...
           0              0           1];

      F = [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];

      UU =  F*G*inv(F);
   */
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
 * @brief compute Rodrigues rotation: Rotates a vector inVec about an axis with an angle
 * (https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula)
 * @param [out] rotatedVec rotated Vector
 * @param [in] inVec input vector
 * @param [in] axis
 * @param [in] angle
 *
 */
__host__ void inline performRodriguesRotation(Real3  & rotatedVec, const Real3 & inVec,const Real3 & axis, const Real & angle){
  const Real cosAlpha = cos(angle);
  const Real sinAlpha = sin(angle);
  const Real3 & crossProduct = computeCrossProduct(inVec,axis);
  const Real & dotProduct = computeDotProduct(axis,inVec);
  Real3 term1{cosAlpha*inVec.x,cosAlpha*inVec.y,cosAlpha*inVec.z};
  Real3 term2{sinAlpha*crossProduct.x,sinAlpha*crossProduct.y,sinAlpha*crossProduct.z};
  Real3 term3{axis.x,axis.y,axis.z};
  scaleVec(term3,term3,(1-cosAlpha)*dotProduct);
  rotatedVec.x = term1.x + term2.x + term3.x;
  rotatedVec.y = term1.y + term2.y + term3.y;
  rotatedVec.z = term1.z + term2.z + term3.z;
}


#endif //CY_RSOXS_ROTATION_H
