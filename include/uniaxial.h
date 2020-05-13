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

#ifndef UNIAXIAL_H
#define UNIAXIAL_H
#include <Datatypes.h>
#include "cudaHeaders.h"
#include "cudaUtils.h"
#include <Input/Input.h>
#include <math.h>
#include <complex.h>

#ifdef EOC
#include <opencv2/opencv.hpp>
#endif

#include <cstring>

#include <omp.h>
#include <Output/writeVTI.h>



/**
 * @brief This function computes the polarization in real space for the uniaxial case.
 * @param [in] material material data for a particular energy level under consideration.
 * @param [in] angle The angle of rotation.
 * @param [in] voxelInput voxel Input data.
 * @param [in] threadID threadID to access the index of the entry.
 * @param [out] polarizationX
 * @param [out] polarizationY
 * @param [out] polarizationZ
 */

__device__ void computePolarizationUniaxial(const Material<NUM_MATERIAL> *material, const Real angle,
                                            const Voxel<NUM_MATERIAL> *voxelInput, const BigUINT threadID,
                                            Complex *polarizationX, Complex *polarizationY, Complex *polarizationZ) {

  Complex pX, pY, pZ, npar[NUM_MATERIAL], nper[NUM_MATERIAL];
  Real4 s1[NUM_MATERIAL];

  pX.x = 0;
  pX.y = 0;
  pY.x = 0;
  pY.y = 0;
  pZ.x = 0;
  pZ.y = 0;

  static constexpr Real OneBy4Pi = static_cast<Real> (1.0 / (4.0 * M_PI));
  Real temp;
  for (int i = 0; i < NUM_MATERIAL; i++) {
    s1[i] = voxelInput[threadID].s1[i];
    temp = cos(angle) * voxelInput[threadID].s1[i].y - sin(angle) * voxelInput[threadID].s1[i].z;
    s1[i].z = sin(angle) * s1[i].y + cos(angle) * s1[i].z;
    s1[i].y = temp;
    npar[i] = material->npara[i];
    nper[i] = material->nperp[i];
  }

  Complex nsum;
  for (int i = 0; i < NUM_MATERIAL; i++) {
    Real phi = s1[i].w + s1[i].x * s1[i].x + s1[i].y * s1[i].y + s1[i].z * s1[i].z;

    nsum.x = npar[i].x + 2 * nper[i].x;
    nsum.y = npar[i].y + 2 * nper[i].y;

    computeComplexSquare(nsum);
    computeComplexSquare(npar[i]);
    computeComplexSquare(nper[i]);

    pX.x += (npar[i].x - nper[i].x) * s1[i].z * s1[i].x;
    pX.y += (npar[i].y - nper[i].y) * s1[i].z * s1[i].x;

    pY.x += (npar[i].x - nper[i].x) * s1[i].z * s1[i].y;
    pY.y += (npar[i].y - nper[i].y) * s1[i].z * s1[i].y;

    pZ.x +=
        (s1[i].w * nsum.x) / 9.0 + npar[i].x * s1[i].z * s1[i].z + nper[i].x * (s1[i].x * s1[i].x + s1[i].y * s1[i].y)
            - phi;
    pZ.y +=
        (s1[i].w * nsum.y) / 9.0 + npar[i].y * s1[i].z * s1[i].z + nper[i].y * (s1[i].x * s1[i].x + s1[i].y * s1[i].y);
  }

  pX.x *= OneBy4Pi;
  pX.y *= OneBy4Pi;

  pY.x *= OneBy4Pi;
  pY.y *= OneBy4Pi;

  pZ.x *= OneBy4Pi;
  pZ.y *= OneBy4Pi;

  polarizationX[threadID] = pX;
  polarizationY[threadID] = pY;
  polarizationZ[threadID] = pZ;

}

__device__ inline BigUINT reshape3Dto1D(UINT i, UINT j, UINT k, uint3 voxel) {
  return ((i) + (j) * voxel.x + (k) * voxel.x * voxel.y);
}

__device__ inline void reshape1Dto3D(BigUINT id, UINT & X, UINT & Y, UINT & Z, uint3 voxel) {
  Z = static_cast<UINT>(floorf(id / (voxel.y * voxel.x * 1.0)));
  Y = static_cast<UINT>(floorf((id - Z * voxel.y * voxel.x) / voxel.x * 1.0));
  X = static_cast<UINT>(id - Y * voxel.x - Z * voxel.y * voxel.x);
}

template<typename T>
__device__ inline void swap(T &var1, T &var2) {
  T temp;
  temp = var1;
  var1 = var2;
  var2 = temp;
}





template<typename T>
__device__ void FFTShift(T *array,  uint3 voxel) {
  uint3 mid;
  mid.x = voxel.x / 2;
  mid.y = voxel.y / 2;
  mid.z = voxel.z / 2;
  UINT threadID = threadIdx.x + blockIdx.x * blockDim.x;
  BigUINT totalSize = (voxel.x * voxel.y * voxel.z);
  if (threadID >= totalSize) {
    return;
  }

  UINT X,Y,Z;
  reshape1Dto3D(threadID,X,Y,Z,voxel);
  if(voxel.z == 1){
    if ((X < mid.x) and (Y < mid.y)){
      BigUINT copyID = reshape3Dto1D(X + mid.x, Y + mid.y, Z + mid.z, voxel);
      BigUINT currID = threadID;
      swap(array[currID], array[copyID]);

      copyID = reshape3Dto1D(mid.x + X, Y, mid.z + Z, voxel);
      currID = reshape3Dto1D(X, mid.y + Y, Z, voxel);
      swap(array[currID], array[copyID]);
    }
    return;
  }
  if ((X < mid.x) and (Y < mid.y) and (Z < mid.z)) {
    BigUINT copyID = reshape3Dto1D(X + mid.x, Y + mid.y, Z + mid.z, voxel);
    BigUINT currID = threadID;
    swap(array[currID], array[copyID]);

    copyID = reshape3Dto1D(mid.x + X, mid.y + Y, Z, voxel);
    currID = reshape3Dto1D(X, Y, mid.z + Z, voxel);
    swap(array[currID], array[copyID]);

    copyID = reshape3Dto1D(mid.x + X, Y, mid.z + Z, voxel);
    currID = reshape3Dto1D(X, mid.y + Y, Z, voxel);
    swap(array[currID], array[copyID]);

    copyID = reshape3Dto1D(X, mid.y + Y, mid.z + Z, voxel);
    currID = reshape3Dto1D(mid.x + X, Y, Z, voxel);
    swap(array[currID], array[copyID]);

    copyID = reshape3Dto1D(mid.x + X, mid.y + Y, mid.z + Z, voxel);
    currID = reshape3Dto1D(X, Y, mid.z + Z, voxel);
    swap(array[currID], array[copyID]);
  }


}


/**
 * @brief This function computes the Scatter3D computation. The current implementation
 * assumes that the Fourier transform of polarization is shifted (in the default order
 * given by Matlab or cuFFT). This is different from the Igor. This offset is taken care
 * on the fly during computation of scatter3D.
 *
 *
 * @param [in] polarizationX  X component of polarization in Fourier space
 * @param [in] polarizationY  Y component of polarization in Fourier space
 * @param [in] polarizationZ  Z component of polarization in Fourier space
 * @param [out] Scatter3D     Scatter 3D result
 * @param [in] elefield       Electric field
 * @param [in] voxelNum       Number of total voxel
 * @param [in] voxel          Number of voxel in each direciton.
 * @param [in] physSize       Physical Size
 */
__global__ void computeScatter3D(Complex *polarizationX,
                                 Complex *polarizationY,
                                 Complex *polarizationZ,
                                 Real *Scatter3D,
                                 const ElectricField elefield,
                                 const BigUINT voxelNum,
                                 const uint3 voxel,
                                 const Real physSize) {



  FFTShift(polarizationX,voxel);
  FFTShift(polarizationY,voxel);
  FFTShift(polarizationZ,voxel);


  UINT threadID = threadIdx.x + blockIdx.x * blockDim.x;
  if (threadID >= voxelNum) {
    return;
  }

  const Complex pX = polarizationX[threadID];
  const Complex pY = polarizationY[threadID];
  const Complex pZ = polarizationZ[threadID];

  UINT X,Y,Z;
  reshape1Dto3D(threadID,X,Y,Z,voxel);

  UINT midX = static_cast<UINT>(voxel.x / 2.0);
  UINT midY = static_cast<UINT>(voxel.y / 2.0);
  UINT midZ = static_cast<UINT>(voxel.z / 2.0);

  Real3 dx;
  dx.x = static_cast<Real>((2 * M_PI / physSize) / ((voxel.x - 1) * 1.0));
  dx.y = static_cast<Real>((2 * M_PI / physSize) / ((voxel.y - 1) * 1.0));

#if ENABLE_2D
  dx.z = 0;
#else
  dx.z = static_cast<Real>((2 * M_PI / physSize) / ((voxel.z - 1) * 1.0));
#endif
  Real3 q;
  q.x = static_cast<Real>((-M_PI / physSize) + X * dx.x);
  q.y = static_cast<Real>((-M_PI / physSize) + Y * dx.y);
#if ENABLE_2D
  q.z = 0;
#else
  q.z = static_cast<Real>((-M_PI / physSize) + Z * dx.z);
#endif
//  if (X <= midX) {
//    q.x = static_cast<Real>((-M_PI / physSize) + (midX - X) * dx.x);
//  } else {
//    q.x = static_cast<Real>((M_PI / physSize) - (X - (midX + 1)) * dx.x);
//  }
//
//  if (Y <= midY) {
//    q.y = static_cast<Real>((-M_PI / physSize) + (midY - Y) * dx.y);
//  } else {
//    q.y = static_cast<Real>((M_PI / physSize) - (Y - (midY + 1)) * dx.y);
//  }
//#if ENABLE_2D
//  q.z = 0;
//#else
//  if (Z <= midZ) {
//    q.z = static_cast<Real>((-M_PI / physSize) + (midZ - Z) * dx.z);
//  } else {
//    q.z = static_cast<Real>((M_PI / physSize) - (Z - (midZ + 1)) * dx.z);
//  }
//#endif

  Complex val;
  Real res;
  // q.z <-> q.x from original Igor code
  val.x = q.z * (2 * elefield.k.z + q.z) * pX.x
      + q.y * (elefield.k.z + q.z) * pY.x + q.x * (elefield.k.z + q.z) * pZ.x;
  val.y = q.z * (2 * elefield.k.z + q.z) * pX.y
      + q.y * (elefield.k.z + q.z) * pY.y + q.x * (elefield.k.z + q.z) * pZ.y;
  res = val.x * val.x + val.y * val.y;

  val.x = q.y * (elefield.k.z + q.z) * pX.x
      + (q.y * q.y - elefield.k.z * elefield.k.z) * pY.x + q.y * q.x * pZ.x;
  val.y = q.y * (elefield.k.z + q.z) * pX.y
      + (q.y * q.y - elefield.k.z * elefield.k.z) * pY.y + q.y * q.x * pZ.y;
  res += val.x * val.x + val.y * val.y;

  val.x = q.x * (elefield.k.z + q.z) * pX.x + (q.y * q.x) * pY.x
      + (q.x * q.x - elefield.k.z * elefield.k.z) * pZ.x;
  val.y = q.x * (elefield.k.z + q.z) * pX.y + (q.y * q.x) * pY.y
      + (q.x * q.x - elefield.k.z * elefield.k.z) * pZ.y;

  res += val.x * val.x + val.y * val.y;
  Scatter3D[threadID] = res;

}

/**
 * This function computes the equivalent id of the given position.
 * The id is calculated keeping in mind the effect of polarization FFT
 * shift done in cufft calculation with respect to Igor.
 * @param [in] pos position of which the index is to be found.
 * @param [in] i corresponding X coordinate
 * @param [in] j correspondi    ng Y coordinate
 * @param [in] start start position
 * @param [in] dx the grid podition in each direction
 * @param [in] voxel voxel dimension
 * @return the equivalent id for the given position.
 */
__host__ __device__ BigUINT computeEquivalentID(const Real3 pos,
                                                const UINT i,
                                                const UINT j,
                                                const Real start,
                                                const Real3 dx,
                                                const uint3 voxel) {

//  uint3 mid;
//  mid.x = voxel.x / 2;
//  mid.y = voxel.y / 2;
//  mid.z = voxel.z / 2;
//
//  uint3 index;
//  if (i <= mid.x) {
//    index.x = mid.x - i;
//  } else {
//    index.x = (voxel.x - 1) - (i - mid.x - 1);
//  }
//
//  if (j <= mid.y) {
//    index.y = mid.y - j;
//  } else {
//    index.y = (voxel.y - 1) - (j - mid.y - 1);
//  }
//#if (ENABLE_2D)
//  index.z = 0;
//#else
//  UINT k = static_cast<UINT >(round((pos.z - start) / (dx.z)));
//  if (k <= mid.z) {
//    index.z = mid.z - k;
//  } else {
//    index.z = (voxel.z - 1) - (k - mid.z - 1);
//  }
//#endif

  UINT k;
  if(voxel.z == 1) {
    k = 0;
  } else{
    k = static_cast<UINT >(round((pos.z - start) / (dx.z)));
  }

  BigUINT id = k * voxel.x * voxel.y + j * voxel.x + i;
  return id;

}

/**
 * @brief This function computes the equivalent Projection of 3D Scatter
 * field on Ewald's sphere on GPU. In this case we assume a constant k in
 * the z - direction. In future, we would rotate the k vector in plain
 * to gain more statistics.
 *
 * @param [out] projection The projection result
 * @param [in] scatter3D Scatter3D result.
 * @param [in] voxel Number of voxel in each direction
 * @param [in] k Electric field k.
 * @param [in] physSize Physical Size.
 */

__global__ void computeEwaldProjectionGPU(Real *projection,
                                          const Real *scatter3D,
                                          const uint3 voxel,
                                          const Real k,
                                          const Real physSize) {
  UINT threadID = threadIdx.x + blockIdx.x * blockDim.x;
  const int totalSize = voxel.x * voxel.y;
  if (threadID >= totalSize) {
    return;
  }
  Real val, start;
  Real3 dx, pos;
  start = -static_cast<Real>(M_PI / physSize);

  dx.x = static_cast<Real>((2.0 * M_PI / physSize) / ((voxel.x - 1) * 1.0));
  dx.y = static_cast<Real>((2.0 * M_PI / physSize) / ((voxel.y - 1) * 1.0));
#if (ENABLE_2D)
  dx.z = 0;
#else
  dx.z = static_cast<Real>((2.0 * M_PI / physSize) / ((voxel.z - 1) * 1.0));
#endif
  UINT Y = static_cast<UINT>(floorf(threadID) / (voxel.x * 1.0));
  UINT X = static_cast<UINT>(threadID - Y * voxel.x);
  pos.y = start + Y * dx.y;
  pos.x = start + X * dx.x;
  val = k * k - pos.x * pos.x - pos.y * pos.y;
  if (val >= 0) {

    pos.z = -k + sqrt(val);
    BigUINT id = computeEquivalentID(pos, X, Y, start, dx, voxel);
    projection[threadID] = scatter3D[id];
  } else {
    projection[threadID] = NAN;
  }

}

#ifdef EOC

/**
 * @brief This function computes the equivalent Projection of 3D Scatter
 * field on Ewald's sphere on GPU. In this case we assume a constant k in
 * the z - direction. In future, we would rotate the k vector in plain
 * to gain more statistics.
 * @param [out] projection The projection result
 * @param [in] scatter3D Scatter3D result.
 * @param [in] voxel Number of voxel in each direction
 * @param [in] k Electric field k.
 */

__host__ void computeEwaldProjectionCPU(Real *projection, const Real *scatter3D, const uint3 voxel, const Real k) {
  Real val, start;
  Real3 dx, pos;
  start = -static_cast<Real>(M_PI / PHYSSIZE);

  dx.x = static_cast<Real>((2.0 * M_PI / PHYSSIZE) / ((voxel.x - 1) * 1.0));
  dx.y = static_cast<Real>((2.0 * M_PI / PHYSSIZE) / ((voxel.y - 1) * 1.0));
#if (ENABLE_2D)
  dx.z = 0;
#else
  dx.z = static_cast<Real>((2.0 * M_PI / PHYSSIZE) / ((voxel.z - 1) * 1.0));
#endif

//#pragma omp parallel for
  for (UINT j = 0; j < voxel.y; j++) {
    pos.y = start + j * dx.y;
    for (UINT i = 0; i < voxel.x; i++) {
      pos.x = start + i * dx.x;
      val = k * k - pos.x * pos.x - pos.y * pos.y;
      if (val >= 0) {

        pos.z = -k + sqrt(val);
        BigUINT id = computeEquivalentID(pos, i, j, start, dx, voxel);
        projection[j * voxel.x + i] = scatter3D[id];
      } else {
        projection[j * voxel.x + i] = NAN;
      }
    }
  }
}

/**
 * This function rotates the image on CPU. The rotation is performed with the library openCV.
 * @param [in,out] projection The rotated projection image
 * @param [in] voxel voxel data
 * @param [in] rotAngle the rotation angle
 */
__host__ void rotateImage(Real *projection, const uint3 voxel, const Real rotAngle){

  const int size[2]{static_cast<int>(voxel.x), static_cast<int>(voxel.y)};

  cv::Mat image(voxel.x,voxel.y,CV_32FC1);
  memcpy(image.data, projection,size[0]*size[1]*sizeof(Real));
  cv::Point2f pc(image.cols/2., image.rows/2.);
  cv::Mat r = cv::getRotationMatrix2D(pc, rotAngle, 1.0);
  cv::Mat rotatedImage;
  warpAffine(image, rotatedImage, r, image.size());
  if (rotatedImage.isContinuous()) {
    std::memcpy(projection,rotatedImage.data,rotatedImage.rows*rotatedImage.cols);
  }
  else{
    throw "Memory not continuous.";
  }
}

#endif
/**
 * @brief This function performs the integration of the projection data.
 * Currently inactive.
 *
 * @param integrate The integrate value.
 * @param angle (min angle , max angle) in radians
 * @param projection Projection input
 * @param voxel 2D voxel size in x and y
 * @param physSize Physical Size
 */
__global__ void radialIntegrate(Real *integrate,
                                const Real2 angle,
                                const Real *projection,
                                const uint2 voxel,
                                const Real physSize) {

  UINT threadID = threadIdx.x + blockIdx.x * blockDim.x;
  const int totalSize = voxel.x * voxel.y;
  if (threadID >= totalSize) {
    return;
  }

  UINT Y = static_cast<UINT>(floorf(threadID) / (voxel.x * 1.0));
  UINT X = static_cast<UINT>(threadID - Y * voxel.x);
  Real2 dx, pos;
  Real start, val;
  start = -static_cast<Real>(M_PI / physSize);
  dx.x = static_cast<Real>((2.0 * M_PI / physSize) / ((voxel.x - 1) * 1.0));
  dx.y = static_cast<Real>((2.0 * M_PI / physSize) / ((voxel.y - 1) * 1.0));

  pos.x = start + Y * dx.x;
  pos.y = start + X * dx.y;

  bool check1 = isnan(projection[threadID]); // if projection is not nan

  val = atan(-pos.y / -pos.x);

  bool check2;
  {

    bool check2_1 = val > angle.x;
    bool check2_2 = val < angle.y;

    check2 = check2_1 and check2_2;
  }
  bool check3;
  {
    bool check3_1 = val > (angle.x - M_PI);
    bool check3_2 = val < (angle.y - M_PI);
    check3 = check3_1 and check3_2;
  }

  bool check = not(check1) and (check2 or check3);

  Real radialDist = sqrt(pos.y * pos.y + pos.x * pos.x);

  if (check) {
    integrate[threadID] = radialDist * projection[threadID];
  } else {
    integrate[threadID] = 0;
  }
}


/**
 * This function dumps the 2D files in parallel.
 * This function is tested with dumping the file for DLEVEL2
 * @param data data to be dumped. The data is in the format [ Energy1 ; Energy2; ..... :Energyn]
 * @param numData number of data
 * @param eachDataSize The size of each data in 2D.
 */

__host__ void dump_files2D(const Real *data, const UINT numData, const UINT *eachDataSize) {
#pragma omp parallel
  {

    UINT threadID = omp_get_thread_num();
    UINT chunkSize = static_cast<UINT>(std::ceil(numData * 1.0 / (omp_get_num_threads() * 1.0)));
    const UINT totalSize = eachDataSize[0] * eachDataSize[1];
    Real *projectionLocal = new Real[totalSize];
    UINT startID = threadID * chunkSize;
    UINT endID = ((threadID + 1) * chunkSize);
    if (startID < numData) {

      for (UINT csize = startID; csize < std::min(endID, numData); csize++) {

        std::memcpy(projectionLocal, &data[csize * totalSize], sizeof(Real) * totalSize);

        std::string dirname = "VTI/";
        std::string fnameFP = dirname + "projectionAverage" + std::to_string(csize);
        VTI::writeDataScalar2DFP(projectionLocal, eachDataSize, fnameFP.c_str(), "projection");

      }
    }
    delete[] projectionLocal;
  }
}
#endif //WRITER_UNIAXIAL_H
