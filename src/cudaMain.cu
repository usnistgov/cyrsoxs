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


#include <cudaMain.h>
#include <Input/Input.h>
#include <bits/ios_base.h>
#include <math.h>
#include <cuda_runtime.h>
#include <Output/writeVTI.h>
#include <uniaxial.h>
#include <cublas_v2.h>
#include <chrono>
#include <ctime>
#include <chrono>
#include <npp.h>
#include <Output/outputUtils.h>

#define START_TIMER(X) timerArrayStart[X] = std::chrono::high_resolution_clock::now();
#define END_TIMER(X) timerArrayEnd[X] = std::chrono::high_resolution_clock::now(); \
                     timings[X] +=  (static_cast<std::chrono::duration<Real>>(timerArrayEnd[X] - timerArrayStart[X])).count();


int warmup() {
  double *d_warmup, *warmup;
  warmup = new double[1000];
  CUDA_CHECK_RETURN(cudaMalloc((void **) &d_warmup, sizeof(double) * 1000));
  gpuErrchk(cudaPeekAtLastError());
  CUDA_CHECK_RETURN(cudaMemcpy(d_warmup, warmup, sizeof(double) * 1000, cudaMemcpyHostToDevice));
  gpuErrchk(cudaPeekAtLastError());
  cudaFree(d_warmup);
  delete[] warmup;
  return EXIT_SUCCESS;
}

__host__ int performFFTShift(Complex *polarization, const UINT &blockSize, const uint3 &vx) {
  FFTIgor<<<blockSize, NUM_THREADS>>>(polarization, vx);
  cudaDeviceSynchronize();
  gpuErrchk(cudaPeekAtLastError());
  return EXIT_SUCCESS;
}

__host__  int performScatter3DComputation(const Complex *d_polarizationX, const Complex *d_polarizationY,
                                          const Complex *d_polarizationZ,
                                          Real *d_scatter3D,
                                          const ElectricField &eleField,
                                          const Real &eAngle,
                                          const Real &kAngle,
                                          const BigUINT &voxelSize,
                                          const uint3 &vx,
                                          const Real &physSize,
                                          const bool &enable2D,
                                          const UINT &blockSize,
                                          const Real3 & kVector) {
  computeScatter3D <<< blockSize, NUM_THREADS >>>(d_polarizationX, d_polarizationY, d_polarizationZ,
                                                  d_scatter3D, eleField, eAngle, kAngle, voxelSize, vx,
                                                  physSize,
                                                  enable2D, kVector);
  cudaDeviceSynchronize();
  gpuErrchk(cudaPeekAtLastError());
  return EXIT_SUCCESS;
}

__host__ int peformEwaldProjectionGPU(Real *d_projection,
                                      const Real *d_scatter,
                                      const Real &k,
                                      const uint3 &vx,
                                      const Real &eAngle,
                                      const Real &kAngle,
                                      const Real &physSize,
                                      const Interpolation::EwaldsInterpolation &interpolation,
                                      const bool &enable2D,
                                      const UINT &blockSize,
                                      const Real3 & kVector) {
  computeEwaldProjectionGPU <<< blockSize, NUM_THREADS >>>(d_projection, d_scatter, vx,
                                                           k, eAngle, kAngle, physSize,
                                                           interpolation,
                                                           enable2D,kVector);
  cudaDeviceSynchronize();
  gpuErrchk(cudaPeekAtLastError());

  return EXIT_SUCCESS;


}

__host__ int peformEwaldProjectionGPU(Real *d_projection,
                                      const Complex *d_polarizationX, const Complex *d_polarizationY,
                                      const Complex *d_polarizationZ,
                                      const Real &k,
                                      const uint3 &vx,
                                      const Real &eAngle,
                                      const Real &kAngle,
                                      const Real &physSize,
                                      const Interpolation::EwaldsInterpolation &interpolation,
                                      const bool &enable2D,
                                      const UINT &blockSize,
                                      const Real3 & kVector) {
  computeEwaldProjectionGPU <<< blockSize, NUM_THREADS >>>(d_projection, d_polarizationX, d_polarizationY,
                                                           d_polarizationZ, vx,
                                                           k, eAngle, kAngle, physSize,
                                                           interpolation,
                                                           enable2D,kVector);
  cudaDeviceSynchronize();
  gpuErrchk(cudaPeekAtLastError());
  return EXIT_SUCCESS;

}

template<ReferenceFrame referenceFrame>
__global__ void computePolarization(Material<NUM_MATERIAL> materialInput,
                                    const Voxel<NUM_MATERIAL> *voxelInput,
                                    const ElectricField elefield,
                                    const Real angle,
                                    const uint3 voxel,
                                    Complex *polarizationX,
                                    Complex *polarizationY,
                                    Complex *polarizationZ,
                                    FFT::FFTWindowing windowing,
                                    const bool enable2D,
                                    const MorphologyType morphologyType,
                                    const Matrix rotationMatrix
) {
  UINT threadID = threadIdx.x + blockIdx.x * blockDim.x;
#ifndef BIAXIAL
  if (morphologyType == MorphologyType::VECTOR_MORPHOLOGY) {
    computePolarizationVectorMorphologyOptimized<referenceFrame>(&materialInput, angle, voxelInput, threadID, polarizationX,
                                                 polarizationY,
                                                 polarizationZ,rotationMatrix);
  } else {
    computePolarizationEulerAngles(&materialInput, angle, voxelInput, threadID, polarizationX, polarizationY,
                                   polarizationZ);
  }
#else
  printf("Kernel not spported\n");
#endif


  if (windowing == FFT::FFTWindowing::HANNING) {
    UINT Z = static_cast<UINT>(threadID / (voxel.y * voxel.x * 1.0));
    UINT Y = static_cast<UINT>((threadID - Z * voxel.y * voxel.x) / (voxel.x * 1.0));
    UINT X = static_cast<UINT>(threadID - Y * voxel.x - Z * voxel.y * voxel.x);
    Real3 hanningWeight;
    hanningWeight.x = static_cast<Real> (0.5 * (1 - cos(2 * M_PI * X / (voxel.x))));
    hanningWeight.y = static_cast<Real> (0.5 * (1 - cos(2 * M_PI * Y / (voxel.y))));
    hanningWeight.z = static_cast<Real>(1.0);
    if (not(enable2D)) {
      hanningWeight.z = static_cast<Real>(0.5 * (1 - cos(2 * M_PI * Z / (voxel.z))));
    }
    Real totalHanningWeight = hanningWeight.x * hanningWeight.y * hanningWeight.z;
    polarizationX[threadID].x *= totalHanningWeight;
    polarizationX[threadID].y *= totalHanningWeight;
    polarizationY[threadID].x *= totalHanningWeight;
    polarizationY[threadID].y *= totalHanningWeight;
    polarizationZ[threadID].x *= totalHanningWeight;
    polarizationZ[threadID].y *= totalHanningWeight;

  }

}

__host__ int computePolarization(const Material<NUM_MATERIAL> &materialInput,
                                 const Voxel<NUM_MATERIAL> *d_voxelInput,
                                 const ElectricField &elefield,
                                 const Real &angle,
                                 const uint3 &vx,
                                 Complex *d_polarizationX,
                                 Complex *d_polarizationY,
                                 Complex *d_polarizationZ,
                                 FFT::FFTWindowing windowing,
                                 const bool &enable2D,
                                 const MorphologyType &morphologyType,
                                 const UINT &blockSize,
                                 const ReferenceFrame & referenceFrame,
                                 const Matrix & rotationMatrix
) {
  if(referenceFrame == ReferenceFrame::MATERIAL) {
    computePolarization<ReferenceFrame::MATERIAL><<< blockSize, NUM_THREADS >>>(materialInput, d_voxelInput, elefield,
                                                                                angle, vx, d_polarizationX,
                                                                                d_polarizationY, d_polarizationZ,
                                                                                windowing,
                                                                                enable2D,
                                                                                morphologyType,rotationMatrix);
  }
  else  {
    computePolarization<ReferenceFrame::LAB><<< blockSize, NUM_THREADS >>>(materialInput, d_voxelInput, elefield,
                                                                                angle, vx, d_polarizationX,
                                                                                d_polarizationY, d_polarizationZ,
                                                                                windowing,
                                                                                enable2D,
                                                                                morphologyType,rotationMatrix);
  }
  cudaDeviceSynchronize();
  gpuErrchk(cudaPeekAtLastError());
  return EXIT_SUCCESS;
}


int cudaMain(const UINT *voxel,
             const InputData &idata,
             const std::vector<Material<NUM_MATERIAL> > &materialInput,
             Real *projectionGPUAveraged,
             const Voxel<NUM_MATERIAL> *voxelInput) {


  if ((static_cast<uint64_t>(voxel[0]) * voxel[1] * voxel[2]) > std::numeric_limits<BigUINT>::max()) {
    std::cout << "Exiting. Compile by Enabling 64 Bit indices\n";
    exit(EXIT_FAILURE);
  }

  const BigUINT numVoxels = voxel[0] * voxel[1] * voxel[2]; /// Voxel size
  const UINT numVoxel2D = voxel[0] * voxel[1];
  const uint3 vx{voxel[0], voxel[1], voxel[2]};
  const UINT
    numAnglesRotation = static_cast<UINT>(std::round((idata.endAngle - idata.startAngle) / idata.incrementAngle + 1));
  const UINT &numEnergyLevel = idata.energies.size();


  int num_gpu;
  cudaGetDeviceCount(&num_gpu);
  std::cout << "Number of CUDA devices:" << num_gpu << "\n";

  if (num_gpu < 1) {
    std::cout << "No GPU found. Exiting" << "\n";
    return (EXIT_FAILURE);
  }


#ifdef PROFILING
  enum TIMERS:UINT{
    MALLOC = 0,
    MEMCOPY_CPU_GPU = 1,
    POLARIZATION = 2,
    FFT = 3,
    SCATTER3D = 4,
    IMAGE_ROTATION=5,
    MEMCOPY_GPU_CPU = 6,
    ENERGY=7,
    MAX = 8
  };
  static const char *timersName[]{"Malloc on CPU + GPU",
                                  "Memcopy CPU -> GPU",
                                  "Polarization",
                                  "FFT",
                                  "Scatter3D + Ewalds",
                                  "Rotation",
                                  "Memcopy GPU -> CPU",
                                  "Total time "};
  static_assert(sizeof(timersName) / sizeof(char*) == TIMERS::MAX,
                "sizes dont match");
  std::array<std::chrono::high_resolution_clock::time_point,TIMERS::MAX> timerArrayStart;
  std::array<std::chrono::high_resolution_clock::time_point,TIMERS::MAX> timerArrayEnd;
  std::array<Real,TIMERS::MAX> timings{};
  timings.fill(0.0);

#endif

#ifdef DUMP_FILES
  createDirectory("Polarize");
  createDirectory("FFT");
  createDirectory("Scatter");

  /** Writing VTI files as a cross check **/

#if (NUM_MATERIAL==2)
  const char * varnameVector[2] = {"material1_s1","material2_s1"};
#elif (NUM_MATERIAL==4)
  const char * varnameVector[4] = {"material1_s","material2_s","material3_s","material4_s"};
#endif
#if (NUM_MATERIAL==2)
  const char * varnameScalar[2] = {"phi0","phi1"};
#elif (NUM_MATERIAL==4)
  const char * varnameScalar[4] = {"phi0","phi1", "phi2", "phi3"};
#endif

  VTI::writeVoxelDataVector(voxelInput, voxel, "S1", varnameVector);
  VTI::writeVoxelDataScalar(voxelInput, voxel, "Phi", varnameScalar);
#endif

  omp_set_num_threads(num_gpu);
#pragma omp parallel
  {


    cudaSetDevice(omp_get_thread_num());
    cudaDeviceProp dprop;
    cudaGetDeviceProperties(&dprop, omp_get_thread_num());

#ifdef PROFILING
    if(warmup() == EXIT_SUCCESS){
      std::cout << "Warmup completed on GPU " << dprop.name << "\n";
    }
    else{
      std::cout << "Warmup failed on GPU " << dprop.name << "\n";
#pragma omp cancel parallel
      exit (EXIT_FAILURE);
    }
#endif

    cufftResult result;
    cufftHandle plan;
#ifdef DOUBLE_PRECISION
    cufftPlan3d(&plan, voxel[2], voxel[1], voxel[0], CUFFT_Z2Z);
#else
    cufftPlan3d(&plan, voxel[2], voxel[1], voxel[0], CUFFT_C2C);
#endif
    cublasHandle_t handle;
    cublasStatus_t stat;
    cublasCreate(&handle);

    NppiSize sizeImage;
    sizeImage.height = voxel[0];
    sizeImage.width = voxel[1];

    NppiRect rect;
    rect.height = voxel[0];
    rect.width = voxel[1];
    rect.x = 0;
    rect.y = 0;


    const UINT ompThreadID = omp_get_thread_num();
    const UINT numEnergyPerGPU = static_cast<UINT>(std::ceil(numEnergyLevel * 1.0 / num_gpu));
    const UINT numStart = (numEnergyPerGPU * ompThreadID);
    UINT numEnd = (numEnergyPerGPU * (ompThreadID + 1));
    numEnd = std::min(numEnd, numEnergyLevel);

    const Real &energyStart = numStart < numEnergyLevel ? idata.energies[numStart] : 0;
    const Real &energyEnd = idata.energies[numEnd - 1];

    if (numStart >= numEnergyLevel) {
      std::cout << "[INFO] [GPU = " << dprop.name << "] -> No computation. Idle\n";
    } else {
      std::cout << "[INFO] [GPU = " << dprop.name << "] : " << energyStart << "eV -> " << energyEnd << "eV\n";
    }


#ifdef PROFILING
    {
      START_TIMER(TIMERS::MALLOC);
    }
#endif
#ifdef DUMP_FILES
    Complex *polarizationZ = new Complex[numVoxels];
    Complex *polarizationX = new Complex[numVoxels];
    Complex *polarizationY = new Complex[numVoxels];
#endif
#if defined(EOC) or defined(DUMP_FILES)
    Real *scatter3D = new Real[numVoxels];
#endif

#ifdef EOC
    Real *projectionCPU = new Real[BATCH * voxel[0] * voxel[1]];
#else

#endif

    Voxel<NUM_MATERIAL> *d_voxelInput;
    mallocGPU(d_voxelInput, numVoxels);

    Complex *d_polarizationZ, *d_polarizationX, *d_polarizationY;
    Real *d_scatter3D;
    UINT *d_mask;
    mallocGPU(d_polarizationX, numVoxels);
    mallocGPU(d_polarizationY, numVoxels);
    mallocGPU(d_polarizationZ, numVoxels);

    if (idata.scatterApproach == ScatterApproach::FULL) {
      mallocGPU(d_scatter3D, numVoxels);
    }
#ifndef EOC
    Real *d_projection, *d_rotProjection, *d_projectionAverage;
    mallocGPU(d_projection, numVoxel2D);
    mallocGPU(d_rotProjection, numVoxel2D);
    if (idata.rotMask) {
      mallocGPU(d_mask, numVoxel2D);
    }
    mallocGPU(d_projectionAverage, numVoxel2D);
#endif


#ifdef PROFILING
    {
      END_TIMER(TIMERS::MALLOC)
      START_TIMER(TIMERS::MEMCOPY_CPU_GPU)
    }
#endif

    hostDeviceExchange(d_voxelInput, voxelInput, numVoxels, cudaMemcpyHostToDevice);

#ifdef PROFILING
    {
      END_TIMER(TIMERS::MEMCOPY_CPU_GPU)
    }
#endif
    const Real3 & kVec = idata.k;
    Matrix rotationMatrixK,rotationMatrix;
    computeRotationMatrixK(kVec,rotationMatrixK);
    Real baseRotAngle;
    computeRotationMatrixBaseConfiguration(kVec,rotationMatrixK,rotationMatrix,baseRotAngle);

    UINT BlockSize  = static_cast<UINT>(ceil(numVoxels * 1.0 / NUM_THREADS));
    UINT BlockSize2 = static_cast<UINT>(ceil(numVoxel2D * 1.0 / NUM_THREADS));

    for (UINT j = numStart; j < numEnd; j++) {

      const Real &energy = (idata.energies[j]);
      std::cout << " [STAT] Energy = " << energy << " starting " << "\n";

      cudaZeroEntries(d_projectionAverage, numVoxel2D);
      if (idata.rotMask) {
        cudaZeroEntries(d_mask, numVoxel2D);
      }

#ifdef  PROFILING
      START_TIMER(TIMERS::ENERGY)
#endif


      ElectricField eleField;
      eleField.e.x = 1;
      eleField.e.y = 0;
      eleField.e.z = 0;
      Real wavelength = static_cast<Real>(1239.84197 / energy);
      eleField.k.x = 0;
      eleField.k.y = 0;
      eleField.k.z = static_cast<Real>(2 * M_PI / wavelength);;
      Real Eangle;
      Matrix ERotationMatrix;
      for (UINT i = 0; i < numAnglesRotation; i++) {
        Eangle = static_cast<Real>((baseRotAngle + i * idata.incrementAngle) * M_PI / 180.0);
        computeRotationMatrix(kVec,rotationMatrixK,ERotationMatrix,Eangle);
#ifdef PROFILING
        {
          START_TIMER(TIMERS::POLARIZATION)
        }
#endif
        computePolarization(materialInput[j], d_voxelInput, eleField, Eangle, vx, d_polarizationX, d_polarizationY,
                            d_polarizationZ, static_cast<FFT::FFTWindowing >(idata.windowingType),
                            idata.if2DComputation(), static_cast<MorphologyType>(idata.morphologyType), BlockSize,
                            ReferenceFrame::MATERIAL,ERotationMatrix);

#ifdef DUMP_FILES

        CUDA_CHECK_RETURN(cudaMemcpy(polarizationX,
                                     d_polarizationX,
                                     sizeof(Complex) * numVoxels,
                                     cudaMemcpyDeviceToHost));
        gpuErrchk(cudaPeekAtLastError());
        CUDA_CHECK_RETURN(cudaMemcpy(polarizationZ,
                                     d_polarizationZ,
                                     sizeof(Complex) * numVoxels,
                                     cudaMemcpyDeviceToHost));
        gpuErrchk(cudaPeekAtLastError());
        CUDA_CHECK_RETURN(cudaMemcpy(polarizationY,
                                     d_polarizationY,
                                     sizeof(Complex) * numVoxels,
                                     cudaMemcpyDeviceToHost));
        gpuErrchk(cudaPeekAtLastError());
        {
          FILE* pX = fopen("polarizeX.dmp", "wb");
          fwrite(polarizationX, sizeof(Complex), numVoxels, pX);
          fclose(pX);
          FILE* pY = fopen("polarizeY.dmp", "wb");
          fwrite(polarizationY, sizeof(Complex), numVoxels, pY);
          fclose(pY);
          FILE* pZ = fopen("polarizeZ.dmp", "wb");
          fwrite(polarizationZ, sizeof(Complex), numVoxels, pZ);
          fclose(pZ);
          std::string dirname = "Polarize/";
          std::string fname = dirname + "polarizationX" + std::to_string(i);
          VTI::writeDataScalar(polarizationX, voxel, fname.c_str(), "polarizeX");
          fname = dirname + "polarizationY" + std::to_string(i);
          VTI::writeDataScalar(polarizationY, voxel, fname.c_str(), "polarizeY");
          fname = dirname + "polarizationZ" + std::to_string(i);
          VTI::writeDataScalar(polarizationZ, voxel, fname.c_str(), "polarizeZ");
        }
#endif

#ifdef PROFILING
        {
          END_TIMER(TIMERS::POLARIZATION)
          START_TIMER(TIMERS::FFT)
        }
#endif
        /** FFT Computation **/
        result = performFFT(d_polarizationX, plan);
        if (result != CUFFT_SUCCESS) {
          std::cout << "CUFFT failed with result " << result << "\n";
#pragma omp cancel parallel
          exit(EXIT_FAILURE);
        }

        result = performFFT(d_polarizationY, plan);
        if (result != CUFFT_SUCCESS) {
          std::cout << "CUFFT failed with result " << result << "\n";
#pragma omp cancel parallel
          exit(EXIT_FAILURE);
        }

        result = performFFT(d_polarizationZ, plan);
        if (result != CUFFT_SUCCESS) {
          std::cout << "CUFFT failed with result " << result << "\n";
          exit(EXIT_FAILURE);
        }

#ifdef DUMP_FILES
        CUDA_CHECK_RETURN(cudaMemcpy(polarizationX,
                                     d_polarizationX,
                                     sizeof(Complex) * numVoxels,
                                     cudaMemcpyDeviceToHost));
        gpuErrchk(cudaPeekAtLastError());
        CUDA_CHECK_RETURN(cudaMemcpy(polarizationY,
                                     d_polarizationY,
                                     sizeof(Complex) * numVoxels,
                                     cudaMemcpyDeviceToHost));
        gpuErrchk(cudaPeekAtLastError());
        CUDA_CHECK_RETURN(cudaMemcpy(polarizationZ,
                                     d_polarizationZ,
                                     sizeof(Complex) * numVoxels,
                                     cudaMemcpyDeviceToHost));
        gpuErrchk(cudaPeekAtLastError());
        {
          FILE* pX = fopen("fftpolarizeXbshift.dmp", "wb");
          fwrite(polarizationX, sizeof(Complex), numVoxels, pX);
          fclose(pX);
          FILE* pY = fopen("fftpolarizeYbshift.dmp", "wb");
          fwrite(polarizationY, sizeof(Complex), numVoxels, pY);
          fclose(pY);
          FILE* pZ = fopen("fftpolarizeZbshift.dmp", "wb");
          fwrite(polarizationZ, sizeof(Complex), numVoxels, pZ);
          fclose(pZ);
          std::string dirname = "FFT/";
          std::string fname = dirname + "polarizationXfftbshift" + std::to_string(i);
          VTI::writeDataScalar(polarizationX, voxel, fname.c_str(), "polarizeXfft");
          fname = dirname + "polarizationYfftbshift" + std::to_string(i);
          VTI::writeDataScalar(polarizationY, voxel, fname.c_str(), "polarizeYfft");
          fname = dirname + "polarizationZfftbshift" + std::to_string(i);
          VTI::writeDataScalar(polarizationZ, voxel, fname.c_str(), "polarizeZfft");
        }
#endif
#pragma message "Make me callbacks"
        performFFTShift(d_polarizationX, BlockSize, vx);
        performFFTShift(d_polarizationY, BlockSize, vx);
        performFFTShift(d_polarizationZ, BlockSize, vx);
#ifdef DUMP_FILES
        CUDA_CHECK_RETURN(cudaMemcpy(polarizationX,
                                     d_polarizationX,
                                     sizeof(Complex) * numVoxels,
                                     cudaMemcpyDeviceToHost));
        gpuErrchk(cudaPeekAtLastError());
        CUDA_CHECK_RETURN(cudaMemcpy(polarizationY,
                                     d_polarizationY,
                                     sizeof(Complex) * numVoxels,
                                     cudaMemcpyDeviceToHost));
        gpuErrchk(cudaPeekAtLastError());
        CUDA_CHECK_RETURN(cudaMemcpy(polarizationZ,
                                     d_polarizationZ,
                                     sizeof(Complex) * numVoxels,
                                     cudaMemcpyDeviceToHost));
        gpuErrchk(cudaPeekAtLastError());
        {
          FILE* pX = fopen("fftpolarizeX.dmp", "wb");
          fwrite(polarizationX, sizeof(Complex), numVoxels, pX);
          fclose(pX);
          FILE* pY = fopen("fftpolarizeY.dmp", "wb");
          fwrite(polarizationY, sizeof(Complex), numVoxels, pY);
          fclose(pY);
          FILE* pZ = fopen("fftpolarizeZ.dmp", "wb");
          fwrite(polarizationZ, sizeof(Complex), numVoxels, pZ);
          fclose(pZ);
          std::string dirname = "FFT/";
          std::string fname = dirname + "polarizationXfft" + std::to_string(i);
          VTI::writeDataScalar(polarizationX, voxel, fname.c_str(), "polarizeXfft");
          fname = dirname + "polarizationYfft" + std::to_string(i);
          VTI::writeDataScalar(polarizationY, voxel, fname.c_str(), "polarizeYfft");
          fname = dirname + "polarizationZfft" + std::to_string(i);
          VTI::writeDataScalar(polarizationZ, voxel, fname.c_str(), "polarizeZfft");
        }
#endif

#ifdef PROFILING
        {
            END_TIMER(TIMERS::FFT)
            START_TIMER(TIMERS::SCATTER3D)
        }
#endif
        cudaZeroEntries(d_rotProjection, numVoxel2D);
        cudaZeroEntries(d_projection, numVoxel2D);

        if (idata.scatterApproach == ScatterApproach::FULL) {

          performScatter3DComputation(d_polarizationX, d_polarizationY, d_polarizationZ, d_scatter3D, eleField, 0.0,
                                      0.0, numVoxels, vx, idata.physSize, idata.if2DComputation(), BlockSize,kVec);

#ifdef DUMP_FILES
          CUDA_CHECK_RETURN(cudaMemcpy(scatter3D, d_scatter3D, sizeof(Real) * numVoxels, cudaMemcpyDeviceToHost));
          gpuErrchk(cudaPeekAtLastError())
          {
            FILE* scatter = fopen("scatter_3D.dmp", "wb");
            fwrite(scatter3D, sizeof(Real), numVoxels, scatter);
            fclose(scatter);
            std::string dirname = "Scatter/";
            std::string fname = dirname + "scatter" + std::to_string(i);
            VTI::writeDataScalar(scatter3D, voxel, fname.c_str(), "scatter3D");
          }

#endif


#ifdef EOC
          CUDA_CHECK_RETURN(cudaMemcpy(scatter3D, d_scatter3D, sizeof(Real) * numVoxels, cudaMemcpyDeviceToHost));
          gpuErrchk(cudaPeekAtLastError());

#ifdef PROFILING
          {

          }
#endif
          computeEwaldProjectionCPU(projectionCPU, scatter3D, vx, eleField.k.x);
#else
          peformEwaldProjectionGPU(d_projection, d_scatter3D, eleField.k.z, vx, 0.0, 0.0, idata.physSize,
                                   static_cast<Interpolation::EwaldsInterpolation>(idata.ewaldsInterpolation),
                                   idata.if2DComputation(), BlockSize2,kVec);
#ifdef DUMP_FILES
          hostDeviceExcange(projectionGPUAveraged,d_projection,voxel[0]*voxel[1],cudaMemcpyDeviceToHost);
          std::string dirname = "Ewald/";
          std::string fname = dirname + "ewlad" + std::to_string(i);
          VTI::writeDataScalar2DFP(projectionGPUAveraged, voxel, fname.c_str(), "ewald");
          FILE* projection = fopen("projection_scatterFull.dmp", "wb");
          fwrite(projectionGPUAveraged, sizeof(Real), numVoxels, projection);
          fclose(projection);
#endif
        } else {
          peformEwaldProjectionGPU(d_projection, d_polarizationX, d_polarizationY, d_polarizationZ, eleField.k.z, vx,
                                   0.0, 0.0, idata.physSize,
                                   static_cast<Interpolation::EwaldsInterpolation>(idata.ewaldsInterpolation),
                                   idata.if2DComputation(), BlockSize2,kVec);
#ifdef DUMP_FILES

          hostDeviceExcange(projectionGPUAveraged,d_projection,voxel[0]*voxel[1],cudaMemcpyDeviceToHost);
          std::string dirname = "Ewald/";
          std::string fname = dirname + "ewlad" + std::to_string(i);
          VTI::writeDataScalar2DFP(projectionGPUAveraged, voxel, fname.c_str(), "ewald");
          FILE* projection = fopen("projection_scatterPartial.dmp", "wb");
          fwrite(projectionGPUAveraged, sizeof(Real), numVoxels, projection);
          fclose(projection);
#endif
        }


        Real _factor;
        _factor = NAN;

        stat = cublasScale(handle, voxel[0] * voxel[1], &_factor, d_rotProjection, 1);
        _factor = 1.0;
        stat = cublasScale(handle, voxel[0] * voxel[1], &_factor, d_projection, 1);

        if (stat != CUBLAS_STATUS_SUCCESS) {
          std::cout << "CUBLAS during scaling failed  with status " << stat << "\n";
          exit(EXIT_FAILURE);
        }


#ifdef PROFILING
        {
          END_TIMER(TIMERS::SCATTER3D)
          START_TIMER(TIMERS::IMAGE_ROTATION)
        }
#endif
        const double alpha = cos(Eangle);
        const double beta = sin(Eangle);

        /**https://docs.opencv.org/2.4/modules/imgproc/doc/geometric_transformations.html?highlight=warpaffine**/
        const double coeffs[2][3]{
          alpha, beta, static_cast<Real>(((1 - alpha) * voxel[0] / 2 - beta * voxel[1] / 2.)),
          -beta, alpha, static_cast<Real>(beta * voxel[0] / 2. + (1 - alpha) * voxel[1] / 2.)
        };


        NppStatus status = warpAffine(d_projection,
                                      sizeImage,
                                      voxel[1] * sizeof(Real),
                                      rect,
                                      d_rotProjection,
                                      voxel[1] * sizeof(Real),
                                      rect,
                                      coeffs,
                                      NPPI_INTER_LINEAR);

        if (status < 0) {
          std::cout << "Image rotation failed with error = " << status << "\n";
          exit(-1);
        }
        if (status != NPP_SUCCESS) {
          std::cout << YLW << "[WARNING] Image rotation warning = " << status << NRM << "\n";
        }

        if (idata.rotMask) {
          computeRotationMask<<< BlockSize2, NUM_THREADS >>>(d_rotProjection, d_mask, vx);
          cudaDeviceSynchronize();
        }

        const Real factor = static_cast<Real>(1.0);
        stat = cublasAXPY(handle, voxel[0] * voxel[1], &factor, d_rotProjection, 1, d_projectionAverage, 1);
        if (stat != CUBLAS_STATUS_SUCCESS) {
          std::cout << "CUBLAS during sum failed  with status " << stat << "\n";
          exit(EXIT_FAILURE);
        }

#ifdef PROFILING
        {
          END_TIMER(TIMERS::IMAGE_ROTATION)
        }
#endif
#endif
      }

      if (idata.rotMask) {
        averageRotation<<<BlockSize2, NUM_THREADS>>>(d_projectionAverage, d_mask, vx);
        cudaDeviceSynchronize();
        gpuErrchk(cudaPeekAtLastError());
      } else {
        /// The averaging out for all angles
        const Real alphaFac = static_cast<Real>(1.0 / numAnglesRotation);
        stat = cublasScale(handle, voxel[0] * voxel[1], &alphaFac, d_projectionAverage, 1);
        if (stat != CUBLAS_STATUS_SUCCESS) {
          std::cout << "CUBLAS during averaging failed  with status " << stat << "\n";
          exit(EXIT_FAILURE);
        }
      }
      //// Rotate Image
      hostDeviceExchange(d_projection,d_rotProjection,numVoxel2D,cudaMemcpyDeviceToDevice);
      const double srcPoints[3][2]{{voxel[0]/2.,voxel[1]/2.},{voxel[0]*0.5,voxel[1]*1.0},{voxel[0]*1.0,voxel[1]*0.5}};
      const double destPoints[3][2]{{voxel[0]/2.,voxel[1]/2.},{voxel[0]*0.5,voxel[1]*1.0},{voxel[0]*1.0,voxel[1]*0.5}};
      double coeffs[2][3];
      computeWarpAffineMatrix(srcPoints,destPoints,coeffs);
      Real _factor;
      _factor = NAN;

      stat = cublasScale(handle, numVoxel2D, &_factor, d_rotProjection, 1);
      NppStatus status = warpAffine(d_projection,
                                    sizeImage,
                                    voxel[1] * sizeof(Real),
                                    rect,
                                    d_rotProjection,
                                    voxel[1] * sizeof(Real),
                                    rect,
                                    coeffs,
                                    NPPI_INTER_LINEAR);

      if (status < 0) {
        std::cout << "Image rotation failed with error = " << status << "\n";
        exit(-1);
      }
      if (status != NPP_SUCCESS) {
        std::cout << YLW << "[WARNING] Image rotation warning = " << status << NRM << "\n";
      }
#ifdef PROFILING
      {
        START_TIMER(TIMERS::MEMCOPY_GPU_CPU)
      }
#endif

      hostDeviceExchange(&projectionGPUAveraged[j * numVoxel2D], d_projectionAverage, numVoxel2D,
                         cudaMemcpyDeviceToHost);
#ifdef PROFILING
      {
        END_TIMER(TIMERS::MEMCOPY_GPU_CPU)
        END_TIMER(TIMERS::ENERGY)
      }
#endif
    }

    /** Freeing bunch of memories not required now **/
    freeCudaMemory(d_polarizationX);
    freeCudaMemory(d_polarizationY);
    freeCudaMemory(d_polarizationZ);
    if (idata.scatterApproach == ScatterApproach::FULL) {
      freeCudaMemory(d_scatter3D);
    }
    freeCudaMemory(d_voxelInput);

#ifndef EOC
    freeCudaMemory(d_projection);
    freeCudaMemory(d_projectionAverage);
    freeCudaMemory(d_rotProjection);
    if (idata.rotMask) {
      freeCudaMemory(d_mask);
    }
#endif
#ifdef DUMP_FILES
    delete[] polarizationX;
    delete[] polarizationY;
    delete[] polarizationZ;

#endif
#if (defined(DUMP_FILES) or defined(EOC))
    delete[] scatter3D;
#endif


    cufftDestroy(plan);
    cublasDestroy(handle);

#ifdef EOC
    delete[] projectionCPU;
#endif
  }


#ifdef PROFILING
  std::cout << "\n\n[INFO] Timings Info\n";
  for(int i = 0; i < TIMERS::MAX; i++){
    std::cout << "[TIMERS] " << std::left << std::setw(20) << timersName[i] << ":" << timings[i] << " s\n";
  }
  std::cout << "\n\n";
#endif


  return (EXIT_SUCCESS);
}

