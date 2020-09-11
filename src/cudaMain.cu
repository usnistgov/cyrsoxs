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
#include <cufft.h>
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


int warmup(){
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





__global__ void computePolarization(Material<NUM_MATERIAL> materialInput,
                                    const Voxel<NUM_MATERIAL> *voxelInput,
                                    const ElectricField elefield,
                                    const Real angle,
                                    const uint3 voxel,
                                    Complex *polarizationX,
                                    Complex *polarizationY,
                                    Complex *polarizationZ,
                                    FFT::FFTWindowing windowing,
                                    const bool enable2D
) {
  UINT threadID = threadIdx.x + blockIdx.x * blockDim.x;
  const UINT voxelNum = voxel.x*voxel.y*voxel.z;

  if (threadID >= voxelNum) {
    return;
  }

#ifndef BIAXIAL
    computePolarizationUniaxial(&materialInput,angle,voxelInput,threadID,polarizationX,polarizationY,polarizationZ);
#else
    printf("Kernel not spported\n");
#endif


if(windowing == FFT::FFTWindowing::HANNING) {
  UINT Z = static_cast<UINT>(floorf(threadID / (voxel.y * voxel.x * 1.0)));
  UINT Y = static_cast<UINT>(floorf((threadID - Z * voxel.y * voxel.x) / (voxel.x * 1.0)));
  UINT X = static_cast<UINT>(threadID - Y * voxel.x - Z * voxel.y * voxel.x);
  Real3 hanningWeight;
  hanningWeight.x = static_cast<Real> (0.5 * (1 - cos(2 * M_PI * X / (voxel.x))));
  hanningWeight.y = static_cast<Real> (0.5 * (1 - cos(2 * M_PI * Y / (voxel.y))));
  hanningWeight.z = static_cast<Real>(1.0);
  if(not(enable2D)){
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


int cudaMain(const UINT *voxel,
             const InputData &idata,
             const std::vector<Material<NUM_MATERIAL> > &materialInput,
             Real * projectionGPUAveraged,
             const Voxel<NUM_MATERIAL> *voxelInput) {

  const BigUINT voxelSize = voxel[0] * voxel[1] * voxel[2]; /// Voxel size
  const UINT
      numAnglesRotation = static_cast<UINT>(std::round((idata.endAngle - idata.startAngle) / idata.incrementAngle + 1));
  const UINT
      numEnergyLevel = static_cast<UINT>(std::round((idata.energyEnd - idata.energyStart) / idata.incrementEnergy + 1));


  int num_gpu;
  cudaGetDeviceCount(&num_gpu);
  std::cout << "Number of CUDA devices:" << num_gpu << "\n";

  if(num_gpu < 1){
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
    EWALDS = 5,
    ROTATION=6,
    MEMCOPY_GPU_CPU = 7,
    ENERGY=8,
    MAX = 9
  };
  static const char *timersName[]{"Malloc on CPU + GPU",
                                  "Memcopy CPU -> GPU",
                                  "Polarization",
                                  "FFT",
                                  "Scatter3D",
                                  "Ewalds",
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


    int device_gpu = -1;
    cudaSetDevice(omp_get_thread_num());
    cudaDeviceProp dprop;
    cudaGetDeviceProperties(&dprop, omp_get_thread_num());
    cudaGetDeviceCount(&device_gpu);
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


    UINT numEnergyPerGPU = static_cast<UINT>(std::ceil(numEnergyLevel*1.0/num_gpu));
    UINT numStart = (numEnergyPerGPU*omp_get_thread_num());
    UINT numEnd   = (numEnergyPerGPU*(omp_get_thread_num() + 1));
    numEnd = std::min(numEnd,numEnergyLevel);

    Real energyStart = idata.energyStart + numStart*idata.incrementEnergy;
    Real energyEnd   = idata.energyStart + (numEnd - 1)*idata.incrementEnergy;

    if(energyStart > idata.energyEnd){
      std::cout << "[INFO] [GPU = " << dprop.name  << "] -> No computation. Idle\n";
    }
    else{
      std::cout << "[INFO] [GPU = " << dprop.name  << "] : " << energyStart << "eV -> " << energyEnd << "eV\n" ;
    }


#ifdef PROFILING
    {
      START_TIMER(TIMERS::MALLOC);
    }
#endif
    uint3 vx;
    vx.x = voxel[0];
    vx.y = voxel[1];
    vx.z = voxel[2];
    Complex *polarizationZ = new Complex[voxelSize];
    Complex *polarizationX = new Complex[voxelSize];
    Complex *polarizationY = new Complex[voxelSize];
    Real *scatter3D = new Real[voxelSize];

#pragma message "Remove me. I am not needed if doing inplace swap"
    Complex * d_temporarySwapVariable;
#ifdef EOC
    Real *projectionCPU = new Real[BATCH * voxel[0] * voxel[1]];
#else

#endif

    Complex *d_polarizationZ, *d_polarizationX, *d_polarizationY;
    Real *d_scatter3D;
    UINT * d_mask;
    CUDA_CHECK_RETURN(cudaMalloc((void **) &d_polarizationZ, sizeof(Complex) * voxelSize));
    gpuErrchk(cudaPeekAtLastError());
    CUDA_CHECK_RETURN(cudaMalloc((void **) &d_polarizationX, sizeof(Complex) * voxelSize));
    gpuErrchk(cudaPeekAtLastError());
    CUDA_CHECK_RETURN(cudaMalloc((void **) &d_polarizationY, sizeof(Complex) * voxelSize));
    gpuErrchk(cudaPeekAtLastError());
    CUDA_CHECK_RETURN(cudaMalloc((void **) &d_scatter3D, sizeof(Real) * voxelSize));
    gpuErrchk(cudaPeekAtLastError());
    CUDA_CHECK_RETURN(cudaMalloc((void **) &d_temporarySwapVariable, sizeof(Complex) * voxelSize));
    gpuErrchk(cudaPeekAtLastError());



#ifndef EOC
    Real *d_projection, *d_rotProjection, *d_projectionAverage;
    CUDA_CHECK_RETURN(cudaMalloc((void **) &d_projection, sizeof(Real) * (voxel[0] * voxel[1])));
    gpuErrchk(cudaPeekAtLastError());
    CUDA_CHECK_RETURN(cudaMalloc((void **) &d_rotProjection, sizeof(Real) * (voxel[0] * voxel[1])));
    gpuErrchk(cudaPeekAtLastError());
    if(idata.rotMask){
      CUDA_CHECK_RETURN(cudaMalloc((void **) &d_mask, sizeof(UINT) * voxel[0]*voxel[1]));
      gpuErrchk(cudaPeekAtLastError());
    }

    CUDA_CHECK_RETURN(cudaMalloc((void **) &d_projectionAverage, sizeof(Real) * (voxel[0] * voxel[1])));
    gpuErrchk(cudaPeekAtLastError());

#endif

    Voxel<NUM_MATERIAL> *d_voxelInput;
    CUDA_CHECK_RETURN(cudaMalloc((void **) &d_voxelInput, sizeof(Voxel<NUM_MATERIAL>) * voxelSize));
    gpuErrchk(cudaPeekAtLastError());

#ifdef PROFILING
    {
      END_TIMER(TIMERS::MALLOC)
      START_TIMER(TIMERS::MEMCOPY_CPU_GPU)
    }
#endif

    CUDA_CHECK_RETURN(cudaMemcpy(d_voxelInput,
                                 voxelInput,
                                 sizeof(Voxel<NUM_MATERIAL>) * voxelSize,
                                 cudaMemcpyHostToDevice));
    gpuErrchk(cudaPeekAtLastError());

#ifdef PROFILING
    {
      END_TIMER(TIMERS::MEMCOPY_CPU_GPU)
    }
#endif

    UINT BlockSize = static_cast<UINT >(ceil(voxelSize * 1.0 / NUM_THREADS));
    UINT BlockSize2 = static_cast<UINT>(ceil(voxel[0] * voxel[1] * 1.0 / NUM_THREADS));

    for (UINT j = numStart; j < numEnd; j++) {

      Real energy = (idata.energyStart + j * idata.incrementEnergy);
      std::cout << " [STAT] Energy = " << energy << " starting "  << "\n";

      CUDA_CHECK_RETURN(cudaMemset(d_projectionAverage, 0, voxel[0] * voxel[1] * sizeof(Real)));
      gpuErrchk(cudaPeekAtLastError());
      if(idata.rotMask) {
        cudaMemset(d_mask, 0, sizeof(UINT) * voxel[0] * voxel[1]);
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
      Real angle;
      for (int i = 0; i < numAnglesRotation; i++) {
        angle = static_cast<Real>((idata.startAngle + i*idata.incrementAngle) * M_PI / 180.0);
#ifdef PROFILING
        {
          START_TIMER(TIMERS::POLARIZATION)
        }
#endif

        computePolarization <<< BlockSize, NUM_THREADS >>> (materialInput[j], d_voxelInput, eleField, angle, vx, d_polarizationX, d_polarizationY, d_polarizationZ,
                static_cast<FFT::FFTWindowing >(idata.windowingType), idata.if2DComputation());

        gpuErrchk(cudaPeekAtLastError());
        cudaDeviceSynchronize();
#ifdef DUMP_FILES

        CUDA_CHECK_RETURN(cudaMemcpy(polarizationX,
                                     d_polarizationX,
                                     sizeof(Complex) * voxelSize,
                                     cudaMemcpyDeviceToHost));
        gpuErrchk(cudaPeekAtLastError());
        CUDA_CHECK_RETURN(cudaMemcpy(polarizationZ,
                                     d_polarizationZ,
                                     sizeof(Complex) * voxelSize,
                                     cudaMemcpyDeviceToHost));
        gpuErrchk(cudaPeekAtLastError());
        CUDA_CHECK_RETURN(cudaMemcpy(polarizationY,
                                     d_polarizationY,
                                     sizeof(Complex) * voxelSize,
                                     cudaMemcpyDeviceToHost));
        gpuErrchk(cudaPeekAtLastError());
        {
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
        result = performFFT(d_polarizationX,plan);
        if (result != CUFFT_SUCCESS) {
          std::cout << "CUFFT failed with result " << result << "\n";
          #pragma omp cancel parallel
          exit(EXIT_FAILURE);
        }

        result = performFFT(d_polarizationY,plan);
        if (result != CUFFT_SUCCESS) {
          std::cout << "CUFFT failed with result " << result << "\n";
          #pragma omp cancel parallel
          exit(EXIT_FAILURE);
        }

        result = performFFT(d_polarizationZ,plan);
        if (result != CUFFT_SUCCESS) {
          std::cout << "CUFFT failed with result " << result << "\n";
          exit(EXIT_FAILURE);
        }

#pragma message "Fix me! Remove me for inplace transofrmation. Not needed and make me as a callback"

        CUDA_CHECK_RETURN(cudaMemcpy(d_temporarySwapVariable,
                                     d_polarizationX,
                                     sizeof(Complex) * voxelSize,
                                     cudaMemcpyDeviceToDevice));
        gpuErrchk(cudaPeekAtLastError());

        FFTIgor<<<BlockSize,NUM_THREADS>>>(d_polarizationX,d_temporarySwapVariable,vx);
        cudaDeviceSynchronize();
        gpuErrchk(cudaPeekAtLastError());

        CUDA_CHECK_RETURN(cudaMemcpy(d_temporarySwapVariable,
                                     d_polarizationY,
                                     sizeof(Complex) * voxelSize,
                                     cudaMemcpyDeviceToDevice));
        gpuErrchk(cudaPeekAtLastError());

        FFTIgor<<<BlockSize,NUM_THREADS>>>(d_polarizationY,d_temporarySwapVariable,vx);
        cudaDeviceSynchronize();
        gpuErrchk(cudaPeekAtLastError());

        CUDA_CHECK_RETURN(cudaMemcpy(d_temporarySwapVariable,
                                     d_polarizationZ,
                                     sizeof(Complex) * voxelSize,
                                     cudaMemcpyDeviceToDevice));
        gpuErrchk(cudaPeekAtLastError());

        FFTIgor<<<BlockSize,NUM_THREADS>>>(d_polarizationZ,d_temporarySwapVariable,vx);
        cudaDeviceSynchronize();
        gpuErrchk(cudaPeekAtLastError());

#ifdef DUMP_FILES
        CUDA_CHECK_RETURN(cudaMemcpy(polarizationX,
                                     d_polarizationX,
                                     sizeof(Complex) * voxelSize,
                                     cudaMemcpyDeviceToHost));
        gpuErrchk(cudaPeekAtLastError());
        CUDA_CHECK_RETURN(cudaMemcpy(polarizationY,
                                     d_polarizationY,
                                     sizeof(Complex) * voxelSize,
                                     cudaMemcpyDeviceToHost));
        gpuErrchk(cudaPeekAtLastError());
        CUDA_CHECK_RETURN(cudaMemcpy(polarizationZ,
                                     d_polarizationZ,
                                     sizeof(Complex) * voxelSize,
                                     cudaMemcpyDeviceToHost));
        gpuErrchk(cudaPeekAtLastError());
        {
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



        /** Scatter 3D computation **/
        computeScatter3D <<< BlockSize, NUM_THREADS >>> (d_polarizationX, d_polarizationY, d_polarizationZ, d_scatter3D, eleField, voxelSize, vx, idata.physSize,
                                                         idata.if2DComputation());
        cudaDeviceSynchronize();
        gpuErrchk(cudaPeekAtLastError());

#ifdef DUMP_FILES
        CUDA_CHECK_RETURN(cudaMemcpy(scatter3D, d_scatter3D, sizeof(Real) * voxelSize, cudaMemcpyDeviceToHost));
        gpuErrchk(cudaPeekAtLastError())
        {
          std::string dirname = "Scatter/";
          std::string fname = dirname + "scatter" + std::to_string(i);
          VTI::writeDataScalar(scatter3D, voxel, fname.c_str(), "scatter3D");
        }

#endif

#ifdef PROFILING
        {
          END_TIMER(TIMERS::SCATTER3D)
          START_TIMER(TIMERS::EWALDS)
        }
#endif

#ifdef EOC
        CUDA_CHECK_RETURN(cudaMemcpy(scatter3D, d_scatter3D, sizeof(Real) * voxelSize, cudaMemcpyDeviceToHost));
        gpuErrchk(cudaPeekAtLastError());

#ifdef PROFILING
        {

        }
#endif
        computeEwaldProjectionCPU(projectionCPU, scatter3D, vx, eleField.k.x);
#else
        cudaMemset(d_rotProjection, 0, voxel[0] * voxel[1] * sizeof(Real));

        computeEwaldProjectionGPU <<< BlockSize2, NUM_THREADS >>> (d_projection,d_rotProjection, d_scatter3D, vx,
            eleField.k.z,idata.physSize, static_cast<Interpolation::EwaldsInterpolation>(idata.ewaldsInterpolation),
                                                                   idata.if2DComputation());
        cudaDeviceSynchronize();
        gpuErrchk(cudaPeekAtLastError());
#ifdef PROFILING
        {
          END_TIMER(TIMERS::EWALDS)
          START_TIMER(TIMERS::ROTATION)
        }
#endif
        const double alpha = cos(angle);
        const double beta = sin(angle);

        /**https://docs.opencv.org/2.4/modules/imgproc/doc/geometric_transformations.html?highlight=warpaffine**/
        const double coeffs[2][3]{
            alpha, beta, static_cast<Real>(((1 - alpha) * voxel[0] / 2 - beta * voxel[1] / 2.)),
            -beta, alpha, static_cast<Real>(beta * voxel[0] / 2. + (1 - alpha) * voxel[1] / 2.)
        };

#ifdef DOUBLE_PRECISION
        NppStatus status = nppiWarpAffine_64f_C1R(d_projection,
                                                  sizeImage,
                                                  voxel[1] * sizeof(Real),
                                                  rect,
                                                  d_rotProjection,
                                                  voxel[1] * sizeof(Real),
                                                  rect,
                                                  coeffs,
                                                  NPPI_INTER_LINEAR);

#else
        NppStatus status = nppiWarpAffine_32f_C1R(d_projection,
                                                  sizeImage,
                                                  voxel[1] * sizeof(Real),
                                                  rect,
                                                  d_rotProjection,
                                                  voxel[1] * sizeof(Real),
                                                  rect,
                                                  coeffs,
                                                  NPPI_INTER_LINEAR);
#endif

        if (status < 0) {
          std::cout << "Image rotation failed with error = " << status << "\n";
          exit(-1);
        }
        if(status != NPP_SUCCESS){
          std::cout << "[WARNING] Image rotation warning = " << status << "\n";
        }

        if(idata.rotMask){
          computeRotationMask<<< BlockSize2, NUM_THREADS >>>(d_rotProjection,d_mask,vx);
          cudaDeviceSynchronize();
        }

        const Real factor = static_cast<Real>(1.0);
#ifdef DOUBLE_PRECISION
        stat = cublasDaxpy(handle, voxel[0] * voxel[1], &factor, d_rotProjection, 1, d_projectionAverage, 1);
#else
        stat = cublasSaxpy(handle, voxel[0] * voxel[1], &factor, d_rotProjection, 1, d_projectionAverage, 1);
#endif
        if (stat != CUBLAS_STATUS_SUCCESS) {
          std::cout << "CUBLAS during sum failed  with status " << stat << "\n";
          exit(EXIT_FAILURE);
        }

#ifdef PROFILING
        {
          END_TIMER(TIMERS::ROTATION)
        }
#endif
#endif
      }

      if(idata.rotMask){
        averageRotation<<<BlockSize2,NUM_THREADS>>>(d_projectionAverage,d_mask,vx);
        cudaDeviceSynchronize();
        gpuErrchk(cudaPeekAtLastError());
      }
      else {
        /// The averaging out for all angles
        const Real alphaFac = static_cast<Real>(1.0 / numAnglesRotation);
#ifdef DOUBLE_PRECISION
        stat = cublasDscal(handle, voxel[0] * voxel[1], &alphaFac, d_projectionAverage, 1);
#else
        stat = cublasSscal(handle, voxel[0] * voxel[1], &alphaFac, d_projectionAverage, 1);
#endif
        if (stat != CUBLAS_STATUS_SUCCESS) {
          std::cout << "CUBLAS during averaging failed  with status " << stat << "\n";
          exit(EXIT_FAILURE);
        }
      }
#ifdef PROFILING
      {
        START_TIMER(TIMERS::MEMCOPY_GPU_CPU)
      }
#endif
      CUDA_CHECK_RETURN(cudaMemcpy(&projectionGPUAveraged[j * voxel[0] * voxel[1]],
                                   d_projectionAverage,
                                   sizeof(Real) * (voxel[0] * voxel[1]),
                                   cudaMemcpyDeviceToHost));
      gpuErrchk(cudaPeekAtLastError());

#ifdef PROFILING
      {
        END_TIMER(TIMERS::MEMCOPY_GPU_CPU)
        END_TIMER(TIMERS::ENERGY)
      }
#endif
    }

    /** Freeing bunch of memories not required now **/
    CUDA_CHECK_RETURN(cudaFree(d_polarizationZ));
    gpuErrchk(cudaPeekAtLastError());
    CUDA_CHECK_RETURN(cudaFree(d_polarizationY));
    gpuErrchk(cudaPeekAtLastError());
    CUDA_CHECK_RETURN(cudaFree(d_polarizationX));
    gpuErrchk(cudaPeekAtLastError());
    CUDA_CHECK_RETURN(cudaFree(d_voxelInput));
    gpuErrchk(cudaPeekAtLastError());
#pragma message "Fix Me! You have not allocated me. No need to delete me"
    CUDA_CHECK_RETURN(cudaFree(d_temporarySwapVariable));
    gpuErrchk(cudaPeekAtLastError());

#ifndef EOC
    CUDA_CHECK_RETURN(cudaFree(d_projection));
    gpuErrchk(cudaPeekAtLastError());
    CUDA_CHECK_RETURN(cudaFree(d_rotProjection));
    gpuErrchk(cudaPeekAtLastError());
    if(idata.rotMask) {
      CUDA_CHECK_RETURN(cudaFree(d_mask));
      gpuErrchk(cudaPeekAtLastError());
    }
#endif

    delete[] polarizationX;
    delete[] polarizationY;
    delete[] polarizationZ;
    delete[] scatter3D;


    cufftDestroy(plan);
    cublasDestroy(handle);

  }



#ifdef PROFILING
  std::cout << "\n\n[INFO] Timings Info\n";
  for(int i = 0; i < TIMERS::MAX; i++){
    std::cout << "[TIMERS] " << std::left << std::setw(20) << timersName[i] << ":" << timings[i] << " s\n";
  }
  std::cout << "\n\n";
#endif

#ifdef EOC
  delete[] projectionCPU;
#else
#endif

  return (EXIT_SUCCESS);
}

