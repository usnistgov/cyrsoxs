
#include <cudaMain.h>
#include <Input/Input.h>
#include <bits/ios_base.h>

#include <cuda_runtime.h>
#include <cufft.h>
#include <Output/writeVTI.h>
#include <uniaxial.h>
#include <cublas_v2.h>
#include <chrono>
#include <ctime>
#include <chrono>
#include <npp.h>



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


void createDirectory(const std::string dirName){

  int ierr = mkdir(dirName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if (ierr != 0 && errno != EEXIST) {
    std::cout  << "Could not create folder for storing results (" <<  strerror(errno) << "\n";
    exit(EXIT_FAILURE);
  }
}

__global__ void computePolarization(Material<NUM_MATERIAL> materialInput,
                                    const Voxel<NUM_MATERIAL> *voxelInput,
                                    const ElectricField elefield,
                                    const Real angle,
                                    const uint3 voxel,
                                    Complex *polarizationX,
                                    Complex *polarizationY,
                                    Complex *polarizationZ
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


#ifdef HANNING
  UINT Z = static_cast<UINT>(floorf(threadID / (voxel.y * voxel.x * 1.0)));
  UINT Y = static_cast<UINT>(floorf((threadID - Z * voxel.y * voxel.x) / voxel.x * 1.0));
  UINT X = static_cast<UINT>(threadID - Y * voxel.x - Z * voxel.y * voxel.x);
  Real3 hanningWeight;
  hanningWeight.x = 0.5*(1 - cos(2*M_PI*X/(voxel.x)));
  hanningWeight.y = 0.5*(1 - cos(2*M_PI*Y/(voxel.y)));
#ifdef ENABLE_2D
  hanningWeight.z = 1.0;
#else
  hanningWeight.z = 0.5*(1 - cos(2*M_PI*Z/(voxel.z)));
#endif
  Real totalHanningWeight =hanningWeight.x * hanningWeight.y * hanningWeight.z;
  polarizationX[threadID].x *= totalHanningWeight;   polarizationX[threadID].y *= totalHanningWeight;
  polarizationY[threadID].x *= totalHanningWeight;   polarizationY[threadID].y *= totalHanningWeight;
  polarizationZ[threadID].x *= totalHanningWeight;   polarizationZ[threadID].y *= totalHanningWeight;

#endif

}

int cudaMain(const UINT *voxel,
             const InputData &idata,
             const std::vector<Material<NUM_MATERIAL> > &materialInput,
             const Voxel<NUM_MATERIAL> *voxelInput) {

  const BigUINT voxelSize = voxel[0] * voxel[1] * voxel[2]; /// Voxel size
  const UINT
      numAnglesRotation = static_cast<UINT>(std::round((idata.endAngle - idata.startAngle) / idata.incrementAngle + 1));
  const UINT
      numEnergyLevel = static_cast<UINT>(std::round((idata.energyEnd - idata.energyStart) / idata.incrementEnergy + 1));


  int num_gpu;
  cudaGetDeviceCount(&num_gpu);
  printf("number of CUDA devices:\t%d\n", num_gpu);

  if(num_gpu < 1){
    std::cout << "No GPU found. Exiting" << "\n";
    return (EXIT_FAILURE);
  }




#ifdef PROFILING
  std::chrono::high_resolution_clock::time_point t0;
  std::chrono::high_resolution_clock::time_point t1;
  std::chrono::high_resolution_clock::time_point t2;
  std::chrono::high_resolution_clock::time_point t3;
  std::chrono::high_resolution_clock::time_point t4;
  std::chrono::high_resolution_clock::time_point t5;
  std::chrono::high_resolution_clock::time_point t6;
  std::chrono::high_resolution_clock::time_point t7;
  std::chrono::high_resolution_clock::time_point t8;
  std::chrono::high_resolution_clock::time_point t9;
  std::chrono::high_resolution_clock::time_point t10;
  std::chrono::high_resolution_clock::time_point tstartEnergy;
  std::chrono::high_resolution_clock::time_point tendEnergy;
  std::chrono::high_resolution_clock::time_point tDumpLevel2Start;
  std::chrono::high_resolution_clock::time_point tDumpLevel2End;

  Real time_span0 = 0;
  Real time_span1 = 0;
  Real time_span2 = 0;
  Real time_span3 = 0;
  Real time_span4 = 0;
  Real time_span5 = 0;
  Real time_span6 = 0;
  Real time_span7 = 0;
  Real time_span8 = 0;
  Real time_spanEnergy = 0;
  Real time_spanDumpLevel2 = 0;

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
#ifdef DLEVEL2
  createDirectory("Projection");
  Real *projectionGPUAveraged = new Real[numEnergyLevel * voxel[0] * voxel[1]];
#endif

  omp_set_num_threads(num_gpu);
#pragma omp parallel
  {


    int device_gpu = -1;
    cudaSetDevice(omp_get_thread_num());
    cudaDeviceProp dprop;
    cudaGetDeviceProperties(&dprop, omp_get_thread_num());
    cudaGetDeviceCount(&device_gpu);
    std::cout << "Number of GPU detected = " << device_gpu << "\n";

    if(warmup() == EXIT_SUCCESS){
      std::cout << "Warmup completed on GPU " << dprop.name << "\n";
    }
    else{
      std::cout << "Warmup failed on GPU " << dprop.name << "\n";
#pragma omp cancel parallel
      exit (EXIT_FAILURE);
    }

    cufftResult result;
    cufftHandle plan;
    cufftPlan3d(&plan, voxel[2], voxel[1], voxel[0], CUFFT_C2C);

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

    std::cout << "[GPU = " << dprop.name  << " " << "]" << " = " << " Energy = [ " << (idata.energyStart + numStart*idata.incrementEnergy)
                << " - > " << (idata.energyStart + numEnd*idata.incrementEnergy) <<"]\n";

#ifdef PROFILING
    {
      t0 = std::chrono::high_resolution_clock::now();
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

#ifdef EOC
    Real *projectionCPU = new Real[BATCH * voxel[0] * voxel[1]];
#else
#ifdef DLEVEL1
    Real *projectionGPU = new Real[numAnglesRotation * voxel[0] * voxel[1]];
#endif

#endif

    Complex *d_polarizationZ, *d_polarizationX, *d_polarizationY;
    Real *d_scatter3D;
    CUDA_CHECK_RETURN(cudaMalloc((void **) &d_polarizationZ, sizeof(Complex) * voxelSize));
    gpuErrchk(cudaPeekAtLastError());
    CUDA_CHECK_RETURN(cudaMalloc((void **) &d_polarizationX, sizeof(Complex) * voxelSize));
    gpuErrchk(cudaPeekAtLastError());
    CUDA_CHECK_RETURN(cudaMalloc((void **) &d_polarizationY, sizeof(Complex) * voxelSize));
    gpuErrchk(cudaPeekAtLastError());
    CUDA_CHECK_RETURN(cudaMalloc((void **) &d_scatter3D, sizeof(Real) * voxelSize));
    gpuErrchk(cudaPeekAtLastError());

#ifndef EOC
    Real *d_projection, *d_rotProjection, *d_projectionAverage;
    CUDA_CHECK_RETURN(cudaMalloc((void **) &d_projection, sizeof(Real) * (voxel[0] * voxel[1])));
    gpuErrchk(cudaPeekAtLastError());
    CUDA_CHECK_RETURN(cudaMalloc((void **) &d_rotProjection, sizeof(Real) * (voxel[0] * voxel[1])));
    gpuErrchk(cudaPeekAtLastError());

#ifdef DLEVEL2
    CUDA_CHECK_RETURN(cudaMalloc((void **) &d_projectionAverage, sizeof(Real) * (voxel[0] * voxel[1])));
    gpuErrchk(cudaPeekAtLastError());
#endif

#endif

    Voxel<NUM_MATERIAL> *d_voxelInput;
    CUDA_CHECK_RETURN(cudaMalloc((void **) &d_voxelInput, sizeof(Voxel<NUM_MATERIAL>) * voxelSize));
    gpuErrchk(cudaPeekAtLastError());

#ifdef PROFILING
    {
      t1 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::milli> time_span = t1 - t0;
      time_span0 += time_span.count();
    }
#endif

    CUDA_CHECK_RETURN(cudaMemcpy(d_voxelInput,
                                 voxelInput,
                                 sizeof(Voxel<NUM_MATERIAL>) * voxelSize,
                                 cudaMemcpyHostToDevice));
    gpuErrchk(cudaPeekAtLastError());

#ifdef PROFILING
    {
      t2 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::milli> time_span = t2 - t1;
      time_span1 += time_span.count();
    }
#endif

    UINT BlockSize = static_cast<UINT >(ceil(voxelSize * 1.0 / NUM_THREADS));

    for (UINT j = numStart; j < numEnd; j++) {

      Real energy = (idata.energyStart + j * idata.incrementEnergy);
      std::cout << "Energy = " << energy << " starting "  << "\n";

#ifdef DLEVEL2
      CUDA_CHECK_RETURN(cudaMemset(d_projectionAverage, 0, voxel[0] * voxel[1] * sizeof(Real)));
      gpuErrchk(cudaPeekAtLastError());
#endif

#ifdef  PROFILING
      tstartEnergy = std::chrono::high_resolution_clock::now();
#endif

      ElectricField eleField;
      eleField.e.x = 0;
      eleField.e.y = 0;
      eleField.e.z = 1;
      Real wavelength = static_cast<Real>(1239.84197 / energy);
      eleField.k.x = 0;
      eleField.k.y = 0;
      eleField.k.z = static_cast<Real>(2 * M_PI / wavelength);;
      Real angle;
      for (int i = 0; i < numAnglesRotation; i++) {
        angle = static_cast<Real>(i * M_PI / 180.0);
#ifdef PROFILING
        {
          t3 = std::chrono::high_resolution_clock::now();
        }
#endif

        computePolarization << < BlockSize, NUM_THREADS >>
            > (materialInput[j], d_voxelInput, eleField, angle, vx, d_polarizationX, d_polarizationY, d_polarizationZ);

        gpuErrchk(cudaPeekAtLastError());
        cudaThreadSynchronize();

#ifdef PROFILING
        {
          t4 = std::chrono::high_resolution_clock::now();
          std::chrono::duration<double, std::milli> time_span = t4 - t3;
          time_span2 += time_span.count();
        }
#endif

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

        /** FFT Computation **/
        result = cufftExecC2C(plan, d_polarizationX, d_polarizationX, CUFFT_FORWARD);
        if (result != CUFFT_SUCCESS) {
          std::cout << "CUFFT failed with result " << result << "\n";
          #pragma omp cancel parallel
          exit(EXIT_FAILURE);
        }

        result = cufftExecC2C(plan, d_polarizationY, d_polarizationY, CUFFT_FORWARD);
        if (result != CUFFT_SUCCESS) {
          std::cout << "CUFFT failed with result " << result << "\n";
          #pragma omp cancel parallel
          exit(EXIT_FAILURE);
        }

        result = cufftExecC2C(plan, d_polarizationZ, d_polarizationZ, CUFFT_FORWARD);
        if (result != CUFFT_SUCCESS) {
          std::cout << "CUFFT failed with result " << result << "\n";
          exit(EXIT_FAILURE);
        }

#ifdef PROFILING
        {
          t5 = std::chrono::high_resolution_clock::now();
          std::chrono::duration<double, std::milli> time_span = t5 - t4;
          time_span3 += time_span.count();
        }
#endif

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

        /** Scatter 3D computation **/
        computeScatter3D << < BlockSize, NUM_THREADS >>
            > (d_polarizationX, d_polarizationY, d_polarizationZ, d_scatter3D, eleField, voxelSize, vx, idata.physSize);
        cudaThreadSynchronize();
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
          t6 = std::chrono::high_resolution_clock::now();
          std::chrono::duration<double, std::milli> time_span = t6 - t5;
          time_span4 += time_span.count();
        }
#endif

#ifdef EOC
        CUDA_CHECK_RETURN(cudaMemcpy(scatter3D, d_scatter3D, sizeof(Real) * voxelSize, cudaMemcpyDeviceToHost));
        gpuErrchk(cudaPeekAtLastError());

#ifdef PROFILING
        {
          t7 = std::chrono::high_resolution_clock::now();
          std::chrono::duration<double, std::milli> time_span =  t7 - t6;
          time_span5 += time_span.count();
        }
#endif
        computeEwaldProjectionCPU(projectionCPU, scatter3D, vx, eleField.k.x);
#else
        UINT BlockSize2 = static_cast<UINT>(ceil(voxel[0] * voxel[1] * 1.0 / NUM_THREADS));
        computeEwaldProjectionGPU << < BlockSize2, NUM_THREADS >> > (d_projection, d_scatter3D, vx, eleField.k.z,idata.physSize);
        cudaThreadSynchronize();
        const double alpha = cos(angle);
        const double beta = sin(angle);

        const double coeffs[2][3]{
            alpha, beta, static_cast<Real>(((1 - alpha) * voxel[0] / 2 - beta * voxel[0] / 2.)),
            -beta, alpha, static_cast<Real>(beta * voxel[0] / 2. + (1 - alpha) * voxel[0] / 2.)
        };

        cudaMemset(d_rotProjection, 0, voxel[0] * voxel[1] * sizeof(Real));

        NppStatus status = nppiWarpAffine_32f_C1R(d_projection,
                                                  sizeImage,
                                                  voxel[0] * sizeof(Real),
                                                  rect,
                                                  d_rotProjection,
                                                  voxel[0] * sizeof(Real),
                                                  rect,
                                                  coeffs,
                                                  NPPI_INTER_LINEAR);
        if (status != NPP_SUCCESS) {
          std::cout << "Image rotation failed with error = " << status << "\n";
          exit(-1);
        }
#ifdef DLEVEL2
        const Real factor = static_cast<Real>(1.0);
        stat = cublasSaxpy(handle, voxel[0] * voxel[1], &factor, d_rotProjection, 1, d_projectionAverage, 1);
        if (stat != CUBLAS_STATUS_SUCCESS) {
          std::cout << "CUBLAS during sum failed  with status " << stat << "\n";
          exit(EXIT_FAILURE);
        }
#endif

#ifdef PROFILING
        {
          t7 = std::chrono::high_resolution_clock::now();
          std::chrono::duration<double, std::milli> time_span = t7 - t6;
          time_span5 += time_span.count();
        }
#endif
#ifdef DLEVEL1
        CUDA_CHECK_RETURN(cudaMemcpy(&projectionGPU[i*voxel[0] * voxel[1]],
                                     d_rotProjection,
                                     sizeof(Real) * (voxel[0] * voxel[1]),
                                     cudaMemcpyDeviceToHost));
        gpuErrchk(cudaPeekAtLastError());
#endif
#endif



#ifdef PROFILING
        {
          t8 = std::chrono::high_resolution_clock::now();
          std::chrono::duration<double, std::milli> time_span = t8 - t7;
          time_span6 += time_span.count();
        }
#endif
      }

#ifdef DLEVEL2
      /// The averaging out for all angles
      const Real alphaFac = static_cast<Real>(1.0 / numAnglesRotation);
      stat = cublasSscal(handle, voxel[0] * voxel[1], &alphaFac, d_projectionAverage, 1);
      if (stat != CUBLAS_STATUS_SUCCESS) {
        std::cout << "CUBLAS during averaging failed  with status " << stat << "\n";
        exit (EXIT_FAILURE);
      }

      CUDA_CHECK_RETURN(cudaMemcpy(&projectionGPUAveraged[j * voxel[0] * voxel[1]],
                                   d_projectionAverage,
                                   sizeof(Real) * (voxel[0] * voxel[1]),
                                   cudaMemcpyDeviceToHost));
      gpuErrchk(cudaPeekAtLastError());
#endif

#ifdef PROFILING
      {
        t9 = std::chrono::high_resolution_clock::now();
      }
#endif

#ifdef DLEVEL1
#pragma omp parallel
      {

        UINT threadID = omp_get_thread_num();
        UINT chunkSize = static_cast<UINT>(std::ceil(numAnglesRotation * 1.0/ omp_get_num_threads()));
        UINT voxelSize2D = voxel[0] * voxel[1];
        Real *projectionLocal = new Real[voxelSize2D];
        UINT startID = threadID * chunkSize;
        UINT endID = ((threadID + 1) * chunkSize);
        if (startID < numAnglesRotation) {

          for (UINT csize = startID; csize < std::min(endID, numAnglesRotation); csize++) {

#ifdef EOC
            Real rotAngle = csize*idata.incrementAngle;
            std::memcpy(projectionLocal,&projectionCPU[threadID], sizeof(Real)*voxel[0]*voxel[1]);
            rotateImage(projectionLocal, vx, rotAngle);
#else
            std::memcpy(projectionLocal, &projectionGPU[csize * voxelSize2D], sizeof(Real) * voxelSize2D);
#endif


            std::string dirname = "Projection/";
  //          std::string fname = dirname + "projection" + std::to_string(j) + "_" + std::to_string(csize);
  //          VTI::writeDataScalar2D(projectionLocal, voxel, fname.c_str(), "projection");
            std::string fnameFP = dirname + "projectionFP" + std::to_string(j) + "_" + std::to_string(csize);
            VTI::writeDataScalar2DFP(projectionLocal, voxel, fnameFP.c_str(), "projection");

          }
        }
        delete[] projectionLocal;
      }
#endif

#ifdef PROFILING
      {
        t10 = std::chrono::high_resolution_clock::now();
        tendEnergy = t10;
        std::chrono::duration<double, std::milli> time_span = t10 - t9;
        time_span8 += time_span.count();
        time_span = tendEnergy - tstartEnergy;
        time_spanEnergy += time_span.count();
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
#ifndef EOC
    CUDA_CHECK_RETURN(cudaFree(d_projection));
    gpuErrchk(cudaPeekAtLastError());
    CUDA_CHECK_RETURN(cudaFree(d_rotProjection));
    gpuErrchk(cudaPeekAtLastError());
#endif

    delete[] polarizationX;
    delete[] polarizationY;
    delete[] polarizationZ;
    delete[] scatter3D;


    cufftDestroy(plan);
    cublasDestroy(handle);

  }

#ifdef DLEVEL2

#ifdef PROFILING
  {
    tDumpLevel2Start = std::chrono::high_resolution_clock::now();
  }
#endif
  omp_set_num_threads(idata.num_threads);
  dump_files2D(projectionGPUAveraged,numEnergyLevel,voxel);
#ifdef PROFILING
  {
    tDumpLevel2End = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> time_span = tDumpLevel2End - tDumpLevel2Start;
    time_spanDumpLevel2 += time_span.count();
  }
#endif
#endif


#ifdef PROFILING
  std::cout << "[TIMER] Malloc on CPU + GPU = " << time_span0 << "\n";
  std::cout << "[TIMER] Memcopy CPU -> GPU = " << time_span1 << "\n";
  std::cout << "[TIMER] Polarization computation (GPU) = " << time_span2 << "\n";
  std::cout << "[TIMER] FFT computation (GPU) = " << time_span3 << "\n";
  std::cout << "[TIMER] computeScatter3D = " << time_span4 << "\n";
#ifdef EOC
  std::cout << "[TIMER] Memcopy GPU -> CPU = " << time_span5 << "\n";
  std::cout << "[TIMER] Ewald Projection (CPU)= " << time_span6 << "\n";
#else
  std::cout << "[TIMER] Ewald Projection (GPU)=  " << time_span5 << "\n";
  std::cout << "[TIMER] Memcopy GPU -> CPU (2D) =" << time_span6 << "\n";
#endif
  std::cout << "[TIMER] Rotation of image (CPU)= " << time_span7 << "\n";
#ifdef DLEVEL1
  std::cout << "[TIMER] File Writing (CPU)= " << time_span8 << "\n";
#endif
#ifdef DLEVEL2
  std::cout << "[TIMER] File Writing (CPU)= " << time_spanDumpLevel2 << "\n";
#endif

  std::cout << "[TIMER] Total time for all Energy calculation = " << time_spanEnergy << "\n";



#endif

#ifdef EOC
  delete[] projectionCPU;
#else
#ifdef DLEVEL1
  delete[] projectionGPU;
#endif
#ifdef DLEVEL2
  delete[] projectionGPUAveraged;
#endif

#endif

  return (EXIT_SUCCESS);
}

