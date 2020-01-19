cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(Cy-RSoXS LANGUAGES CXX CUDA)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 ")

option(HALF_FLOAT "Use 16 bit indices for floating point" OFF)
option(VTI_BINARY "Write VTI in binary with 64 base encoding" ON)
option(ENABLE_2D "Set to ON if only 1 id in Z direction" OFF)
option(DUMP_FILES "Dump files for debugging " OFF)
option(HANNING "Hanning Windowing" OFF)
option(Profiling "Enable Profiling " OFF)
option(EOC "Ewald projection on CPU" OFF)
option(DLEVEL1 "Dump Level 1" OFF)
option(DLEVEL2 "Dump Level 2" ON)
option(BIAXIAL "Biaxial Computation" OFF)
option(BUILD_DOCS "Build Documentation" OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules")

if (HALF_FLOAT)
    add_definitions(-DHALF_FLOAT)
    message(Configuring with Half float)
endif ()


if (ENABLE_2D)
    add_definitions(-DENABLE_2D)
    message(Enabling 2D computation)
endif ()

if (DUMP_FILES)
    add_definitions(-DDUMP_FILES)
    message(Dumping all files.)
endif ()

if (PROFILING)
    add_definitions(-DPROFILING)
    message(Enabling Profiling.)
endif ()

if (EOC)
    add_definitions(-DEOC)
    find_package(OpenCV REQUIRED)
    message("Ewald on CPU")
endif ()

if (VTI_BINARY)
    add_definitions(-DVTI_BINARY)
    message("Binary writing enabled")
else ()
    add_definitions(-DVTI_ASCII)
    message("ASCII writing enabled")
endif ()

if (DLEVEL1)
    add_definitions(-DDLEVEL1)
    message("Dumping current Level 1")
endif ()
if (DLEVEL2)
    add_definitions(-DDLEVEL2)
    message("Dumping current Level 2")
endif ()

if(HANNING)
    add_definitions(-DHANNING)
    message("Switching on Hanning Windowing")
endif()
if (BIAXIAL)
    add_definitions(-DBIAXIAL)
    message("Performing Biaxial computation")
else ()
    message("Performing Uniaxial computation")
endif ()


find_package(CUDA 9 REQUIRED)
find_package(OpenMP REQUIRED)

if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

find_package(Config++)
if (NOT ${CONFIG++_FOUND})
    message(FATAL_ERROR "Libconfig++ could not be located.")
endif ()
find_package(HDF5 COMPONENTS CXX HL REQUIRED)
link_directories(${HDF5_LIBRARY_DIRS})
include_directories(${HDF5_INCLUDE_DIR})
set(LINK_LIBS ${LINK_LIBS} ${HDF5_C_${LIB_TYPE}_LIBRARY})

set(CMAKE_CXX_STANDARD 11)
set(CUDASRC
        src/cudaMain.cu)

set(CUDAINCLUDE
        include/cudaUtils.h
        include/cudaMain.h
        include/cudaHeaders.h
        include/Datatypes.h
        include/uniaxial.h
        include/Output/writeVTI.h
        include/Output/cencode.h
        include/Input/InputData.h
        )


add_library(cuda_base ${CUDAINCLUDE} ${CUDASRC})
find_library(CUDART_LIBRARY cudart ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
target_include_directories(cuda_base PUBLIC include ${CONFIG++_INCLUDE_DIR})
if (EOC)
    target_include_directories(cuda_base PRIVATE ${OpenCV_INCLUDE_DIRS})
endif ()
target_compile_features(cuda_base PUBLIC cxx_std_11)
target_include_directories(cuda_base PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
set_target_properties(cuda_base
        PROPERTIES CUDA_SEPARABLE_COMPILATION ON
        )
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -fopenmp")
#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_30,code=sm_30")
if (NUM_MATERIAL)
    target_compile_definitions(cuda_base PUBLIC NUM_MATERIAL=${NUM_MATERIAL})
endif ()
set(PRS_INC
        include/Input/Input.h
        include/Input/InputData.h)
set(PRS_SRC
        src/main.cpp
        )

add_executable(Cy-RSoXS ${PRS_INC} ${PRS_SRC})

set_target_properties(
        Cy-RSoXS
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
)


if (EOC)
    target_link_libraries(Cy-RSoXS cuda_base ${OpenCV_LIBS} ${HDF5_CXX_LIBRARIES} -lcufft -lcublas -lz "${CUDA_nppc_LIBRARY};${CUDA_nppial_LIBRARY};${CUDA_nppicc_LIBRARY};${CUDA_nppicom_LIBRARY};${CUDA_nppidei_LIBRARY};${CUDA_nppif_LIBRARY};${CUDA_nppig_LIBRARY};${CUDA_nppim_LIBRARY};${CUDA_nppist_LIBRARY};${CUDA_nppisu_LIBRARY};${CUDA_nppitc_LIBRARY};${CUDA_npps_LIBRARY}" ${CONFIG++_LIBRARY})
else ()
    target_link_libraries(Cy-RSoXS cuda_base ${HDF5_CXX_LIBRARIES} -lcufft -lcublas  "${CUDA_nppc_LIBRARY};${CUDA_nppial_LIBRARY};${CUDA_nppicc_LIBRARY};${CUDA_nppicom_LIBRARY};${CUDA_nppidei_LIBRARY};${CUDA_nppif_LIBRARY};${CUDA_nppig_LIBRARY};${CUDA_nppim_LIBRARY};${CUDA_nppist_LIBRARY};${CUDA_nppisu_LIBRARY};${CUDA_nppitc_LIBRARY};${CUDA_npps_LIBRARY}" ${CONFIG++_LIBRARY})
endif ()

if (BUILD_DOCS)
    find_package(Doxygen)
    if (DOXYGEN_FOUND)
        # set input and output files
        set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

        # request to configure the file
        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
        message("Doxygen build started")

        # note the option ALL which allows to build the docs together with the application
        add_custom_target(doc_doxygen ALL
                COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                COMMENT "Generating API documentation with Doxygen"
                VERBATIM)
    else (DOXYGEN_FOUND)
        message("Doxygen need to be installed to generate the doxygen documentation")
    endif (DOXYGEN_FOUND)
endif ()