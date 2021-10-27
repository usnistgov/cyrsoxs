cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
cmake_policy(SET CMP0048 NEW)
project(Cy-RSoXS VERSION 1.0.1.1 DESCRIPTION BetaVersion LANGUAGES CXX CUDA)


option(DOUBLE_PRECISION "Use 64 bit indices for floating point" OFF)
option(VTI_BINARY "Write VTI in binary with 64 base encoding" ON)
option(DUMP_FILES "Dump files for debugging " OFF)
option(Profiling "Enable Profiling " OFF)
option(EOC "Ewald projection on CPU" OFF)
option(BIAXIAL "Biaxial Computation" OFF)
option(BUILD_DOCS "Build Documentation" OFF)
option(PYBIND "Pybind support for Cy-RSOXS " OFF)
option(USE_64_BIT_INDICES, "Use 64 Bit Indices" OFF)
option(ENABLE_TEST, "Enable test" ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules")
if (DOUBLE_PRECISION)
    add_definitions(-DDOUBLE_PRECISION)
    message(Configuring with Double precision)
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

if (USE_64_BIT_INDICES)
    add_definitions(-DUSE_64_BIT_INDICES)
    message("Using 64 Bit indices")
else ()
    message("Using 32 Bit indices")
endif ()

if (BIAXIAL)
    add_definitions(-DBIAXIAL)
    message("Performing Biaxial computation")
else ()
    message("Performing Uniaxial computation")
endif ()

if (PYBIND)
    add_definitions(-DPYBIND)
    message("WARNING: Compiling for PyBind. This will not produce any executable")
endif ()

if (ENABLE_TEST)
    add_definitions(-DENABLE_TEST)
    message("Test enabled")
endif ()

find_package(CUDA 9 REQUIRED)
find_package(OpenMP REQUIRED)

if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

if (PYBIND)
    # find_package(PythonInterp 3.6 REQUIRED)
    # find_package(PythonLibs 3.6 REQUIRED)
    find_package(Python COMPONENTS Interpreter Development REQUIRED)

    set(PYBIND_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/external/pybind11)
    add_subdirectory(${PYBIND_LOCATION})
    include_directories(
            ${PYTHON_INCLUDE_DIRS}
            external/pybind11/include
    )
else ()
    find_package(Config++)
    if (NOT ${CONFIG++_FOUND})
        message(FATAL_ERROR "Libconfig++ could not be located.")
    endif ()
endif ()

find_package(HDF5 COMPONENTS CXX HL REQUIRED)


set(CMAKE_CXX_STANDARD 14)
set(CyRSoXS_SRC
        src/RotationMatrix.cpp
        src/cudaMain.cu)

set(CyRSoXS_INC
        include/cudaUtils.h
        include/cudaMain.h
        include/cudaHeaders.h
        include/Datatypes.h
        include/uniaxial.h
        include/Output/writeVTI.h
        include/Output/cencode.h
        include/Output/outputUtils.h
        include/Input/InputData.h
        include/Input/Input.h
        include/Output/writeH5.h
        include/utils.h
        include/Rotation.h
        include/RotationMatrix.h
        )


find_library(CUDART_LIBRARY cudart ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
include(cmake/version.cmake)
generate_version_h()

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release " FORCE)
endif ()


set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -fopenmp -DNDEBUG --use_fast_math")
if (NUM_MATERIAL)
    add_definitions(-DNUM_MATERIAL=${NUM_MATERIAL})
endif ()


if (PYBIND)

    set(PYBIND_INC
            include/PyClass/RefractiveIndex.h
            include/PyClass/VoxelData.h
            include/PyClass/ScatteringPattern.h
            include/PyClass/Polarization.h
            )
    set(PYBIND_SRC
            src/pymain.cpp
            )


    cuda_add_library(CyRSoXS SHARED
            ${PYBIND_SRC} ${CyRSoXS_SRC} ${CyRSoXS_INC})
    target_include_directories(CyRSoXS PUBLIC ${Python_INCLUDE_DIRS} ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} ${HDF5_INCLUDE_DIR} include)
    target_include_directories(CyRSoXS PRIVATE ${PROJECT_BINARY_DIR}/generated)
    target_link_libraries(CyRSoXS
            ${PYTHON_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES} -lcufft -lcublas  "${CUDA_nppc_LIBRARY};${CUDA_nppial_LIBRARY};${CUDA_nppicc_LIBRARY};${CUDA_nppidei_LIBRARY};${CUDA_nppif_LIBRARY};${CUDA_nppig_LIBRARY};${CUDA_nppim_LIBRARY};${CUDA_nppist_LIBRARY};${CUDA_nppisu_LIBRARY};${CUDA_nppitc_LIBRARY};${CUDA_npps_LIBRARY}" )
    if (EOC)
        target_include_directories(CyRSoXS PRIVATE ${OpenCV_INCLUDE_DIRS})
        target_link_libraries(CyRSoXS ${OpenCV_LIBS})
    endif ()

    set_target_properties(CyRSoXS PROPERTIES
            CUDA_SEPARABLE_COMPILATION ON
            PREFIX ""
            )

else ()
    set(EXE_SRC
            src/main.cpp
            )
    add_executable(CyRSoXS ${CyRSoXS_INC} ${CyRSoXS_SRC} ${EXE_SRC})

    set_target_properties(
            CyRSoXS
            PROPERTIES
            CUDA_SEPARABLE_COMPILATION ON
    )
    target_include_directories(CyRSoXS PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} ${HDF5_INCLUDE_DIR} ${CONFIG++_INCLUDE_DIR} include)
    target_include_directories(CyRSoXS PRIVATE ${PROJECT_BINARY_DIR}/generated)
    target_link_libraries(CyRSoXS
            ${PYTHON_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES} -lcufft -lcublas  "${CUDA_nppc_LIBRARY};${CUDA_nppial_LIBRARY};${CUDA_nppicc_LIBRARY};${CUDA_nppidei_LIBRARY};${CUDA_nppif_LIBRARY};${CUDA_nppig_LIBRARY};${CUDA_nppim_LIBRARY};${CUDA_nppist_LIBRARY};${CUDA_nppisu_LIBRARY};${CUDA_nppitc_LIBRARY};${CUDA_npps_LIBRARY}" ${CONFIG++_LIBRARY})
    if (EOC)
        target_include_directories(CyRSoXS PRIVATE ${OpenCV_INCLUDE_DIRS})
        target_link_libraries(CyRSoXS ${OpenCV_LIBS})
    endif ()
endif ()
set_property(TARGET CyRSoXS PROPERTY CUDA_ARCHITECTURES  52 53 60 61 62 70 72)

if (ENABLE_TEST)
    if (PYBIND)
        message("No tests with pybind with gtest framework")
    else ()
        add_definitions(-DCMAKE_ROOT=\"${CMAKE_SOURCE_DIR}\")
        add_subdirectory(external/googletest)
        set(TEST_SRC
                tests/tests.cpp
                )
        set(TEST_INC
                tests/include/basicFramework.h
                tests/include/testUtils.h
                tests/include/unitTest.h
                tests/include/fullFramework.h
                )
        add_executable(test ${CyRSoXS_INC} ${CyRSoXS_SRC} ${TEST_SRC} ${TEST_INC})
        target_include_directories(test PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} ${HDF5_INCLUDE_DIR} ${CONFIG++_INCLUDE_DIR} include tests/include ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
        target_include_directories(test PRIVATE ${PROJECT_BINARY_DIR}/generated)
        target_link_libraries(test gtest gtest_main
                ${HDF5_CXX_LIBRARIES} -lcufft -lcublas "${CUDA_nppc_LIBRARY};${CUDA_nppial_LIBRARY};${CUDA_nppicc_LIBRARY};${CUDA_nppidei_LIBRARY};${CUDA_nppif_LIBRARY};${CUDA_nppig_LIBRARY};${CUDA_nppim_LIBRARY};${CUDA_nppist_LIBRARY};${CUDA_nppisu_LIBRARY};${CUDA_nppitc_LIBRARY};${CUDA_npps_LIBRARY}" ${CONFIG++_LIBRARY})
    endif ()
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
