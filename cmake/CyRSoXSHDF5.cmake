include(ExternalProject)

if(SKBUILD)
    set(_cyrsoxs_fetch_hdf5_default ON)
else()
    set(_cyrsoxs_fetch_hdf5_default OFF)
endif()

option(CYRSOXS_FETCH_HDF5 "Build HDF5 from source when system HDF5 is not found" ${_cyrsoxs_fetch_hdf5_default})

set(CYRSOXS_HDF5_URL "https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5_1.14.6.tar.gz" CACHE STRING "HDF5 source archive URL when bundling")
set(CYRSOXS_HDF5_URL_HASH "SHA256=09ee1c671a87401a5201c06106650f62badeea5a3b3941e9b1e2e1e08317357f" CACHE STRING "Expected hash for CYRSOXS_HDF5_URL")

set(CYRSOXS_HDF5_EP_TARGET "")
set(CYRSOXS_HDF5_RPATH "")

find_package(HDF5 COMPONENTS CXX HL QUIET)

if(HDF5_FOUND)
    message(STATUS "Using system HDF5 (${HDF5_INCLUDE_DIRS})")
else()
    if(NOT CYRSOXS_FETCH_HDF5)
        find_package(HDF5 COMPONENTS CXX HL REQUIRED)
        message(STATUS "Using system HDF5 (${HDF5_INCLUDE_DIRS})")
    else()
        if(WIN32)
            message(FATAL_ERROR "CYRSOXS_FETCH_HDF5 is not supported on Windows; install HDF5 and set HDF5_ROOT or CMAKE_PREFIX_PATH.")
        endif()
        set(CYRSOXS_HDF5_PREFIX "${CMAKE_BINARY_DIR}/hdf5-install")
        file(MAKE_DIRECTORY "${CYRSOXS_HDF5_PREFIX}")
        set(CYRSOXS_HDF5_EP_TARGET "cyrsoxs_hdf5_ep")
        message(STATUS "HDF5 not found; building static HDF5 into ${CYRSOXS_HDF5_PREFIX} (first build may take several minutes)")
        set(_hdf5_lib_pre "lib")
        set(_hdf5_lib_suf ".a")
        set(_hdf5_lib_dir "${CYRSOXS_HDF5_PREFIX}/lib")
        set(_byproducts
                "${_hdf5_lib_dir}/${_hdf5_lib_pre}hdf5${_hdf5_lib_suf}"
                "${_hdf5_lib_dir}/${_hdf5_lib_pre}hdf5_cpp${_hdf5_lib_suf}"
                "${_hdf5_lib_dir}/${_hdf5_lib_pre}hdf5_hl${_hdf5_lib_suf}"
                "${_hdf5_lib_dir}/${_hdf5_lib_pre}hdf5_hl_cpp${_hdf5_lib_suf}")
        ExternalProject_Add(${CYRSOXS_HDF5_EP_TARGET}
                INSTALL_DIR "${CYRSOXS_HDF5_PREFIX}"
                URL "${CYRSOXS_HDF5_URL}"
                URL_HASH "${CYRSOXS_HDF5_URL_HASH}"
                CMAKE_ARGS
                -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
                -DCMAKE_BUILD_TYPE=Release
                -DCMAKE_POSITION_INDEPENDENT_CODE=ON
                -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                -DBUILD_SHARED_LIBS=OFF
                -DBUILD_TESTING=OFF
                -DHDF5_BUILD_EXAMPLES=OFF
                -DHDF5_BUILD_TOOLS=OFF
                -DHDF5_BUILD_CPP_LIB=ON
                -DHDF5_BUILD_HL_LIB=ON
                -DHDF5_ENABLE_PARALLEL=OFF
                -DHDF5_ENABLE_THREADSAFE=OFF
                -DHDF5_ENABLE_Z_LIB_SUPPORT=OFF
                -DHDF5_ENABLE_SZIP_SUPPORT=OFF
                -DCMAKE_INSTALL_LIBDIR=lib
                BUILD_BYPRODUCTS ${_byproducts})
        set(HDF5_INCLUDE_DIR "${CYRSOXS_HDF5_PREFIX}/include")
        set(HDF5_INCLUDE_DIRS "${CYRSOXS_HDF5_PREFIX}/include")
        set(HDF5_CXX_LIBRARIES
                "${_hdf5_lib_dir}/${_hdf5_lib_pre}hdf5_cpp${_hdf5_lib_suf}"
                "${_hdf5_lib_dir}/${_hdf5_lib_pre}hdf5${_hdf5_lib_suf}")
        set(HDF5_HL_LIBRARIES
                "${_hdf5_lib_dir}/${_hdf5_lib_pre}hdf5_hl_cpp${_hdf5_lib_suf}"
                "${_hdf5_lib_dir}/${_hdf5_lib_pre}hdf5_hl${_hdf5_lib_suf}")
    endif()
endif()

if(HDF5_FOUND AND CYRSOXS_HDF5_EP_TARGET STREQUAL "")
    if(NOT HDF5_INCLUDE_DIR AND HDF5_INCLUDE_DIRS)
        set(HDF5_INCLUDE_DIR "${HDF5_INCLUDE_DIRS}")
    endif()
endif()

if(HDF5_FOUND AND CYRSOXS_HDF5_EP_TARGET STREQUAL "")
    if(DEFINED HDF5_LIBRARY_DIRS AND HDF5_LIBRARY_DIRS)
        set(CYRSOXS_HDF5_RPATH "${HDF5_LIBRARY_DIRS}")
    else()
        foreach(_lib IN LISTS HDF5_CXX_LIBRARIES HDF5_HL_LIBRARIES)
            if(EXISTS "${_lib}")
                get_filename_component(_hdf5_lib_dir "${_lib}" DIRECTORY)
                list(APPEND CYRSOXS_HDF5_RPATH "${_hdf5_lib_dir}")
            endif()
        endforeach()
        if(CYRSOXS_HDF5_RPATH)
            list(REMOVE_DUPLICATES CYRSOXS_HDF5_RPATH)
        endif()
    endif()
endif()

function(cyrsoxs_target_link_hdf5 target_name)
    if(CYRSOXS_HDF5_EP_TARGET AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
        target_link_libraries(${target_name} PRIVATE
                "-Wl,--whole-archive"
                ${HDF5_CXX_LIBRARIES}
                ${HDF5_HL_LIBRARIES}
                "-Wl,--no-whole-archive")
    elseif(CYRSOXS_HDF5_EP_TARGET AND APPLE)
        foreach(_lib IN ITEMS ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES})
            target_link_libraries(${target_name} PRIVATE "-Wl,-force_load,${_lib}")
        endforeach()
    else()
        target_link_libraries(${target_name} PRIVATE ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES})
    endif()
endfunction()

function(cyrsoxs_target_apply_hdf5_rpath target_name)
    if(CYRSOXS_HDF5_RPATH AND NOT CYRSOXS_HDF5_EP_TARGET)
        set_target_properties(${target_name} PROPERTIES
                BUILD_RPATH "${CYRSOXS_HDF5_RPATH}"
                INSTALL_RPATH "${CYRSOXS_HDF5_RPATH}"
                INSTALL_RPATH_USE_LINK_PATH TRUE)
    endif()
endfunction()
