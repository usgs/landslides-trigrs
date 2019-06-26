if ("$ENV{HDF5_DIR}" STREQUAL "")

    find_package(HDF5 REQUIRED COMPONENTS Fortran Fortran_HL)

    set(CMAKE_REQUIRED_INCLUDES ${HDF5_INCLUDE_DIRS} ${HDF5_Fortran_INCLUDE_DIRS})
    set(CMAKE_REQUIRED_LIBRARIES ${HDF5_Fortran_LIBRARIES} ${HDF5_Fortran_HL_LIBRARIES})

    include(CheckFortranSourceCompiles)
    check_fortran_source_compiles("use h5lt; end" hasHDF5 SRC_EXT f90)

    if(NOT hasHDF5)
        message(FATAL_ERROR "HDF5 library not working with ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
    endif()

else()
    set(HDF5_DIR "$ENV{HDF5_DIR}" CACHE INTERNAL "Get the HDF5 install directory")
    message("HDF5_DIR = ${HDF5_DIR}")

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I${HDF5_DIR}/include")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${HDF5_DIR}/lib -lhdf5hl_fortran -lhdf5_fortran")
endif()
