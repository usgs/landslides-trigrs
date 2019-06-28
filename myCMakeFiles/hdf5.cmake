if ("$ENV{HDF5_ROOT}" STREQUAL "")
    message(FATAL_ERROR "Could not find HDF5 install directory.  Please create an environment variable HDF5_ROOT to the install directory.")
else()
    set(HDF5_ROOT "$ENV{HDF5_ROOT}" CACHE INTERNAL "Get the HDF5 install directory")
    message("Using HDF5 installed at ${HDF5_ROOT}")

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I${HDF5_ROOT}/include")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${HDF5_ROOT}/lib -lhdf5hl_fortran -lhdf5_fortran")
endif()
