cmake_minimum_required(VERSION 3.12)
project(landslide_trigrs LANGUAGES Fortran)

# Optional parameters to use
# -DCMAKE_Fortran_COMPILER=/path/to/fortran/compiler
# -DCMAKE_BUILD_TYPE=[DEBUG, RELEASE]  - Build with debug flags or not.
# -DBUILD_SHARED_LIBS=[ON, OFF]  - Build and shared or static library, or link a program using the appropriate flags


# Add our local modules to the module path
if(NOT CMAKE_MODULE_PATH)
  set(CMAKE_MODULE_PATH)
endif()
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../myCMakeFiles")

# GNUInstallDirs is used to install libraries into correct locations
# on all platforms.
include(GNUInstallDirs)


# Setup the Fortran cmake environment
# Sets up the flags for compiling shared libraries or static libraries, and returns the libType parameter.
include("${CMAKE_MODULE_PATH}/FortranEnvironment.cmake")
include("${CMAKE_MODULE_PATH}/hdf5.cmake")

add_subdirectory(${CMAKE_SOURCE_DIR}/TRIGRS/serial)

# target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> # for generated files in build mode
#                                                   $<INSTALL_INTERFACE:include> # for client in install mode
#                                                   $<INSTALL_INTERFACE:lib> # for generated files in install mdoe
#                                            PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}) # for source files


# ================================
# Set up the make install functionality
# ================================
# 'make install' to the correct locations (provided by GNUInstallDirs).
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}-config
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include/${PROJECT_NAME})

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include/${PROJECT_NAME}
    FILES_MATCHING PATTERN "*.mod"
    PATTERN "CMakeFiles" EXCLUDE)

install(EXPORT ${PROJECT_NAME}-config DESTINATION lib/cmake)
install(FILES ../myCMakeFiles/FortranEnvironment.cmake DESTINATION lib/cmake)
install(FILES ../myCMakeFiles/gfortran_flags.cmake DESTINATION lib/cmake)
install(FILES ../myCMakeFiles/intel_flags.cmake DESTINATION lib/cmake)
export(TARGETS ${PROJECT_NAME} FILE ${PROJECT_NAME}-config.cmake)
 

# Thats it! Simple!
