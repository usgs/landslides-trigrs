# =============================
# Link in MPI
# =============================
find_package(MPI REQUIRED)

message(MPI_FOUND)

if (MPI_FOUND)

    project(landslide_triggers_mpi LANGUAGES Fortran)

    # ${PROJECT_NAME} is the name given for the most recent project() call
    add_executable(${PROJECT_NAME} "${CMAKE_CURRENT_LIST_DIR}/trigrs_p.f90"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/flux.f90"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/prpijz.f90"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/svxmdv.f90"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/svijz.f90"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/dzero_brac.f90"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/grids.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/input_vars.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/model_vars.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/dsimps.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/input_file_defs.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/iverson.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/pstpi.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/satfin.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/savage.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/steady.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/trini.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/unsinf.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/ivestp.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/pstpf.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/rnoff.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/satinf.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/smallt.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/svgstp.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/unsfin.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/unsth.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/ssizgrd.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/svlist.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/calerf.f"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/dbsct.f"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/derfc.f"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/irdgrd.f"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/irdswm.f"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/isvgrd.f"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/roots.f"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/srdgrd.f"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/srdswm.f"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/ssvgrd.f"
                                "${CMAKE_CURRENT_LIST_DIR}/../serial/m_writeHDF5.f90"
                                "${CMAKE_CURRENT_LIST_DIR}/partial_p.f90"
                                "${CMAKE_CURRENT_LIST_DIR}/flux_p.f90"
                                "${CMAKE_CURRENT_LIST_DIR}/modules_p.f90"
                                "${CMAKE_CURRENT_LIST_DIR}/trini_p.f90"
                                "${CMAKE_CURRENT_LIST_DIR}/irdgrd_p.f90"
                                "${CMAKE_CURRENT_LIST_DIR}/irdswm_p.f90"
                                "${CMAKE_CURRENT_LIST_DIR}/srdgrd_p.f90"
                                "${CMAKE_CURRENT_LIST_DIR}/srdswm_p.f90"
                                "${CMAKE_CURRENT_LIST_DIR}/pstpi_p.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/satfin_p.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/steady_p.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/unsinf_p.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/ivestp_p.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/pstpf_p.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/rnoff_p.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/satinf_p.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/svgstp_p.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/unsfin_p.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/unsth_p.f95"
                                "${CMAKE_CURRENT_LIST_DIR}/ssizgrd_p.f95")


    include_directories(SYSTEM ${MPI_INCLUDE_PATH})

    include_directories(${MPI_Fortran_ADDITIONAL_INCLUDE_DIRS})

    target_link_libraries(${PROJECT_NAME} ${MPI_Fortran_LIBRARIES})

    if(MPI_COMPILE_FLAGS)
      set_target_properties(${PROJECT_NAME} PROPERTIES
        COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    endif()

    if(MPI_LINK_FLAGS)
      set_target_properties(${PROJECT_NAME} PROPERTIES
        LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
    endif()


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


endif (MPI_FOUND)






