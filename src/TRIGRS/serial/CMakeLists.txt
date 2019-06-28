project(landslide_triggers LANGUAGES Fortran)

# ${PROJECT_NAME} is the name given for the most recent project() call
add_executable(${PROJECT_NAME} "${CMAKE_CURRENT_LIST_DIR}/trigrs.f90"
                               "${CMAKE_CURRENT_LIST_DIR}/flux.f90"
                               "${CMAKE_CURRENT_LIST_DIR}/prpijz.f90"
                               "${CMAKE_CURRENT_LIST_DIR}/svxmdv.f90"
                               "${CMAKE_CURRENT_LIST_DIR}/svijz.f90"
                               "${CMAKE_CURRENT_LIST_DIR}/dzero_brac.f90"
                               "${CMAKE_CURRENT_LIST_DIR}/grids.f95"
                               "${CMAKE_CURRENT_LIST_DIR}/input_vars.f95"
                               "${CMAKE_CURRENT_LIST_DIR}/model_vars.f95"
                               "${CMAKE_CURRENT_LIST_DIR}/dsimps.f95"
                               "${CMAKE_CURRENT_LIST_DIR}/input_file_defs.f95"
                               "${CMAKE_CURRENT_LIST_DIR}/iverson.f95"
                               "${CMAKE_CURRENT_LIST_DIR}/pstpi.f95"
                               "${CMAKE_CURRENT_LIST_DIR}/satfin.f95"
                               "${CMAKE_CURRENT_LIST_DIR}/savage.f95"
                               "${CMAKE_CURRENT_LIST_DIR}/steady.f95"
                               "${CMAKE_CURRENT_LIST_DIR}/trini.f95"
                               "${CMAKE_CURRENT_LIST_DIR}/unsinf.f95"
                               "${CMAKE_CURRENT_LIST_DIR}/ivestp.f95"
                               "${CMAKE_CURRENT_LIST_DIR}/pstpf.f95"
                               "${CMAKE_CURRENT_LIST_DIR}/rnoff.f95"
                               "${CMAKE_CURRENT_LIST_DIR}/satinf.f95"
                               "${CMAKE_CURRENT_LIST_DIR}/smallt.f95"
                               "${CMAKE_CURRENT_LIST_DIR}/svgstp.f95"
                               "${CMAKE_CURRENT_LIST_DIR}/unsfin.f95"
                               "${CMAKE_CURRENT_LIST_DIR}/unsth.f95"
                               "${CMAKE_CURRENT_LIST_DIR}/ssizgrd.f95"
                               "${CMAKE_CURRENT_LIST_DIR}/svlist.f95"
                               "${CMAKE_CURRENT_LIST_DIR}/calerf.f"
                               "${CMAKE_CURRENT_LIST_DIR}/dbsct.f"
                               "${CMAKE_CURRENT_LIST_DIR}/derfc.f"
                               "${CMAKE_CURRENT_LIST_DIR}/irdgrd.f"
                               "${CMAKE_CURRENT_LIST_DIR}/irdswm.f"
                               "${CMAKE_CURRENT_LIST_DIR}/isvgrd.f"
                               "${CMAKE_CURRENT_LIST_DIR}/roots.f"
                               "${CMAKE_CURRENT_LIST_DIR}/srdgrd.f"
                               "${CMAKE_CURRENT_LIST_DIR}/srdswm.f"
                               "${CMAKE_CURRENT_LIST_DIR}/ssvgrd.f"
                               "${CMAKE_CURRENT_LIST_DIR}/m_writeHDF5.f90")                           


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