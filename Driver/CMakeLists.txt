
################################################################################
add_executable(InverseDriver InverseDriver.f90)
target_link_libraries(InverseDriver NTPoly ${MPI_Fortran_LIBRARIES}
                      ${OpenMP_Fortran_LIBRARIES} ${BLAS_LIBRARIES}
                      ${TOOLCHAIN_LIBS})
target_include_directories(InverseDriver PRIVATE ${MPI_INCLUDE_PATH})
install(TARGETS InverseDriver DESTINATION bin)


################################################################################
add_executable(ACCConcept ACCConcept.f90)
target_link_libraries(ACCConcept NTPoly
                      ${OpenMP_Fortran_LIBRARIES} ${BLAS_LIBRARIES}
                      ${TOOLCHAIN_LIBS})
install(TARGETS ACCConcept DESTINATION bin)
