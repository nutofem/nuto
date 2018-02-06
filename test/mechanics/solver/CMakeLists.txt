add_unit_test(EigenSolver)
if(SUITESPARSE_FOUND)
    target_link_libraries(EigenSolver SuiteSparse::UmfPack SuiteSparse::Cholmod)
endif()

if(MUMPS_FOUND)
    add_unit_test(MumpsSolver)
    target_link_libraries(MumpsSolver EigenMumpsSupport)
endif()
