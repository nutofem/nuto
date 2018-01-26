add_unit_test(EigenSolver)
if(SUITESPARSE_FOUND)
    target_link_libraries(EigenSolver SuiteSparse::UmfPack SuiteSparse::Cholmod)
endif()
