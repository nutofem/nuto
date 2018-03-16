if (ARPACK_FOUND)
    add_integrationtest(Arpack)
    target_link_libraries(Arpack Arpack::Arpack)
endif()

add_integrationtest(PolynomialLeastSquaresFitting)
