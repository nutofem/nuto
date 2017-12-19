SET(Trilinos_useInstalledLibs False)

if (${Trilinos_useInstalledLibs})
    SET(Trilinos_PREFIX /home/briedel/install/Trilinos)
    SET(CMAKE_PREFIX_PATH ${Trilinos_PREFIX} ${CMAKE_PREFIX_PATH})
    FIND_PACKAGE(Trilinos REQUIRED)
else()
    SET(Trilinos_INCLUDE_DIRS /usr/include/trilinos)
    SET(Trilinos_TPL_INCLUDE_DIRS /usr/include)

#    SET(Trilinos_Libs "trilinos_epetra;trilinos_tpetra")

##     foreach(trilinosLib ${Trilinos_Libs})
#    foreach(trilinosLib
#            trilinos_epetra)
#        #trilinos_tpetra)

#        FIND_LIBRARY(FOUND_TRILINOS_LIB ${trilinosLib})
#        SET(MyTrilinos_LIBRARIES ${MyTrilinos_LIBRARIES} ${FOUND_TRILINOS_LIB})
#        #UNSET(${FOUND_TRILINOS_LIB} CACHE)
#    endforeach()

#    MESSAGE(${MyTrilinos_LIBRARIES})

#    UNSET(FTL)
#    FIND_LIBRARY(FTL2 NAMES trilinos NAMES_PER_DIR 4)
#    MESSAGE(STATUS ${FTL2})

    SET(MyTrilinos_LIBRARIES
        trilinos_epetra
        trilinos_tpetra
        trilinos_teuchoscomm
        trilinos_teuchoscore
        trilinos_teuchoskokkoscomm
        trilinos_teuchoskokkoscompat
        trilinos_teuchosnumerics
        trilinos_teuchosparameterlist
        trilinos_kokkoscore
        trilinos_amesos2
        trilinos_amesos
        trilinos_aztecoo
        trilinos_ifpack2
        trilinos_ifpack
        trilinos_belos
        trilinos_belosepetra
        trilinos_belostpetra
        trilinos_zoltan
        )
endif()



#INCLUDE_DIRECTORIES(${Trilinos_INCLUDE_DIRS})
#INCLUDE_DIRECTORIES(${Trilinos_TPL_INCLUDE_DIRS})


#add_library(StructureMesh ../TestClasses/StructureMesh.cpp)

foreach(myTestFiles
    StructureMeshTest
    ZoltanTest)
#    LoadBalancing)

    add_executable(${myTestFiles} ${myTestFiles}.cpp)
    #target_link_libraries(${myTestFiles} Mechanics Math Visualize ${Trilinos_LIBRARIES} StructureMesh ConversionTools FetiLib ${MPI_Libraries})
    target_link_libraries(${myTestFiles} Mechanics Math Visualize StructureMesh ConversionTools FetiLib ${MPI_Libraries} ${MyTrilinos_LIBRARIES})
    target_include_directories(${myTestFiles} PRIVATE ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
#    target_link_libraries(${myTestFiles} NuToMechanics NuToVisualize ${MPI_Libraries})

endforeach()
