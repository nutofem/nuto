SET(Trilinos_useInstalledLibs True)

if (${Trilinos_useInstalledLibs})
    SET(Trilinos_PREFIX /home/briedel/install/Trilinos)
    SET(CMAKE_PREFIX_PATH ${Trilinos_PREFIX} ${CMAKE_PREFIX_PATH})
    FIND_PACKAGE(Trilinos REQUIRED)
    INCLUDE_DIRECTORIES(${Trilinos_INCLUDE_DIRS})
    INCLUDE_DIRECTORIES(${Trilinos_TPL_INCLUDE_DIRS})

    #add_library(StructureMesh ../TestClasses/StructureMesh.cpp)

    foreach(myTestFiles
        StructureMeshTest)
    #    LoadBalancing)

        add_executable(${myTestFiles} ${myTestFiles}.cpp)
        target_link_libraries(${myTestFiles} Mechanics Math Visualize ${Trilinos_LIBRARIES} StructureMesh ConversionTools FetiLib ${MPI_Libraries})
    #    target_link_libraries(${myTestFiles} NuToMechanics NuToVisualize ${MPI_Libraries})

    endforeach()
else()
    foreach(myTestFiles
        StructureMeshTest)

    add_executable(${myTestFiles} ${myTestFiles}.cpp)
    target_link_libraries(${myTestFiles} Mechanics Math Visualize StructureMesh ConversionTools FetiLib)

    endforeach()
endif()
