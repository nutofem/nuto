add_subdirectory(TestClasses)
add_subdirectory(TestFiles)

add_executable(MyTests MyTests.cpp)

add_library(Numbering TestClasses/Numbering.cpp)

target_link_libraries(MyTests Numbering)


foreach(myTestFiles
    MyTruss1D2N
    MyHeatedPlate)

    add_library(${myTestFiles} TestFiles/${myTestFiles}.cpp)
    target_link_libraries(${myTestFiles} NuToMechanics NuToMath NuToVisualize)

    target_link_libraries(MyTests ${myTestFiles})

endforeach()

#add_library(MyTruss1D2N TestFiles/MyTruss1D2N.cpp)
#target_link_libraries(MyTruss1D2N NuToMechanics NuToMath NuToVisualize)


