#test files for math module
SET( CMAKE_VERBOSE_MAKEFILE on )

IF(ENABLE_GEOMETRYCONCRETE)
IF(ENABLE_PYTHON_TESTS)
# add all examples

	ADD_TEST(geometryConcrete::CollisionDetection	${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/CollisionDetection.py)	
	ADD_TEST(geometryConcrete::GradingCurve3D		${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/GradingCurve3D.py)		 
	ADD_TEST(geometryConcrete::GradingCurveFileIO	${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/GradingCurveFileIO.py ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})		 



# set PYTHONPATH environment for all examples

	SET(MYTESTS  
	geometryConcrete::CollisionDetection
	geometryConcrete::GradingCurveFileIO
	geometryConcrete::GradingCurve3D
	)   

	SET_PROPERTY(TEST ${MYTESTS}  PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/src")    

IF(ENABLE_MECHANICS)
	ADD_TEST(geometryConcrete::GmshExport			${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/GmshExport.py ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
	SET_PROPERTY(TEST geometryConcrete::GmshExport  PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/src")		 
ENDIF(ENABLE_MECHANICS)



ENDIF(ENABLE_PYTHON_TESTS)




ENDIF(ENABLE_GEOMETRYCONCRETE)