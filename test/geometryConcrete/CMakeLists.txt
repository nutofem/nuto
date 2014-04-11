#test files for math module
SET( CMAKE_VERBOSE_MAKEFILE on )

IF(ENABLE_PYTHON_TESTS)
# add all examples

IF(ENABLE_GEOMETRYCONCRETE)
	ADD_TEST(geometryConcrete::CollisionDetection	${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/CollisionDetection.py)	
	ADD_TEST(geometryConcrete::GradingCurve3D		${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/GradingCurve3D.py)		 
	ADD_TEST(geometryConcrete::GradingCurveFileIO	${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/GradingCurveFileIO.py ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})		 
ENDIF(ENABLE_GEOMETRYCONCRETE)

# set PYTHONPATH environment for all examples

SET(MYTESTS  
geometryConcrete::CollisionDetection
geometryConcrete::GradingCurveFileIO
geometryConcrete::GradingCurve3D
)   

IF(ENABLE_GEOMETRYCONCRETE)
	SET_PROPERTY(TEST ${MYTESTS}  PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/src")    
ENDIF(ENABLE_GEOMETRYCONCRETE)

ENDIF(ENABLE_PYTHON_TESTS)