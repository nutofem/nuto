#test files for math module
SET( CMAKE_VERBOSE_MAKEFILE on )

IF(ENABLE_GEOMETRYCONCRETE)
IF(ENABLE_PYTHON_TESTS)
# add all examples

	ADD_TEST(p::geometryConcrete::CollisionDetection	${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/CollisionDetection.py)	
	ADD_TEST(p::geometryConcrete::GradingCurve3D		${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/GradingCurve3D.py)		 
	ADD_TEST(p::geometryConcrete::GradingCurveFileIO	${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/GradingCurveFileIO.py ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})		 



# set PYTHONPATH environment for all examples

	SET(MYTESTS  
	p::geometryConcrete::CollisionDetection
	p::geometryConcrete::GradingCurveFileIO
	p::geometryConcrete::GradingCurve3D
	)   

	SET_PROPERTY(TEST ${MYTESTS}  PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/src")    


	unset(GMSH_PATH CACHE)
	find_program(GMSH_PATH gmsh)
	# If gmsh is not found, the variables value will be "GMSH_PATH-NOTFOUND", which
	# is recognized by the IF()-statement
	# The link to the mechanics library is required, since the successful import
	# to a structure is tested.

	IF(GMSH_PATH)
		IF(ENABLE_MECHANICS)
		message(STATUS "Gmsh found: ${GMSH_PATH}")
		ADD_TEST(p::geometryConcrete::GmshExport			${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/GmshExport.py ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${GMSH_PATH})
		SET_PROPERTY(TEST p::geometryConcrete::GmshExport  PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/src")
		ENDIF(ENABLE_MECHANICS)		
	ELSE(GMSH_PATH)
		message(STATUS "Gmsh not found. Skipping test geometryConcrete::GmshExport")
	ENDIF(GMSH_PATH)
	

  




ENDIF(ENABLE_PYTHON_TESTS)




ENDIF(ENABLE_GEOMETRYCONCRETE)