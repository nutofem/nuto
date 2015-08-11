import nuto
import sys
import os


lX = 365.
lY = 13.
lZ = 37.

E = 42.
nu = 0.26
deltaL = 0.6174
rho = 3.141529


error = False
errorMsg = ""

# ******************************************************************************************************
# ******************************************************************************************************
# ******************************************************************************************************

def Run(rStructure, rType, rOrder):
	global error, errorMsg
	
	resultFile = os.path.join(pathToResultFiles, rType+rOrder+"_py.vtu")
	
	print "#########################################"
	print "##   Running " + rType+":"+rOrder
	print "#########################################"
	print "##        writing vtu files to "
	print "## " + resultFile
	print "#########################################"
	

	# **********************
	# set constitutive law
	# **********************
	lawId = rStructure.ConstitutiveLawCreate("LINEARELASTICENGINEERINGSTRESS")
        rStructure.ConstitutiveLawSetParameterDouble(lawId,"YoungsModulus", E)
        rStructure.ConstitutiveLawSetParameterDouble(lawId,"PoissonsRatio",nu)
        rStructure.ConstitutiveLawSetParameterDouble(lawId,"Density",rho)
	
	rStructure.ElementTotalSetConstitutiveLaw(lawId)

	# **********************
	# set boundary conditions
	# **********************
	dimension = rStructure.GetDimension()
	directionX = nuto.DoubleFullMatrix(dimension,1)
	directionX.SetValue(0,0,1.)
		
	origin = nuto.DoubleFullVector(dimension)
	nodeGroupOrigin = rStructure.GroupCreate("Nodes");
	rStructure.GroupAddNodeRadiusRange(nodeGroupOrigin, origin, 0, 1.e-5);
	if (rStructure.GroupGetNumMembers(nodeGroupOrigin) != 1):
		errorMsg += "[SetBoundaryConditions:" +rType+":"+rOrder+"] Node at origin (0,0,0) does not exist. \n";
		error = True
		return
            
	#fix origin
	nodeOrigin = rStructure.GroupGetMemberIds(nodeGroupOrigin).GetValue(0);
	for dim in range(1, dimension):
		direction = nuto.DoubleFullMatrix(dimension,1)
		direction.SetValue(dim,0,1.)
		rStructure.ConstraintLinearSetDisplacementNode(nodeOrigin, direction, 0.0);

	# fix x = 0 plane
	nodesX0 = rStructure.GroupCreate("Nodes");
	rStructure.GroupAddNodeCoordinateRange(nodesX0, 0, -1.e-6, 1.e-6);
	rStructure.ConstraintLinearSetDisplacementNodeGroup(nodesX0, directionX, 0.);
	
	# apply displacement on x = lX plane
	nodesXlX = rStructure.GroupCreate("Nodes");
	rStructure.GroupAddNodeCoordinateRange(nodesXlX, 0, lX-1.e-6, lX+1.e-6);
	rStructure.ConstraintLinearSetDisplacementNodeGroup(nodesXlX, directionX, deltaL);
	
	# **********************
	#  SOLVE
	# **********************
	rStructure.NodeBuildGlobalDofs()
	rStructure.CalculateMaximumIndependentSets();

	stiffnessMatrixCSRVector2 = nuto.DoubleSparseMatrixCSRVector2General(0,0)
	dispForceVector = nuto.DoubleFullVector()
	rStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector)
	stiffnessMatrix = nuto.DoubleSparseMatrixCSRGeneral(stiffnessMatrixCSRVector2)

	mySolver = nuto.SparseDirectSolverMUMPS()
	displacementVector = nuto.DoubleFullVector()
	stiffnessMatrix.SetOneBasedIndexing()
	mySolver.Solve(stiffnessMatrix, dispForceVector, displacementVector)

	rStructure.NodeMergeActiveDofValues(displacementVector);
	intForceVector = nuto.DoubleFullVector()
	rStructure.BuildGlobalGradientInternalPotentialVector(intForceVector)
	if ((intForceVector).Abs().Max()>1e-8):
		errorMsg += "[SetBoundaryConditions:" +rType+":"+rOrder+"] residual force vector is not zero. \n"
		error = True
		return

	# **********************
	#  Check Solution
	# **********************

	# check stresses
	analyticStrainX = deltaL / lX
	analyticStressX = analyticStrainX * E
	allElements = rStructure.GroupCreate("Elements")
	allNodes    = rStructure.GroupCreate("Nodes")
	rStructure.GroupAddNodeCoordinateRange(allNodes, 0, -0.1, lX+0.1);
	rStructure.GroupAddElementsFromNodes(allElements, allNodes, True);
	elementIds = rStructure.GroupGetMemberIds(allElements);
	for iElement in range(0, elementIds.GetNumRows()):
		elementId = elementIds.GetValue(iElement)
		stress = nuto.DoubleFullMatrix()
		rStructure.ElementGetEngineeringStress(elementId, stress)
		for iIP in range(0, stress.GetNumColumns()):
			numericStress = stress.GetValue(0,iIP)
			if (abs(numericStress-analyticStressX) > 1.e-8):
				errorMsg += "[CheckSolution:" +rType+":"+rOrder+"] wrong stress calculation. \n"
				error = True
				return
			
	#check reaction forces
	analyticForce = analyticStressX * lY * lZ
	numericForce = 0.
	nodeX0Indices = rStructure.GroupGetMemberIds(nodesX0)
	for iNode in range(0, nodeX0Indices.GetNumRows()):
		nodeId = nodeX0Indices.GetValue(iNode);
		force = nuto.DoubleFullVector()
		rStructure.NodeInternalForce(nodeId, force);
		numericForce += force.GetValue(0);
	
	if (abs(numericForce-numericForce) > 1.e-8):
		errorMsg += "[CheckSolution:" +rType+":"+rOrder+"] wrong reaction force calculation. \n"
		error = True
		return
	
	## **********************
	##  Check Total Mass via Mass matrix
	## **********************
	
	#analyticMass = lX*lY*lZ*rho

	#numActDofs = rStructure.GetNumActiveDofs()
	#numDepDofs = rStructure.GetNumDofs() - numActDofs

	#jj = nuto.DoubleSparseMatrixCSRGeneral(numActDofs,numActDofs)
	#jk = nuto.DoubleSparseMatrixCSRGeneral(numActDofs,numDepDofs)
	#kj = nuto.DoubleSparseMatrixCSRGeneral(numDepDofs,numActDofs)
	#kk = nuto.DoubleSparseMatrixCSRGeneral(numDepDofs,numDepDofs)
	#dummy = nuto.DoubleFullVector()

	#nutoStructureBaseEnumMass = 2

	#rStructure.BuildGlobalCoefficientSubMatricesGeneral(nutoStructureBaseEnumMass, jj, jk, kj, kk)

	#jjFull = nuto.DoubleFullMatrix(jj)
	#jkFull = nuto.DoubleFullMatrix(jk)
	#kjFull = nuto.DoubleFullMatrix(kj)
	#kkFull = nuto.DoubleFullMatrix(kk)
	
	#numericMass = 0.
	#numericMass += jjFull.Sum()
	#numericMass += jkFull.Sum()
	#numericMass += kjFull.Sum()
	#numericMass += kkFull.Sum()

	#numericMass = numericMass / dimension # since the mass is added to nodes in every direction

	#if(abs(numericMass - analyticMass)/numericMass > 1.e-6 ):
		#print "mass analytical : ", analyticMass
		#print "mass numerical  : ", numericMass
		#errorMsg += "[CheckSolution:" +rType+":"+rOrder+"] wrong mass calculation. \n"
		#error = True
		#return
	
	# **********************
	#  Visualize
	# **********************
	
	rStructure.AddVisualizationComponentDisplacements();
	rStructure.AddVisualizationComponentEngineeringStrain();
	rStructure.AddVisualizationComponentEngineeringStress();

	rStructure.ExportVtkDataFileElements(resultFile,True);
	
	
# ******************************************************************************************************
# ******************************************************************************************************
# ******************************************************************************************************


def Run3D(r3DShape, rTypeOrder):
	global error, errorMsg
	myStructure = nuto.Structure(3)
	myStructure.SetShowTime(False)

	numElementsX = 5
	numElementsY = 2
	numElementsZ = 2

	# create nodes
	numNodesX = numElementsX+1
	numNodesY = numElementsY+1
	numNodesZ = numElementsZ+1
		
	deltaX = lX/numElementsX
	deltaY = lY/numElementsY
	deltaZ = lZ/numElementsZ
	
	nodeNum = 0
	for countZ in range(0, numNodesZ):
		for countY in range(0, numNodesY):
			for countX in range(0, numNodesX):
				myStructure.NodeCreate(nodeNum, nuto.DoubleFullVector((countX*deltaX,countY*deltaY, countZ*deltaZ)))
				nodeNum += 1

	myInterpolationType = myStructure.InterpolationTypeCreate(r3DShape);
	myStructure.InterpolationTypeAdd(myInterpolationType, "Coordinates","EQUIDISTANT1");
	myStructure.InterpolationTypeAdd(myInterpolationType, "Displacements", rTypeOrder);

	#create elements
	for countZ in range(0, numElementsZ):
		for countY in range(0, numElementsY):
			for countX in range(0, numElementsX):
				nodes = nuto.IntFullVector(8)
				nodes.SetValue(0, countX  +  countY   *numNodesX +  countZ    * numNodesX * numNodesY)
				nodes.SetValue(1, countX+1+  countY   *numNodesX +  countZ    * numNodesX * numNodesY)
				nodes.SetValue(2, countX+1+ (countY+1)*numNodesX +  countZ    * numNodesX * numNodesY)
				nodes.SetValue(3, countX  + (countY+1)*numNodesX +  countZ    * numNodesX * numNodesY)
				nodes.SetValue(4, countX  +  countY   *numNodesX + (countZ+1) * numNodesX * numNodesY)
				nodes.SetValue(5, countX+1+  countY   *numNodesX + (countZ+1) * numNodesX * numNodesY)
				nodes.SetValue(6, countX+1+ (countY+1)*numNodesX + (countZ+1) * numNodesX * numNodesY)
				nodes.SetValue(7, countX  + (countY+1)*numNodesX + (countZ+1) * numNodesX * numNodesY)
				if (r3DShape == "Brick3D"):
					myStructure.ElementCreate(myInterpolationType, nodes);
				elif (r3DShape == "Tetrahedron3D"):
					myStructure.ElementCreate(myInterpolationType, nuto.IntFullVector((nodes.GetValue(0),nodes.GetValue(1),nodes.GetValue(3),nodes.GetValue(7))));
					myStructure.ElementCreate(myInterpolationType, nuto.IntFullVector((nodes.GetValue(0),nodes.GetValue(1),nodes.GetValue(7),nodes.GetValue(4))));
					myStructure.ElementCreate(myInterpolationType, nuto.IntFullVector((nodes.GetValue(5),nodes.GetValue(4),nodes.GetValue(7),nodes.GetValue(1))));
					myStructure.ElementCreate(myInterpolationType, nuto.IntFullVector((nodes.GetValue(6),nodes.GetValue(5),nodes.GetValue(7),nodes.GetValue(1))));
					myStructure.ElementCreate(myInterpolationType, nuto.IntFullVector((nodes.GetValue(2),nodes.GetValue(7),nodes.GetValue(1),nodes.GetValue(6))));
					myStructure.ElementCreate(myInterpolationType, nuto.IntFullVector((nodes.GetValue(2),nodes.GetValue(3),nodes.GetValue(1),nodes.GetValue(7))));
				else:
					error = True
					errorMsg += "Element shape " + r3DShape +" is invalid. \n"
					return
				

	myStructure.ElementTotalConvertToInterpolationType();

	mySection = myStructure.SectionCreate("Volume");
	myStructure.ElementTotalSetSection(mySection);
	
	Run(myStructure, r3DShape , rTypeOrder);
	

# ******************************************************************************************************
# ******************************************************************************************************
# ******************************************************************************************************

def Run2D(r2DShape, rTypeOrder):
	global error, errorMsg
	myStructure = nuto.Structure(2)
	myStructure.SetShowTime(False)

	numElementsX = 4
	numElementsY = 2

	# create nodes
	numNodesX = numElementsX+1
	numNodesY = numElementsY+1
	deltaX = lX/numElementsX
	deltaY = lY/numElementsY

	nodeNum = 0
	for countY in range(0, numNodesY):
		for countX in range(0, numNodesX):
			myStructure.NodeCreate(nodeNum, nuto.DoubleFullVector((countX*deltaX,countY*deltaY)))
			nodeNum += 1

	myInterpolationType = myStructure.InterpolationTypeCreate(r2DShape);
	myStructure.InterpolationTypeAdd(myInterpolationType, "Coordinates","EQUIDISTANT1");
	myStructure.InterpolationTypeAdd(myInterpolationType, "Displacements", rTypeOrder);

	#create elements
	for countY in range(0, numElementsY):
		for countX in range(0, numElementsX):
			if (r2DShape == "Quad2D"):
				nodes = nuto.IntFullVector(4)
				nodes.SetValue(0, countX  +  countY   *numNodesX)
				nodes.SetValue(1, countX+1+  countY   *numNodesX)
				nodes.SetValue(2, countX+1+ (countY+1)*numNodesX)
				nodes.SetValue(3, countX  + (countY+1)*numNodesX)
				myStructure.ElementCreate(myInterpolationType, nodes);
			elif (r2DShape == "Triangle2D"):
				nodes = nuto.IntFullVector(3)
				nodes.SetValue(0, countX  +  countY   *numNodesX)
				nodes.SetValue(1, countX+1+  countY   *numNodesX)
				nodes.SetValue(2, countX+1+ (countY+1)*numNodesX)
				myStructure.ElementCreate(myInterpolationType, nodes);
				
				nodes = nuto.IntFullVector(3)
				nodes.SetValue(0, countX  +  countY   *numNodesX)
				nodes.SetValue(1, countX+1+ (countY+1)*numNodesX)
				nodes.SetValue(2, countX  + (countY+1)*numNodesX)
				myStructure.ElementCreate(myInterpolationType, nodes);
			else:
				error = True
				errorMsg += "Element shape " + r2DShape +" is invalid. \n"
				return
				

	myStructure.ElementTotalConvertToInterpolationType();

	mySection = myStructure.SectionCreate("Plane_Stress");
	myStructure.SectionSetThickness(mySection, lZ);
	myStructure.ElementTotalSetSection(mySection);
	
	Run(myStructure, r2DShape , rTypeOrder);
	

# ******************************************************************************************************
# ******************************************************************************************************
# ******************************************************************************************************


def Run1D(r1DShape, rTypeOrder):
	global error, errorMsg
	
	myStructure = nuto.Structure(1)
	myStructure.SetShowTime(False)

	numElementsX = 3

	# create nodes
	numNodesX = numElementsX+1
	deltaX = lX/numElementsX


	for countX in range(0, numNodesX):
		coord =  nuto.DoubleFullVector(1)
		coord.SetValue(0,countX*deltaX)
		myStructure.NodeCreate(countX, coord)

	myInterpolationType = myStructure.InterpolationTypeCreate(r1DShape);
	myStructure.InterpolationTypeAdd(myInterpolationType, "Coordinates","EQUIDISTANT1");
	myStructure.InterpolationTypeAdd(myInterpolationType, "Displacements", rTypeOrder);

	#create elements
	for countX in range(0, numElementsX):
		if (r1DShape == "Truss1D"):
			nodes = nuto.IntFullVector(2)
			nodes.SetValue(0, countX  )
			nodes.SetValue(1, countX+1)
			myStructure.ElementCreate(myInterpolationType, nodes);
		else:
			error = True
			errorMsg += "Element shape " + r2DShape +" is invalid. \n"
			return

	myStructure.ElementTotalConvertToInterpolationType();

	mySection = myStructure.SectionCreate("Truss");
	myStructure.SectionSetArea(mySection, lZ*lY);
	myStructure.ElementTotalSetSection(mySection);



	Run(myStructure, r1DShape , rTypeOrder);


# ******************************************************************************************************
# ******************************************************************************************************
# ******************************************************************************************************



#show the results on the screen
printResult = False

#path in the original source directory and current filename at the end
pathToResultFiles = os.path.join(sys.argv[3],"Results"+os.path.basename(sys.argv[0]))

#remove the extension
fileExt = os.path.splitext(sys.argv[0])[1]
pathToResultFiles = pathToResultFiles.replace(fileExt,'')

if not os.path.exists(pathToResultFiles):
    os.makedirs(pathToResultFiles)



# ******************************************************************************************************
# ******************************************************************************************************
# ******************************************************************************************************

Run3D("Brick3D", "Equidistant1")
Run3D("Brick3D", "Equidistant2")

Run3D("Tetrahedron3D", "Equidistant1")
Run3D("Tetrahedron3D", "Equidistant2")

Run2D("Quad2D", "Equidistant1")
Run2D("Quad2D", "Equidistant2")
Run2D("Quad2D", "Lobatto2")
Run2D("Quad2D", "Lobatto3")
Run2D("Quad2D", "Lobatto4")

Run2D("Triangle2D", "Equidistant1")
Run2D("Triangle2D", "Equidistant2")
Run2D("Triangle2D", "Equidistant3")
Run2D("Triangle2D", "Equidistant4")

Run1D("Truss1D", "Equidistant1")
Run1D("Truss1D", "Equidistant2")
Run1D("Truss1D", "Equidistant3")
Run1D("Truss1D", "Equidistant4")
Run1D("Truss1D", "Lobatto2")
Run1D("Truss1D", "Lobatto3")
Run1D("Truss1D", "Lobatto4")

        
if (error):
	print "### \n \n  FAILED \n \n ### "
	print errorMsg
	sys.exit(-1)
else:
	sys.exit(0)
