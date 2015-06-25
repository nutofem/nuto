import nuto
import sys
import os

#show the results on the screen
printResult = False

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% no error in file, modified, if error is detected              %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error = False

import nuto
import sys
import os

printResult = True

# definitions
E = 80000
v = 1./3.

BoundaryDisplacement = 1.0

def RunPatchTest(StressState):
	global error
	#create structure
	myStructure = nuto.Structure(2)
	myStructure.SetShowTime(False)

	# create material law
	myMatLin = myStructure.ConstitutiveLawCreate("LinearElasticEngineeringStress")
	myStructure.ConstitutiveLawSetYoungsModulus(myMatLin, E)
	myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin, v)

	#create section
	mySection = myStructure.SectionCreate("Plane_Stress")
	myStructure.SectionSetThickness(mySection,1)

	numNodes=8
	myStructure.NodesCreate(nuto.DoubleFullMatrix(2,numNodes,(
		 0 ,  0 ,
		10 ,  0 ,
		 2 ,  2 ,
		 8 ,  3 ,
		 4 ,  7 ,
		 8 ,  7 ,
		 0 , 10 ,
		10 , 10	)) )

	elementIncidence = nuto.IntFullMatrix(4,5,(	
		3,2,0,1 ,
		4,6,0,2 ,
		5,4,2,3 ,
		7,5,3,1 ,
		7,6,4,5 ) )

	interpolationType = myStructure.InterpolationTypeCreate("Quad2D")
	myStructure.InterpolationTypeAdd(interpolationType, "Coordinates", "Equidistant1");
	myStructure.InterpolationTypeAdd(interpolationType, "Displacements", "Equidistant1");

	myStructure.ElementsCreate(interpolationType, elementIncidence)
	myStructure.ElementTotalConvertToInterpolationType()
	myStructure.ElementTotalSetConstitutiveLaw(myMatLin)
	myStructure.ElementTotalSetSection(mySection)

	LoadNodesXPos=myStructure.GroupCreate("Nodes")
	LoadNodesXNeg=myStructure.GroupCreate("Nodes")
	LoadNodesYPos=myStructure.GroupCreate("Nodes")
	LoadNodesYNeg=myStructure.GroupCreate("Nodes")

	myStructure.GroupAddNode(LoadNodesXPos,1)
	myStructure.GroupAddNode(LoadNodesXPos,7)
	myStructure.GroupAddNode(LoadNodesYPos,6)
	myStructure.GroupAddNode(LoadNodesYPos,7)
	
	myStructure.GroupAddNode(LoadNodesXNeg,0)
	myStructure.GroupAddNode(LoadNodesXNeg,6)
	myStructure.GroupAddNode(LoadNodesYNeg,0)
	myStructure.GroupAddNode(LoadNodesYNeg,1)

	directionX = nuto.DoubleFullMatrix(2,1,(1,0))
	directionY = nuto.DoubleFullMatrix(2,1,(0,1))

	print "Displacement control with stress state ", StressState

	if StressState == "XX":
		myStructure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesXNeg, directionX, 0.0)
		myStructure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesXPos, directionX, BoundaryDisplacement)
		myStructure.ConstraintLinearSetDisplacementNode(0, directionY, 0)
	elif StressState == "YY":
		myStructure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesYNeg, directionY, 0.0)
		myStructure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesYPos, directionY, BoundaryDisplacement)
		myStructure.ConstraintLinearSetDisplacementNode(0, directionX, 0)
	elif StressState == "XY":
		myStructure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesXNeg, directionX, 0.0)
		myStructure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesXNeg, directionY, 0.0)
		myStructure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesXPos, directionY, BoundaryDisplacement)
		myStructure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesXPos, directionX, 0)
	else:
		print 'Wrong stressstate given'
		error = True;
		sys.exit(-1)

	# start analysis
	# build global dof numbering
	myStructure.NodeBuildGlobalDofs()
	myStructure.CalculateMaximumIndependentSets()

	# build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
	stiffnessMatrixCSRVector2 = nuto.DoubleSparseMatrixCSRVector2General(0,0)
	dispForceVector = nuto.DoubleFullVector()
	myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector)
	stiffnessMatrix = nuto.DoubleSparseMatrixCSRGeneral(stiffnessMatrixCSRVector2)

	# calculate right hand side
	rhsVector = dispForceVector

	# solve
	mySolver = nuto.SparseDirectSolverMUMPS()
	displacementVector = nuto.DoubleFullVector()
	stiffnessMatrix.SetOneBasedIndexing()
	mySolver.Solve(stiffnessMatrix, rhsVector, displacementVector)

	# write displacements to node
	myStructure.NodeMergeActiveDofValues(displacementVector)

	# calculate residual
	intForceVector = nuto.DoubleFullVector()
	myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector)
	if (printResult):
		residualVector = intForceVector*(-1.)
		print "residual: " + str(residualVector.Norm())
		
	if ((intForceVector).Norm()>1e-8):
			print 'Internal and external forces differs.'
			error = True;



	#calculate internal force vector
	Fi = nuto.DoubleFullVector()
	rowIndex = nuto.IntFullVector()
	myStructure.ElementGradientInternalPotential(2,Fi,rowIndex)
	if (printResult):
		print "Internal Force"
		Fi.Info()

	#calculate engineering strain of 2 at all integration points
	#the size the matrix is not important and reallocated within the procedure
	EngineeringStrain = nuto.DoubleFullMatrix(0,0)
	myStructure.ElementGetEngineeringStrain(2, EngineeringStrain)


	#calculate engineering strain of 2 at all integration points
	EngineeringStress = nuto.DoubleFullMatrix(0,0)
	myStructure.ElementGetEngineeringStress(2, EngineeringStress)
	#correct stress
	EngineeringStressCorrect = nuto.DoubleFullMatrix(6,4)
	for i in range(0,4):
		if StressState == "XX":
			sigma = E*BoundaryDisplacement / 10.
			EngineeringStressCorrect.SetValue(0,i, sigma)
		elif StressState == "YY":
			sigma = E*BoundaryDisplacement / 10.
			EngineeringStressCorrect.SetValue(1,i, sigma)
		elif StressState == "XY":
			sigma = E*BoundaryDisplacement / 10. / (2+2*v)
			EngineeringStressCorrect.SetValue(3,i, sigma)



	if (printResult):
		print "EngineeringStressCorrect"
		EngineeringStressCorrect.Info()
		print "EngineeringStress"
		EngineeringStress.Info()

	if ((EngineeringStress-EngineeringStressCorrect).Abs().Max()>1e-4):
			print 'stress is not correct.'
			error = True;
			if (printResult):
				print "(EngineeringStress-EngineeringStressCorrect).Abs().Max()="
				print (EngineeringStress-EngineeringStressCorrect).Abs().Max()


	C1 = 1./E
	C2 = -v/E
	C3 = 2*(1+v)/E

	D = nuto.DoubleFullMatrix(6,6,(
		C1, C2, C2,  0,  0,  0,
		C2, C1, C2,  0,  0,  0,
		C2, C2, C1,  0,  0,  0,
		0,  0,  0, C3,  0,  0,
		0,  0,  0,  0, C3,  0,
		0,  0,  0,  0,  0, C3))

	#correct strain
	EngineeringStrainCorrect = D*EngineeringStressCorrect

	if (printResult):
		print "EngineeringStrainCorrect"
		EngineeringStrainCorrect.Info()
		print "EngineeringStrain"
		EngineeringStrain.Info()

	if ((EngineeringStrain-EngineeringStrainCorrect).Abs().Max()>1e-8):
			print 'strain is not correct.'
			error = True;



RunPatchTest("XX")
RunPatchTest("YY")
RunPatchTest("XY")


if (error):
	sys.exit(-1)
else:
    sys.exit(0)
   
