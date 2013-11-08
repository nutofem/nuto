#!/usr/bin/python

import nuto
import sys
import os

createResult = True
printResult = True

from time import time

# definitions
YoungsModulus = 1.e5
PoissonsRatio = 0.2

Force = 5

#ElementType = "PLANE2D3N"
#IntegragionType = "2D3NGauss3Ip"

ElementType = "PLANE2D4N"
#IntegragionType = "2D4NGauss4Ip"
#IntegragionType2 = IntegragionType
IntegragionType = "2D4NConst9Ip"
IntegragionType2 = "2D4NMod4Ip"

#StressState = "XX"
StressState = "YY"
#StressState = "XY"

EnableDisplacementControl = True
#EnableDisplacementControl = False
BoundaryDisplacement = 1.0
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% no error in file, modified, if error is detected              %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error = False

#create structure
numDim=2
myStructure = nuto.Structure(numDim)

# create material law
myMatLin = myStructure.ConstitutiveLawCreate("LinearElasticEngineeringStress")
myStructure.ConstitutiveLawSetYoungsModulus(myMatLin, YoungsModulus)
myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin, PoissonsRatio)

#create section
mySection = myStructure.SectionCreate("Plane_Strain")
myStructure.SectionSetThickness(mySection,1)

numNodes=8
myStructure.NodesCreate("displacements", nuto.DoubleFullMatrix(2,numNodes,(	 0 ,  0 ,
																			10 ,  0 ,
																			 2 ,  2 ,
																			 8 ,  3 ,
																			 4 ,  7 ,
																			 8 ,  7 ,
																			 0 , 10 ,
																			10 , 10	)) )
																			

if ElementType == "PLANE2D3N":
	numElements=10
	elementIncidence = nuto.IntFullMatrix(3,numElements,(	0,1,3,
															0,2,6,
															0,3,2,
															1,7,3,
															2,4,6,
															2,3,4,
															3,5,4,
															3,7,5,
															4,5,6,
															5,7,6 ) )
elif ElementType == "PLANE2D4N":
	numElements=5
	elementIncidence = nuto.IntFullMatrix(4,numElements,(	3,2,0,1 ,
															4,6,0,2 ,
															5,4,2,3 ,
															7,5,3,1 ,
															7,6,4,5 ) )
else:
	print 'Wrong elementtype given'
	error = True;
	sys.exit(-1)

elements=myStructure.ElementsCreate(ElementType, elementIncidence)
myStructure.ElementTotalSetConstitutiveLaw(myMatLin)
myStructure.ElementTotalSetSection(mySection)

for i in range(0,elements.GetNumRows()):
	if i==2:
		myStructure.ElementSetIntegrationType(i,IntegragionType2,"NOIPDATA")
	else:
		myStructure.ElementSetIntegrationType(i,IntegragionType,"NOIPDATA")

LoadNodesXPos = myStructure.GroupCreate("Nodes")
LoadNodesXNeg = myStructure.GroupCreate("Nodes")
LoadNodesYPos = myStructure.GroupCreate("Nodes")
LoadNodesYNeg = myStructure.GroupCreate("Nodes")
directionX = nuto.DoubleFullMatrix(2,1,(1,0))
directionY = nuto.DoubleFullMatrix(2,1,(0,1))
if StressState == "XX":
	myStructure.GroupAddNode(LoadNodesXPos,1)
	myStructure.GroupAddNode(LoadNodesXPos,7)
	myStructure.ConstraintLinearSetDisplacementNode(0, nuto.DoubleFullMatrix(2,1,(1,0)), 0.0)
	myStructure.ConstraintLinearSetDisplacementNode(0, nuto.DoubleFullMatrix(2,1,(0,1)), 0.0)
	myStructure.ConstraintLinearSetDisplacementNode(6, nuto.DoubleFullMatrix(2,1,(1,0)), 0.0)
elif StressState == "YY":
	myStructure.GroupAddNode(LoadNodesYPos,6)
	myStructure.GroupAddNode(LoadNodesYPos,7)
	myStructure.ConstraintLinearSetDisplacementNode(0, nuto.DoubleFullMatrix(2,1,(1,0)), 0.0)
	myStructure.ConstraintLinearSetDisplacementNode(0, nuto.DoubleFullMatrix(2,1,(0,1)), 0.0)
	myStructure.ConstraintLinearSetDisplacementNode(1, nuto.DoubleFullMatrix(2,1,(0,1)), 0.0)
elif StressState == "XY":
	myStructure.GroupAddNode(LoadNodesXPos,6)
	myStructure.GroupAddNode(LoadNodesXPos,7)
	#~ myStructure.GroupAddNode(LoadNodesXNeg,0)
	myStructure.GroupAddNode(LoadNodesXNeg,1)
	#~ myStructure.GroupAddNode(LoadNodesYPos,1)
	myStructure.GroupAddNode(LoadNodesYPos,7)
	#~ myStructure.GroupAddNode(LoadNodesYNeg,0)
	myStructure.GroupAddNode(LoadNodesYNeg,6)
	myStructure.ConstraintLinearSetDisplacementNode(0, nuto.DoubleFullMatrix(2,1,(1,0)), 0.0)
	myStructure.ConstraintLinearSetDisplacementNode(0, nuto.DoubleFullMatrix(2,1,(0,1)), 0.0)
	myStructure.ConstraintLinearSetDisplacementNode(1, nuto.DoubleFullMatrix(2,1,(0,1)), 0.0)
else:
	print 'Wrong stressstate given'
	error = True;
	sys.exit(-1)


# boundary conditions
# (left border fixed)

if EnableDisplacementControl:
	print "Displacement control"
	# boundary displacments
	myStructure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesXPos, directionX, BoundaryDisplacement)
	myStructure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesXNeg, directionX, -1.0*BoundaryDisplacement)
	myStructure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesYPos, directionY, BoundaryDisplacement)
	myStructure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesYNeg, directionY, -1.0*BoundaryDisplacement)
else:
	print "Load control"
	myStructure.LoadCreateNodeGroupForce(1,LoadNodesXPos, directionX, Force)
	myStructure.LoadCreateNodeGroupForce(1,LoadNodesXNeg, directionX, -1*Force)
	myStructure.LoadCreateNodeGroupForce(1,LoadNodesYPos, directionY, Force)
	myStructure.LoadCreateNodeGroupForce(1,LoadNodesYNeg, directionY, -1*Force)
          

# some Infos
myStructure.Info()
myStructure.NodeInfo(5)
myStructure.ElementInfo(5)
myStructure.GroupInfo(5)

# start analysis
# build global dof numbering
curTime = time()
myStructure.NodeBuildGlobalDofs()
oldTime = curTime
curTime = time()
print "time required for dof numbering: " + str(curTime - oldTime) + " s"

# build maximum independent sets
curTime = time()
myStructure.CalculateMaximumIndependentSets()
oldTime = curTime
curTime = time()
print "time required for build maximum independent sets: " + str(curTime - oldTime) + " s"

# build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
stiffnessMatrix = nuto.DoubleSparseMatrixCSRVector2General()
dispForceVector = nuto.DoubleFullVector()
myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector)
oldTime = curTime
curTime = time()
print "time required for assembling: " + str(curTime - oldTime) + " s"

Ke = nuto.DoubleFullMatrix(0,0)
rowIndex = nuto.IntFullVector(0)
colIndex = nuto.IntFullVector(0)
for i in range( 0 , numElements):
	myStructure.ElementStiffness(i,Ke,rowIndex,colIndex)
	print "Ke"
	Ke.Info()
	
if ElementType == "PLANE2D4Nb":
	myStructure.ElementStiffness(2,Ke,rowIndex,colIndex)
	Stiffness_patchtest_slang = nuto.DoubleFullMatrix(8,8)
	Stiffness_patchtest_slang.ReadFromFile("Stiffness_patchtest_slang.txt",0," ")

	print "Stiffness NuTo"
	Ke.Info()

	print "Stiffness SLang"
	Stiffness_patchtest_slang.Info()

	DiffMat = Ke - Stiffness_patchtest_slang
	print "DiffMat Stiffness"
	DiffMat.Info()

	DiffMatAbsMax=DiffMat.Abs().Max()
	print "DiffMat.Abs().Max()=" 
	print DiffMatAbsMax


# build global external load vector
extForceVector = nuto.DoubleFullVector()
myStructure.BuildGlobalExternalLoadVector(1,extForceVector)
oldTime = curTime
curTime = time()
print "time required for building external load vector: " + str(curTime - oldTime) + " s"

# calculate right hand side
rhsVector = dispForceVector + extForceVector
oldTime = curTime
curTime = time()
print "time required for calculating right-hand-side vector: " + str(curTime - oldTime) + " s"

# solve
mySolver = nuto.SparseDirectSolverMUMPS()
displacementVector = nuto.DoubleFullVector()
stiffnessMatrixCSR = nuto.DoubleSparseMatrixCSRGeneral(stiffnessMatrix)
stiffnessMatrixCSR.SetOneBasedIndexing()
mySolver.Solve(stiffnessMatrixCSR, rhsVector, displacementVector)
oldTime = curTime
curTime = time()
print "time required for solving: " + str(curTime - oldTime) + " s"

# write displacements to node
myStructure.NodeMergeActiveDofValues(displacementVector)
oldTime = curTime
curTime = time()
print "time required for merging dof values: " + str(curTime - oldTime) + " s"

# calculate residual
intForceVector = nuto.DoubleFullVector()
myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector)
oldTime = curTime
curTime = time()
print "time required for building internal forces: " + str(curTime - oldTime) + " s"
residualVector = extForceVector - intForceVector
print "residual: " + str(residualVector.Norm())


# visualize results
myStructure.AddVisualizationComponentDisplacements()
myStructure.AddVisualizationComponentEngineeringStrain()
myStructure.AddVisualizationComponentEngineeringStress()
myStructure.ExportVtkDataFileElements( "Patchtest_" + ElementType + ".vtk")

if (error):
    sys.exit(-1)
else:
    sys.exit(0)
   
