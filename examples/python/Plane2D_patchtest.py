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

ElementType = "PLANE2D3N"
#ElementType = "PLANE2D4N"

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
myMatLin = myStructure.ConstitutiveLawCreate("LinearElastic")
myStructure.ConstitutiveLawSetYoungsModulus(myMatLin, YoungsModulus)
myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin, PoissonsRatio)

#create section
myStructure.SectionCreate("mySection","Plane_Strain")
myStructure.SectionSetThickness("mySection",1)

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

myStructure.ElementsCreate(ElementType, elementIncidence)
myStructure.ElementTotalSetConstitutiveLaw(myMatLin)
myStructure.ElementTotalSetSection("mySection")

myStructure.GroupCreate("LoadNodesXPos","Nodes")
myStructure.GroupCreate("LoadNodesXNeg","Nodes")
myStructure.GroupCreate("LoadNodesYPos","Nodes")
myStructure.GroupCreate("LoadNodesYNeg","Nodes")
directionX = nuto.DoubleFullMatrix(2,1,(1,0))
directionY = nuto.DoubleFullMatrix(2,1,(0,1))
if StressState == "XX":
	myStructure.GroupAddNode("LoadNodesXPos",1)
	myStructure.GroupAddNode("LoadNodesXPos",7)
	myStructure.ConstraintSetDisplacementNode(0, nuto.DoubleFullMatrix(2,1,(1,0)), 0.0)
	myStructure.ConstraintSetDisplacementNode(0, nuto.DoubleFullMatrix(2,1,(0,1)), 0.0)
	myStructure.ConstraintSetDisplacementNode(6, nuto.DoubleFullMatrix(2,1,(1,0)), 0.0)
elif StressState == "YY":
	myStructure.GroupAddNode("LoadNodesYPos",6)
	myStructure.GroupAddNode("LoadNodesYPos",7)
	myStructure.ConstraintSetDisplacementNode(0, nuto.DoubleFullMatrix(2,1,(1,0)), 0.0)
	myStructure.ConstraintSetDisplacementNode(0, nuto.DoubleFullMatrix(2,1,(0,1)), 0.0)
	myStructure.ConstraintSetDisplacementNode(1, nuto.DoubleFullMatrix(2,1,(0,1)), 0.0)
elif StressState == "XY":
	myStructure.GroupAddNode("LoadNodesXPos",6)
	myStructure.GroupAddNode("LoadNodesXPos",7)
	#~ myStructure.GroupAddNode("LoadNodesXNeg",0)
	myStructure.GroupAddNode("LoadNodesXNeg",1)
	#~ myStructure.GroupAddNode("LoadNodesYPos",1)
	myStructure.GroupAddNode("LoadNodesYPos",7)
	#~ myStructure.GroupAddNode("LoadNodesYNeg",0)
	myStructure.GroupAddNode("LoadNodesYNeg",6)
	myStructure.ConstraintSetDisplacementNode(0, nuto.DoubleFullMatrix(2,1,(1,0)), 0.0)
	myStructure.ConstraintSetDisplacementNode(0, nuto.DoubleFullMatrix(2,1,(0,1)), 0.0)
	myStructure.ConstraintSetDisplacementNode(1, nuto.DoubleFullMatrix(2,1,(0,1)), 0.0)
else:
	print 'Wrong stressstate given'
	error = True;
	sys.exit(-1)


# boundary conditions
# (left border fixed)

if EnableDisplacementControl:
	print "Displacement control"
	# boundary displacments
	myStructure.ConstraintSetDisplacementNodeGroup("LoadNodesXPos", directionX, BoundaryDisplacement)
	myStructure.ConstraintSetDisplacementNodeGroup("LoadNodesXNeg", directionX, -1.0*BoundaryDisplacement)
	myStructure.ConstraintSetDisplacementNodeGroup("LoadNodesYPos", directionY, BoundaryDisplacement)
	myStructure.ConstraintSetDisplacementNodeGroup("LoadNodesYNeg", directionY, -1.0*BoundaryDisplacement)
else:
	print "Load control"
	myStructure.LoadCreateNodeGroupForce("LoadNodesXPos", directionX, Force)
	myStructure.LoadCreateNodeGroupForce("LoadNodesXNeg", directionX, -1*Force)
	myStructure.LoadCreateNodeGroupForce("LoadNodesYPos", directionY, Force)
	myStructure.LoadCreateNodeGroupForce("LoadNodesYNeg", directionY, -1*Force)
          

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

# build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
stiffnessMatrix = nuto.DoubleSparseMatrixCSRGeneral()
dispForceVector = nuto.DoubleFullMatrix()
myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector)
oldTime = curTime
curTime = time()
print "time required for assembling: " + str(curTime - oldTime) + " s"

Ke = nuto.DoubleFullMatrix(0,0)
rowIndex = nuto.IntFullMatrix(0,0)
colIndex = nuto.IntFullMatrix(0,0)
for i in range( 0 , numElements):
	myStructure.ElementStiffness(i,Ke,rowIndex,colIndex)
	print "Ke"
	Ke.Info()
	
if ElementType == "PLANE2D4N":
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
extForceVector = nuto.DoubleFullMatrix()
myStructure.BuildGlobalExternalLoadVector(extForceVector)
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
displacementVector = nuto.DoubleFullMatrix()
stiffnessMatrix.SetOneBasedIndexing()
mySolver.Solve(stiffnessMatrix, rhsVector, displacementVector)
oldTime = curTime
curTime = time()
print "time required for solving: " + str(curTime - oldTime) + " s"

# write displacements to node
myStructure.NodeMergeActiveDofValues(displacementVector)
oldTime = curTime
curTime = time()
print "time required for merging dof values: " + str(curTime - oldTime) + " s"

# calculate residual
intForceVector = nuto.DoubleFullMatrix()
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
myStructure.ExportVtkDataFile( "Patchtest_" + ElementType + ".vtk")

if (error):
    sys.exit(-1)
else:
    sys.exit(0)
   
