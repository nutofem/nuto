#!/usr/bin/python

import nuto
import sys
import os

from time import time

# definitions
Width = 5.
Height = 10.
NumElementsX = 13
NumElementsY = 25

YoungsModulus = 20000.
PoissonsRatio = 0.2

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% no error in file, modified, if error is detected              %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error = False

#create structure
numDim=2
myStructure = nuto.Structure(numDim)

myStructure.SetShowTime(False)

# create material law
myMatLin = myStructure.ConstitutiveLawCreate("LinearElastic")
myStructure.ConstitutiveLawSetYoungsModulus(myMatLin, YoungsModulus)
myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin, PoissonsRatio)

#create section
mySection = myStructure.SectionCreate("Plane_Strain")
myStructure.SectionSetThickness(mySection,1)

#create nodes
nodeCoordinates = nuto.DoubleFullMatrix(2,1)
node = 0;
for yCount in range(0,NumElementsY + 1):
	nodeCoordinates.SetValue(1,0, yCount * Height/NumElementsY)
	for xCount in range(0,NumElementsX + 1):
		nodeCoordinates.SetValue(0,0, xCount * Width/NumElementsX)
		#print "node: " + str(node) + " coordinates: " + str(nodeCoordinates.GetValue(0,0)) + "," + str(nodeCoordinates.GetValue(1,0))
		myStructure.NodeCreate(node, "displacements", nodeCoordinates)
		node += 1

#create elements
elementIncidence = nuto.IntFullMatrix(4,1)
element = 0
for yCount in range(0, NumElementsY):
	for xCount in range(0, NumElementsX):
		node1 = yCount * (NumElementsX + 1) + xCount
		elementIncidence.SetValue(0,0,node1)
		elementIncidence.SetValue(1,0,node1 + 1)
		elementIncidence.SetValue(2,0,node1 + NumElementsX + 2)
		elementIncidence.SetValue(3,0,node1 + NumElementsX + 1)
		#print "element: " + str(element) + " incidence: "
		#elementIncidence.Info()
		myStructure.ElementCreate(element, "plane2d4n", elementIncidence, "CONSTITUTIVELAWIPCRACK" , "NOIPDATA")
		myStructure.ElementSetConstitutiveLaw(element,myMatLin)
		myStructure.ElementSetSection(element,mySection)
		element += 1

#Crack initiation
crackPoints=nuto.DoubleFullMatrix(2,2,(	1.6, 5 ,
										6  , 5 ))
print crackPoints
crackPoints.Info(5)
CrackNodes=nuto.IntFullMatrix(2,1)
CrackNodes=myStructure.NodesCreate("coordinates", crackPoints)
myStructure.NodeInfo(5);
print CrackNodes
CrackNodes.Info(5)
crack1= myStructure.CrackCreate(CrackNodes)
print "crack1=", crack1
CrackNodes.Info(5)

myStructure.CrackInfo(5);

#PhantomNode
print "initiate PhantomNode"
myStructure.InitiatePhantomNodeMethod(1000);
print "initiate PhantomNode done"

myStructure.CrackInfo(5);

myStructure.ElementTotalSetConstitutiveLaw(myMatLin)
myStructure.ElementTotalSetSection(mySection)

# visualize results
myStructure.AddVisualizationComponentDisplacements()
myStructure.AddVisualizationComponentElement()
myStructure.AddVisualizationComponentSection()
myStructure.AddVisualizationComponentCracks()
myStructure.AddVisualizationComponentEngineeringStrain()
myStructure.AddVisualizationComponentEngineeringStress()
myStructure.ExportVtkDataFile("PhantomNodeMethodPlane2D4N.vtk")

# boundary conditions
# x-direction
direction = nuto.DoubleFullMatrix(2,1)
direction.SetValue(0,0 , 1)
direction.SetValue(1,0 , 0)
myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0)
# y-direction bottom
direction.SetValue(0,0 , 0)
direction.SetValue(1,0 , 1)
for xCount in range(0,NumElementsX):
	node = xCount
	myStructure.ConstraintLinearSetDisplacementNode(node, direction, 0.0)
# y-direction top
direction.SetValue(0,0 , 0)
direction.SetValue(1,0 , 1)

for xCount in range(0,NumElementsX):
	node = xCount + (NumElementsX + 1) * NumElementsY
	myStructure.ConstraintLinearSetDisplacementNode(node, direction, 0.5)

# start analysis
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
myStructure.ExportVtkDataFile("PhantomNodeMethodPlane2D4N-deformed.vtk");

if (error):
    sys.exit(-1)
else:
    sys.exit(0)
   
##################################################################
