# -*- coding: utf-8 -*-

# load nuto package
import nuto
from time import time

# definitions
YoungsModulus = 20000.
PoissonsRatio = 0.2
Width = 100.
Height = 100.
Length = 1000.
NumElementsX = 10
NumElementsY = 1
NumElementsZ = 1

EnableDisplacementControl = False
BoundaryDisplacement = 0.1
Force = 1.

# create one-dimensional structure
myStructure = nuto.Structure(3)

# create material law
Material1 = myStructure.ConstitutiveLawCreate("LinearElastic")
myStructure.ConstitutiveLawSetYoungsModulus(Material1, YoungsModulus)
myStructure.ConstitutiveLawSetPoissonsRatio(Material1, PoissonsRatio)

# create nodes
nodeCoordinates = nuto.DoubleFullMatrix(3,1)
node = 0
for zCount in range (0, NumElementsZ + 1):
    nodeCoordinates.SetValue(2,0, zCount * Height/NumElementsZ)
    for yCount in range (0, NumElementsY + 1):
        nodeCoordinates.SetValue(1,0, yCount * Width/NumElementsY)
        for xCount in range(0, NumElementsX + 1):
            nodeCoordinates.SetValue(0,0, xCount * Length/NumElementsX)
            #print "node: " + str(node) + " coordinates: " + str(nodeCoordinates.GetValue(0,0)) + "," + str(nodeCoordinates.GetValue(1,0)) + "," + str(nodeCoordinates.GetValue(2,0))
            myStructure.NodeCreate(node, "displacements", nodeCoordinates)
            node += 1

# create elements
elementIncidence = nuto.IntFullMatrix(8,1)
element = 0
for zCount in range (0, NumElementsZ):
    for yCount in range(0, NumElementsY):
        for xCount in range(0, NumElementsX):
            node1 = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1) + xCount
            elementIncidence.SetValue(0,0, node1)
            elementIncidence.SetValue(1,0, node1 + 1)
            elementIncidence.SetValue(2,0, node1 + NumElementsX + 2)
            elementIncidence.SetValue(3,0, node1 + NumElementsX + 1)
            elementIncidence.SetValue(4,0, node1 + (NumElementsX + 1) * (NumElementsY + 1))
            elementIncidence.SetValue(5,0, node1 + (NumElementsX + 1) * (NumElementsY + 1) + 1)
            elementIncidence.SetValue(6,0, node1 + (NumElementsX + 1) * (NumElementsY + 1) + NumElementsX + 2)
            elementIncidence.SetValue(7,0, node1 + (NumElementsX + 1) * (NumElementsY + 1) + NumElementsX + 1)
            #print "element: " + str(element) + " incidence: "
            #elementIncidence.Info()
            myStructure.ElementCreate(element, "Brick8N", elementIncidence)
            myStructure.ElementSetConstitutiveLaw(element,Material1)
            element += 1

# boundary conditions
direction = nuto.DoubleFullMatrix(3,1,(1,0,0))
for zCount in range (0, NumElementsZ + 1):
    for yCount in range (0, NumElementsY + 1):
        node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1)
        #print "node: " + str(node)
        myStructure.ConstraintLinearSetDisplacementNode(node, direction, 0.0)

direction = nuto.DoubleFullMatrix(3,1,(0,0,1))
myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0)
myStructure.ConstraintLinearSetDisplacementNode(NumElementsY * (NumElementsX + 1), direction, 0.0)
#print "node: " + str(NumElementsY * (NumElementsX + 1))

direction = nuto.DoubleFullMatrix(3,1,(0,1,0))
myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0)

if EnableDisplacementControl:
    print "Displacement control"
    # boundary displacments
    direction = nuto.DoubleFullMatrix(3,1,(1,0,0))
    for zCount in range (0, NumElementsZ + 1):
        for yCount in range (0, NumElementsY + 1):
            node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1) + NumElementsX
            #print "node: " + str(node)
            myStructure.ConstraintLinearSetDisplacementNode(node, direction, BoundaryDisplacement)
else:
    #load
    print "Load control"
    direction = nuto.DoubleFullMatrix(3,1,(1,0,0))

    # apply load to nodes
    for zCount in range(0, NumElementsZ + 1):
        if zCount == 0 or zCount == NumElementsZ:
            nodeForce = Force / (4 *NumElementsY * NumElementsZ)
        else:
            nodeForce = Force / (2 *NumElementsY * NumElementsZ)
        node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + NumElementsX
        #print "apply force to node: " + str(node) + " force: " + str(nodeForce)
        myStructure.LoadCreateNodeForce(node, direction, nodeForce)
        for yCount in range(1, NumElementsY):
            node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1) + NumElementsX
            print "apply force to node: " + str(node) + " force: " + str(2 * nodeForce)
            myStructure.LoadCreateNodeForce(node, direction, 2 * nodeForce)
        node = (zCount + 1) * (NumElementsX + 1) * (NumElementsY + 1) - 1
        #print "apply force to node: " + str(node) + " force: " + str(nodeForce)
        myStructure.LoadCreateNodeForce(node, direction, nodeForce)

# start analysis
# build global dof numbering
curTime = time()
myStructure.NodeBuildGlobalDofs()
oldTime = curTime
curTime = time()
print "time required for dof numbering: " + str(curTime - oldTime) + " s"

#build maximum independent sets
myStructure.CalculateMaximumIndependentSets()
oldTime = curTime
curTime = time()
print "time required for calculating maximum independent sets: " + str(curTime - oldTime) + " s"

# build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
stiffnessMatrix = nuto.DoubleSparseMatrixCSRVector2General()
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
myStructure.ExportVtkDataFile("Brick8N.vtk")
