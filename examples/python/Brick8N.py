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

myStructure.SetShowTime(False)

# create material law
Material1 = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress")
myStructure.ConstitutiveLawSetParameterDouble(Material1,"Youngs_Modulus", YoungsModulus)
myStructure.ConstitutiveLawSetParameterDouble(Material1,"Poissons_Ratio", PoissonsRatio)

Section1 = myStructure.SectionCreate("Volume")

# create nodes
nodeCoordinates = nuto.DoubleFullVector(3)
node = 0
for zCount in range (0, NumElementsZ + 1):
    nodeCoordinates.SetValue(2,0, zCount * Height/NumElementsZ)
    for yCount in range (0, NumElementsY + 1):
        nodeCoordinates.SetValue(1,0, yCount * Width/NumElementsY)
        for xCount in range(0, NumElementsX + 1):
            nodeCoordinates.SetValue(0,0, xCount * Length/NumElementsX)
            #print "node: " + str(node) + " coordinates: " + str(nodeCoordinates.GetValue(0,0)) + "," + str(nodeCoordinates.GetValue(1,0)) + "," + str(nodeCoordinates.GetValue(2,0))
            myStructure.NodeCreate(node, nodeCoordinates)
            node += 1
#create interpolation type
myInterpolationType = myStructure.InterpolationTypeCreate("Brick3D")
myStructure.InterpolationTypeAdd(myInterpolationType, "coordinates", "equidistant1")
myStructure.InterpolationTypeAdd(myInterpolationType, "displacements", "equidistant1")

# create elements
elementIncidence = nuto.IntFullVector(8)

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
            myStructure.ElementCreate(myInterpolationType, elementIncidence)

myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10)
myStructure.ElementTotalSetConstitutiveLaw(Material1)
myStructure.ElementTotalSetSection(Section1)



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
myStructure.SetNumLoadCases(1)

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
        myStructure.LoadCreateNodeForce(0,node, direction, nodeForce)
        for yCount in range(1, NumElementsY):
            node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1) + NumElementsX
            print "apply force to node: " + str(node) + " force: " + str(2 * nodeForce)
            myStructure.LoadCreateNodeForce(0,node, direction, 2 * nodeForce)
        node = (zCount + 1) * (NumElementsX + 1) * (NumElementsY + 1) - 1
        #print "apply force to node: " + str(node) + " force: " + str(nodeForce)
        myStructure.LoadCreateNodeForce(0,node, direction, nodeForce)

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

myStructure.SolveGlobalSystemStaticElastic()

# calculate residual
intGradient = myStructure.BuildGlobalInternalGradient()
oldTime = curTime
curTime = time()
print "time required for building internal forces: " + str(curTime - oldTime) + " s"
extGradient = myStructure.BuildGlobalExternalLoadVector(0)
extGradientJ = extGradient.J.Get("Displacements")
intGradientJ = intGradient.J.Get("Displacements")
intGradientK = intGradient.K.Get("Displacements")
# cast FullVector to FullMatrix, python does not get it...
numRows = intGradientK.GetNumRows()
intGradientKcast = nuto.DoubleFullMatrix(numRows,1)
for i in range(numRows):
    intGradientKcast.SetValue(i,0, -intGradientK.GetValue(i))
cmat = myStructure.GetConstraintMatrix()

residual = nuto.DoubleFullVector(intGradientJ + cmat.Get("Displacements", "Displacements").TransMult(intGradientKcast) - extGradientJ)
print "residual: " + str(residual.Norm())

# visualize results
visualizationGroup = myStructure.GroupCreate("Elements");
myStructure.GroupAddElementsTotal(visualizationGroup)

myStructure.AddVisualizationComponent(visualizationGroup, "Displacements");
myStructure.AddVisualizationComponent(visualizationGroup, "EngineeringStrain");
myStructure.AddVisualizationComponent(visualizationGroup, "EngineeringStress");

myStructure.ExportVtkDataFileElements("Brick8N.vtk")
