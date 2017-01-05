# -*- coding: utf-8 -*-

# load nuto package
import nuto
import numpy as np

# definitions
YoungsModulus = 20000.
Area = 100. * 100.
Length = 1000.
NumElements = 10
Force = 1.
EnableDisplacementControl = False
BoundaryDisplacement = 0.1

# create one-dimensional structure
myStructure = nuto.Structure(1)

# create section
Section1 = myStructure.SectionCreate("TRUSS")
myStructure.SectionSetArea(Section1, Area)

# create material law
Material1 = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress")
myStructure.ConstitutiveLawSetParameterDouble(Material1,"Youngs_Modulus", YoungsModulus)

# create nodes
nodeCoordinates = np.zeros(1)
for node in range(0, NumElements + 1):
    print "create node: " + str(node) + " coordinates: " + str(node * Length/NumElements)
    nodeCoordinates[0] = node * Length/NumElements
    myStructure.NodeCreateDOFs(node, "displacements", nodeCoordinates)

#create interpolation type
myInterpolationType = myStructure.InterpolationTypeCreate("Truss1D")
myStructure.InterpolationTypeAdd(myInterpolationType, "coordinates", "equidistant1")
myStructure.InterpolationTypeAdd(myInterpolationType, "displacements", "equidistant1")

# create elements
elementIncidence = nuto.IntVector(2)
for element in range(0, NumElements):
    print "create element: " + str(element) + " nodes: " + str(element) + "," + str(element+1)
    elementIncidence[0] = element
    elementIncidence[1] = element + 1
    myStructure.ElementCreate(element, myInterpolationType, elementIncidence)
    myStructure.ElementSetSection(element,Section1)
    myStructure.ElementSetConstitutiveLaw(element,Material1)

# set boundary conditions and loads
direction = nuto.DoubleFullMatrix(1,1,(1,))
myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0)
myStructure.SetNumLoadCases(1)
if EnableDisplacementControl:
    print "Displacement control"
    myStructure.ConstraintLinearSetDisplacementNode(NumElements, direction, BoundaryDisplacement)
else:
    print "Load control"
    myStructure.LoadCreateNodeForce(0,NumElements, direction, Force)

#build maximum independent sets
myStructure.CalculateMaximumIndependentSets()

# start analysis
# start analysis
myStructure.SolveGlobalSystemStaticElastic(0)
intGradient = myStructure.BuildGlobalInternalGradient()
extGradient = myStructure.BuildGlobalExternalLoadVector(0)
residual = intGradient.J.Get("Displacements") - extGradient.J.Get("Displacements")


print "residual: " + str(np.linalg.norm(residual))

## visualize results
visualizationGroup = myStructure.GroupCreate("Elements")
myStructure.GroupAddElementsTotal(visualizationGroup)

myStructure.AddVisualizationComponent(visualizationGroup, "Displacements")
myStructure.AddVisualizationComponent(visualizationGroup, "EngineeringStrain")
myStructure.AddVisualizationComponent(visualizationGroup, "EngineeringStress")
#myStructure.ExportVtkDataFileElements("Truss1D2N.vtk")
