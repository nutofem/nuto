# -*- coding: utf-8 -*-
import nuto
import numpy as np

myStructure = nuto.Structure(1)

# create section
Section1 = myStructure.SectionCreate("TRUSS")
myStructure.SectionSetArea(Section1, 1)

# create material law
Material1 = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress")
myStructure.ConstitutiveLawSetParameterDouble(Material1,"Youngs_Modulus", 1)

# create nodes
nodeCoordinates = np.zeros(1)
myStructure.NodeCreateDOFs(1, "displacements", nodeCoordinates)
nodeCoordinates[0] = 1.0
myStructure.NodeCreateDOFs(2, "displacements", nodeCoordinates)
nodeCoordinates[0] = 1.0
myStructure.NodeCreateDOFs(3, "displacements", nodeCoordinates)
nodeCoordinates[0] = 2.0
myStructure.NodeCreateDOFs(4, "displacements", nodeCoordinates)

#create interpolation type
myInterpolationType = myStructure.InterpolationTypeCreate("Truss1D")
myStructure.InterpolationTypeAdd(myInterpolationType, "coordinates", "equidistant1")
myStructure.InterpolationTypeAdd(myInterpolationType, "displacements", "equidistant1")

# create elements
elementIncidence = list(range(2))
elementIncidence[0] = 1
elementIncidence[1] = 2
e1 = myStructure.ElementCreate(myInterpolationType, elementIncidence)
myStructure.ElementSetSection(e1,Section1)
myStructure.ElementSetConstitutiveLaw(e1,Material1)
elementIncidence[0] = 3
elementIncidence[1] = 4
e2 = myStructure.ElementCreate(myInterpolationType, elementIncidence)
myStructure.ElementSetSection(e2,Section1)
myStructure.ElementSetConstitutiveLaw(e2,Material1)

# set boundary conditions and loads
direction = np.array([1.0])
myStructure.ConstraintLinearSetDisplacementNode(1, direction, 0.0)
id = myStructure.ConstraintLinearEquationCreate(2, "x_displacement", 1, 0)
myStructure.ConstraintLinearEquationAddTerm(id, 3, "x_displacement", -1)
loadCase = 0
myStructure.SetNumLoadCases(1)
myStructure.LoadCreateNodeForce(loadCase, 4, direction, 1)

#build maximum independent sets
myStructure.CalculateMaximumIndependentSets()

# start analysis
myStructure.SolveGlobalSystemStaticElastic(loadCase)
intGradient = myStructure.BuildGlobalInternalGradient()
extGradient = myStructure.BuildGlobalExternalLoadVector(loadCase)
residual = intGradient.J.Get("Displacements") - extGradient.J.Get("Displacements")


print "residual: " + str(np.linalg.norm(residual))


