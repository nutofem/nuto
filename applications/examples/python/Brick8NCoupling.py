# -*- coding: utf-8 -*-
import nuto
import numpy as np

structure = nuto.Structure(3)

# create material law
Material1 = structure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress")
structure.ConstitutiveLawSetParameterDouble(Material1, "Youngs_Modulus", 20000.)
structure.ConstitutiveLawSetParameterDouble(Material1, "Poissons_Ratio", 0.2)

# create nodes
structure.NodeCreate( 1, np.array([0.0, 0.0, 0.0]))
structure.NodeCreate( 2, np.array([1.0, 0.0, 0.0]))
structure.NodeCreate( 3, np.array([2.0, 0.0, 0.0]))
structure.NodeCreate( 4, np.array([0.0, 1.0, 0.0]))
structure.NodeCreate( 5, np.array([1.0, 1.0, 0.0]))
structure.NodeCreate( 6, np.array([2.0, 1.0, 0.0]))
structure.NodeCreate( 7, np.array([0.0, 2.0, 0.0]))
structure.NodeCreate( 8, np.array([1.0, 2.0, 0.0]))
structure.NodeCreate( 9, np.array([2.0, 2.0, 0.0]))

structure.NodeCreate(10, np.array([0.0, 0.0, 1.0]))
structure.NodeCreate(11, np.array([1.0, 0.0, 1.0]))
structure.NodeCreate(12, np.array([2.0, 0.0, 1.0]))
structure.NodeCreate(13, np.array([0.0, 1.0, 1.0]))
structure.NodeCreate(14, np.array([1.0, 1.0, 1.0]))
structure.NodeCreate(15, np.array([2.0, 1.0, 1.0]))
structure.NodeCreate(16, np.array([0.0, 2.0, 1.0]))
structure.NodeCreate(17, np.array([1.0, 2.0, 1.0]))
structure.NodeCreate(18, np.array([2.0, 2.0, 1.0]))

structure.NodeCreate(19, np.array([0.0, 0.0, 2.0]))
structure.NodeCreate(20, np.array([1.0, 0.0, 2.0]))
structure.NodeCreate(21, np.array([2.0, 0.0, 2.0]))
structure.NodeCreate(22, np.array([0.0, 1.0, 2.0]))
structure.NodeCreate(23, np.array([1.0, 1.0, 2.0]))
structure.NodeCreate(24, np.array([2.0, 1.0, 2.0]))
structure.NodeCreate(25, np.array([0.0, 2.0, 2.0]))
structure.NodeCreate(26, np.array([1.0, 2.0, 2.0]))
structure.NodeCreate(27, np.array([2.0, 2.0, 2.0]))

structure.NodeCreate(28, np.array([4.0, 0.0, 0.0]))
structure.NodeCreate(29, np.array([4.0, 2.0, 0.0]))
structure.NodeCreate(30, np.array([4.0, 0.0, 2.0]))
structure.NodeCreate(31, np.array([4.0, 2.0, 2.0]))

interpolationType = structure.InterpolationTypeCreate("Brick3D")
structure.InterpolationTypeAdd(interpolationType, "Coordinates", "Equidistant1")
structure.InterpolationTypeAdd(interpolationType, "Displacements", "Equidistant1")

# create elements
structure.ElementCreate(interpolationType, [  1,  2,  5,  4, 10, 11, 14, 13])
structure.ElementCreate(interpolationType, [  2,  3,  6,  5, 11, 12, 15, 14])
structure.ElementCreate(interpolationType, [  4,  5,  8,  7, 13, 14, 17, 16])
structure.ElementCreate(interpolationType, [  5,  6,  9,  8, 14, 15, 18, 17])
structure.ElementCreate(interpolationType, [ 10, 11, 14, 13, 19, 20, 23, 22])
structure.ElementCreate(interpolationType, [ 11, 12, 15, 14, 20, 21, 24, 23])
structure.ElementCreate(interpolationType, [ 13, 14, 17, 16, 22, 23, 26, 25])
structure.ElementCreate(interpolationType, [ 14, 15, 18, 17, 23, 24, 27, 26])
structure.ElementCreate(interpolationType, [  3, 28, 29,  9, 21, 30, 31, 27])

structure.ElementTotalSetConstitutiveLaw(Material1)
structure.ElementTotalConvertToInterpolationType()

# boundary conditions
def AddBC(nodeId, directions):
    node = structure.NodeGetNodePtr(nodeId)
    structure.Constraints().Add(nuto.eDof_DISPLACEMENTS, nuto.Component(node, directions))

AddBC(1, [nuto.eDirection_X, nuto.eDirection_Y, nuto.eDirection_Z])
AddBC(4, [nuto.eDirection_X])
AddBC(7, [nuto.eDirection_X, nuto.eDirection_Z])
AddBC(10, [nuto.eDirection_X])
AddBC(13, [nuto.eDirection_X])
AddBC(16, [nuto.eDirection_X])
AddBC(19, [nuto.eDirection_X])
AddBC(22, [nuto.eDirection_X])
AddBC(25, [nuto.eDirection_X])

# coupling conditions
# direction + weight of each node with given id
nodeContributions = [
    (nuto.eDirection_X, {6: -1.0, 3: 0.5, 9: 0.5}),
    (nuto.eDirection_Y, {6: -1.0, 3: 0.5, 9: 0.5}),
    (nuto.eDirection_Z, {6: -1.0, 3: 0.5, 9: 0.5}),
    (nuto.eDirection_X, {12: -1.0, 3: 0.5, 21: 0.5}),
    (nuto.eDirection_Y, {12: -1.0, 3: 0.5, 21: 0.5}),
    (nuto.eDirection_Z, {12: -1.0, 3: 0.5, 21: 0.5}),
    (nuto.eDirection_X, {18: -1.0, 9: 0.5, 27: 0.5}),
    (nuto.eDirection_Y, {18: -1.0, 9: 0.5, 27: 0.5}),
    (nuto.eDirection_Z, {18: -1.0, 9: 0.5, 27: 0.5}),
    (nuto.eDirection_X, {24: -1.0, 21: 0.5, 27: 0.5}),
    (nuto.eDirection_Y, {24: -1.0, 21: 0.5, 27: 0.5}),
    (nuto.eDirection_Z, {24: -1.0, 21: 0.5, 27: 0.5}),
    (nuto.eDirection_X, {15: -1.0, 3: 0.25, 9: 0.25, 21: 0.25, 27: 0.25}),
    (nuto.eDirection_Y, {15: -1.0, 3: 0.25, 9: 0.25, 21: 0.25, 27: 0.25}),
    (nuto.eDirection_Z, {15: -1.0, 3: 0.25, 9: 0.25, 21: 0.25, 27: 0.25}),
    ]

for (direction, nodeContribution) in nodeContributions:
    terms = []
    for nodeId, weight in nodeContribution.items():
        node = structure.NodeGetNodePtr(nodeId)
        terms.append(nuto.Term(node, direction, weight))
    equation = nuto.Equation(terms)
    structure.Constraints().Add(nuto.eDof_DISPLACEMENTS, equation)

# forces
direction = np.array([1.0, 0.0, 0.0])
structure.LoadCreateNodeForce(28, direction, 1)
structure.LoadCreateNodeForce(29, direction, 1)
structure.LoadCreateNodeForce(30, direction, 1)
structure.LoadCreateNodeForce(31, direction, 1)

# start analysis
structure.NodeBuildGlobalDofs()

# Calculate maximum independent sets for parallelization (openmp)
structure.CalculateMaximumIndependentSets()
structure.SolveGlobalSystemStaticElastic()

# calculate residual
print("calculate residual")
# start analysis
intGradient = structure.BuildGlobalInternalGradient()
extGradient = structure.BuildGlobalExternalLoadVector()
extGradientJ = extGradient.J.Get("Displacements")
intGradientJ = intGradient.J.Get("Displacements")
intGradientK = intGradient.K.Get("Displacements")

intGradientKcast = - intGradientK
cmat = structure.GetAssembler().GetConstraintMatrix()

residual = intGradientJ + cmat.Get("Displacements", "Displacements").TransMult(intGradientKcast) - extGradientJ
print("residual: " + str(np.linalg.norm(residual)))

# visualize results
visualizationGroup = structure.GroupCreate("Elements")
structure.GroupAddElementsTotal(visualizationGroup)

structure.AddVisualizationComponent(visualizationGroup, "Displacements")
structure.AddVisualizationComponent(visualizationGroup, "EngineeringStrain")
structure.AddVisualizationComponent(visualizationGroup, "EngineeringStress")
structure.ExportVtkDataFileElements("Brick8NCoupling.vtk")
