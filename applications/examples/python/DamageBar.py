import nuto
import sys
import numpy as np

structure = nuto.Structure(1)
structure.SetNumTimeDerivatives(0)

interpolationType = structure.InterpolationTypeCreate("Truss1D")
structure.InterpolationTypeAdd(interpolationType, "Coordinates", "Equidistant1")
structure.InterpolationTypeAdd(interpolationType, "DISPLACEMENTS", "EQUIDISTANT2")
structure.InterpolationTypeAdd(interpolationType, "NONLOCALEQSTRAIN", "EQUIDISTANT1")

coordinates = np.zeros(1)

length = 100.0
weakened_zone_length = 10.0
area = 10.0
alpha = 0.1

n_elements = 100
delta_l = length / n_elements

totalSection = nuto.SectionTruss.Create(area)
weakenedSection = nuto.SectionTruss.Create((1.0 - alpha)*area)

nodeIDs = list(range(2))
nodeIDs[0] = structure.NodeCreate(coordinates)
for i in range(n_elements):
    coordinates[0] = (i+1)*delta_l
    nodeIDs[1] = structure.NodeCreate(coordinates)
    element = structure.ElementCreate(interpolationType, nodeIDs)
    if (coordinates[0] < (length - weakened_zone_length)/2.0 or coordinates[0] > (length + weakened_zone_length)/2.0):
        structure.ElementSetSection(element, totalSection)
    else:
        structure.ElementSetSection(element, weakenedSection)
    nodeIDs[0] = nodeIDs[1]

structure.ElementTotalConvertToInterpolationType()

damage = structure.ConstitutiveLawCreate("GRADIENT_DAMAGE_ENGINEERING_STRESS")
structure.ConstitutiveLawSetParameterDouble(damage, "DENSITY", 1.0)
structure.ConstitutiveLawSetParameterDouble(damage, "YOUNGS_MODULUS", 20000)
structure.ConstitutiveLawSetParameterDouble(damage, "POISSONS_RATIO", 0.2)
structure.ConstitutiveLawSetParameterDouble(damage, "NONLOCAL_RADIUS", 1)
structure.ConstitutiveLawSetParameterDouble(damage, "TENSILE_STRENGTH", 4.)
structure.ConstitutiveLawSetParameterDouble(damage, "COMPRESSIVE_STRENGTH", 4. * 10)
structure.ConstitutiveLawSetDamageLaw(damage, nuto.DamageLawExponential.Create(4./20000, 4./0.021))

structure.ElementTotalSetConstitutiveLaw(damage)

visualizationGroup = structure.GroupCreate("Elements")
structure.GroupAddElementsTotal(visualizationGroup)

structure.AddVisualizationComponent(visualizationGroup, "Displacements")
structure.AddVisualizationComponent(visualizationGroup, "EngineeringStrain")
structure.AddVisualizationComponent(visualizationGroup, "EngineeringStress")
structure.AddVisualizationComponent(visualizationGroup, "PrincipalEngineeringStress")
structure.AddVisualizationComponent(visualizationGroup, "NonlocalEqStrain")

firstNode = structure.NodeGetAtCoordinate(0)
structure.Constraints().Add(nuto.eDof_DISPLACEMENTS, nuto.Component(firstNode, [nuto.eDirection_X]))

lastNode = structure.NodeGetAtCoordinate(length)
structure.Constraints().Add(nuto.eDof_DISPLACEMENTS, nuto.Component(lastNode, [nuto.eDirection_X], lambda t: 0.05 * t))

newmark = nuto.NewmarkDirect(structure)
newmark.SetTimeStep(0.1)
newmark.SetResultDirectory("damage_bar_results", True)
newmark.SetAutomaticTimeStepping(True)
newmark.AddResultNodeDisplacements("TopDisplacement", nodeIDs[0])
groupID = structure.GroupCreate("Nodes")
structure.GroupAddNode(groupID, nodeIDs[0])
newmark.AddResultGroupNodeForce("TopForce", groupID)
newmark.Solve(1.0)

sol = structure.NodeExtractDofValues(0)

nonlocaleqstrain = sol.J.Get("Nonlocaleqstrain")

# don't plot during testing
if len(sys.argv) != 2:
    import matplotlib.pyplot as plt
    plt.plot(nonlocaleqstrain)
    plt.show()
