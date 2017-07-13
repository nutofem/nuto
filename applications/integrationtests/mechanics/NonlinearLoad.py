#!/usr/bin/env python3
import unittest
import numpy as np
import nuto


def SetConstitutiveLaw(structure):
    damage = structure.ConstitutiveLawCreate("Gradient_Damage_Engineering_Stress")
    structure.ConstitutiveLawSetParameterDouble(damage, "Youngs_Modulus", 25e3)
    structure.ConstitutiveLawSetParameterDouble(damage, "Poissons_Ratio", 0.2)
    structure.ConstitutiveLawSetParameterDouble(damage, "Nonlocal_Radius", 1.0)
    structure.ConstitutiveLawSetParameterDouble(damage, "Tensile_Strength", 4.0)
    structure.ConstitutiveLawSetParameterDouble(damage, "Compressive_Strength", 40.0)
    structure.ConstitutiveLawSetDamageLaw(damage, nuto.DamageLawExponential.Create(4.0/30e3, 4.0/0.021))
    structure.ElementTotalSetConstitutiveLaw(damage)


def SetBCs(structure):
    nodesBottom = structure.GroupGetNodesAtCoordinate(nuto.eDirection_Z, 0)
    nodesTop = structure.GroupGetNodesAtCoordinate(nuto.eDirection_Z, 1)

    structure.Constraints().Add(nuto.eDof_DISPLACEMENTS, nuto.Component(nodesBottom, [nuto.eDirection_Z]))


    nodeZero = structure.NodeGetAtCoordinate(np.r_[0.0, 0.0, 0.0])
    structure.Constraints().Add(nuto.eDof_DISPLACEMENTS, nuto.Component(nodeZero, [nuto.eDirection_X, nuto.eDirection_Y]))

    nodeOne = structure.NodeGetAtCoordinate(np.r_[1.0, 0.0, 0.0])
    structure.Constraints().Add(nuto.eDof_DISPLACEMENTS, nuto.Component(nodeOne, [nuto.eDirection_Y]))


    topNodes = nodesTop.GetMemberIds()
    topNodes = list(topNodes)
    primary = structure.NodeGetNodePtr(topNodes.pop(0))
    for secondary in topNodes:
        secondaryNode = structure.NodeGetNodePtr(secondary) 
        e = nuto.Equation()
        e.AddTerm(nuto.Term(primary, 2, 1))
        e.AddTerm(nuto.Term(secondaryNode, 2, -1))
        structure.Constraints().Add(nuto.eDof_DISPLACEMENTS, e);

    elementsTop = structure.GroupCreate("Elements")
    nodesTopId = structure.GroupGetId(nodesTop)
    structure.GroupAddElementsFromNodes(elementsTop, nodesTopId, False)
    structure.LoadSurfacePressureCreate3D(elementsTop, nodesTopId, 4.0)


class NonlinearLoadTestCase(unittest.TestCase):
    def setUp(self):
        self.structure = nuto.Structure(3)
        self.structure.SetNumTimeDerivatives(1)

        _, interpolation = nuto.MeshGenerator.Grid(self.structure, [1.0, 1.0, 1.0], [1, 1, 1])

        SetConstitutiveLaw(self.structure)

        self.structure.InterpolationTypeAdd(interpolation, "Displacements", "Equidistant2")  
        self.structure.InterpolationTypeAdd(interpolation, "Nonlocaleqstrain", "Equidistant1")
        self.structure.ElementTotalConvertToInterpolationType()

        SetBCs(self.structure)

        force_application = np.array([[0, 0], [1.0, 1.0]])

        newmark = nuto.NewmarkDirect(self.structure)
        newmark.SetTimeDependentLoadCase(0, force_application)
        newmark.SetPerformLineSearch(True)
        newmark.SetTimeStep(1.0)
        newmark.SetResultDirectory("./NonlinearLoadOut", True)
        newmark.Solve(1.0)

    def testStress(self):
        stress = self.structure.ElementGetEngineeringStress(0)
        self.assertTrue(np.allclose(stress[2, :], -4.0))


if __name__ == '__main__':
    unittest.main()
