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
    nodesBottom = structure.GroupCreate("Nodes")
    nodesTop = structure.GroupCreate("Nodes")
    structure.GroupAddNodeCoordinateRange(nodesBottom, 2, -1e-6, 1e-6)
    structure.GroupAddNodeCoordinateRange(nodesTop, 2, 1.0 - 1e-6, 1.0 + 1e-6)

    elementsTop = structure.GroupCreate("Elements")
    structure.GroupAddElementsFromNodes(elementsTop, nodesTop, False)

    structure.ConstraintLinearSetDisplacementNodeGroup(nodesBottom, np.r_[0.0, 0.0, 1.0], 0.0)

    nodeZero = structure.NodeGetIdAtCoordinate(np.r_[0.0, 0.0, 0.0], 1e-5)
    structure.ConstraintLinearSetDisplacementNode(nodeZero, np.r_[1.0, 0.0, 0.0], 0.0)
    structure.ConstraintLinearSetDisplacementNode(nodeZero, np.r_[0.0, 1.0, 0.0], 0.0)

    nodeOne = structure.NodeGetIdAtCoordinate(np.r_[1.0, 0.0, 0.0], 1e-8)
    structure.ConstraintLinearSetDisplacementNode(nodeOne, np.r_[0.0, 1.0, 0.0], 0.0)

    topNodes = structure.GroupGetMemberIds(nodesTop)
    topNodes = list(topNodes)
    primary = topNodes.pop(0)
    for secondary in topNodes:
        constraint = structure.ConstraintLinearEquationCreate(primary, "Z_DISPLACEMENT", 1.0, 0.0)
        structure.ConstraintLinearEquationAddTerm(constraint, secondary, "Z_DISPLACEMENT", -1.0)

    structure.LoadSurfacePressureCreate3D(elementsTop, nodesTop, 4.0)


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
        newmark.SetResultDirectory("blubb")
        newmark.Solve(1.0)

    def testStress(self):
        stress = self.structure.ElementGetEngineeringStress(0)
        self.assertTrue(np.allclose(stress[2, :], -4.0))


if __name__ == '__main__':
    unittest.main()
