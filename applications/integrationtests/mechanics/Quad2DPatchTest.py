#!/usr/bin/env python3
from ctypes import c_int
import unittest
import numpy as np
import nuto

# definitions
E = 80000
v = 1./3.
BoundaryDisplacement = 1.0


def SetupStructure(stressState):
    structure = nuto.Structure(2)
    structure.SetShowTime(False)

    # create material law
    myMatLin = structure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress")
    structure.ConstitutiveLawSetParameterDouble(myMatLin, "Youngs_Modulus", E)
    structure.ConstitutiveLawSetParameterDouble(myMatLin, "Poissons_Ratio", v)

    # create section
    mySection = nuto.SectionPlane.Create(1.0, False)

    structure.NodesCreate(np.array([[0, 10,  2,  8,  4,  8,  0, 10],
                                    [0,  0,  2,  3,  7,  7, 10, 10]], dtype=float))

    elementIncidence = np.array([[3, 4, 5, 7, 7],
                                 [2, 6, 4, 5, 6],
                                 [0, 0, 2, 3, 4],
                                 [1, 2, 3, 1, 5]], dtype=c_int)

    interpolationType = structure.InterpolationTypeCreate("Quad2D")
    structure.InterpolationTypeAdd(interpolationType, "Coordinates", "Equidistant1")
    structure.InterpolationTypeAdd(interpolationType, "Displacements", "Equidistant1")

    structure.ElementsCreate(interpolationType, elementIncidence)
    structure.ElementTotalConvertToInterpolationType()
    structure.ElementTotalSetConstitutiveLaw(myMatLin)
    structure.ElementTotalSetSection(mySection)

    LoadNodesXPos = structure.GroupCreate("Nodes")
    LoadNodesXNeg = structure.GroupCreate("Nodes")
    LoadNodesYPos = structure.GroupCreate("Nodes")
    LoadNodesYNeg = structure.GroupCreate("Nodes")

    structure.GroupAddNode(LoadNodesXPos, 1)
    structure.GroupAddNode(LoadNodesXPos, 7)
    structure.GroupAddNode(LoadNodesYPos, 6)
    structure.GroupAddNode(LoadNodesYPos, 7)

    structure.GroupAddNode(LoadNodesXNeg, 0)
    structure.GroupAddNode(LoadNodesXNeg, 6)
    structure.GroupAddNode(LoadNodesYNeg, 0)
    structure.GroupAddNode(LoadNodesYNeg, 1)

    directionX = np.array([1.0, 0.0])
    directionY = np.array([0.0, 1.0])

    if stressState == "XX":
        structure.Constraints().Add(nuto.eDof_DISPLACEMENTS, nuto.Fix(LoadNodesXNeg, [nuto.eDirection_X]))
        structure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesXPos, directionX, BoundaryDisplacement)
        structure.ConstraintLinearSetDisplacementNode(0, directionY, 0)
    elif stressState == "YY":
        structure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesYNeg, directionY, 0.0)
        structure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesYPos, directionY, BoundaryDisplacement)
        structure.ConstraintLinearSetDisplacementNode(0, directionX, 0)
    elif stressState == "XY":
        structure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesXNeg, directionX, 0.0)
        structure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesXNeg, directionY, 0.0)
        structure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesXPos, directionY, BoundaryDisplacement)
        structure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesXPos, directionX, 0)

    # start analysis
    # build global dof numbering
    structure.NodeBuildGlobalDofs()
    structure.CalculateMaximumIndependentSets()

    structure.SolveGlobalSystemStaticElastic()
    return structure


class PatchTestCase:
    def testInternalGradient(self):
        # calculate residual
        internalGradient = self.structure.BuildGlobalInternalGradient()
        self.assertLess(np.linalg.norm(internalGradient.J.Export()), 1e-8)

    def testStress(self):
        # calculate engineering strain of 2 at all integration points
        EngineeringStress = self.structure.ElementGetEngineeringStress(2)
        # correct stress
        self.assertTrue(np.max(np.abs(EngineeringStress - self.getCorrectStress())) < 1e-4)

    def testStrain(self):
        # calculate engineering strain of 2 at all integration points
        # the size the matrix is not important and reallocated within the procedure
        EngineeringStrain = self.structure.ElementGetEngineeringStrain(2)

        C1 = 1./E
        C2 = -v/E
        C3 = 2*(1+v)/E

        D = np.array([
                [C1, C2, C2,  0,  0,  0],
                [C2, C1, C2,  0,  0,  0],
                [C2, C2, C1,  0,  0,  0],
                [0,  0,  0,  C3,  0,  0],
                [0,  0,  0,   0, C3,  0],
                [0,  0,  0,   0,  0, C3]])

        # correct strain
        EngineeringStrainCorrect = D.dot(self.getCorrectStress())
        self.assertTrue(np.max(np.abs(EngineeringStrain - EngineeringStrainCorrect)) < 1e-8)


class PatchTestCaseXX(PatchTestCase, unittest.TestCase):
    def setUp(self):
        self.structure = SetupStructure("XX")

    def getCorrectStress(self):
        EngineeringStressCorrect = np.zeros((6, 4))
        for i in range(0, 4):
            sigma = E*BoundaryDisplacement / 10.
            EngineeringStressCorrect[0, i] = sigma
        return EngineeringStressCorrect


class PatchTestCaseYY(PatchTestCase, unittest.TestCase):
    def setUp(self):
        self.structure = SetupStructure("YY")

    def getCorrectStress(self):
        EngineeringStressCorrect = np.zeros((6, 4))
        for i in range(0, 4):
            sigma = E*BoundaryDisplacement / 10.
            EngineeringStressCorrect[1, i] = sigma
        return EngineeringStressCorrect


class PatchTestCaseXY(PatchTestCase, unittest.TestCase):
    def setUp(self):
        self.structure = SetupStructure("XY")

    def getCorrectStress(self):
        EngineeringStressCorrect = np.zeros((6, 4))
        for i in range(0, 4):
            sigma = E*BoundaryDisplacement / 10. / (2+2*v)
            EngineeringStressCorrect[5, i] = sigma
        return EngineeringStressCorrect


if __name__ == '__main__':
    unittest.main()
