#!/usr/bin/env python3
import unittest
import nuto
import numpy as np

def GetInterpolationType(interpolationOrder):
    return "Equidistant" + str(interpolationOrder)

def GetIntegrationType(integrationOrder):
    return "INTEGRATIONTYPE1D2NGAUSS" + str(integrationOrder) + "IP"


class MassMatrix:
    def __init__(self, interpolationOrder, integrationOrder):
        self.lx = 42.
        self.ly = 13.
        self.lz = 73.
        self.rho = 0.6174

        self.s = nuto.Structure(1)
        self.s.SetShowTime(False)
        meshInfo = nuto.MeshGenerator.Grid(self.s, [self.lx], [1])
        self.s.InterpolationTypeAdd(meshInfo[1], "Displacements", GetInterpolationType(interpolationOrder)) 
        self.s.InterpolationTypeSetIntegrationType(meshInfo[1], GetIntegrationType(integrationOrder))

        sectionId = self.s.SectionCreate("Truss")
        self.s.SectionSetArea(sectionId, self.ly*self.lz)

        lawId = self.s.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
        self.s.ConstitutiveLawSetParameterDouble(lawId, "Youngs_Modulus", 41)
        self.s.ConstitutiveLawSetParameterDouble(lawId, "Density", self.rho)

        self.s.ElementTotalConvertToInterpolationType()
        self.s.ElementTotalSetSection(sectionId)
        self.s.ElementTotalSetConstitutiveLaw(lawId)

    def GetMassMatrix(self):
        elementGroup = self.s.GroupGetElementsTotal()
        elementId = self.s.GroupGetMemberIds(elementGroup)[0]
        return self.s.ElementBuildHessian2(elementId).Export()


class MassMatrixTest(unittest.TestCase):
    correctMass = 42. * 13. * 73. * 0.6174

    def CheckUnderIntegration(self, interpolationOrder):
        integrationOrderCorrect = interpolationOrder + 1
        integrationOrderUnder = interpolationOrder

        massCorrect = MassMatrix(interpolationOrder, integrationOrderCorrect).GetMassMatrix()
        massUnder= MassMatrix(interpolationOrder, integrationOrderUnder).GetMassMatrix()
        massError = np.abs(massCorrect - massUnder)

        self.assertAlmostEqual(massCorrect.sum(), self.correctMass)
        self.assertAlmostEqual(massUnder.sum(), self.correctMass)
        self.assertNotAlmostEqual(massError.sum(), 0.)
        print("Absolute total mass error order ", interpolationOrder, ": ", massError.sum())


    def test_MassUnderIntegration1(self):
        self.CheckUnderIntegration(1)

    def test_MassUnderIntegration2(self):
        self.CheckUnderIntegration(2)

    def test_MassUnderIntegration3(self):
        self.CheckUnderIntegration(3)

    def test_MassUnderIntegration4(self):
        self.CheckUnderIntegration(4)

if __name__ == '__main__':
    unittest.main()
