#!/usr/bin/env python3
import unittest
import nuto


def OneDimensional(structure):
    _, ipType = nuto.MeshGenerator.Grid(structure, [2.0], [1])
    structure.InterpolationTypeAdd(ipType, "temperature", "equidistant1")

    # set section
    area = 1.0
    section = nuto.SectionTruss.Create(area)
    structure.ElementTotalSetSection(section)

    return structure


def TwoDimensional(structure):
    _, ipType = nuto.MeshGenerator.Grid(structure, [2.0, 2.0], [1, 1])
    structure.InterpolationTypeAdd(ipType, "temperature", "equidistant1")

    # create section
    thickness = 1.0
    section = nuto.SectionPlane.Create(thickness, True)
    structure.ElementTotalSetSection(section)

    return structure


def ThreeDimensional(structure):
    _, ipType = nuto.MeshGenerator.Grid(structure, [2.0, 2.0, 2.0], [1, 1, 1])
    structure.InterpolationTypeAdd(ipType, "temperature", "equidistant1")

    return structure


def SetupStructure(dimension):
    assert(dimension in [1, 2, 3])
    structure = nuto.Structure(dimension)

    if dimension == 1:
        structure = OneDimensional(structure)
    elif dimension == 2:
        structure = TwoDimensional(structure)
    else:
        structure = ThreeDimensional(structure)

    # create material law
    conductivity = 1.0
    material = structure.ConstitutiveLawCreate("Heat_Conduction")
    structure.ConstitutiveLawSetParameterDouble(material, "Thermal_Conductivity", conductivity)

    structure.ElementTotalConvertToInterpolationType()
    structure.ElementTotalSetConstitutiveLaw(material)

    return structure


class HeatConductionXd:
    def testHessian(self):
        delta = 1e-16
        rel_tolerance = 1e-15
        self.assertTrue(self.structure.CheckHessian0(delta, rel_tolerance, True))


class HeatConduction1D(HeatConductionXd, unittest.TestCase):
    def setUp(self):
        self.structure = SetupStructure(1)


class HeatConduction2D(HeatConductionXd, unittest.TestCase):
    def setUp(self):
        self.structure = SetupStructure(2)


class HeatConduction3D(HeatConductionXd, unittest.TestCase):
    def setUp(self):
        self.structure = SetupStructure(3)


if __name__ == '__main__':
    unittest.main()
