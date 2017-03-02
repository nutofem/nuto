#!/usr/bin/env python3
import unittest
import nuto
import numpy as np
import matplotlib.pyplot as plt

L = 20.
A = 12.
E = 6174.
F = 42.
rho = 0.0001

u = F*L/(E*A)


class TestResultForceDirection(unittest.TestCase):

    def DefineTestStructure1D(self, s):
        s.SetShowTime(False)
        s.SetNumTimeDerivatives(2)
        meshInfo = nuto.MeshGenerator.Grid(s, [L], [1])
        s.InterpolationTypeAdd(meshInfo[1], "Displacements", "Equidistant2")
    
        sectionId = s.SectionCreate("truss")
        s.SectionSetArea(sectionId, A)
        lawId = s.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress")
        s.ConstitutiveLawSetParameterDouble(lawId, "Youngs_Modulus", E)
        s.ConstitutiveLawSetParameterDouble(lawId, "Density", rho)
    
    
        s.ElementTotalConvertToInterpolationType()
        s.ElementTotalSetSection(sectionId)
        s.ElementTotalSetConstitutiveLaw(lawId)
    
        nodeLeft = s.NodeGetIdAtCoordinate(np.array([0.0]), 1.e-10)
        s.ConstraintLinearSetDisplacementNode(nodeLeft, np.array([1.0]), 0)
        nodeRight = s.NodeGetIdAtCoordinate(np.array([L]), 1.e-10)
        nodeRightGroup = s.GroupCreate("nodes")
        s.GroupAddNode(nodeRightGroup, nodeRight)
    
        constraintId = s.ConstraintLinearSetDisplacementNode(nodeRight, np.array([1.0]), 0)
    
        return nodeRightGroup, constraintId
    
    def RunTest(self, s, TI, name):
        BCGroup, constraintId = self.DefineTestStructure1D(s)
        
        TI.AddTimeDependentConstraint(constraintId, np.array([[0.,0.], [1.,u]]))
        TI.SetResultDirectory(".", False)
        TI.AddResultGroupNodeForce("Force" +name, BCGroup)
        
        TI.Solve(1.)
        
        force = np.loadtxt("Force" + name + ".dat")[-1]
        self.assertAlmostEqual(F, force, places = 3)
    
    def test_DisplacementControlNewmark(self):
        s = nuto.Structure(1)
        TI = nuto.NewmarkDirect(s)
        TI.SetTimeStep(1.)
        self.RunTest(s, TI, "LoadNewmark");

    def test_DisplacementControlRK4(self):
        s = nuto.Structure(1)
        TI = nuto.RungeKutta4(s)
        TI.SetTimeStep(0.001)
        self.RunTest(s, TI, "LoadRK4")

if __name__ == '__main__':
    unittest.main()
