#!/usr/bin/env python3
import unittest
import nuto
import numpy as np

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
    
        section = nuto.SectionTruss.Create(A)
        lawId = s.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress")
        s.ConstitutiveLawSetParameterDouble(lawId, "Youngs_Modulus", E)
        s.ConstitutiveLawSetParameterDouble(lawId, "Density", rho)
    
    
        s.ElementTotalConvertToInterpolationType()
        s.ElementTotalSetSection(section)
        s.ElementTotalSetConstitutiveLaw(lawId)
    
        nodeLeft = s.NodeGetAtCoordinate(np.array([0.0]), 1.e-10)
        s.Constraints().Add(nuto.eDof_DISPLACEMENTS, nuto.Value(nodeLeft))
        nodeRight = s.NodeGetIdAtCoordinate(np.array([L]), 1.e-10)
        nodeRightGroup = s.GroupCreate("nodes")
        s.GroupAddNode(nodeRightGroup, nodeRight)
    
        nodeRight = s.NodeGetNodePtr(nodeRight)
        s.Constraints().Add(nuto.eDof_DISPLACEMENTS, nuto.Value(nodeRight, lambda t: u*t))
    
        return nodeRightGroup
    
    def RunTest(self, s, TI, name):
        BCGroup = self.DefineTestStructure1D(s)
        
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
