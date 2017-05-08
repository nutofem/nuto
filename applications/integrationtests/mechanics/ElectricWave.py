#!/usr/bin/env python3
import unittest
import nuto
import sys
import numpy as np

class Wave1D(unittest.TestCase):
    def setUp(self):
        self.s = nuto.Structure(1)
        self.s.SetNumTimeDerivatives(2)
        self.s.SetShowTime(False)
        self.s.SetVerboseLevel(7) # sounds about right 

        mesh_info = nuto.MeshGenerator_Grid(self.s, [2.], [25])
        self.s.InterpolationTypeAdd(mesh_info[1], "ElectricPotential", "Lobatto4")
        self.s.ElementTotalConvertToInterpolationType()

        law_id = self.s.ConstitutiveLawCreate("Linear_Dielectric")
        dielectric_tensor = np.array([1., 0., 0., 0., 1., 0., 0., 0., 1.])
        self.s.ConstitutiveLawSetParameterFullVectorDouble(law_id, "Dielectric_Tensor", -dielectric_tensor)
        self.s.ElementTotalSetConstitutiveLaw(law_id)
        self.s.ElementTotalSetSection(nuto.SectionTruss_Create(12.))
        self.s.AddVisualizationComponent(self.s.GroupGetElementsTotal(), "ElectricPotential")

        n = self.s.NodeGetAtCoordinate(0)
        self.s.Constraints().Add(nuto.eDof_ELECTRICPOTENTIAL, nuto.Value(n, self.load)) 

        self.rk4 = nuto.RungeKutta4(self.s)
        self.rk4.SetTimeStep(0.01)
        self.rk4.SetShowTime(False)

    def load(self, time):
        return np.sin(time * 2 * np.pi )

    def test_time_step(self):
        self.assertLess(self.rk4.GetTimeStep(), self.rk4.CalculateCriticalTimeStep())

    def test_check_solution(self):
        self.rk4.SetResultDirectory("./ElectricWave/", True)
        self.rk4.Solve(2 - 1.e-6) # without 1.e-6, it would end at t = 2.0049999
        nodes = self.s.GroupCreate("Nodes")
        self.s.GroupAddNodeCoordinateRange(nodes, 0, 0, 1.)
        for node_id in self.s.GroupGetMemberIds(nodes):
            node = self.s.NodeGetNodePtr(node_id)
            x = node.Get(nuto.eDof_COORDINATES)[0][0]
            potential = node.Get(nuto.eDof_ELECTRICPOTENTIAL)[0][0]
            analytic = -np.sin(2 * np.pi * x)
            self.assertAlmostEqual(potential, analytic, 2) # correct up to 3 decimal places

if __name__ == '__main__':
    unittest.main()
