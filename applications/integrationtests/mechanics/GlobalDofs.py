#!/usr/bin/env python3
import unittest
import numpy as np
import nuto


class GlobalDofTestCase(unittest.TestCase):
    def setUp(self):
        # create structure
        self.structure = nuto.Structure(1)

        # create nodes
        node1 = self.structure.NodeCreate(np.array([0.]))
        node2 = self.structure.NodeCreate(np.array([5.]))
        node3 = self.structure.NodeCreate(np.array([10.]))

        interpolationType = self.structure.InterpolationTypeCreate("Truss1D")
        self.structure.InterpolationTypeAdd(interpolationType, "coordinates", "equidistant2")
        self.structure.InterpolationTypeAdd(interpolationType, "displacements", "equidistant2")

        self.structure.ElementCreate(interpolationType, [node1, node2, node3])
        self.structure.ElementTotalConvertToInterpolationType()

        # create group of nodes
        nodeGroup = self.structure.GroupCreate("Nodes")
        self.structure.GroupAddNode(nodeGroup, node1)
        self.structure.GroupAddNode(nodeGroup, node3)

        # add constraints for a single node
        self.structure.ConstraintLinearSetDisplacementNode(node2, np.array([1.0, 1.0, -1.0]), 0.5)

        # add constraints for a group of nodes
        self.structure.ConstraintLinearSetDisplacementNodeGroup(nodeGroup, np.array([1.0, 0.0, 0.0]), 2.0)

        self.structure.NodeBuildGlobalDofs()

    def testNumConstraints(self):
        numConstraints = self.structure.ConstraintGetNumLinearConstraints("Displacements")
        self.assertEqual(numConstraints, 3)

    def testGlobalDofs(self):
        numberGlobalDofs = self.structure.GetNumDofs("Displacements")
        self.assertEqual(numberGlobalDofs, 3)

    def testConstraints(self):
        # build constraint matrix and rhs
        rhs = self.structure.ConstraintGetRHSBeforeGaussElimination().Export()
        rhs = rhs.squeeze()
        constraintMatrixFull = self.structure.ConstraintGetConstraintMatrixBeforeGaussElimination().ExportToFullMatrix()

        # correct constraint matrix
        constraintMatrixFullCorrect = np.eye(3)

        # correct rhs
        rhsCorrect = np.r_[0.5, 2.0, 2.0]

        self.assertTrue(np.allclose(constraintMatrixFull, constraintMatrixFullCorrect))
        self.assertTrue(np.allclose(rhs, rhsCorrect))


if __name__ == '__main__':
    unittest.main()
