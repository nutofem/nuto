import numpy as np
import nuto
import unittest


class GlobalDofTestCase(unittest.TestCase):
    def setUp(self):
        # create structure
        self.structure = nuto.Structure(1)

        # create nodes
        myNode1 = self.structure.NodeCreate(np.array([0.]))
        myNode2 = self.structure.NodeCreate(np.array([5.]))
        myNode3 = self.structure.NodeCreate(np.array([10.]))

        myInterpolationType = self.structure.InterpolationTypeCreate("Truss1D")
        self.structure.InterpolationTypeAdd(myInterpolationType, "coordinates", "equidistant2")
        self.structure.InterpolationTypeAdd(myInterpolationType, "displacements", "equidistant2")

        self.structure.ElementCreate(myInterpolationType, [myNode1, myNode2, myNode3])
        self.structure.ElementTotalConvertToInterpolationType()

        # create group of nodes
        myNodeGroup = self.structure.GroupCreate("Nodes")
        self.structure.GroupAddNode(myNodeGroup, myNode1)
        self.structure.GroupAddNode(myNodeGroup, myNode3)

        # add constraints for a single node
        self.structure.ConstraintLinearSetDisplacementNode(myNode2, np.array([1.0, 1.0, -1.0]), 0.5)

        # add constraints for a group of nodes
        self.structure.ConstraintLinearSetDisplacementNodeGroup(myNodeGroup, np.array([1.0, 0.0, 0.0]), 2.0)

        self.structure.NodeBuildGlobalDofs()


class CheckConstraints(GlobalDofTestCase):
    def runTest(self):
        numConstraints = self.structure.ConstraintGetNumLinearConstraints("Displacements")
        self.assertEqual(numConstraints, 3)


class CheckGlobalDofs(GlobalDofTestCase):
    def runTest(self):
        numberGlobalDofs = self.structure.GetNumDofs("Displacements")
        self.assertEqual(numberGlobalDofs, 3)


class CheckConstraints(GlobalDofTestCase):
    def runTest(self):
        # build constraint matrix and rhs
        rhs = self.structure.ConstraintGetRHSBeforeGaussElimination().Export()
        rhs = rhs.squeeze()
        constraintMatrixFull = self.structure.ConstraintGetConstraintMatrixBeforeGaussElimination().ExportToFullMatrix()

        # correct constraint matrix
        constraintMatrixFullCorrect = np.eye(3)

        # correct rhs
        rhsCorrect = np.r_[0.5, 2.0, 2.0]

        self.assertTrue(np.max(np.abs(constraintMatrixFull - constraintMatrixFullCorrect)) < 1e-8)
        self.assertTrue(np.max(np.abs(rhs - rhsCorrect)) < 1e-8)


if __name__ == '__main__':
    unittest.main()
