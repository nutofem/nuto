#!/usr/bin/env python3
import os
import unittest
import numpy as np
import nuto

# path in the original source directory and current filename at the end
path = os.path.dirname(os.path.abspath(__file__))
pathToResultFiles = os.path.join(path, "results")


class Tetrahedron10NTestCase(unittest.TestCase):
    def setUp(self):
        # create structure
        self.structure = nuto.Structure(3)

        # create nodes
        node1 = self.structure.NodeCreate(np.array([0.0, 0.0, 0.0]))
        node2 = self.structure.NodeCreate(np.array([1.0, 0.0, 0.0]))
        node3 = self.structure.NodeCreate(np.array([0.0, 1.0, 0.0]))
        node4 = self.structure.NodeCreate(np.array([0.0, 0.0, 1.0]))
        node5 = self.structure.NodeCreate(np.array([0.5, 0.0, 0.0]))
        node6 = self.structure.NodeCreate(np.array([0.5, 0.5, 0.0]))
        node7 = self.structure.NodeCreate(np.array([0.0, 0.5, 0.0]))
        node8 = self.structure.NodeCreate(np.array([0.0, 0.0, 0.5]))
        node9 = self.structure.NodeCreate(np.array([0.0, 0.5, 0.5]))
        node10 = self.structure.NodeCreate(np.array([0.5, 0.0, 0.5]))

        interpolationType = self.structure.InterpolationTypeCreate("TETRAHEDRON3D")
        self.structure.InterpolationTypeAdd(interpolationType, "Coordinates", "EQUIDISTANT2")
        self.structure.InterpolationTypeAdd(interpolationType, "Displacements", "EQUIDISTANT2")

        # create element
        self.element = self.structure.ElementCreate(interpolationType,
                [node1, node2, node3, node4, node5, node6, node7, node8, node9, node10])
        self.structure.ElementTotalConvertToInterpolationType()

        # create constitutive law
        material = self.structure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress")
        self.structure.ConstitutiveLawSetParameterDouble(material, "Youngs_Modulus", 10)
        self.structure.ConstitutiveLawSetParameterDouble(material, "Poissons_Ratio", 0.25)

        # assign constitutive law
        self.structure.ElementSetConstitutiveLaw(self.element, material)

        # make group of boundary nodes
        groupBoundaryNodes = self.structure.GroupCreate("Nodes")
        self.structure.GroupAddNode(groupBoundaryNodes, node1)
        self.structure.GroupAddNode(groupBoundaryNodes, node3)
        self.structure.GroupAddNode(groupBoundaryNodes, node4)
        self.structure.GroupAddNode(groupBoundaryNodes, node7)
        self.structure.GroupAddNode(groupBoundaryNodes, node8)
        self.structure.GroupAddNode(groupBoundaryNodes, node9)

        # make group of boundary elements (in this case it is just one)
        groupBoundaryElements = self.structure.GroupCreate("Elements")
        self.structure.GroupAddElementsFromNodes(groupBoundaryElements, groupBoundaryNodes, False)

        # create surface loads (0 - pressure on X, 1-const-direction Y)
        self.structure.LoadSurfacePressureCreate3D(0, groupBoundaryElements, groupBoundaryNodes, 2.)
        self.structure.LoadSurfaceConstDirectionCreate3D(1, groupBoundaryElements, groupBoundaryNodes, np.array([0., 5., 0.]))

        # set displacements of right node
        self.structure.NodeSetDisplacements(node2, np.array([0.2, 0.2, 0.2]))
        self.structure.NodeSetDisplacements(node3, np.array([0.2, 0.2, 0.2]))
        self.structure.NodeSetDisplacements(node6, np.array([0.2, 0.2, 0.2]))
        self.structure.NodeSetDisplacements(node7, np.array([0.2, 0.2, 0.2]))

    def testStiffness(self):
        # calculate element stiffness matrix
        Ke = self.structure.ElementBuildHessian0(self.element).Export()

        # correct stiffness matrix
        resultFile = os.path.join(pathToResultFiles, "Tetrahedron10NStiffness.txt")
        KeCorrect = np.loadtxt(resultFile, skiprows=1)
        self.assertTrue(np.allclose(Ke, KeCorrect))
        self.assertTrue(self.structure.CheckHessian0(1e-7, 1e-8))

    def testInternalGradient(self):
        # calculate internal force vector (this is only due to the prescribed
        # displacements, not in equilibrium with external forces
        Fi = self.structure.ElementBuildInternalGradient(self.element).Get("Displacements")
        Fi = Fi.squeeze()

        # correct resforce vector
        resultFile = os.path.join(pathToResultFiles, "Tetrahedron10NInternalforce.txt")
        FiCorrect = np.loadtxt(resultFile, skiprows=1)
        self.assertTrue(np.allclose(Fi, FiCorrect))

    def testStrain(self):
        # calculate engineering strain of myelement1 at all integration points
        # the size the matrix is not important and reallocated within the procedure
        engineeringStrain = self.structure.ElementGetEngineeringStrain(self.element)

        # correct strain
        engineeringStrainCorrect = np.array([
            [-0.08944272, 0.37888544, -0.11055728,  0.26832816, -0.20000000,  0.28944272],
            [ 0.26832816, 0.37888544, -0.11055728,  0.26832816,  0.15777088,  0.64721360],
            [-0.08944272, 0.02111456, -0.46832816, -0.44721360, -0.55777088, -0.06832816],
            [-0.08944272, 0.02111456, -0.11055728, -0.08944272, -0.20000000, -0.06832816]
            ]).transpose()

        self.assertTrue(np.allclose(engineeringStrain, engineeringStrainCorrect))

    def testStress(self):
        # calculate engineering strain of myelement1 at all integration points
        engineeringStress = self.structure.ElementGetEngineeringStress(self.element)

        # correct stress
        engineeringStressCorrect = np.array([
            [ 0.00000000,  3.74662528, -0.16891648,  1.07331264, -0.80000000,  1.15777088],
            [ 4.29325056,  5.17770880,  1.26216704,  1.07331264,  0.63108352,  2.58885440],
            [-2.86216704, -1.97770880, -5.89325056, -1.78885440, -2.23108352, -0.27331264],
            [-1.43108352, -0.54662528, -1.60000000, -0.35777088, -0.80000000, -0.27331264]
            ]).transpose()

        self.assertTrue(np.allclose(engineeringStress, engineeringStressCorrect))

    def testExternalForce(self):
        # calculate external force vector for the first load case (pressure)
        Fe = self.structure.BuildGlobalExternalLoadVector(0)
        sumX = np.sum(Fe.J.Get("Displacements"))
        self.assertAlmostEqual(sumX, 1.0)

        # calculate external force vector for the second load cases (constDirection)
        Fe = self.structure.BuildGlobalExternalLoadVector(1)

        # correct external force for pressure load vector
        # (sum up the load in x direction eveything else should be zero)
        sumY = np.sum(Fe.J.Get("Displacements"))
        self.assertAlmostEqual(sumY, 2.5)


if __name__ == '__main__':
    unittest.main()
