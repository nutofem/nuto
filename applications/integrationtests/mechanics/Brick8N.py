#!/usr/bin/env python3
import os
import unittest
import numpy as np
import nuto

# path in the original source directory and current filename at the end
path = os.path.dirname(os.path.abspath(__file__))
pathToResultFiles = os.path.join(path, "results")


class Brick8NTestCase(unittest.TestCase):
    def setUp(self):
        # create structure
        self.structure = nuto.Structure(3)

        # create nodes
        node1 = self.structure.NodeCreate(np.array([-1., -1., -1.]))
        node2 = self.structure.NodeCreate(np.array([+1., -1., -1.]))
        node3 = self.structure.NodeCreate(np.array([+1., +1., -1.]))
        node4 = self.structure.NodeCreate(np.array([-1., +1., -1.]))
        node5 = self.structure.NodeCreate(np.array([-1., -1., +1.]))
        node6 = self.structure.NodeCreate(np.array([+1., -1., +1.]))
        node7 = self.structure.NodeCreate(np.array([+1., +1., +1.]))
        node8 = self.structure.NodeCreate(np.array([-1., +1., +1.]))

        # create interpolation type
        interpolationType = self.structure.InterpolationTypeCreate("Brick3D")
        self.structure.InterpolationTypeAdd(interpolationType, "coordinates", "equidistant1")
        self.structure.InterpolationTypeAdd(interpolationType, "displacements", "equidistant1")

        # create element
        nodeIds = [node1, node2, node3, node4, node5, node6, node7, node8]
        self.element = self.structure.ElementCreate(interpolationType, nodeIds)
        self.structure.ElementTotalConvertToInterpolationType()

        # create constitutive law
        material = self.structure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress")
        self.structure.ConstitutiveLawSetParameterDouble(material, "Youngs_Modulus", 10)
        self.structure.ConstitutiveLawSetParameterDouble(material, "Poissons_Ratio", 0.25)

        # create section
        section = self.structure.SectionCreate("Volume")

        # assign constitutive law
        self.structure.ElementSetConstitutiveLaw(self.element, material)
        self.structure.ElementSetSection(self.element, section)

        # make group of boundary nodes
        groupBoundaryNodes = self.structure.GroupCreate("Nodes")
        self.structure.GroupAddNode(groupBoundaryNodes, node1)
        self.structure.GroupAddNode(groupBoundaryNodes, node4)
        self.structure.GroupAddNode(groupBoundaryNodes, node5)
        self.structure.GroupAddNode(groupBoundaryNodes, node8)

        # make group of boundary elements (in this case it is just one
        groupBoundaryElements = self.structure.GroupCreate("Elements")
        self.structure.GroupAddElementsFromNodes(groupBoundaryElements, groupBoundaryNodes, False)

        # create surface loads (0 - pressure on X, 1-const-direction Y)
        self.structure.LoadSurfacePressureCreate3D(0, groupBoundaryElements, groupBoundaryNodes, 2.)
        self.structure.LoadSurfaceConstDirectionCreate3D(1, groupBoundaryElements, groupBoundaryNodes, np.array((0.0, 5.0, 0.0)))

        # set displacements of right node
        self.structure.NodeSetDisplacements(node2, np.array([0.2, 0.2, 0.2]))
        self.structure.NodeSetDisplacements(node3, np.array([0.2, 0.2, 0.2]))
        self.structure.NodeSetDisplacements(node6, np.array([0.2, 0.2, 0.2]))
        self.structure.NodeSetDisplacements(node7, np.array([0.2, 0.2, 0.2]))

    def testStiffness(self):
        # calculate element stiffness matrix
        Ke = self.structure.ElementBuildHessian0(self.element).Get("Displacements", "Displacements")

        # correct stiffness matrix
        resultFile = os.path.join(pathToResultFiles, "Brick8NStiffness.txt")
        KeCorrect = np.loadtxt(resultFile, skiprows=1)
        self.assertTrue(np.allclose(Ke, KeCorrect))
        self.assertTrue(self.structure.CheckHessian0(1e-7, 1e-8))

    def testInternalGradient(self):
        # calculate internal force vector (this is only due to the prescribed
        # displacements, not in equilibrium with external forces
        Fi = self.structure.ElementBuildInternalGradient(self.element).Get("Displacements")
        Fi = Fi.squeeze()

        # correct resforce vector
        resultFile = os.path.join(pathToResultFiles, "Brick8NInternalforce.txt")
        FiCorrect = np.loadtxt(resultFile, skiprows=1)
        self.assertTrue(np.allclose(Fi, FiCorrect))

    def testStrain(self):
        # calculate engineering strain of myelement1 at all integration points
        # the size the matrix is not important and reallocated within the procedure
        engineeringStrain = self.structure.ElementGetEngineeringStrain(self.element)

        # correct strain
        engineeringStrainCorrect = np.array([
            [0.1, 0, 0, 0, 0.1, 0.1],
            [0.1, 0, 0, 0, 0.1, 0.1],
            [0.1, 0, 0, 0, 0.1, 0.1],
            [0.1, 0, 0, 0, 0.1, 0.1],
            [0.1, 0, 0, 0, 0.1, 0.1],
            [0.1, 0, 0, 0, 0.1, 0.1],
            [0.1, 0, 0, 0, 0.1, 0.1],
            [0.1, 0, 0, 0, 0.1, 0.1]]).transpose()

        self.assertTrue(np.allclose(engineeringStrain, engineeringStrainCorrect))

    def testStress(self):
        # calculate engineering strain of myelement1 at all integration points
        engineeringStress = self.structure.ElementGetEngineeringStress(self.element)
        # correct stress
        engineeringStressCorrect = np.array([
            [1.2, 0.4, 0.4, 0.0, 0.4, 0.4],
            [1.2, 0.4, 0.4, 0.0, 0.4, 0.4],
            [1.2, 0.4, 0.4, 0.0, 0.4, 0.4],
            [1.2, 0.4, 0.4, 0.0, 0.4, 0.4],
            [1.2, 0.4, 0.4, 0.0, 0.4, 0.4],
            [1.2, 0.4, 0.4, 0.0, 0.4, 0.4],
            [1.2, 0.4, 0.4, 0.0, 0.4, 0.4],
            [1.2, 0.4, 0.4, 0.0, 0.4, 0.4]]).transpose()

        self.assertTrue(np.allclose(engineeringStress, engineeringStressCorrect))

    def testExternalForce(self):
        # calculate external force vector for the first load case (pressure)
        Fe = self.structure.BuildGlobalExternalLoadVector(0)

        # correct external force for pressure load vector (sum up the load in
        # x direction; everything else should be zero)
        sumX = np.sum(Fe.J.Get("Displacements"))
        self.assertAlmostEqual(sumX, 8.0)

        # calculate external force vector for the second load cases (constDirection)
        Fe = self.structure.BuildGlobalExternalLoadVector(1)

        # correct external force for pressure load vector (sum up the load in
        # y direction; everything else should be zero)
        sumY = np.sum(Fe.J.Get("Displacements"))
        self.assertAlmostEqual(sumY, 20.0)


if __name__ == '__main__':
    unittest.main()
