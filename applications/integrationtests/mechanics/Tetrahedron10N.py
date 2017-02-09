# -*- coding: utf-8 -*-
import nuto
import os
import numpy as np
import unittest

# path in the original source directory and current filename at the end
path = os.path.dirname(os.path.abspath(__file__))
pathToResultFiles = os.path.join(path, "results")


class Tetrahedron10NTestCase(unittest.TestCase):
    def setUp(self):
        # create structure
        self.structure = nuto.Structure(3)

        # create nodes
        myNode1 = self.structure.NodeCreate(np.array([0.0, 0.0, 0.0]))
        myNode2 = self.structure.NodeCreate(np.array([1.0, 0.0, 0.0]))
        myNode3 = self.structure.NodeCreate(np.array([0.0, 1.0, 0.0]))
        myNode4 = self.structure.NodeCreate(np.array([0.0, 0.0, 1.0]))
        myNode5 = self.structure.NodeCreate(np.array([0.5, 0.0, 0.0]))
        myNode6 = self.structure.NodeCreate(np.array([0.5, 0.5, 0.0]))
        myNode7 = self.structure.NodeCreate(np.array([0.0, 0.5, 0.0]))
        myNode8 = self.structure.NodeCreate(np.array([0.0, 0.0, 0.5]))
        myNode9 = self.structure.NodeCreate(np.array([0.0, 0.5, 0.5]))
        myNode10 = self.structure.NodeCreate(np.array([0.5, 0.0, 0.5]))

        # create section
        mySection = self.structure.SectionCreate("Volume")

        myInterpolationType = self.structure.InterpolationTypeCreate("TETRAHEDRON3D")
        self.structure.InterpolationTypeAdd(myInterpolationType, "Coordinates", "EQUIDISTANT2")
        self.structure.InterpolationTypeAdd(myInterpolationType, "Displacements", "EQUIDISTANT2")

        # create element
        self.element = self.structure.ElementCreate(myInterpolationType,
                [myNode1, myNode2, myNode3, myNode4, myNode5, myNode6, myNode7, myNode8, myNode9, myNode10])
        self.structure.ElementTotalConvertToInterpolationType()

        # create constitutive law
        myMatLin = self.structure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress")
        self.structure.ConstitutiveLawSetParameterDouble(myMatLin, "Youngs_Modulus", 10)
        self.structure.ConstitutiveLawSetParameterDouble(myMatLin, "Poissons_Ratio", 0.25)

        # assign constitutive law
        self.structure.ElementSetConstitutiveLaw(self.element, myMatLin)
        self.structure.ElementSetSection(self.element, mySection)

        # make group of boundary nodes
        groupBoundaryNodes = self.structure.GroupCreate("Nodes")
        self.structure.GroupAddNode(groupBoundaryNodes, myNode1)
        self.structure.GroupAddNode(groupBoundaryNodes, myNode3)
        self.structure.GroupAddNode(groupBoundaryNodes, myNode4)
        self.structure.GroupAddNode(groupBoundaryNodes, myNode7)
        self.structure.GroupAddNode(groupBoundaryNodes, myNode8)
        self.structure.GroupAddNode(groupBoundaryNodes, myNode9)

        # make group of boundary elements (in this case it is just one)
        groupBoundaryElements = self.structure.GroupCreate("Elements")
        self.structure.GroupAddElementsFromNodes(groupBoundaryElements, groupBoundaryNodes, False)

        # create surface loads (0 - pressure on X, 1-const-direction Y)
        self.structure.LoadSurfacePressureCreate3D(0, groupBoundaryElements, groupBoundaryNodes, 2.)
        self.structure.LoadSurfaceConstDirectionCreate3D(1, groupBoundaryElements, groupBoundaryNodes, np.array([0., 5., 0.]))

        # set displacements of right node
        self.structure.NodeSetDisplacements(myNode2, np.array([0.2, 0.2, 0.2]))
        self.structure.NodeSetDisplacements(myNode3, np.array([0.2, 0.2, 0.2]))
        self.structure.NodeSetDisplacements(myNode6, np.array([0.2, 0.2, 0.2]))
        self.structure.NodeSetDisplacements(myNode7, np.array([0.2, 0.2, 0.2]))


class CheckStiffness(Tetrahedron10NTestCase):
    def runTest(self):
        # calculate element stiffness matrix
        Ke = self.structure.ElementBuildHessian0(self.element).Export()

        # correct stiffness matrix
        resultFile = os.path.join(pathToResultFiles, "Tetrahedron10NStiffness.txt")
        KeCorrect = np.loadtxt(resultFile, skiprows=1)
        self.assertTrue(np.max(np.abs(Ke-KeCorrect)) < 1e-8)
        self.assertTrue(self.structure.CheckHessian0(1e-7, 1e-8))


class CheckInternalGradient(Tetrahedron10NTestCase):
    def runTest(self):
        # calculate internal force vector (this is only due to the prescribed
        # displacements, not in equilibrium with external forces
        Fi = self.structure.ElementBuildInternalGradient(self.element).Get("Displacements")
        Fi = Fi.squeeze()

        # correct resforce vector
        resultFile = os.path.join(pathToResultFiles, "Tetrahedron10NInternalforce.txt")
        FiCorrect = np.loadtxt(resultFile, skiprows=1)
        self.assertTrue(np.max(np.abs(Fi-FiCorrect)) < 1e-8)


class CheckStrain(Tetrahedron10NTestCase):
    def runTest(self):
        # calculate engineering strain of myelement1 at all integration points
        # the size the matrix is not important and reallocated within the procedure
        EngineeringStrain = self.structure.ElementGetEngineeringStrain(self.element)

        # correct strain
        EngineeringStrainCorrect = np.array([
            [-0.08944272, 0.37888544, -0.11055728,  0.26832816, -0.20000000, 0.28944272],
            [ 0.26832816, 0.37888544, -0.11055728,  0.26832816,  0.15777088, 0.64721360],
            [-0.08944272, 0.02111456, -0.46832816, -0.44721360, -0.55777088,-0.06832816],
            [-0.08944272, 0.02111456, -0.11055728, -0.08944272, -0.20000000,-0.06832816]
            ]).transpose()

        self.assertTrue(np.max(np.abs(EngineeringStrain - EngineeringStrainCorrect)) < 1e-5)


class CheckStress(Tetrahedron10NTestCase):
    def runTest(self):
        # calculate engineering strain of myelement1 at all integration points
        EngineeringStress = self.structure.ElementGetEngineeringStress(self.element)

        # correct stress
        EngineeringStressCorrect = np.array([
            [ 0.00000000, 3.74662528 , -0.16891648 , 1.07331264 ,-0.80000000, 1.15777088], 
            [ 4.29325056, 5.17770880 ,  1.26216704 , 1.07331264 , 0.63108352, 2.58885440],
            [-2.86216704, -1.97770880, -5.89325056 ,-1.78885440 ,-2.23108352,-0.27331264],
            [-1.43108352, -0.54662528, -1.60000000 ,-0.35777088 ,-0.80000000,-0.27331264]
            ]).transpose()

        self.assertTrue(np.max(np.abs(EngineeringStress - EngineeringStressCorrect)) < 1e-5)


class CheckExternalForce(Tetrahedron10NTestCase):
    def runTest(self):
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
