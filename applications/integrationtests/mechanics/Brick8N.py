import nuto
import os
import numpy as np
import unittest

# path in the original source directory and current filename at the end
path = os.path.dirname(os.path.abspath(__file__))
pathToResultFiles = os.path.join(path, "results")


class Brick8NTestCase(unittest.TestCase):
    def setUp(self):
        # create structure
        self.structure = nuto.Structure(3)

        # create nodes
        myNode1 = self.structure.NodeCreate(np.array([-1., -1., -1.]))
        myNode2 = self.structure.NodeCreate(np.array([+1., -1., -1.]))
        myNode3 = self.structure.NodeCreate(np.array([+1., +1., -1.]))
        myNode4 = self.structure.NodeCreate(np.array([-1., +1., -1.]))
        myNode5 = self.structure.NodeCreate(np.array([-1., -1., +1.]))
        myNode6 = self.structure.NodeCreate(np.array([+1., -1., +1.]))
        myNode7 = self.structure.NodeCreate(np.array([+1., +1., +1.]))
        myNode8 = self.structure.NodeCreate(np.array([-1., +1., +1.]))

        # create interpolation type
        myInterpolationType = self.structure.InterpolationTypeCreate("Brick3D")
        self.structure.InterpolationTypeAdd(myInterpolationType, "coordinates", "equidistant1")
        self.structure.InterpolationTypeAdd(myInterpolationType, "displacements", "equidistant1")

        # create element
        nodeIds = [myNode1, myNode2, myNode3, myNode4, myNode5, myNode6, myNode7, myNode8]
        self.element = self.structure.ElementCreate(myInterpolationType, nodeIds)
        self.structure.ElementTotalConvertToInterpolationType()

        # create constitutive law
        myMatLin = self.structure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress")
        self.structure.ConstitutiveLawSetParameterDouble(myMatLin, "Youngs_Modulus", 10)
        self.structure.ConstitutiveLawSetParameterDouble(myMatLin, "Poissons_Ratio", 0.25)

        # create section
        mySection = self.structure.SectionCreate("Volume")

        # assign constitutive law
        self.structure.ElementSetConstitutiveLaw(self.element, myMatLin)
        self.structure.ElementSetSection(self.element, mySection)

        # make group of boundary nodes
        groupBoundaryNodes = self.structure.GroupCreate("Nodes")
        self.structure.GroupAddNode(groupBoundaryNodes, myNode1)
        self.structure.GroupAddNode(groupBoundaryNodes, myNode4)
        self.structure.GroupAddNode(groupBoundaryNodes, myNode5)
        self.structure.GroupAddNode(groupBoundaryNodes, myNode8)

        # make group of boundary elements (in this case it is just one
        groupBoundaryElements = self.structure.GroupCreate("Elements")
        self.structure.GroupAddElementsFromNodes(groupBoundaryElements, groupBoundaryNodes, False)

        # create surface loads (0 - pressure on X, 1-const-direction Y)
        self.structure.LoadSurfacePressureCreate3D(0, groupBoundaryElements, groupBoundaryNodes, 2.)
        self.structure.LoadSurfaceConstDirectionCreate3D(1, groupBoundaryElements, groupBoundaryNodes, np.array((0.0, 5.0, 0.0)))

        # set displacements of right node
        self.structure.NodeSetDisplacements(myNode2, np.array([0.2, 0.2, 0.2]))
        self.structure.NodeSetDisplacements(myNode3, np.array([0.2, 0.2, 0.2]))
        self.structure.NodeSetDisplacements(myNode6, np.array([0.2, 0.2, 0.2]))
        self.structure.NodeSetDisplacements(myNode7, np.array([0.2, 0.2, 0.2]))


class CheckStiffness(Brick8NTestCase):
    def runTest(self):
        # calculate element stiffness matrix
        Ke = self.structure.ElementBuildHessian0(self.element).Get("Displacements", "Displacements")

        # correct stiffness matrix
        resultFile = os.path.join(pathToResultFiles, "Brick8NStiffness.txt")
        KeCorrect = np.loadtxt(resultFile, skiprows=1)
        self.assertTrue(np.max(np.abs(Ke-KeCorrect)) < 1e-8)
        self.assertTrue(self.structure.CheckHessian0(1e-7, 1e-8))


class CheckInternalGradient(Brick8NTestCase):
    def runTest(self):
        # calculate internal force vector (this is only due to the prescribed
        # displacements, not in equilibrium with external forces
        Fi = self.structure.ElementBuildInternalGradient(self.element).Get("Displacements")
        Fi = Fi.squeeze()

        # correct resforce vector
        resultFile = os.path.join(pathToResultFiles, "Brick8NInternalforce.txt")
        FiCorrect = np.loadtxt(resultFile, skiprows=1)
        self.assertTrue(np.max(np.abs(Fi-FiCorrect)) < 1e-8)


class CheckStrain(Brick8NTestCase):
    def runTest(self):
        # calculate engineering strain of myelement1 at all integration points
        # the size the matrix is not important and reallocated within the procedure
        EngineeringStrain = self.structure.ElementGetEngineeringStrain(self.element)

        # correct strain
        EngineeringStrainCorrect = np.array([
            [0.1, 0, 0, 0, 0.1, 0.1],
            [0.1, 0, 0, 0, 0.1, 0.1],
            [0.1, 0, 0, 0, 0.1, 0.1],
            [0.1, 0, 0, 0, 0.1, 0.1],
            [0.1, 0, 0, 0, 0.1, 0.1],
            [0.1, 0, 0, 0, 0.1, 0.1],
            [0.1, 0, 0, 0, 0.1, 0.1],
            [0.1, 0, 0, 0, 0.1, 0.1]]).transpose()

        self.assertTrue(np.max(np.abs(EngineeringStrain-EngineeringStrainCorrect)) < 1e-8)


class CheckStress(Brick8NTestCase):
    def runTest(self):
        # calculate engineering strain of myelement1 at all integration points
        EngineeringStress = self.structure.ElementGetEngineeringStress(self.element)
        # correct stress
        EngineeringStressCorrect = np.array([
            [1.2, 0.4, 0.4, 0.0, 0.4, 0.4],
            [1.2, 0.4, 0.4, 0.0, 0.4, 0.4],
            [1.2, 0.4, 0.4, 0.0, 0.4, 0.4],
            [1.2, 0.4, 0.4, 0.0, 0.4, 0.4],
            [1.2, 0.4, 0.4, 0.0, 0.4, 0.4],
            [1.2, 0.4, 0.4, 0.0, 0.4, 0.4],
            [1.2, 0.4, 0.4, 0.0, 0.4, 0.4],
            [1.2, 0.4, 0.4, 0.0, 0.4, 0.4]]).transpose()

        self.assertTrue(np.max(np.abs(EngineeringStress-EngineeringStressCorrect)) < 1e-8)


class CheckExternalForce(Brick8NTestCase):
    def runTest(self):
        # calculate external force vector for the first load case (pressure)
        Fe = self.structure.BuildGlobalExternalLoadVector(0)

        # correct external force for pressure load vector (sum up the load in
        # x direction eveything else should be zero
        sumX = np.sum(Fe.J.Get("Displacements"))
        self.assertAlmostEqual(sumX, 8.0)

        # calculate external force vector for the second load cases (constDirection)
        Fe = self.structure.BuildGlobalExternalLoadVector(1)

        # correct external force for pressure load vector (sum up the load in
        # x direction eveything else should be zero
        sumY = np.sum(Fe.J.Get("Displacements"))
        self.assertAlmostEqual(sumY, 20.0)


if __name__ == '__main__':
    unittest.main()
