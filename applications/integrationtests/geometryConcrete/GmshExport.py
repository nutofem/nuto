#!/usr/bin/env python3
import math
import unittest
import sys
import os
import subprocess
import numpy as np
import nuto


def mesh(geoFile, dimension):
    """ Mesh a given geoFile using gmsh and return the filename of the mesh file."""
    name, _ = os.path.splitext(geoFile)
    mshFile = name + ".msh"
    dimensionFlag = "-{0}".format(dimension)
    subprocess.check_call([gmsh, geoFile, dimensionFlag, "-order", "2", "-o", mshFile], stdout=subprocess.DEVNULL)
    return mshFile


class CylinderMesh(unittest.TestCase):
    def testVolumes(self):
        r = 10.0
        h = 20.0
        boundingBox = np.array([[-r, r], [-r, r], [0.0, h]])
        # the magic number "2" here means cylinder specimen
        # the interface is not my fault
        cylinder = nuto.Specimen(boundingBox, 2)
        r_sphere = 5.0
        sphereArray = np.array([[0.0, 0.0, h/2.0, r_sphere]])
        sphere = nuto.ParticleHandler(sphereArray, 0, 0, 0)
        geoFile = os.path.join(resultDir, "cylinder.geo")
        sphere.ExportParticlesToGmsh3D(geoFile, cylinder, r/4.0)
        mshFile = mesh(geoFile, 3)
        structure = nuto.Structure(3)
        groups = structure.ImportFromGmsh(mshFile)

        cylinderVolume = math.pi * r**2 * h
        sphereVolume = 4.0/3.0 * math.pi * r_sphere**3
        meshedVolumeMatrix = structure.ElementGroupGetVolume(groups[0][0])
        meshedVolumeSphere = structure.ElementGroupGetVolume(groups[1][0])
        self.assertAlmostEqual(meshedVolumeMatrix, cylinderVolume - sphereVolume, delta=1)
        self.assertAlmostEqual(meshedVolumeSphere, sphereVolume, delta=1)


class BoxMesh(unittest.TestCase):
    def setUp(self):
        self.a = a = 10.0
        boundingBox = np.array([[-a, a], [-a, a], [-a, a]])
        # the magic number "0" here means box specimen
        # the interface is still not my fault
        self.box = nuto.Specimen(boundingBox, 0)

        # create a sphere at origin with radius 5
        self.r = 5.0
        sphereMatrix = np.array([[0., 0., 0., self.r]])
        self.spheres = nuto.ParticleHandler(sphereMatrix, 0, 0, 0)

    def test3DVolume(self):
        meshSize = 5.0
        geoFile = resultDir + "/3D.geo"
        self.spheres.ExportParticlesToGmsh3D(geoFile, self.box, meshSize)
        mshFile = mesh(geoFile, 3)
        structure = nuto.Structure(3)
        groups = structure.ImportFromGmsh(mshFile)

        boxVolume = (2*self.a)**3
        sphereVolume = 4.0/3.0 * math.pi * self.r**3
        meshedVolumeMatrix = structure.ElementGroupGetVolume(groups[0][0])
        meshedVolumeSphere = structure.ElementGroupGetVolume(groups[1][0])
        self.assertAlmostEqual(meshedVolumeMatrix, boxVolume - sphereVolume, delta=4)
        self.assertAlmostEqual(self.box.GetVolume(), boxVolume)
        self.assertAlmostEqual(meshedVolumeSphere, sphereVolume, delta=4)
        self.assertAlmostEqual(self.spheres.GetVolume(), sphereVolume)

    def test2DSomething(self):
        meshSize = 5.0
        geoFile = resultDir + "/2D.geo"
        zPlane = 0.0
        minRadius = 1.0
        self.spheres.ExportParticlesToGmsh2D(geoFile, self.box, meshSize, zPlane, minRadius)
        mshFile = mesh(geoFile, 2)
        structure = nuto.Structure(2)
        groups = structure.ImportFromGmsh(mshFile)

        rectangleArea = (2*self.a)**2
        circleArea = math.pi * self.r**2
        meshedVolumeMatrix = structure.ElementGroupGetVolume(groups[0][0])
        meshedVolumeCircle = structure.ElementGroupGetVolume(groups[1][0])
        self.assertAlmostEqual(meshedVolumeMatrix, rectangleArea - circleArea, delta=1)
        self.assertAlmostEqual(meshedVolumeCircle, circleArea, delta=1)


if __name__ == "__main__":
    # Get path to build dir, and create result dir if it doesn't exist
    resultDir = os.path.join(os.getcwd(), "GmshExport")
    if not os.path.exists(resultDir):
        os.makedirs(resultDir)

    gmsh = sys.argv[1]
    sys.argv[1:] = []

    unittest.main()
