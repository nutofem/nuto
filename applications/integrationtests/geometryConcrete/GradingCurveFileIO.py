#!/usr/bin/env python3
import os
import numpy as np
import nuto

""" Test Input File reading

For 2D and 3D do:
  1) Read information from an input file
  2) Perform the take-and-place algorithm
  3) Perform the EDMD algorithm

If everything runs without SEGFAULTS or exceptions, the test passes.
"""


def RunSimulationFromInputFile(inputFile, visuFileName):
    reader = nuto.InputReader(inputFile)
    reader.ReadFile()

    specimen = nuto.Specimen(reader.GetBoundingBox(), reader.GetTypeOfSpecimen())

    creator = nuto.ParticleCreator(specimen, reader.GetShrinkage())
    spheresBoundary = np.zeros((0, 4))
    seed = 6174
    spheresMatrix = creator.CreateSpheresInSpecimen(
        reader.GetVolumeFraction(),
        reader.GetGradingCurve(),
        0.0,
        reader.GetAbsoluteDistance(),
        seed,
        spheresBoundary)

    # particle handler
    spheres = nuto.ParticleHandler(
        spheresMatrix,
        reader.GetRandomVelocityRange(),
        reader.GetRelativeGrowthRate(),
        reader.GetAbsoluteGrowthRate())
    spheres.SetVisualizationFileName(visuFileName)
    minDist = spheres.GetAbsoluteMininimalDistance(specimen)
    print(minDist)

    # sub box handler
    subBoxes = nuto.SubBoxHandler(spheres, specimen, 10)
    subBoxes.VisualizeBorders(resultDir + "/borders.vtu")

    collisions = nuto.CollisionHandler(spheres, subBoxes, resultDir)

    wallTimeMax = (1. / (1. - reader.GetShrinkage()) - 1.) / reader.GetRelativeGrowthRate()

    collisions.Simulate(
        reader.GetNumEventsMax(),
        reader.GetTimeMax(),
        wallTimeMax,
        reader.GetTimePrintOut(),
        reader.GetInitialTimeBarrier())

    spheres.ExportParticlesToVTU3D(resultDir, 0, 0, False)
    spheres.ExportParticlesToVTU3D(resultDir, 1, 1, True)

    reader.Close()


if __name__ == "__main__":
    cwd = os.getcwd()
    resultDir = os.path.join(cwd, "GradingCurveFileIO")

    if not os.path.exists(resultDir):
        os.makedirs(resultDir)

    RunSimulationFromInputFile("input3D_test.dat", "3dSpheres")
    RunSimulationFromInputFile("input2D_test.dat", "2dSpheres")
