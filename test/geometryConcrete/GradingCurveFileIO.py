# -*- coding: utf-8 -*-

import nuto
import math
import sys
import os



# ======================================================
# ==                                                  ==
# ==  FOR 2D and 3D DO:                               ==
# ==    1) Read information from an input file        ==
# ==    2) Perform the take-and-place algorithm       ==
# ==    3) Perform the EDMD algorithm                 ==
# ==                                                  ==
# ==  If everything runs without SEGFAULTS or         ==
# ==  exceptions, the test passes.                    ==
# ==                                                  ==
# ======================================================

inputFileName3D = "input3D_test.dat"
IN_seed = 6174


# ======================================================
# ==       Function for running an input file         ==
# ======================================================
def RunSimulationFromInputFile (rInputFile, rWorkDir, rVisuFileName):
  IN = nuto.InputReader(rInputFile)
  IN.ReadFile()

  specimen = nuto.Specimen(IN.GetBoundingBox(), IN.GetTypeOfSpecimen())
  
  creator = nuto.ParticleCreator(specimen, IN.GetShrinkage())
  spheresBoundary = nuto.DoubleFullMatrix(0,4,[])
  spheresMatrix = creator.CreateSpheresInSpecimen(
    IN.GetVolumeFraction(),
    IN.GetGradingCurve(),
    0., 
    IN.GetAbsoluteDistance(),
    IN_seed,
    spheresBoundary)
  
  # particle handler
  spheres = nuto.ParticleHandler(
    spheresMatrix, 
    IN.GetRandomVelocityRange(), 
    IN.GetRelativeGrowthRate(), 
    IN.GetAbsoluteGrowthRate())
  spheres.SetVisualizationFileName(rVisuFileName)
  minDist = spheres.GetAbsoluteMininimalDistance(specimen)
  print minDist
  
  # sub box handler
  subBoxes = nuto.SubBoxHandler(spheres, specimen, 10)
  subBoxes.VisualizeBorders(rWorkDir + "/borders.vtu")

  collisions = nuto.CollisionHandler(spheres, subBoxes, rWorkDir)

  wallTimeMax = (1. / (1. - IN.GetShrinkage()) - 1.) / IN.GetRelativeGrowthRate();

  collisions.Simulate(
    IN.GetNumEventsMax(), 
    IN.GetTimeMax(), 
    wallTimeMax , 
    IN.GetTimePrintOut(), 
    IN.GetInitialTimeBarrier())

  spheres.ExportParticlesToVTU3D(rWorkDir, 0, 0, False);
  spheres.ExportParticlesToVTU3D(rWorkDir, 1, 1, True);

  IN.Close()  
  del collisions
  del subBoxes
  del spheres
  del spheresBoundary
  del spheresMatrix
  del creator

    
# ======================================================
# ==       Build Pathes for workdir and input         ==
# ======================================================

# Get path to work dir as cmake_current_binary_dir
pathToWorkDir = sys.argv[1] 

# Split the current filename (basename) from its extension (most likely "GradingCurve") ... 
testName = fileExt = os.path.splitext(os.path.basename(sys.argv[0]))[0]

# ... and join it with the path to work dir
workDir = os.path.join(pathToWorkDir, testName)

# get path to input file as cmake_current_source_dir
inputDir = sys.argv[2]

# create the work dir 
if not os.path.exists(workDir):
    os.makedirs(workDir)

print
print "   |--->  Read input files from"
print "   |---> ",inputDir
print
print "   |--->  Print result files to"
print "   |---> ",workDir
print

inputFile3D = os.path.join(inputDir, inputFileName3D)

# ======================================================
# ==           Run Examples                           ==
# ======================================================
RunSimulationFromInputFile(inputFile3D, workDir, "3dSpheres")

