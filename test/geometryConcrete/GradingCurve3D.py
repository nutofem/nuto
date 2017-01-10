# -*- coding: utf-8 -*-

import nuto
import math
import sys
import os
import numpy as np

# This test is designed to run without error.
# If the c++ code exits with an unhandled exeption,
# the test fails automatically.


IN_shrinkage = 0.02
IN_boxType = 0
IN_boundingBox = np.array([[0., 40.], [0., 40.], [0., 40.]])
IN_volumeFraction = 0.5
IN_gradingCurve = np.array([[8., 16., 0.40],
                            [4.,  8., 0.24],
                            [2.,  4., 0.15]])
IN_relativeDistance = 0.0
IN_absoluteDistance = 0.0
IN_seed = 6174
IN_spheresBoundary = np.zeros((0,4))

IN_velocityRange = 0.1
IN_relativeGrowthRate = 1.
IN_absoluteGrowthRate = 0.

IN_InitialTimeBarrier = 0.1


specimen = nuto.Specimen(IN_boundingBox, IN_boxType)

creator = nuto.ParticleCreator(specimen, IN_shrinkage)

spheresMatrix = creator.CreateSpheresInSpecimen(
  IN_volumeFraction, 
  IN_gradingCurve, 
  IN_relativeDistance, 
  IN_absoluteDistance, 
  IN_seed, 
  IN_spheresBoundary)


# particle handler
spheres = nuto.ParticleHandler(spheresMatrix, IN_velocityRange, IN_relativeGrowthRate, IN_absoluteGrowthRate)

# sub box handler
subBoxes = nuto.SubBoxHandler(spheres, specimen, 12)


collisions = nuto.CollisionHandler(spheres,subBoxes, "")

wallTimeMax = (1. / (1. - IN_shrinkage) - 1.) / IN_relativeGrowthRate;

collisions.Simulate(10000000, 60., wallTimeMax , 1., IN_InitialTimeBarrier)










