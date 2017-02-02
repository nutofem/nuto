# -*- coding: utf-8 -*-

import nuto
import sys
import os
import numpy as np

# Get path to work dir as cmake_current_binary_dir
pathToWorkDir = sys.argv[1]

# Split the current filename (basename) from its extension ...
testName = fileExt = os.path.splitext(os.path.basename(sys.argv[0]))[0]

# ... and join it with the path to work dir
workDir = os.path.join(pathToWorkDir, testName)

# create the work dir
if not os.path.exists(workDir):
    os.makedirs(workDir)

gmsh = sys.argv[3]


print "   |--->  Gmsh found at"
print "   |---> ",gmsh
print "   |--->  Print result files to"
print "   |---> ",workDir
print

# create a box specimen
specimen = nuto.Specimen(np.array([[-10., 10.], [-10., 10.], [-5., 5.]]), 0)

# create some spheres
sphereMatrix = np.array([[-4., -4., 0., 2.],
                         [-4.,  4., 0., 2.],
                         [ 4.,  4., 0., 2.],
                         [ 4., -4., 0., 2.]])
spheres = nuto.ParticleHandler(sphereMatrix, 0,0,0)


meshSize = 4

# define files

geoFile2D = workDir+"/2D.geo"
mshFile2D = workDir+"/2D.msh"
errFile2D = workDir+"/2D.err"

geoFile3D = workDir+"/3D.geo"
mshFile3D = workDir+"/3D.msh"
errFile3D = workDir+"/3D.err"

# run the export

spheres.ExportParticlesToGmsh2D(geoFile2D,specimen, meshSize,0,1)
spheres.ExportParticlesToGmsh3D(geoFile3D,specimen, meshSize)

# do the meshing
# redirect the stderr ("2>") to the error file (errors and warnings)
# -v 1 sets the verbose level to stderr only (as far as I know)
os.system(gmsh + ' ' + geoFile2D + ' -2 -order 1 -v 1 -o '+mshFile2D +' 2> ' + errFile2D)
os.system(gmsh + " " + geoFile3D + ' -3 -order 1 -v 1 -o '+mshFile3D +' 2> ' + errFile3D)

# check the emptyness of the error files
if (os.stat(errFile2D).st_size >0):
	print "Errors or warnings occured while meshing 2D. See error file for details:"
	print errFile2D
	sys.exit(-1)

if (os.stat(errFile3D).st_size >0):
	print "Errors or warnings occured while meshing 3D. See error file for details:"
	print errFile3D
	sys.exit(-1)

# check the mesh proper import of the mesh

structure = nuto.Structure(2)
structure.ImportFromGmsh(mshFile2D);
if structure.GetNumElements() == 0:
	print "Empty Structure"
	sys.exit(-1)

structure = nuto.Structure(3)
structure.ImportFromGmsh(mshFile3D);
if structure.GetNumElements() == 0:
	print "Empty Structure"
	sys.exit(-1)
