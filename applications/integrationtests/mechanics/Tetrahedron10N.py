# -*- coding: utf-8 -*-
import nuto
import sys
import os
import numpy as np

# call of the test file, e.g.
# /usr/local/bin/python ~/develop/nuto/test/mechanics/Tetrahedron10N.py Linux x86_64 ~/develop/nuto/test/mechanics

# if set to true, the result will be generated (for later use in the test routine)
# otherwise, the current result will be compared to the stored result
createResult = False

# show the results on the screen
printResult = True

# system name and processor
system = sys.argv[1]+sys.argv[2]

# path in the original source directory and current filename at the end
pathToResultFiles = os.path.join(sys.argv[3], "results", system, os.path.basename(sys.argv[0]))

# remove the extension
fileExt = os.path.splitext(sys.argv[0])[1]
pathToResultFiles = pathToResultFiles.replace(fileExt, '')

# no error in file, modified, if error is detected
error = False

# create structure
myStructure = nuto.Structure(3)

# create nodes
myNode1  = myStructure.NodeCreate(np.array([0.0, 0.0, 0.0]))
myNode2  = myStructure.NodeCreate(np.array([1.0, 0.0, 0.0]))
myNode3  = myStructure.NodeCreate(np.array([0.0, 1.0, 0.0]))
myNode4  = myStructure.NodeCreate(np.array([0.0, 0.0, 1.0]))
myNode5  = myStructure.NodeCreate(np.array([0.5, 0.0, 0.0]))
myNode6  = myStructure.NodeCreate(np.array([0.5, 0.5, 0.0]))
myNode7  = myStructure.NodeCreate(np.array([0.0, 0.5, 0.0]))
myNode8  = myStructure.NodeCreate(np.array([0.0, 0.0, 0.5]))
myNode9  = myStructure.NodeCreate(np.array([0.0, 0.5, 0.5]))
myNode10 = myStructure.NodeCreate(np.array([0.5, 0.0, 0.5]))

# create section
mySection = myStructure.SectionCreate("Volume")

myInterpolationType = myStructure.InterpolationTypeCreate("TETRAHEDRON3D")
myStructure.InterpolationTypeAdd(myInterpolationType, "Coordinates","EQUIDISTANT2")
myStructure.InterpolationTypeAdd(myInterpolationType, "Displacements", "EQUIDISTANT2")

# create element
myElement1 = myStructure.ElementCreate(myInterpolationType,
        [myNode1, myNode2, myNode3, myNode4, myNode5, myNode6, myNode7, myNode8, myNode9, myNode10])
myStructure.ElementTotalConvertToInterpolationType()

# create constitutive law
myMatLin = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress")
myStructure.ConstitutiveLawSetParameterDouble(myMatLin, "Youngs_Modulus", 10)
myStructure.ConstitutiveLawSetParameterDouble(myMatLin, "Poissons_Ratio", 0.25)

# assign constitutive law
myStructure.ElementSetConstitutiveLaw(myElement1, myMatLin)
myStructure.ElementSetSection(myElement1, mySection)

# make group of boundary nodes
groupBoundaryNodes = myStructure.GroupCreate("Nodes")
myStructure.GroupAddNode(groupBoundaryNodes, myNode1)
myStructure.GroupAddNode(groupBoundaryNodes, myNode3)
myStructure.GroupAddNode(groupBoundaryNodes, myNode4)
myStructure.GroupAddNode(groupBoundaryNodes, myNode7)
myStructure.GroupAddNode(groupBoundaryNodes, myNode8)
myStructure.GroupAddNode(groupBoundaryNodes, myNode9)

# make group of boundary elements (in this case it is just one)
groupBoundaryElements = myStructure.GroupCreate("Elements")
myStructure.GroupAddElementsFromNodes(groupBoundaryElements, groupBoundaryNodes, False)

# create surface loads (0 - pressure on X, 1-const-direction Y)
myStructure.LoadSurfacePressureCreate3D(0, groupBoundaryElements, groupBoundaryNodes, 2.)
myStructure.LoadSurfaceConstDirectionCreate3D(1, groupBoundaryElements, groupBoundaryNodes, np.array([0., 5., 0.]))

# set displacements of right node
myStructure.NodeSetDisplacements(myNode2, np.array([0.2, 0.2, 0.2]))
myStructure.NodeSetDisplacements(myNode3, np.array([0.2, 0.2, 0.2]))
myStructure.NodeSetDisplacements(myNode6, np.array([0.2, 0.2, 0.2]))
myStructure.NodeSetDisplacements(myNode7, np.array([0.2, 0.2, 0.2]))

# calculate element stiffness matrix
Ke = myStructure.ElementBuildHessian0(myElement1).Export()
if (printResult):
    print "Ke"
    print Ke

# correct stiffness matrix
if createResult:
    print pathToResultFiles+"Stiffness.txt"
    np.savetxt(pathToResultFiles+"Stiffness.txt", Ke, header="#Correct result")
else:
    KeCorrect = np.loadtxt(pathToResultFiles+"Stiffness.txt", skiprows=1)
    if (printResult):
        print "KeCorrect"
        print KeCorrect
    if (np.max(np.abs(Ke-KeCorrect)) > 1e-8):
        print '[' + system, sys.argv[0] + '] : stiffness is not correct.'
        error = True

# calculate internal force vector
Fi = myStructure.ElementBuildInternalGradient(myElement1).Export()
Fi = Fi.squeeze()

# correct resforce vector
if createResult:
    print pathToResultFiles+"Internalforce.txt"
    np.savetxt(pathToResultFiles+"Internalforce.txt", Fi, header="#Correct result")
else:
    FiCorrect = np.loadtxt(pathToResultFiles+"Internalforce.txt", skiprows=1)
    if (printResult):
        print "FiCorrect"
        print FiCorrect
    if (np.max(np.abs(Fi-FiCorrect)) > 1e-8):
        print '[' + system, sys.argv[0] + '] : internal force is not correct.'
        error = True


# check stiffness with internal force vector
prevDisp = np.zeros((3,1))

delta = 1e-4
rows, cols = Ke.shape
KeApprox = np.zeros((rows, cols))
curColumn = 0
for theNode in range(0, cols/myStructure.GetDimension()):
    for theDof in range(0, myStructure.GetDimension()):
        myStructure.NodeGetDisplacements(theNode, prevDisp)
        prevDisp[theDof] += delta
        myStructure.NodeSetDisplacements(theNode, prevDisp)
        Fi_new = myStructure.ElementBuildInternalGradient(myElement1).Export()
        Fi_new = Fi_new.squeeze()
        prevDisp[theDof] -= delta
        myStructure.NodeSetDisplacements(theNode, prevDisp)
        KeApprox[:, curColumn] = (Fi_new-Fi)*(1./delta)
        curColumn += 1

if (printResult):
    print "KeApprox with central differences"
    print KeApprox


# check stiffness with internal force vector
if (np.max(np.abs(KeApprox-Ke)) > 1e-8):
    print '[' + system, sys.argv[0] + '] : stiffness matrix via central differences and resforces not correct.'
    error = True

# calculate engineering strain of myelement1 at all integration points
# the size the matrix is not important and reallocated within the procedure
EngineeringStrain = myStructure.ElementGetEngineeringStrain(myElement1)

# correct strain
EngineeringStrainCorrect = np.array([
    [-0.08944272, 0.37888544, -0.11055728,  0.26832816, -0.20000000, 0.28944272],
    [ 0.26832816, 0.37888544, -0.11055728,  0.26832816,  0.15777088, 0.64721360],
    [-0.08944272, 0.02111456, -0.46832816, -0.44721360, -0.55777088,-0.06832816],
    [-0.08944272, 0.02111456, -0.11055728, -0.08944272, -0.20000000,-0.06832816]
    ]).transpose()

if (printResult):
    print "EngineeringStrainCorrect"
    print EngineeringStrainCorrect
    print "EngineeringStrain"
    print EngineeringStrain

if (np.max(np.abs(EngineeringStrain - EngineeringStrainCorrect)) > 1e-5):
    print '[' + system, sys.argv[0] + '] : strain is not correct.'
    error = True

# calculate engineering strain of myelement1 at all integration points
EngineeringStress = myStructure.ElementGetEngineeringStress(myElement1)

# correct stress
EngineeringStressCorrect = np.array([
    [ 0.00000000, 3.74662528 , -0.16891648 , 1.07331264 ,-0.80000000, 1.15777088], 
    [ 4.29325056, 5.17770880 ,  1.26216704 , 1.07331264 , 0.63108352, 2.58885440],
    [-2.86216704, -1.97770880, -5.89325056 ,-1.78885440 ,-2.23108352,-0.27331264],
    [-1.43108352, -0.54662528, -1.60000000 ,-0.35777088 ,-0.80000000,-0.27331264]
    ]).transpose()

if (printResult):
    print "EngineeringStressCorrect"
    print EngineeringStressCorrect
    print "EngineeringStress"
    print EngineeringStress

if (np.max(np.abs(EngineeringStress - EngineeringStressCorrect)) > 1e-5):
    print '[' + system, sys.argv[0] + '] : stress is not correct.'
    error = True

# calculate external force vector for the first load case (pressure)
Fe = myStructure.BuildGlobalExternalLoadVector(0)

if (printResult):
    print "Fe.J for pressure load"
    print Fe.J.Get("Displacements")

# correct external force for pressure load vector
# (sum up the load in x direction eveything else should be zero)
sumX = np.sum(Fe.J.Get("Displacements"))
if (abs(sumX-1.) > 1e-8):
    print '[' + system, sys.argv[0] + '] : pressure load is not correct.'
    error = True

# calculate external force vector for the second load cases (constDirection)
Fe = myStructure.BuildGlobalExternalLoadVector(1)

if (printResult):
    print "Fe.J const direction load"
    print Fe.J.Get("Displacements")

# correct external force for pressure load vector
# (sum up the load in x direction eveything else should be zero)
sumY = np.sum(Fe.J.Get("Displacements"))
if (abs(sumY - 2.5) > 1e-8):
    print '[' + system, sys.argv[0] + '] : const direction load is not correct.'
    error = True

# visualize results
visualizationGroup = myStructure.GroupCreate("Elements")
myStructure.GroupAddElementsTotal(visualizationGroup)

myStructure.AddVisualizationComponent(visualizationGroup, "Displacements")
myStructure.AddVisualizationComponent(visualizationGroup, "EngineeringStrain")
myStructure.AddVisualizationComponent(visualizationGroup, "EngineeringStress")
myStructure.ExportVtkDataFileElements("Tetrahedron10N.vtk")

if (error):
    sys.exit(-1)
else:
    sys.exit(0)
