import nuto
import sys
import os
import numpy as np

# show the results on the screen
printResult = False

# no error in file, modified, if error is detected
error = False

# definitions
E = 80000
v = 1./3.

BoundaryDisplacement = 1.0


def RunPatchTest(StressState):
    global error
    # create structure
    myStructure = nuto.Structure(2)
    myStructure.SetShowTime(False)

    # create material law
    myMatLin = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress")
    myStructure.ConstitutiveLawSetParameterDouble(myMatLin, "Youngs_Modulus", E)
    myStructure.ConstitutiveLawSetParameterDouble(myMatLin, "Poissons_Ratio", v)

    # create section
    mySection = myStructure.SectionCreate("Plane_Stress")
    myStructure.SectionSetThickness(mySection, 1)

    numNodes = 8
    myStructure.NodesCreate(np.array([[  0, 10,  2,  8,  4,  8,  0, 10],
                                      [  0,  0,  2,  3,  7,  7, 10, 10]], dtype=float))

    elementIncidence = nuto.IntFullMatrix(4,5,(	
            3,2,0,1 ,
            4,6,0,2 ,
            5,4,2,3 ,
            7,5,3,1 ,
            7,6,4,5 ) )

    interpolationType = myStructure.InterpolationTypeCreate("Quad2D")
    myStructure.InterpolationTypeAdd(interpolationType, "Coordinates", "Equidistant1")
    myStructure.InterpolationTypeAdd(interpolationType, "Displacements", "Equidistant1")

    myStructure.ElementsCreate(interpolationType, elementIncidence)
    myStructure.ElementTotalConvertToInterpolationType()
    myStructure.ElementTotalSetConstitutiveLaw(myMatLin)
    myStructure.ElementTotalSetSection(mySection)

    LoadNodesXPos = myStructure.GroupCreate("Nodes")
    LoadNodesXNeg = myStructure.GroupCreate("Nodes")
    LoadNodesYPos = myStructure.GroupCreate("Nodes")
    LoadNodesYNeg = myStructure.GroupCreate("Nodes")

    myStructure.GroupAddNode(LoadNodesXPos, 1)
    myStructure.GroupAddNode(LoadNodesXPos, 7)
    myStructure.GroupAddNode(LoadNodesYPos, 6)
    myStructure.GroupAddNode(LoadNodesYPos, 7)

    myStructure.GroupAddNode(LoadNodesXNeg, 0)
    myStructure.GroupAddNode(LoadNodesXNeg, 6)
    myStructure.GroupAddNode(LoadNodesYNeg, 0)
    myStructure.GroupAddNode(LoadNodesYNeg, 1)

    directionX = np.array([1.0, 0.0])
    directionY = np.array([0.0, 1.0])

    print "Displacement control with stress state ", StressState

    if StressState == "XX":
        myStructure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesXNeg, directionX, 0.0)
        myStructure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesXPos, directionX, BoundaryDisplacement)
        myStructure.ConstraintLinearSetDisplacementNode(0, directionY, 0)
    elif StressState == "YY":
        myStructure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesYNeg, directionY, 0.0)
        myStructure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesYPos, directionY, BoundaryDisplacement)
        myStructure.ConstraintLinearSetDisplacementNode(0, directionX, 0)
    elif StressState == "XY":
        myStructure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesXNeg, directionX, 0.0)
        myStructure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesXNeg, directionY, 0.0)
        myStructure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesXPos, directionY, BoundaryDisplacement)
        myStructure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesXPos, directionX, 0)
    else:
        print 'Wrong stressstate given'
        error = True
        sys.exit(-1)

    # start analysis
    # build global dof numbering
    myStructure.NodeBuildGlobalDofs()
    myStructure.CalculateMaximumIndependentSets()

    myStructure.SolveGlobalSystemStaticElastic()

    # calculate residual
    internalGradient = myStructure.BuildGlobalInternalGradient()
    if (printResult):
        residualVector = internalGradient.J.Export()*(-1.)
        print "residual: " + str(np.linalg.norm(residualVector))

    if (np.linalg.norm(internalGradient.J.Export()) > 1e-8):
        print 'Internal and external forces differs.'
        error = True

    # calculate internal force vector
    Fi = myStructure.ElementBuildInternalGradient(2)
    if (printResult):
        print "Internal Force"
        print Fi

    # calculate engineering strain of 2 at all integration points
    # the size the matrix is not important and reallocated within the procedure
    EngineeringStrain = myStructure.ElementGetEngineeringStrain(2)

    # calculate engineering strain of 2 at all integration points
    EngineeringStress = myStructure.ElementGetEngineeringStress(2)
    # correct stress
    EngineeringStressCorrect = np.zeros((6, 4))
    for i in range(0, 4):
        if StressState == "XX":
            sigma = E*BoundaryDisplacement / 10.
            EngineeringStressCorrect[0, i] = sigma
        elif StressState == "YY":
            sigma = E*BoundaryDisplacement / 10.
            EngineeringStressCorrect[1, i] = sigma
        elif StressState == "XY":
            sigma = E*BoundaryDisplacement / 10. / (2+2*v)
            EngineeringStressCorrect[5, i] = sigma

    if (printResult):
        print "EngineeringStressCorrect"
        print EngineeringStressCorrect
        print "EngineeringStress"
        print EngineeringStress

    if (np.max(np.abs(EngineeringStress - EngineeringStressCorrect)) > 1e-4):
        print 'stress is not correct.'
        error = True
        if (printResult):
            print "np.max(np.abs(EngineeringStress - EngineeringStressCorrect))"
            print np.max(np.abs(EngineeringStress - EngineeringStressCorrect))

    C1 = 1./E
    C2 = -v/E
    C3 = 2*(1+v)/E

    D = np.array([
            [C1, C2, C2,  0,  0,  0],
            [C2, C1, C2,  0,  0,  0],
            [C2, C2, C1,  0,  0,  0],
            [0,  0,  0,  C3,  0,  0],
            [0,  0,  0,   0, C3,  0],
            [0,  0,  0,   0,  0, C3]])

    # correct strain
    EngineeringStrainCorrect = D.dot(EngineeringStressCorrect)

    if (printResult):
        print "EngineeringStrainCorrect"
        print EngineeringStrainCorrect
        print "EngineeringStrain"
        print EngineeringStrain

    if (np.max(np.abs(EngineeringStrain - EngineeringStrainCorrect)) > 1e-8):
        print 'strain is not correct.'
        error = True


RunPatchTest("XX")
RunPatchTest("YY")
RunPatchTest("XY")

if (error):
    sys.exit(-1)
else:
    sys.exit(0)
