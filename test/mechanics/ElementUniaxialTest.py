import nuto
import sys
import os
import numpy as np

lX = 365.
lY = 13.
lZ = 37.

E = 42.
nu = 0.26
deltaL = 0.6174
rho = 3.141529

error = False
errorMsg = ""


def Run(rStructure, rType, rOrder):
    global error, errorMsg

    resultFile = os.path.join(pathToResultFiles, rType+rOrder+"_py.vtu")

    print "#########################################"
    print "##   Running " + rType + ":" + rOrder
    print "#########################################"
    print "##        writing vtu files to "
    print "## " + resultFile
    print "#########################################"

    # set constitutive law
    lawId = rStructure.ConstitutiveLawCreate("LINEAR_ELASTIC_ENGINEERING_STRESS")
    rStructure.ConstitutiveLawSetParameterDouble(lawId, "Youngs_Modulus", E)
    rStructure.ConstitutiveLawSetParameterDouble(lawId, "Poissons_Ratio", nu)
    rStructure.ConstitutiveLawSetParameterDouble(lawId, "Density", rho)

    rStructure.ElementTotalSetConstitutiveLaw(lawId)

    # set boundary conditions
    dimension = rStructure.GetDimension()
    directionX = np.zeros((dimension, 1))
    directionX[0, 0] = 1.0

    origin = np.zeros((dimension, 1))
    nodeGroupOrigin = rStructure.GroupCreate("Nodes")
    rStructure.GroupAddNodeRadiusRange(nodeGroupOrigin, origin, 0, 1.e-5)
    if (rStructure.GroupGetNumMembers(nodeGroupOrigin) != 1):
        errorMsg += "[SetBoundaryConditions:" + rType + ":" + rOrder + "] Node at origin (0,0,0) does not exist. \n"
        error = True
        return

    # fix origin
    nodeOrigin = rStructure.GroupGetMemberIds(nodeGroupOrigin)[0]
    for dim in range(1, dimension):
        direction = np.zeros((dimension, 1))
        direction[dim, 0] = 1.0
        rStructure.ConstraintLinearSetDisplacementNode(nodeOrigin, direction, 0.0)

    # fix x = 0 plane
    nodesX0 = rStructure.GroupCreate("Nodes")
    rStructure.GroupAddNodeCoordinateRange(nodesX0, 0, -1.e-6, 1.e-6)
    rStructure.ConstraintLinearSetDisplacementNodeGroup(nodesX0, directionX, 0.)

    # apply displacement on x = lX plane
    nodesXlX = rStructure.GroupCreate("Nodes")
    rStructure.GroupAddNodeCoordinateRange(nodesXlX, 0, lX-1.e-6, lX+1.e-6)
    rStructure.ConstraintLinearSetDisplacementNodeGroup(nodesXlX, directionX, deltaL)

    # solve
    rStructure.CalculateMaximumIndependentSets()
    rStructure.SolveGlobalSystemStaticElastic()

    internalGradient = rStructure.BuildGlobalInternalGradient()
    if (np.max(np.abs(internalGradient.J.Export())) > 1e-8):
        internalGradient.J.Export()
        errorMsg += "[SetBoundaryConditions:" + rType + ":" + rOrder + "] residual force vector is not zero. \n"
        error = True
        return

    # check stresses
    analyticStrainX = deltaL / lX
    analyticStressX = analyticStrainX * E
    allElements = rStructure.GroupCreate("Elements")
    allNodes = rStructure.GroupCreate("Nodes")
    rStructure.GroupAddNodeCoordinateRange(allNodes, 0, -0.1, lX+0.1)
    rStructure.GroupAddElementsFromNodes(allElements, allNodes, True)
    elementIds = rStructure.GroupGetMemberIds(allElements)
    for elementId in elementIds:
        stress = rStructure.ElementGetEngineeringStress(elementId)
        for iIP in range(0, stress.shape[1]):
            numericStress = stress[0, iIP]
            if (abs(numericStress-analyticStressX) > 1.e-6):
                errorMsg += "[CheckSolution:" + rType + ":" + rOrder + "] wrong stress calculation. \n"
                error = True
                return

    # check reaction forces
    analyticForce = analyticStressX * lY * lZ
    numericForce = 0.
    nodeX0Indices = rStructure.GroupGetMemberIds(nodesX0)
    for nodeId in nodeX0Indices:
        force = np.zeros((dimension, 1))
        rStructure.NodeInternalForce(nodeId, force)
        numericForce += force[0]

    if (abs(numericForce-numericForce) > 1.e-8):
        errorMsg += "[CheckSolution:" + rType + ":" + rOrder + "] wrong reaction force calculation. \n"
        error = True
        return

    # Check Total Mass via Mass matrix
    analyticMass = lX*lY*lZ*rho

    hessian2 = rStructure.BuildGlobalHessian2()

    numericMass = 0.
    numericMass += np.sum(hessian2.JJ.ExportToFullMatrix())
    numericMass += np.sum(hessian2.JK.ExportToFullMatrix())
    numericMass += np.sum(hessian2.KJ.ExportToFullMatrix())
    numericMass += np.sum(hessian2.KK.ExportToFullMatrix())

    # since the mass is added to nodes in every direction
    numericMass = numericMass / dimension

    if(abs(numericMass - analyticMass)/numericMass > 1.e-6):
        print "mass analytical : ", analyticMass
        print "mass numerical  : ", numericMass
        errorMsg += "[CheckSolution:" + rType + ":" + rOrder + "] wrong mass calculation. \n"
        error = True
        return

    hessian2 = rStructure.BuildGlobalHessian2Lumped()

    numericMass = 0.
    numericMass += np.sum(hessian2.JJ.ExportToFullMatrix())
    numericMass += np.sum(hessian2.JK.ExportToFullMatrix())
    numericMass += np.sum(hessian2.KJ.ExportToFullMatrix())
    numericMass += np.sum(hessian2.KK.ExportToFullMatrix())

    # since the mass is added to nodes in every direction
    numericMass = numericMass / dimension

    if(abs(numericMass - analyticMass)/numericMass > 1.e-6):
        print "mass analytical : ", analyticMass
        print "mass numerical  : ", numericMass
        errorMsg += "[CheckSolution:" + rType + ":" + rOrder + "] wrong lumped mass calculation. \n"
        error = True
        return

    #  Visualize
    visualizationGroup = rStructure.GroupCreate("Elements")
    rStructure.GroupAddElementsTotal(visualizationGroup)

    rStructure.AddVisualizationComponent(visualizationGroup, "Displacements")
    rStructure.AddVisualizationComponent(visualizationGroup, "EngineeringStrain")
    rStructure.AddVisualizationComponent(visualizationGroup, "EngineeringStress")

    rStructure.ExportVtkDataFileElements(resultFile, True)


def Run3D(r3DShape, rTypeOrder):
    global error, errorMsg
    myStructure = nuto.Structure(3)
    myStructure.SetShowTime(False)

    numElementsX = 2
    numElementsY = 1
    numElementsZ = 1

    # create nodes
    numNodesX = numElementsX+1
    numNodesY = numElementsY+1
    numNodesZ = numElementsZ+1

    deltaX = lX/numElementsX
    deltaY = lY/numElementsY
    deltaZ = lZ/numElementsZ

    nodeNum = 0
    for countZ in range(0, numNodesZ):
        for countY in range(0, numNodesY):
            for countX in range(0, numNodesX):
                myStructure.NodeCreate(nodeNum, np.array([countX*deltaX, countY*deltaY, countZ*deltaZ]))
                nodeNum += 1

    myInterpolationType = myStructure.InterpolationTypeCreate(r3DShape);
    myStructure.InterpolationTypeAdd(myInterpolationType, "Coordinates", "EQUIDISTANT1")
    myStructure.InterpolationTypeAdd(myInterpolationType, "Displacements", rTypeOrder)

    # create elements
    nodes = list(range(8))
    for countZ in range(0, numElementsZ):
        for countY in range(0, numElementsY):
            for countX in range(0, numElementsX):
                nodes[0] = countX  +  countY   *numNodesX +  countZ    * numNodesX * numNodesY
                nodes[1] = countX+1+  countY   *numNodesX +  countZ    * numNodesX * numNodesY
                nodes[2] = countX+1+ (countY+1)*numNodesX +  countZ    * numNodesX * numNodesY
                nodes[3] = countX  + (countY+1)*numNodesX +  countZ    * numNodesX * numNodesY
                nodes[4] = countX  +  countY   *numNodesX + (countZ+1) * numNodesX * numNodesY
                nodes[5] = countX+1+  countY   *numNodesX + (countZ+1) * numNodesX * numNodesY
                nodes[6] = countX+1+ (countY+1)*numNodesX + (countZ+1) * numNodesX * numNodesY
                nodes[7] = countX  + (countY+1)*numNodesX + (countZ+1) * numNodesX * numNodesY
                if (r3DShape == "Brick3D"):
                    myStructure.ElementCreate(myInterpolationType, nodes)
                elif (r3DShape == "Tetrahedron3D"):
                    myStructure.ElementCreate(myInterpolationType, [nodes[0], nodes[1], nodes[3], nodes[7]])
                    myStructure.ElementCreate(myInterpolationType, [nodes[0], nodes[1], nodes[7], nodes[4]])
                    myStructure.ElementCreate(myInterpolationType, [nodes[5], nodes[4], nodes[7], nodes[1]])
                    myStructure.ElementCreate(myInterpolationType, [nodes[6], nodes[5], nodes[7], nodes[1]])
                    myStructure.ElementCreate(myInterpolationType, [nodes[2], nodes[7], nodes[1], nodes[6]])
                    myStructure.ElementCreate(myInterpolationType, [nodes[2], nodes[3], nodes[1], nodes[7]])
                elif (r3DShape == "Prism3D"):
                    myStructure.ElementCreate(myInterpolationType, [nodes[0], nodes[1], nodes[2], nodes[4], nodes[5], nodes[6]])
                    myStructure.ElementCreate(myInterpolationType, [nodes[0], nodes[2], nodes[3], nodes[4], nodes[6], nodes[7]])
                else:
                    error = True
                    errorMsg += "Element shape " + r3DShape + " is invalid. \n"
                    return


    myStructure.ElementTotalConvertToInterpolationType()

    mySection = myStructure.SectionCreate("Volume")
    myStructure.ElementTotalSetSection(mySection)

    Run(myStructure, r3DShape , rTypeOrder)


def Run2D(r2DShape, rTypeOrder):
    global error, errorMsg
    myStructure = nuto.Structure(2)
    myStructure.SetShowTime(False)

    numElementsX = 4
    numElementsY = 2

    # create nodes
    numNodesX = numElementsX+1
    numNodesY = numElementsY+1
    deltaX = lX/numElementsX
    deltaY = lY/numElementsY

    nodeNum = 0
    for countY in range(0, numNodesY):
        for countX in range(0, numNodesX):
            myStructure.NodeCreate(nodeNum, np.array([countX*deltaX,countY*deltaY]))
            nodeNum += 1

    myInterpolationType = myStructure.InterpolationTypeCreate(r2DShape)
    myStructure.InterpolationTypeAdd(myInterpolationType, "Coordinates", "EQUIDISTANT1")
    myStructure.InterpolationTypeAdd(myInterpolationType, "Displacements", rTypeOrder)

    # create elements
    nodes = nuto.IntVector()
    for countY in range(0, numElementsY):
        for countX in range(0, numElementsX):
            if (r2DShape == "Quad2D"):
                nodes.resize(4)
                nodes[0] = countX  +  countY   *numNodesX
                nodes[1] = countX+1+  countY   *numNodesX
                nodes[2] = countX+1+ (countY+1)*numNodesX
                nodes[3] = countX  + (countY+1)*numNodesX
                myStructure.ElementCreate(myInterpolationType, nodes)
            elif (r2DShape == "Triangle2D"):
                nodes.resize(3)
                nodes[0] = countX  +  countY   *numNodesX
                nodes[1] = countX+1+  countY   *numNodesX
                nodes[2] = countX+1+ (countY+1)*numNodesX
                myStructure.ElementCreate(myInterpolationType, nodes)

                nodes.resize(3)
                nodes[0] = countX  +  countY   *numNodesX
                nodes[1] = countX+1+ (countY+1)*numNodesX
                nodes[2] = countX  + (countY+1)*numNodesX
                myStructure.ElementCreate(myInterpolationType, nodes)
            else:
                error = True
                errorMsg += "Element shape " + r2DShape +" is invalid. \n"
                return

    myStructure.ElementTotalConvertToInterpolationType()

    mySection = myStructure.SectionCreate("Plane_Stress")
    myStructure.SectionSetThickness(mySection, lZ)
    myStructure.ElementTotalSetSection(mySection)

    Run(myStructure, r2DShape, rTypeOrder)


def Run1D(r1DShape, rTypeOrder):
    global error, errorMsg

    myStructure = nuto.Structure(1)
    myStructure.SetShowTime(False)

    numElementsX = 3

    # create nodes
    numNodesX = numElementsX+1
    deltaX = lX/numElementsX

    for countX in range(0, numNodesX):
        myStructure.NodeCreate(countX, np.array([countX*deltaX]))

    myInterpolationType = myStructure.InterpolationTypeCreate(r1DShape)
    myStructure.InterpolationTypeAdd(myInterpolationType, "Coordinates", "EQUIDISTANT1")
    myStructure.InterpolationTypeAdd(myInterpolationType, "Displacements", rTypeOrder)

    # create elements
    nodes = list(range(2))
    for countX in range(0, numElementsX):
        if (r1DShape == "Truss1D"):
            nodes[0] = countX
            nodes[1] = countX+1
            myStructure.ElementCreate(myInterpolationType, nodes)
        else:
            error = True
            errorMsg += "Element shape " + r2DShape + " is invalid. \n"
            return

    myStructure.ElementTotalConvertToInterpolationType()

    mySection = myStructure.SectionCreate("Truss")
    myStructure.SectionSetArea(mySection, lZ*lY)
    myStructure.ElementTotalSetSection(mySection)

    Run(myStructure, r1DShape, rTypeOrder)

# show the results on the screen
printResult = False

# path in the original source directory and current filename at the end
pathToResultFiles = os.path.join(sys.argv[4], "Results" + os.path.basename(sys.argv[0]))

# remove the extension
fileExt = os.path.splitext(sys.argv[0])[1]
pathToResultFiles = pathToResultFiles.replace(fileExt, '')

if not os.path.exists(pathToResultFiles):
    os.makedirs(pathToResultFiles)

Run3D("Brick3D", "Equidistant1")
Run3D("Brick3D", "Equidistant2")

Run3D("Tetrahedron3D", "Equidistant1")
Run3D("Tetrahedron3D", "Equidistant2")

Run3D("Prism3D", "Equidistant1")
Run3D("Prism3D", "Equidistant2")

Run2D("Quad2D", "Equidistant1")
Run2D("Quad2D", "Equidistant2")
Run2D("Quad2D", "Lobatto2")
Run2D("Quad2D", "Lobatto3")
Run2D("Quad2D", "Lobatto4")

Run2D("Triangle2D", "Equidistant1")
Run2D("Triangle2D", "Equidistant2")
Run2D("Triangle2D", "Equidistant3")
Run2D("Triangle2D", "Equidistant4")

Run1D("Truss1D", "Equidistant1")
Run1D("Truss1D", "Equidistant2")
Run1D("Truss1D", "Equidistant3")
Run1D("Truss1D", "Equidistant4")
Run1D("Truss1D", "Lobatto2")
Run1D("Truss1D", "Lobatto3")
Run1D("Truss1D", "Lobatto4")

if (error):
    print "### \n \n  FAILED \n \n ### "
    print errorMsg
    sys.exit(-1)
else:
    sys.exit(0)
