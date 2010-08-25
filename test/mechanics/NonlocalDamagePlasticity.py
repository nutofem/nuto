import nuto
import sys
import os

#call of the test file, e.g.
#/usr/local/bin/python ~/develop/nuto/test/mechanics/NonlocalDamagePlasticity.py Linux x86_64 ~/develop/nuto/test/mechanics

#if set to true, the result will be generated (for later use in the test routine)
#otherwise, the current result will be compared to the stored result
#createResult = True
createResult = False

#show the results on the screen
printResult = True

#system name and processor
system = sys.argv[1]+sys.argv[2]

#path in the original source directory and current filename at the end
pathToResultFiles = os.path.join(sys.argv[3],"results",system,os.path.basename(sys.argv[0]))

#remove the extension
fileExt = os.path.splitext(sys.argv[0])[1]
pathToResultFiles = pathToResultFiles.replace(fileExt,'')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% no error in file, modified, if error is detected              %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error = False

#create structure
myStructure = nuto.Structure(2)

#create nodes
node1 = myStructure.NodeCreate("displacements",nuto.DoubleFullMatrix(2,1,(0,0)))
node2 = myStructure.NodeCreate("displacements",nuto.DoubleFullMatrix(2,1,(1,0)))
node3 = myStructure.NodeCreate("displacements",nuto.DoubleFullMatrix(2,1,(2,0)))
node4 = myStructure.NodeCreate("displacements",nuto.DoubleFullMatrix(2,1,(0,1)))
node5 = myStructure.NodeCreate("displacements",nuto.DoubleFullMatrix(2,1,(1,1)))
node6 = myStructure.NodeCreate("displacements",nuto.DoubleFullMatrix(2,1,(2,1)))
node7 = myStructure.NodeCreate("displacements",nuto.DoubleFullMatrix(2,1,(0,2)))
node8 = myStructure.NodeCreate("displacements",nuto.DoubleFullMatrix(2,1,(1,2)))
node9 = myStructure.NodeCreate("displacements",nuto.DoubleFullMatrix(2,1,(2,2)))

#create element
element1 = myStructure.ElementCreate("PLANE2D4N",nuto.IntFullMatrix(4,1,(node1,node2,node5,node4)),"ConstitutiveLawIpNonlocal","StaticDataNonlocal")
myStructure.ElementSetIntegrationType(element1,"2D4NGauss1Ip","StaticDataNonlocal")
element2 = myStructure.ElementCreate("PLANE2D4N",nuto.IntFullMatrix(4,1,(node2,node3,node6,node5)),"ConstitutiveLawIpNonlocal","StaticDataNonlocal")
myStructure.ElementSetIntegrationType(element2,"2D4NGauss4Ip","StaticDataNonlocal")
element3 = myStructure.ElementCreate("PLANE2D4N",nuto.IntFullMatrix(4,1,(node4,node5,node8,node7)),"ConstitutiveLawIpNonlocal","StaticDataNonlocal")
myStructure.ElementSetIntegrationType(element3,"2D4NGauss1Ip","StaticDataNonlocal")
element4 = myStructure.ElementCreate("PLANE2D4N",nuto.IntFullMatrix(4,1,(node5,node6,node9,node8)),"ConstitutiveLawIpNonlocal","StaticDataNonlocal")
myStructure.ElementSetIntegrationType(element4,"2D4NGauss4Ip","StaticDataNonlocal")

#create constitutive law
myMatDamage = myStructure.ConstitutiveLawCreate("NonlocalDamagePlasticity")
myStructure.ConstitutiveLawSetYoungsModulus(myMatDamage,9)
myStructure.ConstitutiveLawSetPoissonsRatio(myMatDamage,0.25)
myStructure.ConstitutiveLawSetNonlocalRadius(myMatDamage,0.7)
myStructure.ConstitutiveLawSetTensileStrength(myMatDamage,2)
myStructure.ConstitutiveLawSetCompressiveStrength(myMatDamage,20)
myStructure.ConstitutiveLawSetBiaxialCompressiveStrength(myMatDamage,25)
myStructure.ConstitutiveLawSetFractureEnergy(myMatDamage,0.2)

myMatLin = myStructure.ConstitutiveLawCreate("LinearElastic")
myStructure.ConstitutiveLawSetYoungsModulus(myMatLin,10)
myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin,0.25)

#create section
mySection = myStructure.SectionCreate("Plane_Strain")
myStructure.SectionSetThickness(mySection,0.5)

##assign constitutive law 
myStructure.ElementTotalSetSection(mySection)
myStructure.ElementTotalSetConstitutiveLaw(myMatDamage)

#Build nonlocal elements
myStructure.BuildNonlocalData(myMatDamage)

#visualize results (nonlocal weights)
myStructure.AddVisualizationComponentNonlocalWeights(element1,0)
myStructure.AddVisualizationComponentNonlocalWeights(element2,0)
myStructure.AddVisualizationComponentNonlocalWeights(element2,1)
myStructure.AddVisualizationComponentNonlocalWeights(element2,2)
myStructure.AddVisualizationComponentNonlocalWeights(element2,3)

#calculate linear elastic matrix
stiffnessMatrix = nuto.DoubleSparseMatrixCSRGeneral(0,0)
dispForceVector = nuto.DoubleFullMatrix(0,0)

myStructure.ElementTotalUpdateTmpStaticData()
myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector)
stiffnessMatrix.RemoveZeroEntries(0,1e-14)

fullStiffnessMatrixElastic = nuto.DoubleFullMatrix(stiffnessMatrix)
if printResult:
    print "stiffnessMatrix elastic"
    fullStiffnessMatrixElastic.Info()

displacements = nuto.DoubleFullMatrix(0,0) 
dependentDofs = nuto.DoubleFullMatrix(0,0)
intForce	  = nuto.DoubleFullMatrix(0,0)
intForce2	  = nuto.DoubleFullMatrix(0,0)

#check the stiffness three times
#loadstep 0 : uniform plastic loading
#loadstep 1 : unloading to zero
#loadstep 2 : nonuniform loading, some elements unloading
for theLoadStep in range(0,3):
    #apply displacements
    if theLoadStep==0:
        rightDisp = 0.5
    elif theLoadStep==1:
        rightDisp = 0.0
    else:
        rightDisp = 0.6

    matrixRightDisp = nuto.DoubleFullMatrix(2,1)
    matrixRightDisp.SetValue(0,0,rightDisp)
    matrixRightDisp.SetValue(1,0,0.)

    myStructure.NodeSetDisplacements(node3,matrixRightDisp)
    myStructure.NodeSetDisplacements(node6,matrixRightDisp)
    myStructure.NodeSetDisplacements(node9,matrixRightDisp)

    matrixCenterDisp = nuto.DoubleFullMatrix(2,1)
    if theLoadStep!=2:
        matrixCenterDisp.SetValue(0,0,0.5*rightDisp)
    else:
        matrixCenterDisp.SetValue(0,0,0.4*rightDisp)
    matrixCenterDisp.SetValue(1,0,0.)

    myStructure.NodeSetDisplacements(node2,matrixCenterDisp)
    myStructure.NodeSetDisplacements(node5,matrixCenterDisp)
    myStructure.NodeSetDisplacements(node8,matrixCenterDisp)

    matrixLeftDisp = nuto.DoubleFullMatrix(2,1)
    matrixLeftDisp.SetValue(0,0,0.0)
    matrixLeftDisp.SetValue(1,0,0.)

    myStructure.NodeSetDisplacements(node1,matrixLeftDisp)
    myStructure.NodeSetDisplacements(node4,matrixLeftDisp)
    myStructure.NodeSetDisplacements(node7,matrixLeftDisp)

    #nuto.FullMatrix<double>matrixLeftDispNode1(2,1)
    #matrixLeftDispNode1.SetValue(0,0,0.0)
    #matrixLeftDispNode1.SetValue(1,0,0.)
    #myStructure.NodeSetDisplacements(node1,matrixLeftDispNode1)

    myStructure.ElementTotalUpdateTmpStaticData()

    myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector)
    stiffnessMatrix.RemoveZeroEntries(0,1e-14)

    fullStiffnessMatrix = nuto.DoubleFullMatrix(stiffnessMatrix)

    if printResult:
        print "fullStiffnessMatrix analytic"
        fullStiffnessMatrix.Info()

    myStructure.NodeExtractDofValues(displacements,dependentDofs)
    myStructure.BuildGlobalGradientInternalPotentialVector(intForce)

    delta = 1e-8
    stiffnessMatrixCD = nuto.DoubleFullMatrix(displacements.GetNumRows(),displacements.GetNumRows())

    #check with central differences
    for count2 in range(0,displacements.GetNumRows()):
        displacements.AddValue(count2,0,delta)
        myStructure.NodeMergeActiveDofValues(displacements)
        myStructure.ElementTotalUpdateTmpStaticData()
        myStructure.BuildGlobalGradientInternalPotentialVector(intForce2)
        stiffnessMatrixCD.SetColumn(count2,(intForce2-intForce)*(1/delta))
        displacements.AddValue(count2,0,-delta)

    if printResult:
        print "stiffness CD"
        stiffnessMatrixCD.Info()
    maxerror=(fullStiffnessMatrix-stiffnessMatrixCD).Abs().Max()
    if (theLoadStep==0):
        if printResult:
            print "max difference in stiffness matrix for uniform plastic loading ", maxerror 
        if (maxerror[0]>1e-6):
            print '[' + system,sys.argv[0] + '] : stiffness for uniform plastic loading is not correct.'
            error = True;
    elif (theLoadStep==1):
        if printResult:
            print "max difference in stiffness matrix for unloading " , maxerror 
        if (maxerror[0]>1e-6):
            print '[' + system,sys.argv[0] + '] : stiffness after unloading is not correct.'
            error = True;
        omega = fullStiffnessMatrix.GetValue(0,0)/fullStiffnessMatrixElastic.GetValue(0,0)
        maxerror2 = (fullStiffnessMatrix - fullStiffnessMatrixElastic*omega).Abs().Max()
        if printResult:
            print "difference in stiffness matrix for unloading and scaled elastic matrix"
            (fullStiffnessMatrixElastic*omega-fullStiffnessMatrix).Info()
        if printResult:
            print "max difference in stiffness matrix for unloading and scaled elastic matrix " , maxerror2 
        if (maxerror2[0]>1e-6):
            print '[' + system,sys.argv[0] + '] : stiffness matrix for unloading and scaled elastic matrix are not identical.'
            error = True;
    else:
        if printResult:
            print "max difference in stiffness matrix for nonuniform plastic loading/unloading " , maxerror
        if (maxerror[0]>1e-6):
            print '[' + system,sys.argv[0] + '] :stiffness matrix for nonuniform plastic loading/unloading is not correct.'
            error = True;

    #update the structure, and then recalculate stiffness
    myStructure.ElementTotalUpdateStaticData()


myStructure.AddVisualizationComponentDisplacements()
myStructure.AddVisualizationComponentEngineeringStrain()
myStructure.AddVisualizationComponentEngineeringStress()
myStructure.AddVisualizationComponentDamage()
myStructure.AddVisualizationComponentEngineeringPlasticStrain()
myStructure.ExportVtkDataFile("NonlocalDamagePlasticityModel.vtk")

if (error):
    sys.exit(-1)
else:
    sys.exit(0)
