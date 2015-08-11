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
node1 = myStructure.NodeCreate(nuto.DoubleFullVector((0,0)))
node2 = myStructure.NodeCreate(nuto.DoubleFullVector((1,0)))
node3 = myStructure.NodeCreate(nuto.DoubleFullVector((2,0)))
node4 = myStructure.NodeCreate(nuto.DoubleFullVector((0,1)))
node5 = myStructure.NodeCreate(nuto.DoubleFullVector((1,1)))
node6 = myStructure.NodeCreate(nuto.DoubleFullVector((2,1)))
node7 = myStructure.NodeCreate(nuto.DoubleFullVector((0,2)))
node8 = myStructure.NodeCreate(nuto.DoubleFullVector((1,2)))
node9 = myStructure.NodeCreate(nuto.DoubleFullVector((2,2)))

it1IP = myStructure.InterpolationTypeCreate("Quad2D")
myStructure.InterpolationTypeAdd(it1IP, "coordinates", "equidistant1")
myStructure.InterpolationTypeAdd(it1IP, "displacements", "equidistant1")

it4IP = myStructure.InterpolationTypeCreate("Quad2D")
myStructure.InterpolationTypeAdd(it4IP, "coordinates", "equidistant1")
myStructure.InterpolationTypeAdd(it4IP, "displacements", "equidistant1")


#create element
element1 = myStructure.ElementCreate(it1IP,nuto.IntFullVector((node1,node2,node5,node4)),"ConstitutiveLawIpNonlocal","StaticDataNonlocal")
element2 = myStructure.ElementCreate(it4IP,nuto.IntFullVector((node2,node3,node6,node5)),"ConstitutiveLawIpNonlocal","StaticDataNonlocal")
element3 = myStructure.ElementCreate(it1IP,nuto.IntFullVector((node4,node5,node8,node7)),"ConstitutiveLawIpNonlocal","StaticDataNonlocal")
element4 = myStructure.ElementCreate(it4IP,nuto.IntFullVector((node5,node6,node9,node8)),"ConstitutiveLawIpNonlocal","StaticDataNonlocal")

myStructure.InterpolationTypeSetIntegrationType(it1IP, "2D4NGauss1Ip","StaticDataNonlocal")
myStructure.InterpolationTypeSetIntegrationType(it1IP, "2D4NGauss4Ip","StaticDataNonlocal")

myStructure.ElementTotalConvertToInterpolationType()

#create constitutive law
myMatDamage = myStructure.ConstitutiveLawCreate("NonlocalDamagePlasticityEngineeringStress")
myStructure.ConstitutiveLawSetParameterDouble(myMatDamage,"YoungsModulus",9)
myStructure.ConstitutiveLawSetParameterDouble(myMatDamage,"PoissonsRatio",0.25)
myStructure.ConstitutiveLawSetParameterDouble(myMatDamage,"NonlocalRadius",0.7)
myStructure.ConstitutiveLawSetParameterDouble(myMatDamage,"TensileStrength",2)
myStructure.ConstitutiveLawSetParameterDouble(myMatDamage,"CompressiveStrength",20)
myStructure.ConstitutiveLawSetParameterDouble(myMatDamage,"BiaxialCompressiveStrength",25)
myStructure.ConstitutiveLawSetParameterDouble(myMatDamage,"FractureEnergy",0.2)

myMatLin = myStructure.ConstitutiveLawCreate("LinearElasticEngineeringStress")
myStructure.ConstitutiveLawSetParameterDouble(myMatLin,"YoungsModulus",9)
myStructure.ConstitutiveLawSetParameterDouble(myMatLin,"PoissonsRatio",0.25)

#create section
mySection = myStructure.SectionCreate("Plane_Strain")
myStructure.SectionSetThickness(mySection,0.5)

##assign constitutive law 
myStructure.ElementTotalSetSection(mySection)
myStructure.ElementTotalSetConstitutiveLaw(myMatDamage)
#myStructure.ElementTotalSetConstitutiveLaw(myMatLin)

#Build nonlocal elements
myStructure.BuildNonlocalData(myMatDamage)

#Calculate maximum independent sets for parallelization (openmp)
myStructure.CalculateMaximumIndependentSets();

#visualize results (nonlocal weights)
myStructure.AddVisualizationComponentNonlocalWeights(element1,0)
myStructure.AddVisualizationComponentNonlocalWeights(element2,0)
myStructure.AddVisualizationComponentNonlocalWeights(element2,1)
myStructure.AddVisualizationComponentNonlocalWeights(element2,2)
myStructure.AddVisualizationComponentNonlocalWeights(element2,3)

#calculate linear elastic matrix
stiffnessMatrix = nuto.DoubleSparseMatrixCSRVector2General(0,0)
dispForceVector = nuto.DoubleFullVector(0)

myStructure.ElementTotalUpdateTmpStaticData()
myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector)
stiffnessMatrix.RemoveZeroEntries(0,1e-14)

fullStiffnessMatrixElastic = nuto.DoubleFullMatrix(stiffnessMatrix)
if printResult:
    print "stiffnessMatrix elastic"
    fullStiffnessMatrixElastic.Info()

displacements = nuto.DoubleFullVector(0) 
dependentDofs = nuto.DoubleFullVector(0)
intForce	  = nuto.DoubleFullVector(0)
intForce2	  = nuto.DoubleFullVector(0)

#check the stiffness three times
#loadstep 0 : uniform plastic loading
#loadstep 1 : unloading to zero
#loadstep 2 : nonuniform loading, some elements unloading
for theLoadStep in range(0,1):
    #apply displacements
    if theLoadStep==0:
        rightDisp = 0.5
    elif theLoadStep==1:
        rightDisp = 0.0
    else:
        rightDisp = 0.6

    print "test1"
    matrixRightDisp = nuto.DoubleFullVector(2)
    matrixRightDisp.SetValue(0,0,rightDisp)
    matrixRightDisp.SetValue(1,0,0.)

    print "test2"
    myStructure.NodeSetDisplacements(node3,matrixRightDisp)
    myStructure.NodeSetDisplacements(node6,matrixRightDisp)
    myStructure.NodeSetDisplacements(node9,matrixRightDisp)

    print "test3"
    matrixCenterDisp = nuto.DoubleFullVector(2)
    if theLoadStep!=2:
        matrixCenterDisp.SetValue(0,0,0.5*rightDisp)
    else:
        matrixCenterDisp.SetValue(0,0,0.4*rightDisp)
    matrixCenterDisp.SetValue(1,0,0.)

    myStructure.NodeSetDisplacements(node2,matrixCenterDisp)
    myStructure.NodeSetDisplacements(node5,matrixCenterDisp)
    myStructure.NodeSetDisplacements(node8,matrixCenterDisp)

    print "test4"
    matrixLeftDisp = nuto.DoubleFullVector(2)
    matrixLeftDisp.SetValue(0,0,0.0)
    matrixLeftDisp.SetValue(1,0,0.)

    myStructure.NodeSetDisplacements(node1,matrixLeftDisp)
    myStructure.NodeSetDisplacements(node4,matrixLeftDisp)
    myStructure.NodeSetDisplacements(node7,matrixLeftDisp)

    #nuto.FullMatrix<double>matrixLeftDispNode1(2,1)
    #matrixLeftDispNode1.SetValue(0,0,0.0)
    #matrixLeftDispNode1.SetValue(1,0,0.)
    #myStructure.NodeSetDisplacements(node1,matrixLeftDispNode1)

    print "test5"
    myStructure.ElementTotalUpdateTmpStaticData()

    print "test6"
    myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector)
    stiffnessMatrix.RemoveZeroEntries(0,1e-14)
    print "test7"

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
        print "fullStiffnessMatrix"
        fullStiffnessMatrix.Info()
        print "stiffness CD"
        stiffnessMatrixCD.Info()
    maxerror=(fullStiffnessMatrix-stiffnessMatrixCD).Abs().Max()
    if (theLoadStep==0):
        if printResult:
            print "max difference in stiffness matrix for uniform plastic loading ", maxerror 
        if (maxerror>1e-6):
            print '[' + system,sys.argv[0] + '] : stiffness for uniform plastic loading is not correct.'
            error = True;
    elif (theLoadStep==1):
        if printResult:
            print "max difference in stiffness matrix for unloading " , maxerror 
        if (maxerror>1e-6):
            print '[' + system,sys.argv[0] + '] : stiffness after unloading is not correct.'
            error = True;
        omega = fullStiffnessMatrix.GetValue(0,0)/fullStiffnessMatrixElastic.GetValue(0,0)
        maxerror2 = (fullStiffnessMatrix - fullStiffnessMatrixElastic*omega).Abs().Max()
        if printResult:
            print "difference in stiffness matrix for unloading and scaled elastic matrix"
            (fullStiffnessMatrixElastic*omega-fullStiffnessMatrix).Info()
        if printResult:
            print "max difference in stiffness matrix for unloading and scaled elastic matrix " , maxerror2 
        if (maxerror2>1e-6):
            print '[' + system,sys.argv[0] + '] : stiffness matrix for unloading and scaled elastic matrix are not identical.'
            error = True;
    else:
        if printResult:
            print "max difference in stiffness matrix for nonuniform plastic loading/unloading " , maxerror
        if (maxerror>1e-6):
            print '[' + system,sys.argv[0] + '] :stiffness matrix for nonuniform plastic loading/unloading is not correct.'
            error = True;

    #update the structure, and then recalculate stiffness
    myStructure.ElementTotalUpdateStaticData()

myStructure.AddVisualizationComponentDisplacements()
myStructure.AddVisualizationComponentEngineeringStrain()
myStructure.AddVisualizationComponentEngineeringStress()
myStructure.AddVisualizationComponentDamage()
myStructure.AddVisualizationComponentEngineeringPlasticStrain()
myStructure.ExportVtkDataFileElements("NonlocalDamagePlasticityModel.vtk")

if (error):
    sys.exit(-1)
else:
    sys.exit(0)
