import nuto
import sys
import os

#call of the test file, e.g.
#/usr/local/bin/python ~/develop/nuto/test/mechanics/Truss1D2N.py Linux x86_64 ~/develop/nuto/test/mechanics

#if set to true, the result will be generated (for later use in the test routine)
#otherwise, the current result will be compared to the stored result
#createResult = True
createResult = False

#show the results on the screen
printResult = False

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
myStructure = nuto.Structure(1)

#create nodes
myNode1 = myStructure.NodeCreate("displacements",nuto.DoubleFullMatrix(1,1,(1,)))
myNode2 = myStructure.NodeCreate("displacements",nuto.DoubleFullMatrix(1,1,(6,)))

#create element
myElement1 = myStructure.ElementCreate("Truss1D2N",nuto.IntFullMatrix(2,1,(myNode1,myNode2)))

#create constitutive law
myMatLin = myStructure.ConstitutiveLawCreate("LinearElastic")
myStructure.ConstitutiveLawSetYoungsModulus(myMatLin,10)
myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin,0.1)
myStructure.ConstitutiveLawSetDensity(myMatLin,0.3)

#create section
mySection1 = myStructure.SectionCreate("TRUSS")
myStructure.SectionSetArea(mySection1,0.01)

#assign constitutive law and section to elements
#myStructure.ElementSetIntegrationType("myElement1","1D2NGauss2Ip")
#variable integration type, for the 2-node truss element it makes no sense, 
#but it demonstrates the idea
myStructure.ElementSetIntegrationType(myElement1,"1D2NConst3Ip","NOIPDATA")
myStructure.ElementSetConstitutiveLaw(myElement1,myMatLin)
myStructure.ElementSetSection(myElement1,mySection1)

#set displacements of right node
myStructure.NodeSetDisplacements(myNode2,nuto.DoubleFullMatrix(1,1,(0.5,)))

#calculate element stiffness matrix
Ke = nuto.DoubleFullMatrix(0,0)
rowIndex = nuto.IntFullMatrix(0,0)
colIndex = nuto.IntFullMatrix(0,0)
myStructure.ElementStiffness(myElement1,Ke,rowIndex,colIndex)

#correct stiffness matrix
KeCorrect = nuto.DoubleFullMatrix(2,2,(0.02, -0.02, -0.02, 0.02))

if (printResult):
    print "KeCorrect"
    KeCorrect.Info()
    print "Ke"
    Ke.Info()

if ((Ke-KeCorrect).Abs().Max()[0]>1e-8):
        print '[' + system,sys.argv[0] + '] : stiffness is not correct.'
        error = True;



#calculate internal force vector
Fi = nuto.DoubleFullMatrix(0,0)
rowIndex = nuto.IntFullMatrix(0,0)
myStructure.ElementGradientInternalPotential(myElement1,Fi,rowIndex)

#correct internal force vector
FiCorrect = nuto.DoubleFullMatrix(2,1,(-0.01, 0.01))
if (printResult):
    print "FiCorrect"
    FiCorrect.Info()
    print "Fi"
    Fi.Info()

if ((Fi-FiCorrect).Abs().Max()[0]>1e-8):
        print '[' + system,sys.argv[0] + '] : resforce is not correct.'
        error = True;


#check stiffness with internal force vector
prevDisp = nuto.DoubleFullMatrix(0,0) 

delta=1e-4;
Fi1 = nuto.DoubleFullMatrix(0,0)
myStructure.NodeGetDisplacements(myNode1,prevDisp)
prevDisp.AddValue(0,0,delta)
myStructure.NodeSetDisplacements(myNode1,prevDisp)
myStructure.ElementGradientInternalPotential(myElement1,Fi1,rowIndex)
prevDisp.AddValue(0,0,-delta)
myStructure.NodeSetDisplacements(myNode1,prevDisp)
#if (printResult):
#    print "Fi1"
#    Fi1.Info()

Fi2 = nuto.DoubleFullMatrix(0,0)
myStructure.NodeGetDisplacements(myNode2,prevDisp)
prevDisp.AddValue(0,0,delta)
myStructure.NodeSetDisplacements(myNode2,prevDisp)
myStructure.ElementGradientInternalPotential(myElement1,Fi2,rowIndex)
prevDisp.AddValue(0,0,-delta)
myStructure.NodeSetDisplacements(myNode2,prevDisp)
#if (printResult):
#    print "Fi2"
#    Fi2.Info()

KeApprox = nuto.DoubleFullMatrix(2,2)
KeApprox.SetBlock ( 0, 0, (Fi1-Fi)*(1./delta))
KeApprox.SetBlock ( 0, 1, (Fi2-Fi)*(1./delta))
if (printResult):
    print "KeApprox"
    KeApprox.Info()


#correct approximated stiffness via central differences of the internal force vector
if ((KeApprox-KeCorrect).Abs().Max()[0]>1e-8):
        print '[' + system,sys.argv[0] + '] : stiffness matrix via central differences and resforces not correct.'
        error = True;



#calculate engineering strain of myelement1 at all integration points
EngineeringStrain = nuto.DoubleFullMatrix(6,3)
myStructure.ElementGetEngineeringStrain(myElement1, EngineeringStrain)

#correct strain
EngineeringStrainCorrect = nuto.DoubleFullMatrix(6,3,(0.1,-0.01,-0.01,0,0,0,
0.1,-0.01,-0.01,0,0,0,
0.1,-0.01,-0.01,0,0,0))

if (printResult):
    print "EngineeringStrainCorrect"
    EngineeringStrainCorrect.Info()
    print "EngineeringStrain"
    EngineeringStrain.Info()

if ((EngineeringStrain-EngineeringStrainCorrect).Abs().Max()[0]>1e-8):
        print '[' + system,sys.argv[0] + '] : strain is not correct.'
        error = True;

#calculate engineering strain of myelement1 at all integration points
EngineeringStress = nuto.DoubleFullMatrix(6,3)
myStructure.ElementGetEngineeringStress(myElement1, EngineeringStress)
#correct stress
EngineeringStressCorrect = nuto.DoubleFullMatrix(6,3,(1,0,0,0,0,0,
1,0,0,0,0,0,
1,0,0,0,0,0))

if (printResult):
    print "EngineeringStressCorrect"
    EngineeringStressCorrect.Info()
    print "EngineeringStress"
    EngineeringStress.Info()

if ((EngineeringStress-EngineeringStressCorrect).Abs().Max()[0]>1e-8):
        print '[' + system,sys.argv[0] + '] : stress is not correct.'
        error = True;

# calculate mass matrix of myelement1 (requires higher integration order)
myStructure.ElementSetIntegrationType(myElement1,"1D2NGauss2Ip","NOIPDATA")

#calculate element stiffness matrix
Me = nuto.DoubleFullMatrix(0,0)
rowIndex = nuto.IntFullMatrix(0,0)
colIndex = nuto.IntFullMatrix(0,0)
myStructure.ElementCoefficientMatrix_2(myElement1,Me,rowIndex,colIndex)

#correct stiffness matrix
MeCorrect = nuto.DoubleFullMatrix(2,2,(5e-3,2.5e-3,2.5e-3,5e-3))

if (printResult):
    print "MeCorrect"
    MeCorrect.Info()
    print "Me"
    Me.Info()
    
if ((Me-MeCorrect).Abs().Max()[0]>1e-8):
    print '[' + system,sys.argv[0] + '] : mass is not correct.'
    error = True;

if (error):
    sys.exit(-1)
else:
    sys.exit(0)
