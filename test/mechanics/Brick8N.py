import nuto
import sys
import os

#call of the test file, e.g.
#/usr/local/bin/python ~/develop/nuto/test/mechanics/Brick8N.py Linux x86_64 ~/develop/nuto/test/mechanics

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
myStructure = nuto.Structure(3)

#create nodes
myNode1 = myStructure.NodeCreate("displacements",nuto.DoubleFullMatrix(3,1,(-1,-1,-1)))
myNode2 = myStructure.NodeCreate("displacements",nuto.DoubleFullMatrix(3,1,(+1,-1,+1)))
myNode3 = myStructure.NodeCreate("displacements",nuto.DoubleFullMatrix(3,1,(-1,+1,-1)))
myNode4 = myStructure.NodeCreate("displacements",nuto.DoubleFullMatrix(3,1,(+1,+1,-1)))
myNode5 = myStructure.NodeCreate("displacements",nuto.DoubleFullMatrix(3,1,(-1,-1,+1)))
myNode6 = myStructure.NodeCreate("displacements",nuto.DoubleFullMatrix(3,1,(+1,-1,+1)))
myNode7 = myStructure.NodeCreate("displacements",nuto.DoubleFullMatrix(3,1,(-1,+1,+1)))
myNode8 = myStructure.NodeCreate("displacements",nuto.DoubleFullMatrix(3,1,(+1,+1,+1)))

#create element
myElement1 = myStructure.ElementCreate("Brick8N",nuto.IntFullMatrix(8,1,(myNode1,myNode2,myNode3,myNode4,myNode5,myNode6,myNode7,myNode8)))

#create constitutive law
myMatLin = myStructure.ConstitutiveLawCreate("LinearElastic")
myStructure.ConstitutiveLawSetYoungsModulus(myMatLin,10)
myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin,0.25)

#assign constitutive law 
#myStructure.ElementSetIntegrationType(myElement1,"3D8NGauss1Ip")
myStructure.ElementSetConstitutiveLaw(myElement1,myMatLin)

#set displacements of right node
myStructure.NodeSetDisplacements(myNode2,nuto.DoubleFullMatrix(3,1,(0.2,0.2,0.2)))
myStructure.NodeSetDisplacements(myNode4,nuto.DoubleFullMatrix(3,1,(0.2,0.2,0.2)))
myStructure.NodeSetDisplacements(myNode6,nuto.DoubleFullMatrix(3,1,(0.2,0.2,0.2)))
myStructure.NodeSetDisplacements(myNode8,nuto.DoubleFullMatrix(3,1,(0.2,0.2,0.2)))

#calculate element stiffness matrix
Ke = nuto.DoubleFullMatrix(0,0)
rowIndex = nuto.IntFullMatrix(0,0)
colIndex = nuto.IntFullMatrix(0,0)
myStructure.ElementStiffness(myElement1,Ke,rowIndex,colIndex)
if (printResult):
    print "Ke"
    Ke.Info()

#correct stiffness matrix
if createResult:
    print pathToResultFiles+"Stiffness.txt"
    Ke.WriteToFile(pathToResultFiles+"Stiffness.txt"," ","#Correct result","  ")
else:
    KeCorrect = nuto.DoubleFullMatrix(24,24)
    KeCorrect.ReadFromFile(pathToResultFiles+"Stiffness.txt",1," ")
    if (printResult):
        print "KeCorrect"
        KeCorrect.Info()
    if ((Ke-KeCorrect).Abs().Max()[0]>1e-8):
        print '[' + system,sys.argv[0] + '] : stiffness is not correct.'
        error = True;

#calculate internal force vector
Fi = nuto.DoubleFullMatrix(0,0)
rowIndex = nuto.IntFullMatrix(0,0)
myStructure.ElementGradientInternalPotential(myElement1,Fi,rowIndex)

#correct resforce vector
if createResult:
    print pathToResultFiles+"Internalforce.txt"
    Fi.WriteToFile(pathToResultFiles+"Internalforce.txt"," ","#Correct result","  ")
else:
    FiCorrect = nuto.DoubleFullMatrix(24,1)
    FiCorrect.ReadFromFile(pathToResultFiles+"Internalforce.txt",1," ")
    if (printResult):
        print "FiCorrect"
        FiCorrect.Info()
    if ((Fi-FiCorrect).Abs().Max()[0]>1e-8):
        print '[' + system,sys.argv[0] + '] : internal force is not correct.'
        error = True;

#check stiffness with internal force vector
prevDisp = nuto.DoubleFullMatrix(0,0) 

delta=1e-4;
KeApprox = nuto.DoubleFullMatrix(Ke.GetNumRows(),Ke.GetNumColumns())
curColumn=0
for theNode in range(0, Ke.GetNumColumns()/myStructure.GetDimension()):
   for theDof in range(0,myStructure.GetDimension()):
      Fi_new = nuto.DoubleFullMatrix(0,0)
      myStructure.NodeGetDisplacements(theNode,prevDisp)
      prevDisp.AddValue(theDof,0,delta)
      myStructure.NodeSetDisplacements(theNode,prevDisp)
      myStructure.ElementGradientInternalPotential(myElement1,Fi_new,rowIndex)
      prevDisp.AddValue(theDof,0,-delta)
      myStructure.NodeSetDisplacements(theNode,prevDisp)
      KeApprox.SetBlock ( 0, curColumn, (Fi_new-Fi)*(1./delta))
      curColumn+=1

if (printResult):
    print "KeApprox with central differences"
    KeApprox.Info()


#check stiffness with internal force vector
if ((KeApprox-Ke).Abs().Max()[0]>1e-8):
        print '[' + system,sys.argv[0] + '] : stiffness matrix via central differences and resforces not correct.'
        error = True;

#calculate engineering strain of myelement1 at all integration points
#the size the matrix is not important and reallocated within the procedure
EngineeringStrain = nuto.DoubleFullMatrix(6,1)
myStructure.ElementGetEngineeringStrain(myElement1, EngineeringStrain)

#correct strain
EngineeringStrainCorrect = nuto.DoubleFullMatrix(6,8,(
0.1,0,0,0.1,0,0.1,
0.1,0,0,0.1,0,0.1,
0.1,0,0,0.1,0,0.1,
0.1,0,0,0.1,0,0.1,
0.1,0,0,0.1,0,0.1,
0.1,0,0,0.1,0,0.1,
0.1,0,0,0.1,0,0.1,
0.1,0,0,0.1,0,0.1
))

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
EngineeringStressCorrect = nuto.DoubleFullMatrix(6,8,(
1.2,0.4,0.4,0.4,0,0.4,
1.2,0.4,0.4,0.4,0,0.4,
1.2,0.4,0.4,0.4,0,0.4,
1.2,0.4,0.4,0.4,0,0.4,
1.2,0.4,0.4,0.4,0,0.4,
1.2,0.4,0.4,0.4,0,0.4,
1.2,0.4,0.4,0.4,0,0.4,
1.2,0.4,0.4,0.4,0,0.4
))

if (printResult):
    print "EngineeringStressCorrect"
    EngineeringStressCorrect.Info()
    print "EngineeringStress"
    EngineeringStress.Info()

if ((EngineeringStress-EngineeringStressCorrect).Abs().Max()[0]>1e-8):
        print '[' + system,sys.argv[0] + '] : stress is not correct.'
        error = True;

if (error):
    sys.exit(-1)
else:
    sys.exit(0)
