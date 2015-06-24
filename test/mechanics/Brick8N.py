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
myNode1 = myStructure.NodeCreate(nuto.DoubleFullVector((-1.,-1.,-1.)))
myNode2 = myStructure.NodeCreate(nuto.DoubleFullVector((+1.,-1.,-1.)))
myNode3 = myStructure.NodeCreate(nuto.DoubleFullVector((+1.,+1.,-1.)))
myNode4 = myStructure.NodeCreate(nuto.DoubleFullVector((-1.,+1.,-1.)))
myNode5 = myStructure.NodeCreate(nuto.DoubleFullVector((-1.,-1.,+1.)))
myNode6 = myStructure.NodeCreate(nuto.DoubleFullVector((+1.,-1.,+1.)))
myNode7 = myStructure.NodeCreate(nuto.DoubleFullVector((+1.,+1.,+1.)))
myNode8 = myStructure.NodeCreate(nuto.DoubleFullVector((-1.,+1.,+1.)))

#create interpolation type
myInterpolationType = myStructure.InterpolationTypeCreate("Brick3D");
myStructure.InterpolationTypeAdd(myInterpolationType, "coordinates", "equidistant1");
myStructure.InterpolationTypeAdd(myInterpolationType, "displacements", "equidistant1");

#create element
myElement1 = myStructure.ElementCreate(myInterpolationType,nuto.IntFullVector((myNode1,myNode2,myNode3,myNode4,myNode5,myNode6,myNode7,myNode8)))
myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);

#create constitutive law
myMatLin = myStructure.ConstitutiveLawCreate("LinearElasticEngineeringStress")
myStructure.ConstitutiveLawSetYoungsModulus(myMatLin,10)
myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin,0.25)

#create section
mySection = myStructure.SectionCreate("Volume")

#assign constitutive law 
#myStructure.ElementSetIntegrationType(myElement1,"3D8NGauss1Ip")
myStructure.ElementSetConstitutiveLaw(myElement1,myMatLin)
myStructure.ElementSetSection(myElement1,mySection)

#make group of boundary nodes
groupBoundaryNodes = myStructure.GroupCreate("Nodes")
myStructure.GroupAddNode(groupBoundaryNodes,myNode1)
myStructure.GroupAddNode(groupBoundaryNodes,myNode4)
myStructure.GroupAddNode(groupBoundaryNodes,myNode5)
myStructure.GroupAddNode(groupBoundaryNodes,myNode8)

#make group of boundary elements (in this case it is just one
groupBoundaryElements = myStructure.GroupCreate("Elements")
myStructure.GroupAddElementsFromNodes(groupBoundaryElements,groupBoundaryNodes,False);

#create surface loads (0 - pressure on X, 1-const-direction Y)
myStructure.LoadSurfacePressureCreate3D(0, groupBoundaryElements, groupBoundaryNodes, 2.);
myStructure.LoadSurfaceConstDirectionCreate3D(1, groupBoundaryElements, groupBoundaryNodes, nuto.DoubleFullVector((0.,5.,0.)))

#set displacements of right node
myStructure.NodeSetDisplacements(myNode2,nuto.DoubleFullVector((0.2,0.2,0.2)))
myStructure.NodeSetDisplacements(myNode3,nuto.DoubleFullVector((0.2,0.2,0.2)))
myStructure.NodeSetDisplacements(myNode6,nuto.DoubleFullVector((0.2,0.2,0.2)))
myStructure.NodeSetDisplacements(myNode7,nuto.DoubleFullVector((0.2,0.2,0.2)))

#calculate element stiffness matrix
Ke = nuto.DoubleFullMatrix(0,0)
rowIndex = nuto.IntFullVector(0)
colIndex = nuto.IntFullVector(0)
print "test1"
myStructure.ElementStiffness(myElement1,Ke,rowIndex,colIndex)
print "test2"
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
    if ((Ke-KeCorrect).Abs().Max()>1e-8):
        print '[' + system,sys.argv[0] + '] : stiffness is not correct.'
        error = True;

#calculate internal force vector (this is only due to the prescribed displacements, not in equilibrium with external forces
Fi = nuto.DoubleFullVector(0)
rowIndex = nuto.IntFullVector(0)
myStructure.ElementGradientInternalPotential(myElement1,Fi,rowIndex)

if (printResult):
    print "Fi"
    Fi.Info()

#correct resforce vector
if createResult:
    print pathToResultFiles+"Internalforce.txt"
    Fi.WriteToFile(pathToResultFiles+"Internalforce.txt"," ","#Correct result","  ")
else:
    FiCorrect = nuto.DoubleFullVector(24)
    FiCorrect.ReadFromFile(pathToResultFiles+"Internalforce.txt",1," ")
    if (printResult):
        print "FiCorrect"
        FiCorrect.Info()
    if ((Fi-FiCorrect).Abs().Max()>1e-8):
        print '[' + system,sys.argv[0] + '] : internal force is not correct.'
        error = True;

#check stiffness with internal force vector
prevDisp = nuto.DoubleFullVector(0) 

delta=1e-4;
KeApprox = nuto.DoubleFullMatrix(Ke.GetNumRows(),Ke.GetNumColumns())
curColumn=0
for theNode in range(0, Ke.GetNumColumns()/myStructure.GetDimension()):
   for theDof in range(0,myStructure.GetDimension()):
      Fi_new = nuto.DoubleFullVector(0)
      myStructure.NodeGetDisplacements(theNode,prevDisp)
      prevDisp.AddValue(theDof,0,delta)
      myStructure.NodeSetDisplacements(theNode,prevDisp)
      myStructure.ElementGradientInternalPotential(myElement1,Fi_new,rowIndex)
      prevDisp.AddValue(theDof,0,-delta)
      myStructure.NodeSetDisplacements(theNode,prevDisp)
      KeApprox.SetColumn(curColumn, (Fi_new-Fi)*(1./delta))
      curColumn+=1

if (printResult):
    print "KeApprox with central differences"
    KeApprox.Info()


#check stiffness with internal force vector
if ((KeApprox-Ke).Abs().Max()>1e-8):
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

if ((EngineeringStrain-EngineeringStrainCorrect).Abs().Max()>1e-8):
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

if ((EngineeringStress-EngineeringStressCorrect).Abs().Max()>1e-8):
        print '[' + system,sys.argv[0] + '] : stress is not correct.'
        error = True;

        
#calculate external force vector for the first load case (pressure)
Fe1 = nuto.DoubleFullVector(0)
Fe2 = nuto.DoubleFullVector(0) # no displacements are constrained so this vector is empty
myStructure.BuildGlobalExternalLoadVector(0,Fe1,Fe2)

if (printResult):
    print "Fe1 for pressure load"
    Fe1.Info()

#correct external force for pressure load vector (sum up the load in x direction eveything else should be zero
sumX = Fe1.Sum()
if (abs(sumX-8.)>1e-8):
        print '[' + system,sys.argv[0] + '] : pressure load is not correct.'
        error = True;
        
#calculate external force vector for the second load cases (constDirection)
myStructure.BuildGlobalExternalLoadVector(1,Fe1,Fe2)

if (printResult):
    print "Fe1 const direction load"
    Fe1.Info()

#correct external force for pressure load vector (sum up the load in x direction eveything else should be zero
sumY = Fe1.Sum()
if (abs(sumY-20.)>1e-8):
        print '[' + system,sys.argv[0] + '] : const direction load is not correct.'
        error = True; 
        
if (error):
    sys.exit(-1)
else:
    sys.exit(0)
