import nuto
import sys
import os

#call of the test file, e.g.
#/usr/local/bin/python ~/develop/nuto_phantom_nodes/nuto/test/mechanics/Plane2D4N.py Linux x86_64 ~/develop/nuto_phantom_nodes/nuto/test/mechanics

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

import nuto
import sys
import os

createResult = True
printResult = True

from time import time

# definitions
YoungsModulus = 1.e5
PoissonsRatio = 0.2

Force = 5

#ElementType = "PLANE2D3N"
ElementType = "PLANE2D4N"

StressState = "XX"
#StressState = "YY"
#StressState = "XY"

EnableDisplacementControl = True
#EnableDisplacementControl = False
BoundaryDisplacement = 1.0

#create structure
numDim=2
myStructure = nuto.Structure(numDim)

# create material law
myMatLin = myStructure.ConstitutiveLawCreate("LinearElastic")
myStructure.ConstitutiveLawSetYoungsModulus(myMatLin, YoungsModulus)
myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin, PoissonsRatio)

#create section
mySection = myStructure.SectionCreate("Plane_Strain")
myStructure.SectionSetThickness(mySection,1)

numNodes=8
myStructure.NodesCreate("displacements", nuto.DoubleFullMatrix(2,numNodes,(	 0 ,  0 ,
																			10 ,  0 ,
																			 2 ,  2 ,
																			 8 ,  3 ,
																			 4 ,  7 ,
																			 8 ,  7 ,
																			 0 , 10 ,
																			10 , 10	)) )
																			
myStructure.SetShowTime(True)
if ElementType == "PLANE2D3N":
	elementIncidence = nuto.IntFullMatrix(3,1)
elif ElementType == "PLANE2D4N":
	numElements=5
	elementIncidence = nuto.IntFullMatrix(4,numElements,(	3,2,0,1 ,
															4,6,0,2 ,
															5,4,2,3 ,
															7,5,3,1 ,
															7,6,4,5 ) )
else:
	print 'Wrong elementtype given'
	error = True;
	sys.exit(-1)

myStructure.ElementsCreate(ElementType, elementIncidence)
myStructure.ElementTotalSetConstitutiveLaw(myMatLin)
myStructure.ElementTotalSetSection(mySection)

LoadNodesXPos=myStructure.GroupCreate("Nodes")
LoadNodesXNeg=myStructure.GroupCreate("Nodes")
LoadNodesYPos=myStructure.GroupCreate("Nodes")
LoadNodesYNeg=myStructure.GroupCreate("Nodes")
directionX = nuto.DoubleFullMatrix(2,1,(1,0))
directionY = nuto.DoubleFullMatrix(2,1,(0,1))
if StressState == "XX":
	myStructure.GroupAddNode(LoadNodesXPos,1)
	myStructure.GroupAddNode(LoadNodesXPos,7)
	myStructure.ConstraintLinearSetDisplacementNode(0, nuto.DoubleFullMatrix(2,1,(1,0)), 0.0)
	myStructure.ConstraintLinearSetDisplacementNode(0, nuto.DoubleFullMatrix(2,1,(0,1)), 0.0)
	myStructure.ConstraintLinearSetDisplacementNode(6, nuto.DoubleFullMatrix(2,1,(1,0)), 0.0)
elif StressState == "YY":
	myStructure.GroupAddNode(LoadNodesYPos,6)
	myStructure.GroupAddNode(LoadNodesYPos,7)
	myStructure.ConstraintLinearSetDisplacementNode(0, nuto.DoubleFullMatrix(2,1,(1,0)), 0.0)
	myStructure.ConstraintLinearSetDisplacementNode(0, nuto.DoubleFullMatrix(2,1,(0,1)), 0.0)
	myStructure.ConstraintLinearSetDisplacementNode(1, nuto.DoubleFullMatrix(2,1,(0,1)), 0.0)
elif StressState == "XY":
	myStructure.GroupAddNode(LoadNodesXPos,6)
	myStructure.GroupAddNode(LoadNodesXPos,7)
	#~ myStructure.GroupAddNode("LoadNodesXNeg",0)
	myStructure.GroupAddNode(LoadNodesXNeg,1)
	#~ myStructure.GroupAddNode("LoadNodesYPos",1)
	myStructure.GroupAddNode(LoadNodesYPos,7)
	#~ myStructure.GroupAddNode("LoadNodesYNeg",0)
	myStructure.GroupAddNode(LoadNodesYNeg,6)
	myStructure.ConstraintLinearSetDisplacementNode(0, nuto.DoubleFullMatrix(2,1,(1,0)), 0.0)
	myStructure.ConstraintLinearSetDisplacementNode(0, nuto.DoubleFullMatrix(2,1,(0,1)), 0.0)
	myStructure.ConstraintLinearSetDisplacementNode(1, nuto.DoubleFullMatrix(2,1,(0,1)), 0.0)
else:
	print 'Wrong stressstate given'
	error = True;
	sys.exit(-1)


# boundary conditions
# (left border fixed)

if EnableDisplacementControl:
	print "Displacement control"
	# boundary displacments
	myStructure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesXPos, directionX, BoundaryDisplacement)
	myStructure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesXNeg, directionX, -1.0*BoundaryDisplacement)
	myStructure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesYPos, directionY, BoundaryDisplacement)
	myStructure.ConstraintLinearSetDisplacementNodeGroup(LoadNodesYNeg, directionY, -1.0*BoundaryDisplacement)
else:
	print "Load control"
	myStructure.LoadCreateNodeGroupForce(LoadNodesXPos, directionX, Force)
	myStructure.LoadCreateNodeGroupForce(LoadNodesXNeg, directionX, -1*Force)
	myStructure.LoadCreateNodeGroupForce(LoadNodesYPos, directionY, Force)
	myStructure.LoadCreateNodeGroupForce(LoadNodesYNeg, directionY, -1*Force)
          
# start analysis
# build global dof numbering
myStructure.NodeBuildGlobalDofs()

# build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
stiffnessMatrix = nuto.DoubleSparseMatrixCSRGeneral()
dispForceVector = nuto.DoubleFullMatrix()
myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector)

#calculate element stiffness matrix
Ke = nuto.DoubleFullMatrix(0,0)
rowIndex = nuto.IntFullMatrix(0,0)
colIndex = nuto.IntFullMatrix(0,0)
myStructure.ElementStiffness(2,Ke,rowIndex,colIndex)
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
   if (printResult):
       print "(Ke-KeCorrect).Abs().Max()=" 
       print (Ke-KeCorrect).Abs().Max()

# build global external load vector
extForceVector = nuto.DoubleFullMatrix()
myStructure.BuildGlobalExternalLoadVector(extForceVector)

# calculate right hand side
rhsVector = dispForceVector + extForceVector

# solve
mySolver = nuto.SparseDirectSolverMUMPS()
displacementVector = nuto.DoubleFullMatrix()
stiffnessMatrix.SetOneBasedIndexing()
mySolver.Solve(stiffnessMatrix, rhsVector, displacementVector)

# write displacements to node
myStructure.NodeMergeActiveDofValues(displacementVector)

# calculate residual
intForceVector = nuto.DoubleFullMatrix()
myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector)
if (printResult):
	residualVector = extForceVector - intForceVector
	print "residual: " + str(residualVector.Norm())
	
if ((extForceVector - intForceVector).Norm()>1e-8):
        print '[' + system,sys.argv[0] + '] : Internal and external forces differs.'
        error = True;



#calculate internal force vector
Fi = nuto.DoubleFullMatrix(0,0)
rowIndex = nuto.IntFullMatrix(0,0)
myStructure.ElementGradientInternalPotential(2,Fi,rowIndex)
if (printResult):
    print "Internal Force"
    Fi.Info()

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


#calculate engineering strain of 2 at all integration points
#the size the matrix is not important and reallocated within the procedure
EngineeringStrain = nuto.DoubleFullMatrix(0,0)
myStructure.ElementGetEngineeringStrain(2, EngineeringStrain)

#correct strain
EngineeringStrainCorrect = nuto.DoubleFullMatrix(4,6,(
 0.1   ,  0.1   ,  0.1   ,  0.1  ,
-0.025 , -0.025 , -0.025 , -0.025,
 0.0   ,  0.0   ,  0.0   ,  0.0  ,
 0.0   ,  0.0   ,  0.0   ,  0.0  ,
 0.0   ,  0.0   ,  0.0   ,  0.0  ,
 0.0   ,  0.0   ,  0.0   ,  0.0  )).Trans()

if (printResult):
    print "EngineeringStrainCorrect"
    EngineeringStrainCorrect.Info()
    print "EngineeringStrain"
    EngineeringStrain.Info()

if ((EngineeringStrain-EngineeringStrainCorrect).Abs().Max()[0]>1e-8):
        print '[' + system,sys.argv[0] + '] : strain is not correct.'
        error = True;

#calculate engineering strain of 2 at all integration points
EngineeringStress = nuto.DoubleFullMatrix(0,0)
myStructure.ElementGetEngineeringStress(2, EngineeringStress)
#correct stress
EngineeringStressCorrect = nuto.DoubleFullMatrix(4,6,(
10416.6666667 ,  10416.6666667 , 10416.6666667 , 10416.6666667 ,
    0.0000000 ,      0.0000000 ,     0.0000000 ,     0.0000000 , 
 2083.3333333 ,   2083.3333333 ,  2083.3333333 ,  2083.3333333 ,
    0.0000000 ,      0.0000000 ,     0.0000000 ,     0.0000000 , 
    0.0000000 ,      0.0000000 ,     0.0000000 ,     0.0000000 , 
    0.0000000 ,      0.0000000 ,     0.0000000 ,     0.0000000 )).Trans()

if (printResult):
    print "EngineeringStressCorrect"
    EngineeringStressCorrect.Info()
    print "EngineeringStress"
    EngineeringStress.Info()

if ((EngineeringStress-EngineeringStressCorrect).Abs().Max()[0]>1e-4):
        print '[' + system,sys.argv[0] + '] : stress is not correct.'
        error = True;
        if (printResult):
			print "(EngineeringStress-EngineeringStressCorrect).Abs().Max()="
			print (EngineeringStress-EngineeringStressCorrect).Abs().Max()

if (error):
    sys.exit(-1)
else:
    sys.exit(0)
   
