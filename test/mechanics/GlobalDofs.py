import nuto
import sys
import os
from math import sqrt

#call of the test file, e.g.
#/usr/local/bin/python ~/develop/nuto/test/mechanics/GlobalDofs.py Linux x86_64 ~/develop/nuto/test/mechanics

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
myStructure = nuto.Structure(1)

c = nuto.DoubleFullVector(1)
#create nodes
c.SetValue(0,0) 
myNode1 = myStructure.NodeCreate(c)
c.SetValue(0,5) 
myNode2 = myStructure.NodeCreate(c)
c.SetValue(0,10) 
myNode3 = myStructure.NodeCreate(c)

myInterpolationType = myStructure.InterpolationTypeCreate("Truss1D");
myStructure.InterpolationTypeAdd(myInterpolationType, "coordinates", "equidistant2");
myStructure.InterpolationTypeAdd(myInterpolationType, "displacements", "equidistant2");

myElement1 =  myStructure.ElementCreate(myInterpolationType,nuto.IntVector([myNode1,myNode2,myNode3]))
myStructure.ElementTotalConvertToInterpolationType();

#create group of nodes
myNodeGroup = myStructure.GroupCreate("Nodes")
myStructure.GroupAddNode(myNodeGroup,myNode1)
myStructure.GroupAddNode(myNodeGroup,myNode3)


#create constitutive law
myMatLin = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress")
myStructure.ConstitutiveLawSetParameterDouble(myMatLin,"Youngs_Modulus",10)
myStructure.ConstitutiveLawSetParameterDouble(myMatLin,"Poissons_Ratio",0.1)

#add constraints for a single node
Constraint1 = myStructure.ConstraintLinearSetDisplacementNode(myNode2,nuto.DoubleFullMatrix(3,1,(1,1,-1)),0.5)

#add constraints for a group of nodes
Constraint2 = myStructure.ConstraintLinearSetDisplacementNodeGroup(myNodeGroup,nuto.DoubleFullMatrix(3,1,(1,0,0)),2)
numConstraints = myStructure.ConstraintGetNumLinearConstraints("Displacements")

if (printResult):
    print "Number of constraints : " + str(numConstraints) 
if (numConstraints!=3):
        print '[' + system,sys.argv[0] + '] : number of constraints is not correct.'
        error = True;

#number global dofs of the nodes
myStructure.NodeBuildGlobalDofs()

numberGlobalDofs = myStructure.GetNumDofs("Displacements")
if (printResult):
    print "Number of global dofs: " + str(numberGlobalDofs) 
if (numberGlobalDofs!=3):
        print '[' + system,sys.argv[0] + '] : number of global dofs is not correct.'
        error = True;

#build constraint matrix and rhs
rhs = myStructure.ConstraintGetRHSBeforeGaussElimination().Export()
constraintMatrixFull = myStructure.ConstraintGetConstraintMatrixBeforeGaussElimination().ExportToFullMatrix()


#correct constraint matrix
constraintMatrixFullCorrect = nuto.DoubleFullMatrix(3,3)
constraintMatrixFullCorrect.SetValue(0,0,1)
constraintMatrixFullCorrect.SetValue(1,1,1)
constraintMatrixFullCorrect.SetValue(2,2,1)

#correct rhs
rhsCorrect = nuto.DoubleFullVector((0.5,2,2))

if (printResult):
    print "constraintMatrixCorrect"
    constraintMatrixFullCorrect.Info()
    print "constraintMatrix"
    constraintMatrixFull.Info()
    print "rhsCorrect"
    rhsCorrect.Info()
    print "rhs"
    rhs.Info()

if ((constraintMatrixFull-constraintMatrixFullCorrect).Abs().Max()>1e-8):
        print '[' + system,sys.argv[0] + '] : constraint matrix is not correct.'
        error = True;

if ((rhs-rhsCorrect).Abs().Max()>1e-8):
        print '[' + system,sys.argv[0] + '] : right hand side is not correct.'
        error = True;

if (error):
    sys.exit(-1)
else:
    sys.exit(0)
