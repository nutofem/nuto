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
myStructure = nuto.Structure(3)

#create nodes
myNode1 = myStructure.NodeCreate("displacements",nuto.DoubleFullMatrix(3,1,(0,0,0)))
myNode2 = myStructure.NodeCreate("displacements rotations",nuto.DoubleFullMatrix(3,1,(5,0,0)))
myNode3 = myStructure.NodeCreate("displacements rotations",nuto.DoubleFullMatrix(3,1,(10,0,0)))

#create group of nodes
myStructure.GroupCreate("myNodeGroup","Nodes")
myStructure.GroupAddNode("myNodeGroup",myNode1)
myStructure.GroupAddNode("myNodeGroup",myNode3)

#create constitutive law
myStructure.ConstitutiveLawCreate("myMatLin","LinearElastic")
myStructure.ConstitutiveLawSetYoungsModulus("myMatLin",10)
myStructure.ConstitutiveLawSetPoissonsRatio("myMatLin",0.1)

#add constraints for a single node
Constraint1 = myStructure.ConstraintSetDisplacementNode(myNode2,nuto.DoubleFullMatrix(3,1,(1,1,-1)),0.5)

#add constraints for a group of nodes
Constraint2 = myStructure.ConstraintSetDisplacementNodeGroup("myNodeGroup",nuto.DoubleFullMatrix(3,1,(1,0,0)),2)
numConstraints = myStructure.ConstraintGetNumConstraintEquations()
if (printResult):
    print "Number of constraints : " + str(numConstraints) 
if (numConstraints!=3):
        print '[' + system,sys.argv[0] + '] : number of constraints is not correct.'
        error = True;

#number global dofs of the nodes
myStructure.NodeNumberGlobalDofs()
numberGlobalDofs = myStructure.NodeGetNumberGlobalDofs()
if (printResult):
    print "Number of global dofs: " + str(numberGlobalDofs) 
if (numberGlobalDofs!=15):
        print '[' + system,sys.argv[0] + '] : number of global dofs is not correct.'
        error = True;

#build constraint matrix and rhs
constraintMatrixSparse = nuto.DoubleSparseMatrixCSRGeneral(numConstraints,numberGlobalDofs)
rhs = nuto.DoubleFullMatrix(numConstraints,1)
myStructure.ConstraintGetConstraintMatrix(constraintMatrixSparse,rhs)
constraintMatrixFull = nuto.DoubleFullMatrix(constraintMatrixSparse)

#correct constraint matrix
constraintMatrixFullCorrect = nuto.DoubleFullMatrix(3,15)
constraintMatrixFullCorrect.SetValue(0,3,1./sqrt(3))
constraintMatrixFullCorrect.SetValue(0,4,1./sqrt(3))
constraintMatrixFullCorrect.SetValue(0,5,-1./sqrt(3))
constraintMatrixFullCorrect.SetValue(1,0,1)
constraintMatrixFullCorrect.SetValue(2,9,1)

#correct rhs
rhsCorrect = nuto.DoubleFullMatrix(3,1,(0.5,2,2))

if (printResult):
    print "constraintMatrixCorrect"
    constraintMatrixFullCorrect.Info()
    print "constraintMatrix"
    constraintMatrixFull.Info()
    print "rhsCorrect"
    rhsCorrect.Info()
    print "rhs"
    rhs.Info()

if ((constraintMatrixFull-constraintMatrixFullCorrect).Abs().Max()[0]>1e-8):
        print '[' + system,sys.argv[0] + '] : constraint matrix is not correct.'
        error = True;

if ((rhs-rhsCorrect).Abs().Max()[0]>1e-8):
        print '[' + system,sys.argv[0] + '] : right hand side is not correct.'
        error = True;

if (error):
    sys.exit(-1)
else:
    sys.exit(0)
