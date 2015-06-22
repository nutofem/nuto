import nuto
import sys
import os

#call of the test file, e.g.
#/usr/local/bin/python ~/develop/nuto/test/math/Matrix.py Linux x86_64 ~/develop/nuto/test/math

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

myStructure = nuto.Structure(2)
#coordinates = nuto.FullMatrixDouble(2,1,(1,1))
node1 = myStructure.NodeCreate(nuto.DoubleFullVector((1,1)))
node2 = myStructure.NodeCreate(nuto.DoubleFullVector((1,0)))
node3 = myStructure.NodeCreate(nuto.DoubleFullVector((0,0)))
node4 = myStructure.NodeCreate(nuto.DoubleFullVector((0,1)))

NodeGroup1 = myStructure.GroupCreate("Nodes")
NodeGroup2 = myStructure.GroupCreate("Nodes")

myStructure.GroupAddNode(NodeGroup1,node1)
myStructure.GroupAddNode(NodeGroup1,node2)
myStructure.GroupAddNode(NodeGroup1,node3)
myStructure.GroupAddNode(NodeGroup2,node3)
myStructure.GroupAddNode(NodeGroup2,node4)

# union
NodeGroupUnion = myStructure.GroupUnion(NodeGroup1,NodeGroup2)
if (myStructure.GroupGetNumMembers(NodeGroupUnion)!=4):
    error = True

# difference
NodeGroupDifference = myStructure.GroupDifference(NodeGroup1,NodeGroup2)
if (myStructure.GroupGetNumMembers(NodeGroupDifference)!=2):
    error = True

# intersection
NodeGroupIntersection = myStructure.GroupIntersection(NodeGroup1,NodeGroup2)
if (myStructure.GroupGetNumMembers(NodeGroupIntersection)!=1):
    error = True

# symmetric difference
NodeGroupSymmetricDifference = myStructure.GroupSymmetricDifference(NodeGroup1,NodeGroup2)
if (myStructure.GroupGetNumMembers(NodeGroupSymmetricDifference)!=3):
    error = True

if (printResult):
    myStructure.SetVerboseLevel(10)
    myStructure.Info()

if (error):
    sys.exit(-1)
else:
    sys.exit(0)
