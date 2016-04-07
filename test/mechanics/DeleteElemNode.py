import nuto
import sys
import os

#call of the test file, e.g.
#/usr/local/bin/python ~/develop/nuto/test/mechanics/DeleteElemNode.py Linux x86_64 ~/develop/nuto/test/mechanics

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


# definitions
YoungsModulus = 20000.
PoissonsRatio = 0.2
Width = 100.
Height = 100.
Length = 200.
NumElementsX = 4
NumElementsY = 2
NumElementsZ = 2

EnableDisplacementControl = False
BoundaryDisplacement = 0.1
Force = 1.

# create three-dimensional structure
myStructure = nuto.Structure(3)

# create material law
myMatLin = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress")
myStructure.ConstitutiveLawSetParameterDouble(myMatLin,"Youngs_Modulus", YoungsModulus)
myStructure.ConstitutiveLawSetParameterDouble(myMatLin,"Poissons_Ratio",PoissonsRatio)

# create section
mySection = myStructure.SectionCreate("Volume")

# create nodes
nodeCoordinates = nuto.DoubleFullVector(3)
#create group of nodes at right boundary
NodeGroup1 = myStructure.GroupCreate("Nodes")
node = 0
for zCount in range (0, NumElementsZ + 1):
    nodeCoordinates.SetValue(2,0, zCount * Height/NumElementsZ)
    for yCount in range (0, NumElementsY + 1):
        nodeCoordinates.SetValue(1,0, yCount * Width/NumElementsY)
        for xCount in range(0, NumElementsX + 1):
            nodeCoordinates.SetValue(0,0, xCount * Length/NumElementsX)
            #print "node: " + str(node) + " coordinates: " + str(nodeCoordinates.GetValue(0,0)) + "," + str(nodeCoordinates.GetValue(1,0)) + "," + str(nodeCoordinates.GetValue(2,0))
            myStructure.NodeCreate(node,nodeCoordinates)
            myStructure.GroupAddNode(NodeGroup1,node)
            node += 1

interpolationType = myStructure.InterpolationTypeCreate("Brick3D");
myStructure.InterpolationTypeAdd(interpolationType, "Coordinates", "Equidistant1");
myStructure.InterpolationTypeAdd(interpolationType, "Displacements", "Equidistant1");


# create elements
elementIncidence = nuto.IntFullVector(8)
element = 0
ElementGroup1 = myStructure.GroupCreate("Elements")
for zCount in range (0, NumElementsZ):
    for yCount in range(0, NumElementsY):
        for xCount in range(0, NumElementsX):
            node1 = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1) + xCount
            elementIncidence.SetValue(0,0, node1)
            elementIncidence.SetValue(1,0, node1 + 1)
            elementIncidence.SetValue(2,0, node1 + NumElementsX + 2)
            elementIncidence.SetValue(3,0, node1 + NumElementsX + 1)
            elementIncidence.SetValue(4,0, node1 + (NumElementsX + 1) * (NumElementsY + 1))
            elementIncidence.SetValue(5,0, node1 + (NumElementsX + 1) * (NumElementsY + 1) + 1)
            elementIncidence.SetValue(6,0, node1 + (NumElementsX + 1) * (NumElementsY + 1) + NumElementsX + 2)
            elementIncidence.SetValue(7,0, node1 + (NumElementsX + 1) * (NumElementsY + 1) + NumElementsX + 1)
            #print "element: " + str(element) + " incidence: "
            #elementIncidence.Info()
            myStructure.ElementCreate(element, interpolationType, elementIncidence,"CONSTITUTIVELAWIP","NOIPDATA")
            myStructure.ElementSetConstitutiveLaw(element,myMatLin)
            myStructure.ElementSetSection(element,mySection)
            myStructure.GroupAddElement(ElementGroup1,element)
            element += 1
            
ElementGroup2 = myStructure.GroupCreate("Elements")
myStructure.GroupAddElement(ElementGroup2,1)
myStructure.GroupAddElement(ElementGroup2,2)
myStructure.GroupAddElement(ElementGroup2,3)

NodeGroup2 = myStructure.GroupCreate("Nodes")
myStructure.GroupAddNode(NodeGroup2,1)
myStructure.GroupAddNode(NodeGroup2,2)
myStructure.GroupAddNode(NodeGroup2,3)

myStructure.ElementDelete(0)
myStructure.NodeDelete(0)

myStructure.ElementDelete(NumElementsX)
myStructure.NodeDelete(NumElementsX+1)
myStructure.NodeDelete(2*(NumElementsX+1))

numNodes=myStructure.GetNumNodes()
numElements=myStructure.GetNumElements()

myStructure.ExportVtkDataFileElements("ElementDelete.vtk")

if(printResult):
    myStructure.ElementInfo(5)
    myStructure.NodeInfo(5)
    myStructure.GroupInfo(3)

if(numElements != element-2):
    print '[' + system,sys.argv[0] + '] : number of elements is not correct.'
    error = True

if(numNodes != node-3):
    print '[' + system,sys.argv[0] + '] : number of nodes is not correct.'
    error = True

if (error):
    sys.exit(-1)
else:
    sys.exit(0)
