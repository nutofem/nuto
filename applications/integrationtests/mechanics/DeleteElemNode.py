import nuto
import numpy as np
import unittest


class DeleteElemNode(unittest.TestCase):
    def runTest(self):
        Width = 100.
        Height = 100.
        Length = 200.
        NumElementsX = 4
        NumElementsY = 2
        NumElementsZ = 2

        # create three-dimensional structure
        myStructure = nuto.Structure(3)

        # create nodes
        nodeCoordinates = np.zeros((3, 1))
        # create group of nodes at right boundary
        NodeGroup1 = myStructure.GroupCreate("Nodes")
        node = 0
        for zCount in range(0, NumElementsZ + 1):
            nodeCoordinates[2] = zCount * Height/NumElementsZ
            for yCount in range(0, NumElementsY + 1):
                nodeCoordinates[1] = yCount * Width/NumElementsY
                for xCount in range(0, NumElementsX + 1):
                    nodeCoordinates[0] = xCount * Length/NumElementsX
                    myStructure.NodeCreate(node, nodeCoordinates)
                    myStructure.GroupAddNode(NodeGroup1, node)
                    node += 1

        interpolationType = myStructure.InterpolationTypeCreate("Brick3D")
        myStructure.InterpolationTypeAdd(interpolationType, "Coordinates", "Equidistant1")
        myStructure.InterpolationTypeAdd(interpolationType, "Displacements", "Equidistant1")

        # create elements
        elementIncidence = list(range(8))
        element = 0
        ElementGroup1 = myStructure.GroupCreate("Elements")
        for zCount in range(0, NumElementsZ):
            for yCount in range(0, NumElementsY):
                for xCount in range(0, NumElementsX):
                    node1 = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1) + xCount
                    elementIncidence[0] = node1
                    elementIncidence[1] = node1 + 1
                    elementIncidence[2] = node1 + NumElementsX + 2
                    elementIncidence[3] = node1 + NumElementsX + 1
                    elementIncidence[4] = node1 + (NumElementsX + 1) * (NumElementsY + 1)
                    elementIncidence[5] = node1 + (NumElementsX + 1) * (NumElementsY + 1) + 1
                    elementIncidence[6] = node1 + (NumElementsX + 1) * (NumElementsY + 1) + NumElementsX + 2
                    elementIncidence[7] = node1 + (NumElementsX + 1) * (NumElementsY + 1) + NumElementsX + 1
                    myStructure.ElementCreate(element, interpolationType, elementIncidence)
                    myStructure.GroupAddElement(ElementGroup1,element)
                    element += 1

        ElementGroup2 = myStructure.GroupCreate("Elements")
        myStructure.GroupAddElement(ElementGroup2, 1)
        myStructure.GroupAddElement(ElementGroup2, 2)
        myStructure.GroupAddElement(ElementGroup2, 3)

        NodeGroup2 = myStructure.GroupCreate("Nodes")
        myStructure.GroupAddNode(NodeGroup2, 1)
        myStructure.GroupAddNode(NodeGroup2, 2)
        myStructure.GroupAddNode(NodeGroup2, 3)

        myStructure.ElementDelete(0)
        myStructure.NodeDelete(0)

        myStructure.ElementDelete(NumElementsX)
        myStructure.NodeDelete(NumElementsX+1)
        myStructure.NodeDelete(2*(NumElementsX+1))

        numNodes = myStructure.GetNumNodes()
        numElements = myStructure.GetNumElements()

        self.assertEqual(numElements, element-2)
        self.assertEqual(numNodes, node-3)


if __name__ == '__main__':
    unittest.main()
