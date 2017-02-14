#!/usr/bin/env python3
import unittest
import nuto


class DeleteElemNode(unittest.TestCase):
    def runTest(self):
        width = 100.
        height = 100.
        length = 200.
        numElementsX = 4
        numElementsY = 2
        numElementsZ = 2

        # create three-dimensional structure
        structure = nuto.Structure(3)

        nuto.MeshGenerator.Grid(structure, [length, width, height],
                                [numElementsX, numElementsY, numElementsZ])

        initialElements = structure.GetNumElements()
        initialNodes = structure.GetNumNodes()

        structure.ElementDelete(0)
        structure.NodeDelete(0)

        structure.ElementDelete(numElementsX)
        structure.NodeDelete(numElementsX + 1)
        structure.NodeDelete(2 * (numElementsX + 1))

        numNodes = structure.GetNumNodes()
        numElements = structure.GetNumElements()

        self.assertEqual(numElements, initialElements - 2)
        self.assertEqual(numNodes, initialNodes - 3)


if __name__ == '__main__':
    unittest.main()
