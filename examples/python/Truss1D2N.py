# -*- coding: utf-8 -*-

# load nuto package
import nuto

# definitions
YoungsModulus = 20000.
Area = 100. * 100.
Length = 1000.
NumElements = 10
Force = 1.
EnableDisplacementControl = False
BoundaryDisplacement = 0.1

# create one-dimensional structure
myStructure = nuto.Structure(1)

# create section
Section1 = myStructure.SectionCreate("TRUSS")
myStructure.SectionSetArea(Section1, Area)

# create material law
Material1 = myStructure.ConstitutiveLawCreate("LinearElastic")
myStructure.ConstitutiveLawSetYoungsModulus(Material1, YoungsModulus)

# create nodes
nodeCoordinates = nuto.DoubleFullMatrix(1,1)
for node in range(0, NumElements + 1):
    print "create node: " + str(node) + " coordinates: " + str(node * Length/NumElements)
    nodeCoordinates.SetValue(0, 0, node * Length/NumElements)
    myStructure.NodeCreate(node, "displacements", nodeCoordinates)

# create elements
elementIncidence = nuto.IntFullMatrix(2,1)
for element in range(0, NumElements):
    print "create element: " + str(element) + " nodes: " + str(element) + "," + str(element+1)
    elementIncidence.SetValue(0, 0, element)
    elementIncidence.SetValue(1, 0, element + 1)
    myStructure.ElementCreate(element, "Truss1D2N", elementIncidence)
    myStructure.ElementSetSection(element,Section1)
    myStructure.ElementSetConstitutiveLaw(element,Material1)

# set boundary conditions and loads
direction = nuto.DoubleFullMatrix(1,1,(1,))
myStructure.ConstraintSetDisplacementNode(0, direction, 0.0)
if EnableDisplacementControl:
    print "Displacement control"
    myStructure.ConstraintSetDisplacementNode(NumElements, direction, BoundaryDisplacement)
else:
    print "Load control"
    myStructure.LoadCreateNodeForce(NumElements, direction, Force)

# start analysis
# build global dof numbering
myStructure.NodeBuildGlobalDofs()

# build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
stiffnessMatrix = nuto.DoubleSparseMatrixCSRGeneral()
dispForceVector = nuto.DoubleFullMatrix()
myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector)

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
residualVector = extForceVector - intForceVector
print "residual: " + str(residualVector.Norm())

# visualize results
myStructure.AddVisualizationComponentDisplacements()
myStructure.AddVisualizationComponentEngineeringStrain()
myStructure.AddVisualizationComponentEngineeringStress()
myStructure.ExportVtkDataFile("Truss1D2N.vtk")
