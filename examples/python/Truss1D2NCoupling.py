# -*- coding: utf-8 -*-
import nuto

myStructure = nuto.Structure(1)

# create section
Section1 = myStructure.SectionCreate("TRUSS")
myStructure.SectionSetArea(Section1, 1)

# create material law
Material1 = myStructure.ConstitutiveLawCreate("LinearElastic")
myStructure.ConstitutiveLawSetYoungsModulus(Material1, 1)

# create nodes
nodeCoordinates = nuto.DoubleFullMatrix(1,1)
nodeCoordinates.SetValue(0,0,0)
myStructure.NodeCreate(1, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,1)
myStructure.NodeCreate(2, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,1)
myStructure.NodeCreate(3, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,2)
myStructure.NodeCreate(4, "displacements", nodeCoordinates)

# create elements
elementIncidence = nuto.IntFullMatrix(2,1)
elementIncidence.SetValue(0,0,1)
elementIncidence.SetValue(1,0,2)
myStructure.ElementCreate(1, "Truss1D2N", elementIncidence)
myStructure.ElementSetSection(1,Section1)
myStructure.ElementSetConstitutiveLaw(1,Material1)
elementIncidence.SetValue(0,0,3)
elementIncidence.SetValue(1,0,4)
myStructure.ElementCreate(2, "Truss1D2N", elementIncidence)
myStructure.ElementSetSection(2,Section1)
myStructure.ElementSetConstitutiveLaw(2,Material1)

# set boundary conditions and loads
direction = nuto.DoubleFullMatrix(1,1,(1,))
myStructure.ConstraintSetDisplacementNode(1, direction, 0.0)
id = myStructure.ConstraintEquationCreate(2, "x_displacement", 1, 0)
myStructure.ConstraintEquationAddTerm(id, 3, "x_displacement", -1)
myStructure.LoadCreateNodeForce(4, direction, 1)

# start analysis
myStructure.NodeBuildGlobalDofs()

# build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
print "build stiffness matrix"
stiffnessMatrix = nuto.DoubleSparseMatrixCSRSymmetric()
dispForceVector = nuto.DoubleFullMatrix()
myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector)
stiffnessMatrix.RemoveZeroEntries(0,1e-14)

# build global external load vector
print "build external force vector"
extForceVector = nuto.DoubleFullMatrix()
myStructure.BuildGlobalExternalLoadVector(extForceVector)

# calculate right hand side
print "build right-hand-side vector"
rhsVector = dispForceVector + extForceVector

# solve
print "solve"
mySolver = nuto.SparseDirectSolverMUMPS()
displacementVector = nuto.DoubleFullMatrix()
stiffnessMatrix.SetOneBasedIndexing()
mySolver.Solve(stiffnessMatrix, rhsVector, displacementVector)

# write displacements to node
print "merge displacements"
myStructure.NodeMergeActiveDofValues(displacementVector)

# calculate residual
print "calculate residual"
intForceVector = nuto.DoubleFullMatrix()
myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector)
residualVector = extForceVector - intForceVector
print "residual: " + str(residualVector.Norm())


