# -*- coding: utf-8 -*-
import nuto

myStructure = nuto.Structure(3)

# create material law
Material1 = myStructure.ConstitutiveLawCreate("LinearElastic")
myStructure.ConstitutiveLawSetYoungsModulus(Material1, 20000.)
myStructure.ConstitutiveLawSetPoissonsRatio(Material1, 0.2)

# create nodes
nodeCoordinates = nuto.DoubleFullMatrix(3,1)
nodeCoordinates.SetValue(0,0,0)
nodeCoordinates.SetValue(1,0,0)
nodeCoordinates.SetValue(2,0,0)
myStructure.NodeCreate(1, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,1)
nodeCoordinates.SetValue(1,0,0)
nodeCoordinates.SetValue(2,0,0)
myStructure.NodeCreate(2, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,2)
nodeCoordinates.SetValue(1,0,0)
nodeCoordinates.SetValue(2,0,0)
myStructure.NodeCreate(3, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,0)
nodeCoordinates.SetValue(1,0,1)
nodeCoordinates.SetValue(2,0,0)
myStructure.NodeCreate(4, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,1)
nodeCoordinates.SetValue(1,0,1)
nodeCoordinates.SetValue(2,0,0)
myStructure.NodeCreate(5, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,2)
nodeCoordinates.SetValue(1,0,1)
nodeCoordinates.SetValue(2,0,0)
myStructure.NodeCreate(6, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,0)
nodeCoordinates.SetValue(1,0,2)
nodeCoordinates.SetValue(2,0,0)
myStructure.NodeCreate(7, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,1)
nodeCoordinates.SetValue(1,0,2)
nodeCoordinates.SetValue(2,0,0)
myStructure.NodeCreate(8, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,2)
nodeCoordinates.SetValue(1,0,2)
nodeCoordinates.SetValue(2,0,0)
myStructure.NodeCreate(9, "displacements", nodeCoordinates)

nodeCoordinates.SetValue(0,0,0)
nodeCoordinates.SetValue(1,0,0)
nodeCoordinates.SetValue(2,0,1)
myStructure.NodeCreate(10, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,1)
nodeCoordinates.SetValue(1,0,0)
nodeCoordinates.SetValue(2,0,1)
myStructure.NodeCreate(11, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,2)
nodeCoordinates.SetValue(1,0,0)
nodeCoordinates.SetValue(2,0,1)
myStructure.NodeCreate(12, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,0)
nodeCoordinates.SetValue(1,0,1)
nodeCoordinates.SetValue(2,0,1)
myStructure.NodeCreate(13, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,1)
nodeCoordinates.SetValue(1,0,1)
nodeCoordinates.SetValue(2,0,1)
myStructure.NodeCreate(14, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,2)
nodeCoordinates.SetValue(1,0,1)
nodeCoordinates.SetValue(2,0,1)
myStructure.NodeCreate(15, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,0)
nodeCoordinates.SetValue(1,0,2)
nodeCoordinates.SetValue(2,0,1)
myStructure.NodeCreate(16, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,1)
nodeCoordinates.SetValue(1,0,2)
nodeCoordinates.SetValue(2,0,1)
myStructure.NodeCreate(17, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,2)
nodeCoordinates.SetValue(1,0,2)
nodeCoordinates.SetValue(2,0,1)
myStructure.NodeCreate(18, "displacements", nodeCoordinates)

nodeCoordinates.SetValue(0,0,0)
nodeCoordinates.SetValue(1,0,0)
nodeCoordinates.SetValue(2,0,2)
myStructure.NodeCreate(19, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,1)
nodeCoordinates.SetValue(1,0,0)
nodeCoordinates.SetValue(2,0,2)
myStructure.NodeCreate(20, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,2)
nodeCoordinates.SetValue(1,0,0)
nodeCoordinates.SetValue(2,0,2)
myStructure.NodeCreate(21, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,0)
nodeCoordinates.SetValue(1,0,1)
nodeCoordinates.SetValue(2,0,2)
myStructure.NodeCreate(22, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,1)
nodeCoordinates.SetValue(1,0,1)
nodeCoordinates.SetValue(2,0,2)
myStructure.NodeCreate(23, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,2)
nodeCoordinates.SetValue(1,0,1)
nodeCoordinates.SetValue(2,0,2)
myStructure.NodeCreate(24, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,0)
nodeCoordinates.SetValue(1,0,2)
nodeCoordinates.SetValue(2,0,2)
myStructure.NodeCreate(25, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,1)
nodeCoordinates.SetValue(1,0,2)
nodeCoordinates.SetValue(2,0,2)
myStructure.NodeCreate(26, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,2)
nodeCoordinates.SetValue(1,0,2)
nodeCoordinates.SetValue(2,0,2)
myStructure.NodeCreate(27, "displacements", nodeCoordinates)

nodeCoordinates.SetValue(0,0,4)
nodeCoordinates.SetValue(1,0,0)
nodeCoordinates.SetValue(2,0,0)
myStructure.NodeCreate(28, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,4)
nodeCoordinates.SetValue(1,0,2)
nodeCoordinates.SetValue(2,0,0)
myStructure.NodeCreate(29, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,4)
nodeCoordinates.SetValue(1,0,0)
nodeCoordinates.SetValue(2,0,2)
myStructure.NodeCreate(30, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,4)
nodeCoordinates.SetValue(1,0,2)
nodeCoordinates.SetValue(2,0,2)
myStructure.NodeCreate(31, "displacements", nodeCoordinates)

#create elements
elementIncidence = nuto.IntFullMatrix(8,1)
elementIncidence.SetValue(0,0,1)
elementIncidence.SetValue(1,0,2)
elementIncidence.SetValue(2,0,5)
elementIncidence.SetValue(3,0,4)
elementIncidence.SetValue(4,0,10)
elementIncidence.SetValue(5,0,11)
elementIncidence.SetValue(6,0,14)
elementIncidence.SetValue(7,0,13)
myStructure.ElementCreate(1, "Brick8N", elementIncidence)
myStructure.ElementSetConstitutiveLaw(1,Material1)
elementIncidence.SetValue(0,0,2)
elementIncidence.SetValue(1,0,3)
elementIncidence.SetValue(2,0,6)
elementIncidence.SetValue(3,0,5)
elementIncidence.SetValue(4,0,11)
elementIncidence.SetValue(5,0,12)
elementIncidence.SetValue(6,0,15)
elementIncidence.SetValue(7,0,14)
myStructure.ElementCreate(2, "Brick8N", elementIncidence)
myStructure.ElementSetConstitutiveLaw(2,Material1)
elementIncidence.SetValue(0,0,4)
elementIncidence.SetValue(1,0,5)
elementIncidence.SetValue(2,0,8)
elementIncidence.SetValue(3,0,7)
elementIncidence.SetValue(4,0,13)
elementIncidence.SetValue(5,0,14)
elementIncidence.SetValue(6,0,17)
elementIncidence.SetValue(7,0,16)
myStructure.ElementCreate(3, "Brick8N", elementIncidence)
myStructure.ElementSetConstitutiveLaw(3,Material1)
elementIncidence.SetValue(0,0,5)
elementIncidence.SetValue(1,0,6)
elementIncidence.SetValue(2,0,9)
elementIncidence.SetValue(3,0,8)
elementIncidence.SetValue(4,0,14)
elementIncidence.SetValue(5,0,15)
elementIncidence.SetValue(6,0,18)
elementIncidence.SetValue(7,0,17)
myStructure.ElementCreate(4, "Brick8N", elementIncidence)
myStructure.ElementSetConstitutiveLaw(4,Material1)
elementIncidence.SetValue(0,0,10)
elementIncidence.SetValue(1,0,11)
elementIncidence.SetValue(2,0,14)
elementIncidence.SetValue(3,0,13)
elementIncidence.SetValue(4,0,19)
elementIncidence.SetValue(5,0,20)
elementIncidence.SetValue(6,0,23)
elementIncidence.SetValue(7,0,22)
myStructure.ElementCreate(5, "Brick8N", elementIncidence)
myStructure.ElementSetConstitutiveLaw(5,Material1)
elementIncidence.SetValue(0,0,11)
elementIncidence.SetValue(1,0,12)
elementIncidence.SetValue(2,0,15)
elementIncidence.SetValue(3,0,14)
elementIncidence.SetValue(4,0,20)
elementIncidence.SetValue(5,0,21)
elementIncidence.SetValue(6,0,24)
elementIncidence.SetValue(7,0,23)
myStructure.ElementCreate(6, "Brick8N", elementIncidence)
myStructure.ElementSetConstitutiveLaw(6,Material1)
elementIncidence.SetValue(0,0,13)
elementIncidence.SetValue(1,0,14)
elementIncidence.SetValue(2,0,17)
elementIncidence.SetValue(3,0,16)
elementIncidence.SetValue(4,0,22)
elementIncidence.SetValue(5,0,23)
elementIncidence.SetValue(6,0,26)
elementIncidence.SetValue(7,0,25)
myStructure.ElementCreate(7, "Brick8N", elementIncidence)
myStructure.ElementSetConstitutiveLaw(7,Material1)
elementIncidence.SetValue(0,0,14)
elementIncidence.SetValue(1,0,15)
elementIncidence.SetValue(2,0,18)
elementIncidence.SetValue(3,0,17)
elementIncidence.SetValue(4,0,23)
elementIncidence.SetValue(5,0,24)
elementIncidence.SetValue(6,0,27)
elementIncidence.SetValue(7,0,26)
myStructure.ElementCreate(8, "Brick8N", elementIncidence)
myStructure.ElementSetConstitutiveLaw(8,Material1)
elementIncidence.SetValue(0,0,3)
elementIncidence.SetValue(1,0,28)
elementIncidence.SetValue(2,0,29)
elementIncidence.SetValue(3,0,9)
elementIncidence.SetValue(4,0,21)
elementIncidence.SetValue(5,0,30)
elementIncidence.SetValue(6,0,31)
elementIncidence.SetValue(7,0,27)
myStructure.ElementCreate(9, "Brick8N", elementIncidence)
myStructure.ElementSetConstitutiveLaw(9,Material1)

#boundary conditions
myStructure.ConstraintEquationCreate(1, "x_displacement", 1, 0)
myStructure.ConstraintEquationCreate(1, "y_displacement", 1, 0)
myStructure.ConstraintEquationCreate(1, "z_displacement", 1, 0)
myStructure.ConstraintEquationCreate(4, "x_displacement", 1, 0)
myStructure.ConstraintEquationCreate(7, "x_displacement", 1, 0)
myStructure.ConstraintEquationCreate(7, "z_displacement", 1, 0)
myStructure.ConstraintEquationCreate(10, "x_displacement", 1, 0)
myStructure.ConstraintEquationCreate(13, "x_displacement", 1, 0)
myStructure.ConstraintEquationCreate(16, "x_displacement", 1, 0)
myStructure.ConstraintEquationCreate(19, "x_displacement", 1, 0)
myStructure.ConstraintEquationCreate(22, "x_displacement", 1, 0)
myStructure.ConstraintEquationCreate(25, "x_displacement", 1, 0)

# coupling conditions
id = myStructure.ConstraintEquationCreate(6, "x_displacement", -1, 0)
myStructure.ConstraintEquationAddTerm(id, 3, "x_displacement", 0.5)
myStructure.ConstraintEquationAddTerm(id, 9, "x_displacement", 0.5)
id = myStructure.ConstraintEquationCreate(6, "y_displacement", -1, 0)
myStructure.ConstraintEquationAddTerm(id, 3, "y_displacement", 0.5)
myStructure.ConstraintEquationAddTerm(id, 9, "y_displacement", 0.5)
id = myStructure.ConstraintEquationCreate(6, "z_displacement", -1, 0)
myStructure.ConstraintEquationAddTerm(id, 3, "z_displacement", 0.5)
myStructure.ConstraintEquationAddTerm(id, 9, "z_displacement", 0.5)
id = myStructure.ConstraintEquationCreate(12, "x_displacement", -1, 0)
myStructure.ConstraintEquationAddTerm(id, 3, "x_displacement", 0.5)
myStructure.ConstraintEquationAddTerm(id, 21, "x_displacement", 0.5)
id = myStructure.ConstraintEquationCreate(12, "y_displacement", -1, 0)
myStructure.ConstraintEquationAddTerm(id, 3, "y_displacement", 0.5)
myStructure.ConstraintEquationAddTerm(id, 21, "y_displacement", 0.5)
id = myStructure.ConstraintEquationCreate(12, "z_displacement", -1, 0)
myStructure.ConstraintEquationAddTerm(id, 3, "z_displacement", 0.5)
myStructure.ConstraintEquationAddTerm(id, 21, "z_displacement", 0.5)
id = myStructure.ConstraintEquationCreate(18, "x_displacement", -1, 0)
myStructure.ConstraintEquationAddTerm(id, 9, "x_displacement", 0.5)
myStructure.ConstraintEquationAddTerm(id, 27, "x_displacement", 0.5)
id = myStructure.ConstraintEquationCreate(18, "y_displacement", -1, 0)
myStructure.ConstraintEquationAddTerm(id, 9, "y_displacement", 0.5)
myStructure.ConstraintEquationAddTerm(id, 27, "y_displacement", 0.5)
id = myStructure.ConstraintEquationCreate(18, "z_displacement", -1, 0)
myStructure.ConstraintEquationAddTerm(id, 9, "z_displacement", 0.5)
myStructure.ConstraintEquationAddTerm(id, 27, "z_displacement", 0.5)
id = myStructure.ConstraintEquationCreate(24, "x_displacement", -1, 0)
myStructure.ConstraintEquationAddTerm(id, 21, "x_displacement", 0.5)
myStructure.ConstraintEquationAddTerm(id, 27, "x_displacement", 0.5)
id = myStructure.ConstraintEquationCreate(24, "y_displacement", -1, 0)
myStructure.ConstraintEquationAddTerm(id, 21, "y_displacement", 0.5)
myStructure.ConstraintEquationAddTerm(id, 27, "y_displacement", 0.5)
id = myStructure.ConstraintEquationCreate(24, "z_displacement", -1, 0)
myStructure.ConstraintEquationAddTerm(id, 21, "z_displacement", 0.5)
myStructure.ConstraintEquationAddTerm(id, 27, "z_displacement", 0.5)
id = myStructure.ConstraintEquationCreate(15, "x_displacement", -1, 0)
myStructure.ConstraintEquationAddTerm(id, 3, "x_displacement", 0.25)
myStructure.ConstraintEquationAddTerm(id, 9, "x_displacement", 0.25)
myStructure.ConstraintEquationAddTerm(id, 21, "x_displacement", 0.25)
myStructure.ConstraintEquationAddTerm(id, 27, "x_displacement", 0.25)
id = myStructure.ConstraintEquationCreate(15, "y_displacement", -1, 0)
myStructure.ConstraintEquationAddTerm(id, 3, "y_displacement", 0.25)
myStructure.ConstraintEquationAddTerm(id, 9, "y_displacement", 0.25)
myStructure.ConstraintEquationAddTerm(id, 21, "y_displacement", 0.25)
myStructure.ConstraintEquationAddTerm(id, 27, "y_displacement", 0.25)
id = myStructure.ConstraintEquationCreate(15, "z_displacement", -1, 0)
myStructure.ConstraintEquationAddTerm(id, 3, "z_displacement", 0.25)
myStructure.ConstraintEquationAddTerm(id, 9, "z_displacement", 0.25)
myStructure.ConstraintEquationAddTerm(id, 21, "z_displacement", 0.25)
myStructure.ConstraintEquationAddTerm(id, 27, "z_displacement", 0.25)

# forces
direction = nuto.DoubleFullMatrix(3,1,(1,0,0))
myStructure.LoadCreateNodeForce(28, direction, 1)
myStructure.LoadCreateNodeForce(29, direction, 1)
myStructure.LoadCreateNodeForce(30, direction, 1)
myStructure.LoadCreateNodeForce(31, direction, 1)

# start analysis
myStructure.NodeBuildGlobalDofs()

# build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
print "build stiffness matrix"
#stiffnessMatrix = nuto.DoubleSparseMatrixCSRGeneral()
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

# visualize results
myStructure.AddVisualizationComponentDisplacements()
myStructure.AddVisualizationComponentEngineeringStrain()
myStructure.AddVisualizationComponentEngineeringStress()
myStructure.ExportVtkDataFile("Brick8NCoupling.vtk")
