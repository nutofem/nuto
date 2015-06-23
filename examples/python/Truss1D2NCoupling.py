# -*- coding: utf-8 -*-
import nuto

myStructure = nuto.Structure(1)

# create section
Section1 = myStructure.SectionCreate("TRUSS")
myStructure.SectionSetArea(Section1, 1)

# create material law
Material1 = myStructure.ConstitutiveLawCreate("LinearElasticEngineeringStress")
myStructure.ConstitutiveLawSetYoungsModulus(Material1, 1)

# create nodes
nodeCoordinates = nuto.DoubleFullVector(1)
nodeCoordinates.SetValue(0,0,0)
myStructure.NodeCreateDOFs(1, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,1)
myStructure.NodeCreateDOFs(2, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,1)
myStructure.NodeCreateDOFs(3, "displacements", nodeCoordinates)
nodeCoordinates.SetValue(0,0,2)
myStructure.NodeCreateDOFs(4, "displacements", nodeCoordinates)

#create interpolation type
myInterpolationType = myStructure.InterpolationTypeCreate("Truss1D")
myStructure.InterpolationTypeAdd(myInterpolationType, "coordinates", "equidistant1")
myStructure.InterpolationTypeAdd(myInterpolationType, "displacements", "equidistant1")

# create elements
elementIncidence = nuto.IntFullVector(2)
elementIncidence.SetValue(0,0,1)
elementIncidence.SetValue(1,0,2)
e1 = myStructure.ElementCreate(myInterpolationType, elementIncidence)
myStructure.ElementSetSection(e1,Section1)
myStructure.ElementSetConstitutiveLaw(e1,Material1)
elementIncidence.SetValue(0,0,3)
elementIncidence.SetValue(1,0,4)
e2 = myStructure.ElementCreate(myInterpolationType, elementIncidence)
myStructure.ElementSetSection(e2,Section1)
myStructure.ElementSetConstitutiveLaw(e2,Material1)

# set boundary conditions and loads
direction = nuto.DoubleFullMatrix(1,1,(1,))
myStructure.ConstraintLinearSetDisplacementNode(1, direction, 0.0)
id = myStructure.ConstraintLinearEquationCreate(2, "x_displacement", 1, 0)
myStructure.ConstraintLinearEquationAddTerm(id, 3, "x_displacement", -1)
loadCase = 0;
myStructure.SetNumLoadCases(1)
myStructure.LoadCreateNodeForce(loadCase, 4, direction, 1)

#build maximum independent sets
myStructure.CalculateMaximumIndependentSets()

# start analysis
myStructure.NodeBuildGlobalDofs()

# build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
print "build stiffness matrix"
stiffnessMatrix = nuto.DoubleSparseMatrixCSRVector2General()
dispForceVector = nuto.DoubleFullVector()
myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector)
stiffnessMatrix.RemoveZeroEntries(0,1e-14)

# build global external load vector
print "build external force vector"
extForceVector = nuto.DoubleFullVector()
myStructure.BuildGlobalExternalLoadVector(loadCase, extForceVector)

# calculate right hand side
print "build right-hand-side vector"
rhsVector = dispForceVector + extForceVector

# solve
print "solve"
mySolver = nuto.SparseDirectSolverMUMPS()
displacementVector = nuto.DoubleFullVector()
stiffnessMatrixCSR = nuto.DoubleSparseMatrixCSRGeneral(stiffnessMatrix)
stiffnessMatrixCSR.SetOneBasedIndexing()
mySolver.Solve(stiffnessMatrixCSR, rhsVector, displacementVector)

# write displacements to node
print "merge displacements"
myStructure.NodeMergeActiveDofValues(displacementVector)

# calculate residual
print "calculate residual"
intForceVector = nuto.DoubleFullVector()
myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector)
residualVector = extForceVector - intForceVector
print "residual: " + str(residualVector.Norm())


