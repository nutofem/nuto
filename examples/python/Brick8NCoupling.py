# -*- coding: utf-8 -*-
import nuto

myStructure = nuto.Structure(3)

# create material law
Material1 = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress")
myStructure.ConstitutiveLawSetParameterDouble(Material1,"Youngs_Modulus", 20000.)
myStructure.ConstitutiveLawSetParameterDouble(Material1,"Poissons_Ratio", 0.2)

# create section
Section1 = myStructure.SectionCreate("Volume")

# create nodes
myStructure.NodeCreate( 1, nuto.DoubleFullVector(( 0, 0, 0)))
myStructure.NodeCreate( 2, nuto.DoubleFullVector(( 1, 0, 0)))
myStructure.NodeCreate( 3, nuto.DoubleFullVector(( 2, 0, 0)))
myStructure.NodeCreate( 4, nuto.DoubleFullVector(( 0, 1, 0)))
myStructure.NodeCreate( 5, nuto.DoubleFullVector(( 1, 1, 0)))
myStructure.NodeCreate( 6, nuto.DoubleFullVector(( 2, 1, 0)))
myStructure.NodeCreate( 7, nuto.DoubleFullVector(( 0, 2, 0)))
myStructure.NodeCreate( 8, nuto.DoubleFullVector(( 1, 2, 0)))
myStructure.NodeCreate( 9, nuto.DoubleFullVector(( 2, 2, 0)))


myStructure.NodeCreate(10, nuto.DoubleFullVector(( 0, 0, 1)))
myStructure.NodeCreate(11, nuto.DoubleFullVector(( 1, 0, 1)))
myStructure.NodeCreate(12, nuto.DoubleFullVector(( 2, 0, 1)))
myStructure.NodeCreate(13, nuto.DoubleFullVector(( 0, 1, 1)))
myStructure.NodeCreate(14, nuto.DoubleFullVector(( 1, 1, 1)))
myStructure.NodeCreate(15, nuto.DoubleFullVector(( 2, 1, 1)))
myStructure.NodeCreate(16, nuto.DoubleFullVector(( 0, 2, 1)))
myStructure.NodeCreate(17, nuto.DoubleFullVector(( 1, 2, 1)))
myStructure.NodeCreate(18, nuto.DoubleFullVector(( 2, 2, 1)))

myStructure.NodeCreate(19, nuto.DoubleFullVector(( 0, 0, 2)))
myStructure.NodeCreate(20, nuto.DoubleFullVector(( 1, 0, 2)))
myStructure.NodeCreate(21, nuto.DoubleFullVector(( 2, 0, 2)))
myStructure.NodeCreate(22, nuto.DoubleFullVector(( 0, 1, 2)))
myStructure.NodeCreate(23, nuto.DoubleFullVector(( 1, 1, 2)))
myStructure.NodeCreate(24, nuto.DoubleFullVector(( 2, 1, 2)))
myStructure.NodeCreate(25, nuto.DoubleFullVector(( 0, 2, 2)))
myStructure.NodeCreate(26, nuto.DoubleFullVector(( 1, 2, 2)))
myStructure.NodeCreate(27, nuto.DoubleFullVector(( 2, 2, 2)))

myStructure.NodeCreate(28, nuto.DoubleFullVector(( 4, 0, 0)))
myStructure.NodeCreate(29, nuto.DoubleFullVector(( 4, 2, 0)))
myStructure.NodeCreate(30, nuto.DoubleFullVector(( 4, 0, 2)))
myStructure.NodeCreate(31, nuto.DoubleFullVector(( 4, 2, 2)))

interpolationType = myStructure.InterpolationTypeCreate("Brick3D")
myStructure.InterpolationTypeAdd(interpolationType, "Coordinates", "Equidistant1")
myStructure.InterpolationTypeAdd(interpolationType, "Displacements", "Equidistant1")

#create elements
elementIncidence = nuto.IntVector(8)

myStructure.ElementCreate(interpolationType, nuto.IntVector([  1,  2,  5,  4, 10, 11, 14, 13]))
myStructure.ElementCreate(interpolationType, nuto.IntVector([  2,  3,  6,  5, 11, 12, 15, 14]))
myStructure.ElementCreate(interpolationType, nuto.IntVector([  4,  5,  8,  7, 13, 14, 17, 16]))
myStructure.ElementCreate(interpolationType, nuto.IntVector([  5,  6,  9,  8, 14, 15, 18, 17]))
myStructure.ElementCreate(interpolationType, nuto.IntVector([ 10, 11, 14, 13, 19, 20, 23, 22]))
myStructure.ElementCreate(interpolationType, nuto.IntVector([ 11, 12, 15, 14, 20, 21, 24, 23]))
myStructure.ElementCreate(interpolationType, nuto.IntVector([ 13, 14, 17, 16, 22, 23, 26, 25]))
myStructure.ElementCreate(interpolationType, nuto.IntVector([ 14, 15, 18, 17, 23, 24, 27, 26]))
myStructure.ElementCreate(interpolationType, nuto.IntVector([  3, 28, 29,  9, 21, 30, 31, 27]))


myStructure.ElementTotalSetConstitutiveLaw(Material1)
myStructure.ElementTotalSetSection(Section1)
myStructure.ElementTotalConvertToInterpolationType()

#boundary conditions
myStructure.ConstraintLinearEquationCreate( 1, "x_displacement", 1, 0)
myStructure.ConstraintLinearEquationCreate( 1, "y_displacement", 1, 0)
myStructure.ConstraintLinearEquationCreate( 1, "z_displacement", 1, 0)
myStructure.ConstraintLinearEquationCreate( 4, "x_displacement", 1, 0)
myStructure.ConstraintLinearEquationCreate( 7, "x_displacement", 1, 0)
myStructure.ConstraintLinearEquationCreate( 7, "z_displacement", 1, 0)
myStructure.ConstraintLinearEquationCreate(10, "x_displacement", 1, 0)
myStructure.ConstraintLinearEquationCreate(13, "x_displacement", 1, 0)
myStructure.ConstraintLinearEquationCreate(16, "x_displacement", 1, 0)
myStructure.ConstraintLinearEquationCreate(19, "x_displacement", 1, 0)
myStructure.ConstraintLinearEquationCreate(22, "x_displacement", 1, 0)
myStructure.ConstraintLinearEquationCreate(25, "x_displacement", 1, 0)

# coupling conditions
id = myStructure.ConstraintLinearEquationCreate( 6, "x_displacement", -1, 0)
myStructure.ConstraintLinearEquationAddTerm(id,  3, "x_displacement", 0.5)
myStructure.ConstraintLinearEquationAddTerm(id,  9, "x_displacement", 0.5)

id = myStructure.ConstraintLinearEquationCreate( 6, "y_displacement", -1, 0)
myStructure.ConstraintLinearEquationAddTerm(id,  3, "y_displacement", 0.5)
myStructure.ConstraintLinearEquationAddTerm(id,  9, "y_displacement", 0.5)

id = myStructure.ConstraintLinearEquationCreate( 6, "z_displacement", -1, 0)
myStructure.ConstraintLinearEquationAddTerm(id,  3, "z_displacement", 0.5)
myStructure.ConstraintLinearEquationAddTerm(id,  9, "z_displacement", 0.5)

id = myStructure.ConstraintLinearEquationCreate(12, "x_displacement", -1, 0)
myStructure.ConstraintLinearEquationAddTerm(id,  3, "x_displacement", 0.5)
myStructure.ConstraintLinearEquationAddTerm(id, 21, "x_displacement", 0.5)

id = myStructure.ConstraintLinearEquationCreate(12, "y_displacement", -1, 0)
myStructure.ConstraintLinearEquationAddTerm(id,  3, "y_displacement", 0.5)
myStructure.ConstraintLinearEquationAddTerm(id, 21, "y_displacement", 0.5)

id = myStructure.ConstraintLinearEquationCreate(12, "z_displacement", -1, 0)
myStructure.ConstraintLinearEquationAddTerm(id,  3, "z_displacement", 0.5)
myStructure.ConstraintLinearEquationAddTerm(id, 21, "z_displacement", 0.5)

id = myStructure.ConstraintLinearEquationCreate(18, "x_displacement", -1, 0)
myStructure.ConstraintLinearEquationAddTerm(id,  9, "x_displacement", 0.5)
myStructure.ConstraintLinearEquationAddTerm(id, 27, "x_displacement", 0.5)

id = myStructure.ConstraintLinearEquationCreate(18, "y_displacement", -1, 0)
myStructure.ConstraintLinearEquationAddTerm(id,  9, "y_displacement", 0.5)
myStructure.ConstraintLinearEquationAddTerm(id, 27, "y_displacement", 0.5)

id = myStructure.ConstraintLinearEquationCreate(18, "z_displacement", -1, 0)
myStructure.ConstraintLinearEquationAddTerm(id,  9, "z_displacement", 0.5)
myStructure.ConstraintLinearEquationAddTerm(id, 27, "z_displacement", 0.5)

id = myStructure.ConstraintLinearEquationCreate(24, "x_displacement", -1, 0)
myStructure.ConstraintLinearEquationAddTerm(id, 21, "x_displacement", 0.5)
myStructure.ConstraintLinearEquationAddTerm(id, 27, "x_displacement", 0.5)

id = myStructure.ConstraintLinearEquationCreate(24, "y_displacement", -1, 0)
myStructure.ConstraintLinearEquationAddTerm(id, 21, "y_displacement", 0.5)
myStructure.ConstraintLinearEquationAddTerm(id, 27, "y_displacement", 0.5)

id = myStructure.ConstraintLinearEquationCreate(24, "z_displacement", -1, 0)
myStructure.ConstraintLinearEquationAddTerm(id, 21, "z_displacement", 0.5)
myStructure.ConstraintLinearEquationAddTerm(id, 27, "z_displacement", 0.5)

id = myStructure.ConstraintLinearEquationCreate(15, "x_displacement", -1, 0)
myStructure.ConstraintLinearEquationAddTerm(id,  3, "x_displacement", 0.25)
myStructure.ConstraintLinearEquationAddTerm(id,  9, "x_displacement", 0.25)
myStructure.ConstraintLinearEquationAddTerm(id, 21, "x_displacement", 0.25)
myStructure.ConstraintLinearEquationAddTerm(id, 27, "x_displacement", 0.25)

id = myStructure.ConstraintLinearEquationCreate(15, "y_displacement", -1, 0)
myStructure.ConstraintLinearEquationAddTerm(id,  3, "y_displacement", 0.25)
myStructure.ConstraintLinearEquationAddTerm(id,  9, "y_displacement", 0.25)
myStructure.ConstraintLinearEquationAddTerm(id, 21, "y_displacement", 0.25)
myStructure.ConstraintLinearEquationAddTerm(id, 27, "y_displacement", 0.25)

id = myStructure.ConstraintLinearEquationCreate(15, "z_displacement", -1, 0)
myStructure.ConstraintLinearEquationAddTerm(id,  3, "z_displacement", 0.25)
myStructure.ConstraintLinearEquationAddTerm(id,  9, "z_displacement", 0.25)
myStructure.ConstraintLinearEquationAddTerm(id, 21, "z_displacement", 0.25)
myStructure.ConstraintLinearEquationAddTerm(id, 27, "z_displacement", 0.25)

# forces
myStructure.SetNumLoadCases(1)
direction = nuto.DoubleFullMatrix(3,1,(1,0,0))
myStructure.LoadCreateNodeForce(0,28, direction, 1)
myStructure.LoadCreateNodeForce(0,29, direction, 1)
myStructure.LoadCreateNodeForce(0,30, direction, 1)
myStructure.LoadCreateNodeForce(0,31, direction, 1)

# start analysis
myStructure.NodeBuildGlobalDofs()

#Calculate maximum independent sets for parallelization (openmp)
myStructure.CalculateMaximumIndependentSets();
myStructure.SolveGlobalSystemStaticElastic(0)

# calculate residual
print "calculate residual"
# start analysis
intGradient = myStructure.BuildGlobalInternalGradient()
extGradient = myStructure.BuildGlobalExternalLoadVector(0)
extGradientJ = extGradient.J.Get("Displacements")
intGradientJ = intGradient.J.Get("Displacements")
intGradientK = intGradient.K.Get("Displacements")
# cast FullVector to FullMatrix, python does not get it...
numRows = intGradientK.GetNumRows()
intGradientKcast = nuto.DoubleFullMatrix(numRows,1)
for i in range(numRows):
    intGradientKcast.SetValue(i,0, -intGradientK.GetValue(i))
cmat = myStructure.GetConstraintMatrix()

residual = nuto.DoubleFullVector(intGradientJ + cmat.Get("Displacements", "Displacements").TransMult(intGradientKcast) - extGradientJ)
print "residual: " + str(residual.Norm())

# visualize results
visualizationGroup = myStructure.GroupCreate("Elements");
myStructure.GroupAddElementsTotal(visualizationGroup)

myStructure.AddVisualizationComponent(visualizationGroup, "Displacements");
myStructure.AddVisualizationComponent(visualizationGroup, "EngineeringStrain");
myStructure.AddVisualizationComponent(visualizationGroup, "EngineeringStress");
myStructure.ExportVtkDataFileElements("Brick8NCoupling.vtk")
