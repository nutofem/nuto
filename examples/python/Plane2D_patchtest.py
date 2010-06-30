#!/usr/bin/python

import nuto
import sys
import os

createResult = True
printResult = True

from time import time

# definitions
YoungsModulus = 1.
PoissonsRatio = 0.0

Force = 1.

#ElementType = "PLANE2D3N"
ElementType = "PLANE2D4N"

EnableDisplacementControl = False
BoundaryDisplacement = 1
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% no error in file, modified, if error is detected              %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error = False

#create structure
myStructure = nuto.Structure(2)

# create material law
myMatLin = myStructure.ConstitutiveLawCreate("LinearElastic")
myStructure.ConstitutiveLawSetYoungsModulus(myMatLin, YoungsModulus)
myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin, PoissonsRatio)

#create section
myStructure.SectionCreate("mySection","Plane_Strain")
myStructure.SectionSetThickness("mySection",1)


myStructure.NodesCreate("displacements", nuto.DoubleFullMatrix(2,8,(	 0 ,  0 ,
																		10 ,  0 ,
																		 0 , 10 ,
																		10 , 10 ,
																		 2 ,  2 ,
																		 8 ,  3 ,
																		 4 ,  7 ,
																		 8 ,  7	)) )



if ElementType == "PLANE2D3N":
	elementIncidence = nuto.IntFullMatrix(3,1)
elif ElementType == "PLANE2D4N":
	elementIncidence = nuto.IntFullMatrix(4,5,( 0 , 1 , 5 , 4 ,
												0 , 4 , 6 , 2 ,
												4 , 5 , 7 , 6 ,
												5 , 1 , 3 , 7 ,
												6 , 7 , 3 , 2 ) )
else:
	print 'Wrong elementtype given'
	error = True;
	sys.exit(-1)

myStructure.ElementsCreate(ElementType, elementIncidence)
myStructure.ElementTotalSetConstitutiveLaw(myMatLin)
myStructure.ElementTotalSetSection("mySection")

myStructure.Info()
myStructure.NodeInfo(5)
myStructure.ElementInfo(5)

# boundary conditions
# (left border fixed)
direction = nuto.DoubleFullMatrix(2,1,(1,0))
myStructure.ConstraintSetDisplacementNode(0, direction, 0.0)
myStructure.ConstraintSetDisplacementNode(2, direction, 0.0)

direction = nuto.DoubleFullMatrix(2,1,(0,1))
myStructure.ConstraintSetDisplacementNode(0, direction, 0.0)

direction = nuto.DoubleFullMatrix(2,1,(1,0))
if EnableDisplacementControl:
	print "Displacement control"
	# boundary displacments
	myStructure.ConstraintSetDisplacementNode(1, direction, BoundaryDisplacement)
	myStructure.ConstraintSetDisplacementNode(3, direction, BoundaryDisplacement)
else:
	print "Load control"
	myStructure.LoadCreateNodeForce(1, direction, Force)
	myStructure.LoadCreateNodeForce(3, direction, Force)
          
# start analysis
# build global dof numbering
curTime = time()
myStructure.NodeBuildGlobalDofs()
oldTime = curTime
curTime = time()
print "time required for dof numbering: " + str(curTime - oldTime) + " s"

# build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
stiffnessMatrix = nuto.DoubleSparseMatrixCSRGeneral()
dispForceVector = nuto.DoubleFullMatrix()
myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector)
oldTime = curTime
curTime = time()
print "time required for assembling: " + str(curTime - oldTime) + " s"

# build global external load vector
extForceVector = nuto.DoubleFullMatrix()
myStructure.BuildGlobalExternalLoadVector(extForceVector)
oldTime = curTime
curTime = time()
print "time required for building external load vector: " + str(curTime - oldTime) + " s"

# calculate right hand side
rhsVector = dispForceVector + extForceVector
oldTime = curTime
curTime = time()
print "time required for calculating right-hand-side vector: " + str(curTime - oldTime) + " s"

# solve
mySolver = nuto.SparseDirectSolverMUMPS()
displacementVector = nuto.DoubleFullMatrix()
stiffnessMatrix.SetOneBasedIndexing()
mySolver.Solve(stiffnessMatrix, rhsVector, displacementVector)
oldTime = curTime
curTime = time()
print "time required for solving: " + str(curTime - oldTime) + " s"

# write displacements to node
myStructure.NodeMergeActiveDofValues(displacementVector)
oldTime = curTime
curTime = time()
print "time required for merging dof values: " + str(curTime - oldTime) + " s"

# calculate residual
intForceVector = nuto.DoubleFullMatrix()
myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector)
oldTime = curTime
curTime = time()
print "time required for building internal forces: " + str(curTime - oldTime) + " s"
residualVector = extForceVector - intForceVector
print "residual: " + str(residualVector.Norm())



# visualize results
myStructure.AddVisualizationComponentDisplacements()
myStructure.AddVisualizationComponentEngineeringStrain()
myStructure.AddVisualizationComponentEngineeringStress()
myStructure.ExportVtkDataFile( ElementType + ".vtk")

if (error):
    sys.exit(-1)
else:
    sys.exit(0)
   
