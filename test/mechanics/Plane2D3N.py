import nuto
import sys
import os

#call of the test file, e.g.
#/usr/local/bin/python ~/develop/nuto/test/mechanics/Plane2D3N.py Linux x86_64 ~/develop/nuto/test/mechanics

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

#create structure
myStructure = nuto.Structure(2)

#create nodes
myStructure.NodesCreate("displacements", nuto.DoubleFullMatrix(2,8,(	 
 0 ,  0 ,
10 ,  0 ,
 2 ,  2 ,
 8 ,  3 ,
 4 ,  7 ,
 8 ,  7 ,
 0 , 10 ,
10 , 10	)))

#create element
elementIncidence = nuto.IntFullMatrix(3,10,(	
0,1,3,
0,2,6,
0,3,2,
1,7,3,
2,4,6,
2,3,4,
3,5,4,
3,7,5,
4,5,6,
5,7,6 ))
myStructure.ElementsCreate("PLANE2D3N", elementIncidence)

#create constitutive law
myMatLin = myStructure.ConstitutiveLawCreate("LinearElastic")
myStructure.ConstitutiveLawSetYoungsModulus(myMatLin,10)
myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin,0.2)

#create section
mySection = myStructure.SectionCreate("Plane_Strain")
myStructure.SectionSetThickness(mySection,1)

#assign constitutive law 
myStructure.ElementTotalSetConstitutiveLaw(myMatLin)
myStructure.ElementTotalSetSection(mySection)

#set displacements of right node
myStructure.ConstraintLinearSetDisplacementNode(0, nuto.DoubleFullMatrix(2,1,(1,0)), 0.0)
myStructure.ConstraintLinearSetDisplacementNode(0, nuto.DoubleFullMatrix(2,1,(0,1)), 0.0)
myStructure.ConstraintLinearSetDisplacementNode(6, nuto.DoubleFullMatrix(2,1,(1,0)), 0.0)
myStructure.ConstraintLinearSetDisplacementNode(1, nuto.DoubleFullMatrix(2,1,(1,0)), 1.0)
myStructure.ConstraintLinearSetDisplacementNode(7, nuto.DoubleFullMatrix(2,1,(1,0)), 1.0)

myStructure.SetVerboseLevel(10)
myStructure.Info()
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
if ((residualVector).Abs().Max()[0]>1e-8):
        print '[' + system,sys.argv[0] + '] : residual force vector is not zero.'
        error = True;

#calculate engineering strain of myelement1 at all integration points
#the size the matrix is not important and reallocated within the procedure
EngineeringStrain = nuto.DoubleFullMatrix(6,1)
#correct strain
EngineeringStrainCorrect = nuto.DoubleFullMatrix(6,1,(
0.1,-0.025,0,0.0,0,0
))
EngineeringStress = nuto.DoubleFullMatrix(6,1)
EngineeringStressCorrect = nuto.DoubleFullMatrix(6,1,(
1.0416666666666667,0.0,0.2083333333333333,0.0,0.0,0.0
))

for element in range(0,8):
 	myStructure.ElementGetEngineeringStrain(element, EngineeringStrain)

	if (printResult):
	    print "EngineeringStrainCorrect"
	    EngineeringStrainCorrect.Info()
	    print "EngineeringStrain"
	    EngineeringStrain.Info()

	if ((EngineeringStrain-EngineeringStrainCorrect).Abs().Max()[0]>1e-8):
	        if (not printResult):
	            print "EngineeringStrainCorrect"
	            EngineeringStrainCorrect.Info()
	            print "EngineeringStrain"
	            EngineeringStrain.Info()
        	print '[' + system,sys.argv[0] + '] : strain is not correct.'
        	error = True;

	#calculate engineering strain at all integration points
	myStructure.ElementGetEngineeringStress(element, EngineeringStress)
	#correct stress

	if (printResult):
	    print "EngineeringStressCorrect"
	    EngineeringStressCorrect.Info()
	    print "EngineeringStress"
	    EngineeringStress.Info()

	if ((EngineeringStress-EngineeringStressCorrect).Abs().Max()[0]>1e-8):
	        if (not printResult):
                    print "EngineeringStressCorrect"
	            EngineeringStressCorrect.Info()
	            print "EngineeringStress"
	            EngineeringStress.Info()
        	print '[' + system,sys.argv[0] + '] : stress is not correct.'
        	error = True;


# visualize results
myStructure.AddVisualizationComponentDisplacements()
myStructure.AddVisualizationComponentEngineeringStrain()
myStructure.AddVisualizationComponentEngineeringStress()
myStructure.ExportVtkDataFile( "Plane2D3N.vtk")

if (error):
    sys.exit(-1)
else:
    sys.exit(0)
