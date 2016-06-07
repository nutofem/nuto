import nuto
import sys
import os
#import numpy
#import Gnuplot

printResult = False

#create structure
myStructure = nuto.Structure(3)

#create nodes
myNode1 = myStructure.NodeCreateDOFs("displacements",nuto.DoubleFullVector((+1,+1,+1)))
myNode2 = myStructure.NodeCreateDOFs("displacements",nuto.DoubleFullVector(( 0,+1,+1)))
myNode3 = myStructure.NodeCreateDOFs("displacements",nuto.DoubleFullVector(( 0, 0,+1)))
myNode4 = myStructure.NodeCreateDOFs("displacements",nuto.DoubleFullVector((+1, 0,+1)))
myNode5 = myStructure.NodeCreateDOFs("displacements",nuto.DoubleFullVector((+1,+1, 0)))
myNode6 = myStructure.NodeCreateDOFs("displacements",nuto.DoubleFullVector(( 0,+1, 0)))
myNode7 = myStructure.NodeCreateDOFs("displacements",nuto.DoubleFullVector(( 0, 0, 0)))
myNode8 = myStructure.NodeCreateDOFs("displacements",nuto.DoubleFullVector((+1, 0, 0)))

#create interpolation type
myInterpolationType = myStructure.InterpolationTypeCreate("Brick3D");
myStructure.InterpolationTypeAdd(myInterpolationType, "coordinates", "equidistant1");
myStructure.InterpolationTypeAdd(myInterpolationType, "displacements", "equidistant1");

#create element
myElement1 = myStructure.ElementCreate(myInterpolationType,nuto.IntFullVector((myNode5,myNode6,myNode7,myNode8,myNode1,myNode2,myNode3,myNode4)),"ConstitutiveLawIp","StaticData")

#create section
mySection = myStructure.SectionCreate("Volume")

#create constitutive law
myMat = myStructure.ConstitutiveLawCreate("Mises_Plasticity_Engineering_Stress")
#myStructure.ConstitutiveLawCreate("myMat","LinearElastic")

myStructure.ConstitutiveLawSetParameterDouble(myMat,"Youngs_Modulus",100)
myStructure.ConstitutiveLawSetParameterDouble(myMat,"Poissons_Ratio",0.0)
myStructure.ConstitutiveLawSetParameterDouble(myMat,"Initial_Yield_Strength",100)
myStructure.ConstitutiveLawAddYieldStrength(myMat,0.25,150)
myStructure.ConstitutiveLawAddYieldStrength(myMat,0.3,150)
#myStructure.ConstitutiveLawSetInitialHardeningModulus(myMat,0)
#myStructure.ConstitutiveLawAddHardeningModulus(myMat,0.5,0)
#myStructure.ConstitutiveLawAddHardeningModulus(myMat,0.2,1)
myStructure.ConstitutiveLawInfo(10)

#assign constitutive law 
myStructure.ElementSetConstitutiveLaw(myElement1,myMat)
myStructure.ElementSetSection(myElement1,mySection)

#apply constraints of left boundary
direction = nuto.DoubleFullMatrix(3,1,(1,0,0))
myStructure.ConstraintLinearSetDisplacementNode(myNode2, direction, 0)
myStructure.ConstraintLinearSetDisplacementNode(myNode3, direction, 0)
myStructure.ConstraintLinearSetDisplacementNode(myNode6, direction, 0)
myStructure.ConstraintLinearSetDisplacementNode(myNode7, direction, 0)
direction = nuto.DoubleFullMatrix(3,1,(0,1,0))
myStructure.ConstraintLinearSetDisplacementNode(myNode3, direction, 0)
myStructure.ConstraintLinearSetDisplacementNode(myNode7, direction, 0)
direction = nuto.DoubleFullMatrix(3,1,(0,0,1))
myStructure.ConstraintLinearSetDisplacementNode(myNode6, direction, 0)
myStructure.ConstraintLinearSetDisplacementNode(myNode7, direction, 0)

#create group of nodes at right boundary
NodeGroupRightBoundary = myStructure.GroupCreate("Nodes")
myStructure.GroupAddNode(NodeGroupRightBoundary,myNode1)
myStructure.GroupAddNode(NodeGroupRightBoundary,myNode4)
myStructure.GroupAddNode(NodeGroupRightBoundary,myNode5)
myStructure.GroupAddNode(NodeGroupRightBoundary,myNode8)

#apply displacement at right boundary
direction = nuto.DoubleFullMatrix(3,1,(1,0,0))
constraint_right_side = myStructure.ConstraintLinearSetDisplacementNodeGroup(NodeGroupRightBoundary, direction, 0)
   
#number dofs and perform gauss elimination of the constraint matrix
myStructure.NodeBuildGlobalDofs()

#initialize gnuplot
#g = Gnuplot.Gnuplot(debug=1)
#g.title('stress strain curve')   # (optional)
#g("set style line 1 lt 1 lw 1.5 pt 5 ps 0.5 lc rgb 'red'")
#g("set style line 2 lt 1 lw 1.5 pt 7 ps 0.5 lc rgb 'blue'")
#g("set yrange [0:170]")


#loop over boundaryDisplacement
max_disp = 2.
num_steps = 10
stiffnessMatrix = nuto.DoubleSparseMatrixCSRVector2General()
dispForceVector = nuto.DoubleFullVector()
extForceVector = nuto.DoubleFullVector()
mySolver = nuto.SparseDirectSolverMUMPS()
displacementVector = nuto.DoubleFullVector()
deltaDisplacementVector = nuto.DoubleFullVector()
stiffnessMatrix.SetOneBasedIndexing()
intForceVector = nuto.DoubleFullVector()
boundaryForceVector = nuto.DoubleFullVector()
numActiveDofs = 0
plotMatrixLoadDisp = nuto.DoubleFullMatrix(2,1)
plotVectorLoadDisp = nuto.DoubleFullVector(2)

visualizationGroup = myStructure.GroupCreate("Elements");
myStructure.GroupAddElementsTotal(visualizationGroup)

myStructure.AddVisualizationComponent(visualizationGroup, "EngineeringPlasticStrain");
myStructure.AddVisualizationComponent(visualizationGroup, "Displacements");
myStructure.AddVisualizationComponent(visualizationGroup, "EngineeringStrain");
myStructure.AddVisualizationComponent(visualizationGroup, "EngineeringStress");


myStructure.CalculateMaximumIndependentSets()
for i in range(0, num_steps):
    boundaryDisplacement = max_disp*(i+1)/num_steps
    print "boundary displacement :" + str(boundaryDisplacement)
    myStructure.ConstraintSetRHS(constraint_right_side,boundaryDisplacement)
    
    if (numActiveDofs==0):
        numActiveDofs = myStructure.GetNumActiveDofs("Displacements")
        intForceVector.Resize(numActiveDofs)
        displacementVector.Resize(numActiveDofs)

    #do while loop does not exist in python
    iteration = 0
    while True:
        #increment iteration
        iteration+=1
    
        # build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
        myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector)

        # build global external load vector
        myStructure.BuildGlobalExternalLoadVector(1,extForceVector)

        # calculate right hand side
        rhsVector = dispForceVector + extForceVector - intForceVector

        # solve
        stiffnessMatrixCSR = nuto.DoubleSparseMatrixCSRGeneral(stiffnessMatrix)
        mySolver.Solve(stiffnessMatrixCSR, rhsVector, deltaDisplacementVector)

        #add delta
        displacementVector+=deltaDisplacementVector
        
        # write displacements to node
        myStructure.NodeMergeActiveDofValues(displacementVector)

        # calculate internal force vector
        myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector)

        # calculate residual
        residualVector = extForceVector - intForceVector
        print "residual: " + str(residualVector.Norm())

        #calculate engineering strain of myelement at all integration points
        EngineeringStrain = nuto.DoubleFullMatrix(6,3)
        myStructure.ElementGetEngineeringStrain(myElement1, EngineeringStrain)
        print "strain in element 1"
        EngineeringStrain.Info()

        #calculate engineering plastic strain of myelement at all integration points
        EngineeringPlasticStrain = nuto.DoubleFullMatrix(6,3)
        myStructure.ElementGetEngineeringPlasticStrain(myElement1, EngineeringPlasticStrain)
        print "plastic strain in element 1"
        EngineeringPlasticStrain.Info()

        #calculate engineering stress of myelement at all integration points
        EngineeringStress = nuto.DoubleFullMatrix(6,3)
        myStructure.ElementGetEngineeringStress(myElement1, EngineeringStress)
        print "stress in element 1"
        EngineeringStress.Info()
        
        # check for convergence
        if (residualVector.Norm()<1e-6):
            break
        else:
            print "iteration " + str(iteration)
            if printResult:
	        _ = raw_input('More than one iteration required, press enter to continue...') 

    #update static data
    myStructure.ElementTotalUpdateStaticData()
    
    #boundary force
    myStructure.NodeGroupInternalForce(NodeGroupRightBoundary, boundaryForceVector)
    print "boundary force vector"
    boundaryForceVector.Info()
    print ""
    
    #add load displacement curve and write to file
    plotVectorLoadDisp.SetValue(0,0,boundaryDisplacement)
    plotVectorLoadDisp.SetValue(1,0,boundaryForceVector.GetValue(0,0))
    plotMatrixLoadDisp.AppendColumns(plotVectorLoadDisp)
    plotMatrixLoadDisp.Trans().WriteToFile( "MisesLoadDisp.dat", " ", "#load-disp for mises plasticity", "" )
    
    #g("plot 'MisesLoadDisp.dat' using 1:2 title'training set' with lines linestyle 1")

    # visualize results
    myStructure.ExportVtkDataFileElements("MisesPlasticity.vtk")
    
if printResult:
    _ = raw_input('press enter to continue...') 
sys.exit(0)
