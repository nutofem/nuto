# -*- coding: utf-8 -*-
# $Id: ConstraintsLagrange1D.py 299 2010-07-21 10:38:18Z eckardt4 $
# call with python ~/develop/nuto/test/mechanics/ConstraintsLagrange1D.py Linux x86_64 ~/develop/nuto/test/mechanics

import sys
import nuto
import os
import math

PRINTRESULT = True
MAXNUMNEWTONITERATIONS = 20

#if set to true, the result will be generated (for later use in the test routine)
#otherwise, the current result will be compared to the stored result
createResult = False

#show the results on the screen
printResult = False

#system name and processor
system = sys.argv[1]+sys.argv[2]

#path in the original source directory and current filename at the and
pathToResultFiles = os.path.join(sys.argv[3],"results",system,os.path.basename(sys.argv[0]))

#remove the extension
fileExt = os.path.splitext(sys.argv[0])[1]
pathToResultFiles = pathToResultFiles.replace(fileExt,'')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% no error in file, modified, if error is detected              %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error = False

#create structure
myStructure = nuto.Structure(1);

lx = 1.;

#2 nodes 1 element grid
#create nodes
Coordinates = nuto.DoubleFullMatrix(1,1);
Coordinates.SetValue(0,0,0.0);
node1 = myStructure.NodeCreate("displacements",Coordinates);

Coordinates.SetValue(0,0,lx);
node2 = myStructure.NodeCreate("displacements",Coordinates);

Incidence = nuto.IntFullMatrix(2,1);
Incidence.SetValue(0,0,node1);
Incidence.SetValue(1,0,node2);
myStructure.ElementCreate("TRUSS1D2N",Incidence);

#create constitutive law
myMatLin = myStructure.ConstitutiveLawCreate("LinearElasticEngineeringStress");
myStructure.ConstitutiveLawSetYoungsModulus(myMatLin,1);
myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin,0);

#create section
mySection = myStructure.SectionCreate("TRUSS");
area=1.;
myStructure.SectionSetArea(mySection,area);

myStructure.ElementTotalSetSection(mySection);
myStructure.ElementTotalSetConstitutiveLaw(myMatLin);

#Create groups to apply the boundary conditions
#left boundary
GrpNodesLeftBoundary = myStructure.GroupCreate("Nodes");
direction = 0; #either 0,1,2
min = 0.;
max = 0.;
myStructure.GroupAddNodeCoordinateRange(GrpNodesLeftBoundary,direction,min,max);

GrpNodesRightBoundary = myStructure.GroupCreate("Nodes");
direction = 0; #either 0,1,2
min = lx;
max = lx;
myStructure.GroupAddNodeCoordinateRange(GrpNodesRightBoundary,direction,min,max);

#fix bottom left corner node
DirectionX = nuto.DoubleFullMatrix(1,1);
DirectionX.SetValue(0,0,1.0);

constraintLHS = myStructure.ConstraintLagrangeSetDisplacementNodeGroup(GrpNodesLeftBoundary,DirectionX, "EQUAL",0.0);
myStructure.ConstraintLagrangeSetPenaltyStiffness(constraintLHS,1.);
constraintRHS = myStructure.ConstraintLagrangeSetDisplacementNodeGroup(GrpNodesRightBoundary,DirectionX, "SMALLER",0.5);
myStructure.ConstraintLagrangeSetPenaltyStiffness(constraintRHS,1.);

#ifdef ENABLE_VISUALIZE
myStructure.AddVisualizationComponentSection();
myStructure.AddVisualizationComponentConstitutive();
myStructure.AddVisualizationComponentDisplacements();
myStructure.AddVisualizationComponentEngineeringStrain();
myStructure.AddVisualizationComponentEngineeringStress();
myStructure.AddVisualizationComponentDamage();
myStructure.AddVisualizationComponentEngineeringPlasticStrain();
myStructure.AddVisualizationComponentPrincipalEngineeringStress();
myStructure.ElementTotalUpdateTmpStaticData();
myStructure.ExportVtkDataFileElements("ConstraintsLagrange1D.vtk");
#endif

# init some result data
PlotData = nuto.DoubleFullMatrix(1,7);

# start analysis
maxDisp = 1.;
deltaDispFactor = 0.2;
maxDeltaDispFactor = 0.2;
curDispFactor = 0.2;

#calculate maximum independent sets for openmp parallelization
myStructure.CalculateMaximumIndependentSets();

#update conre mat
myStructure.NodeBuildGlobalDofs();

#update tmpstatic data with zero displanuto.DoubleFullMatrix(cements
myStructure.ElementTotalUpdateTmpStaticData();

#init some auxiliary variables
stiffnessMatrixCSRVector2 = nuto.DoubleSparseMatrixCSRVector2General(0,0);
dispForceVector = nuto.DoubleFullMatrix(0,0);
intForceVector = nuto.DoubleFullMatrix(0,0);
extForceVector = nuto.DoubleFullMatrix(0,0);
rhsVector = nuto.DoubleFullMatrix(0,0);

#allocate solver
mySolver = nuto.SparseDirectSolverMUMPS();
if (PRINTRESULT):
   mySolver.SetShowTime(True);

#calculate stiffness
myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);

#apply initial displacement
curDisp = maxDisp*curDispFactor;
myStructure.ConstraintSetRHS(constraintLHS,curDisp);

#update conre mat (not necessary any more)
#myStructure.NodeBuildGlobalDofs();

#update displacements of all nodes according to the new conre mat
displacementsActiveDOFsCheck = nuto.DoubleFullMatrix(0,0);
displacementsDependentDOFsCheck = nuto.DoubleFullMatrix(0,0);
myStructure.NodeExtractDofValues(0,displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
myStructure.NodeMergeActiveDofValues(0,displacementsActiveDOFsCheck);
myStructure.ElementTotalUpdateTmpStaticData();

#build global external load vector and RHS vector
myStructure.BuildGlobalExternalLoadVector(extForceVector);
#extForceVector.Info(10,13);
myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector);
#intForceVector.Info(10,13);
rhsVector = extForceVector + dispForceVector - intForceVector;

#calculate absolute tolerance for matrix entries to be not considered as zero
maxValue = 0.2/4.2; minValue = 0.; ToleranceZeroStiffness = 0.;
maxValue = stiffnessMatrixCSRVector2.Max();
minValue = stiffnessMatrixCSRVector2.Min();
#std::cout << "min and max " << minValue << " , " << maxValue << std::endl;

if (abs(maxValue)>abs(minValue)):
    ToleranceZeroStiffness = (1e-14) * abs(maxValue);
else:
    ToleranceZeroStiffness = (1e-14) * abs(minValue);
myStructure.SetToleranceStiffnessEntries(ToleranceZeroStiffness);
numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
if (PRINTRESULT):
   print "stiffnessMatrix: num zero removed ", numRemoved, ", numEntries ", numEntries;


#repeat until max displacement is reached
convergenceStatusLoadSteps = False;
loadstep = 1;
displacementsActiveDOFsLastConverged = nuto.DoubleFullMatrix();
displacementsDependentDOFsLastConverged = nuto.DoubleFullMatrix();
while (not convergenceStatusLoadSteps):
    normResidual = 1.;
    maxResidual = 1.;
    numNewtonIterations = 0;
    normRHS = 1.;
    alpha = 1.;
    convergenceStatus = 0;
    #0 - not converged, continue Newton iteration
    #1 - converged
    #2 - stop iteration, decrease load step
    while(convergenceStatus==0):
	numNewtonIterations=numNewtonIterations+1;

        if (numNewtonIterations>MAXNUMNEWTONITERATIONS and alpha<0.5):
            if (PRINTRESULT):
                print "numNewtonIterations (" , numNewtonIterations , ") > MAXNUMNEWTONITERATIONS (" , MAXNUMNEWTONITERATIONS , ")";
            convergenceStatus = 2; #decrease load step
            break;

        normRHS = rhsVector.Norm();

        # solve
        deltaDisplacementsActiveDOFs = nuto.DoubleFullMatrix();
        oldDisplacementsActiveDOFs = nuto.DoubleFullMatrix();
        displacementsActiveDOFs = nuto.DoubleFullMatrix();
        displacementsDependentDOFs = nuto.DoubleFullMatrix();
        if (PRINTRESULT):
            stiffnessMatrixCSRVector2Full = nuto.DoubleFullMatrix(stiffnessMatrixCSRVector2);
            stiffnessMatrixCSRVector2Full.Info(20,10);

        stiffnessMatrixCSR = nuto.DoubleSparseMatrixCSRGeneral(stiffnessMatrixCSRVector2);
        stiffnessMatrixCSR.SetOneBasedIndexing();
        mySolver.Solve(stiffnessMatrixCSR, rhsVector, deltaDisplacementsActiveDOFs);

        # write displacements to node
        myStructure.NodeExtractDofValues(0,oldDisplacementsActiveDOFs, displacementsDependentDOFs);

        #perform a linesearch
        alpha = 1.;
        while True:
            #add new displacement state
            displacementsActiveDOFs = oldDisplacementsActiveDOFs + deltaDisplacementsActiveDOFs*alpha;
            myStructure.NodeMergeActiveDofValues(0,displacementsActiveDOFs);
            myStructure.ElementTotalUpdateTmpStaticData();

            # calculate residual
            myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector);
            rhsVector = extForceVector - intForceVector;
            normResidual = rhsVector.Norm();
            if (PRINTRESULT):
		        print "alpha " , alpha , ", normResidual " , normResidual , ", normResidualInit ", normRHS , ", normRHS*(1-0.5*alpha) " , normRHS*(1-0.5*alpha);
            alpha*=0.5;
            if(not(alpha>1e-3 and normResidual>normRHS*(1-0.5*alpha) and normResidual>1e-5)):
		        break;
        if (normResidual>normRHS*(1-0.5*alpha)  and normResidual>1e-5):
            convergenceStatus=2;
            break;

        maxResidual = rhsVector.Abs().Max();
 
        if (PRINTRESULT):
            print "\n", "Newton iteration " , numNewtonIterations , ", final alpha " , 2*alpha , ", normResidual " , normResidual , ", maxResidual " , maxResidual;

        #check convergence
	if (normResidual<1e-5 or maxResidual<1e-5):
            if (PRINTRESULT):
                print "Convergence after " , numNewtonIterations , " Newton iterations, curdispFactor " , curDispFactor , ", deltaDispFactor ", deltaDispFactor , "\n\n";
            convergenceStatus=1;
            break;

        #convergence status == 0 (continue Newton iteration)
        #build new stiffness matrix
        myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
        numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
        if (PRINTRESULT):
            numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
            print "stiffnessMatrix: num zero removed " , numRemoved , ", numEntries " , numEntries;

    if (convergenceStatus==1):
        myStructure.ElementTotalUpdateStaticData();

        # visualize results
        #ifdef ENABLE_VISUALIZE
        ss = "ConstraintsLagrange1D" + str(loadstep) + ".vtk";
        myStructure.ExportVtkDataFile(ss);
        #endif

        #store result/plot data
	SinglePlotData = nuto.DoubleFullMatrix(1,7);

        #displacements
        dispNode = nuto.DoubleFullMatrix();
        myStructure.NodeGetDisplacements(node1,dispNode);
        SinglePlotData.SetValue(0,0,dispNode.GetValue(0,0));
        myStructure.NodeGetDisplacements(node2,dispNode);
        SinglePlotData.SetValue(0,1,dispNode.GetValue(0,0));

        #boundary force
        SupportingForce = nuto.DoubleFullMatrix();
        myStructure.NodeGroupInternalForce(GrpNodesLeftBoundary,SupportingForce);
        SinglePlotData.SetValue(0,2,SupportingForce.GetValue(0,0));
        myStructure.NodeGroupInternalForce(GrpNodesRightBoundary,SupportingForce);
        SinglePlotData.SetValue(0,3,SupportingForce.GetValue(0,0));

        #lagrange multiplier
        lagrangeMultiplier = nuto.DoubleFullMatrix();
        myStructure.ConstraintLagrangeGetMultiplier(constraintLHS,lagrangeMultiplier);
        SinglePlotData.SetValue(0,4,lagrangeMultiplier.GetValue(0,0));
        myStructure.ConstraintLagrangeGetMultiplier(constraintRHS,lagrangeMultiplier);
        SinglePlotData.SetValue(0,5,lagrangeMultiplier.GetValue(0,0));

        #number of Newton iterations
        SinglePlotData.SetValue(0,6,numNewtonIterations);

        PlotData.AppendRows(SinglePlotData);
        PlotData.WriteToFile("ConstraintsLagrange1DLoadDisp.txt"," ","#disp left and right, boundary force left and right, lagrange multiplier left and right","  ");
        if (PRINTRESULT):
            print "disp left and right, boundary force left and right, lagrange multiplier left and right";
            SinglePlotData.Trans().Info();

        myStructure.NodeExtractDofValues(0,displacementsActiveDOFsLastConverged,displacementsDependentDOFsLastConverged);
        if (curDispFactor==1):
            convergenceStatusLoadSteps=True;
        else:
          #eventually increase load step
            if (numNewtonIterations<MAXNUMNEWTONITERATIONS/3):
                deltaDispFactor*=1.5;
            if (deltaDispFactor>maxDeltaDispFactor):
                deltaDispFactor = maxDeltaDispFactor;

        #increase displacement
        curDispFactor+=deltaDispFactor;
        if (curDispFactor>1):
          deltaDispFactor -= curDispFactor -1.;
          curDispFactor=1;

        curDisp = maxDisp*curDispFactor;

        #old stiffness matrix is used in first step of next load increment in order to prevent spurious problems at the boundary
        #std::cout << "press enter to next load increment, delta disp factor " << deltaDispFactor << " max delta disp factor " <<  maxDeltaDispFactor << std::endl << std::endl;
        #char cDummy[100]="";
        #std::cin.getline(cDummy, 100);
        loadstep=loadstep+1;
    else:
        #calculate stiffness of previous loadstep (used as initial stiffness in the next load step)
        #this is done within the loop in order to ensure, that for the first step the stiffness matrix of the previous step is used
        #otherwise, the additional boundary displacements will result in an artificial localization in elements at the boundary
        curDispFactor-=deltaDispFactor;
        curDisp = maxDisp*curDispFactor;

        myStructure.ConstraintSetRHS(constraintLHS,curDisp);

        # build global dof numbering
        myStructure.NodeBuildGlobalDofs();

        #set previous converged displacements
        myStructure.NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);
        myStructure.ElementTotalUpdateTmpStaticData();

        #decrease load step
        deltaDispFactor*=0.5;
        curDispFactor+=deltaDispFactor;
        curDisp = maxDisp*curDispFactor;

        #check for minimum delta (this mostly indicates an error in the software
        if (deltaDispFactor<MIN_DELTA_STRAIN_FACTOR):
            deltaDispFactor = 0;
            #throw nuto::MechanicsException("Example ConcurrentMultiscale : No convergence, delta strain factor < 1e-7");

        print "press enter to reduce load increment";

    if (not convergenceStatusLoadSteps):
        #update new displacement of RHS
        myStructure.ConstraintSetRHS(constraintLHS,curDisp);

        # build global dof numbering and conre mat
        myStructure.NodeBuildGlobalDofs();

        #update stiffness in order to calculate new dispForceVector
        myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
        numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
        numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
	if (PRINTRESULT):
            print "stiffnessMatrix: num zero removed " , numRemoved , ", numEntries " , numEntries;

        #update displacements of all nodes according to the new conre mat
        displacementsActiveDOFsCheck = nuto.DoubleFullMatrix();
        displacementsDependentDOFsCheck = nuto.DoubleFullMatrix();
	myStructure.NodeExtractDofValues(0,displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
        myStructure.NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
        myStructure.ElementTotalUpdateTmpStaticData();

        # calculate initial residual for next load step
        myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector);

        #update rhs vector for next Newton iteration
        rhsVector = dispForceVector + extForceVector - intForceVector;

PlotDataRef = nuto.DoubleFullMatrix (6,7);
PlotDataRef.SetValue(0,0,0.0);
PlotDataRef.SetValue(1,0,0.2);
PlotDataRef.SetValue(2,0,0.4);
PlotDataRef.SetValue(3,0,0.6);
PlotDataRef.SetValue(4,0,0.8);
PlotDataRef.SetValue(5,0,1.0);

PlotDataRef.SetValue(0,0,0.0);
PlotDataRef.SetValue(1,1,0.2);
PlotDataRef.SetValue(2,1,0.4);
PlotDataRef.SetValue(3,1,0.5);
PlotDataRef.SetValue(4,1,0.5);
PlotDataRef.SetValue(5,1,0.5);

PlotDataRef.SetValue(0,2,0.0);
PlotDataRef.SetValue(1,2,0.0);
PlotDataRef.SetValue(2,2,0.0);
PlotDataRef.SetValue(3,2,0.1);
PlotDataRef.SetValue(4,2,0.3);
PlotDataRef.SetValue(5,2,0.5);

PlotDataRef.SetValue(0,3,0.0);
PlotDataRef.SetValue(1,3,0.0);
PlotDataRef.SetValue(2,3,0.0);
PlotDataRef.SetValue(3,3,-0.1);
PlotDataRef.SetValue(4,3,-0.3);
PlotDataRef.SetValue(5,3,-0.5);

PlotDataRef.SetValue(0,4,0.0);
PlotDataRef.SetValue(1,4,0.0);
PlotDataRef.SetValue(2,4,0.0);
PlotDataRef.SetValue(3,4,-0.1);
PlotDataRef.SetValue(4,4,-0.3);
PlotDataRef.SetValue(5,4,-0.5);

PlotDataRef.SetValue(0,5,0.0);
PlotDataRef.SetValue(1,5,0.0);
PlotDataRef.SetValue(2,5,0.0);
PlotDataRef.SetValue(3,5,0.1);
PlotDataRef.SetValue(4,5,0.3);
PlotDataRef.SetValue(5,5,0.5);

PlotDataRef.SetValue(0,6,0);
PlotDataRef.SetValue(1,6,1);
PlotDataRef.SetValue(2,6,1);
PlotDataRef.SetValue(3,6,2);
PlotDataRef.SetValue(4,6,1);
PlotDataRef.SetValue(5,6,1);

if ((PlotDataRef-PlotData).Abs().Max()[0]>1e-4):
    print "final results stored in load disp file as well";
    PlotData.Info();
    print "reference results";
    PlotDataRef.Info();
    print "delta results";
    (PlotDataRef-PlotData).Info(10,7);
    print "[ConstraintLagrange1D] result is not correct." , (PlotDataRef-PlotData).Abs().Max()[0];
else:
    print "[ConstraintLagrange1D] nice, result is correct";

if (error):
    sys.exit(-1)
else:
    sys.exit(0)
