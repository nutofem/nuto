#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/MechanicsException.h"

#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"

#define MAXNUMNEWTONITERATIONS 20
//#define MAXNORMRHS 100
#define PRINTRESULT true
#define MIN_DELTA_STRAIN_FACTOR 1e-7

int main()
{
try
{
    //
    double lX(100);
    double lY(100);

    //create structure
    NuTo::Structure myStructureCoarseScale(2);
    myStructureCoarseScale.SetShowTime(true);

    NuTo::FullMatrix<int> createdGroupIds;
    myStructureCoarseScale.ImportFromGmsh("/home/unger3/develop/nuto/examples/c++/ConcurrentMultiscale.msh","displacements", "ConstitutiveLawIp", "StaticData",createdGroupIds);
    //myStructureCoarseScale.ImportFromGmsh("ConcurrentMultiscale.msh","displacements", "ConstitutiveLawIpNonlocal", "StaticData",createdGroupIds);

    //create constitutive law from fine scale model
    int myMatCoarseScale = myStructureCoarseScale.ConstitutiveLawCreate("Multiscale");

    //create section
    double thickness(1);
    int mySectionCoarseScale = myStructureCoarseScale.SectionCreate("Plane_Strain");
    myStructureCoarseScale.SectionSetThickness(mySectionCoarseScale,thickness);

    myStructureCoarseScale.ElementTotalSetSection(mySectionCoarseScale);
	myStructureCoarseScale.ElementTotalSetConstitutiveLaw(myMatCoarseScale);

	//set fine scale model as ip data for all integration points
	//myStructureCoarseScale.ElementTotalSetFineScaleModel("myStructureFineScale.bin");
    myStructureCoarseScale.ElementTotalSetFineScaleModel("/home/unger3/develop/nuto_build/examples/c++/myStructureFineScale.bin");

	//Create groups to apply the boundary conditions
	//left boundary
	int GrpNodesLeftBoundary = myStructureCoarseScale.GroupCreate("Nodes");
	int direction = 0; //either 0,1,2
	double min(0.);
	double max(0.);
	myStructureCoarseScale.GroupAddNodeCoordinateRange(GrpNodesLeftBoundary,direction,min,max);

    //right boundary
    int GrpNodesRightBoundary = myStructureCoarseScale.GroupCreate("Nodes");
    direction = 0; //either 0,1,2
    min=lX;
    max=lX;
    myStructureCoarseScale.GroupAddNodeCoordinateRange(GrpNodesRightBoundary,direction,min,max);

    //bottom boundary
    int GrpNodesBottomBoundary = myStructureCoarseScale.GroupCreate("Nodes");
    direction=1;
    min=0;
    max=0;
    myStructureCoarseScale.GroupAddNodeCoordinateRange(GrpNodesBottomBoundary,direction,min,max);

    //left lower node
    int GrpNodesBottomLeftNodeBoundary = myStructureCoarseScale.GroupIntersection(GrpNodesBottomBoundary,GrpNodesLeftBoundary);

    //fix bottom left corner node
	NuTo::FullMatrix<double> DirectionX(2,1);
	DirectionX.SetValue(0,0,1.0);
	DirectionX.SetValue(1,0,0.0);

	NuTo::FullMatrix<double>DirectionY(2,1);
	DirectionY.SetValue(0,0,0.0);
	DirectionY.SetValue(1,0,1.0);

    myStructureCoarseScale.ConstraintLinearSetDisplacementNodeGroup(GrpNodesLeftBoundary ,DirectionX, 0);
    int ConstraintRHS = myStructureCoarseScale.ConstraintLinearSetDisplacementNodeGroup(GrpNodesRightBoundary,DirectionX, 0);
    myStructureCoarseScale.ConstraintLinearSetDisplacementNodeGroup(GrpNodesBottomLeftNodeBoundary,DirectionY, 0);
/*
#ifdef ENABLE_VISUALIZE
    myStructureCoarseScale.AddVisualizationComponentSection();
    myStructureCoarseScale.AddVisualizationComponentConstitutive();
    myStructureCoarseScale.AddVisualizationComponentDisplacements();
    myStructureCoarseScale.AddVisualizationComponentEngineeringStrain();
    myStructureCoarseScale.AddVisualizationComponentEngineeringStress();
    myStructureCoarseScale.AddVisualizationComponentDamage();
    myStructureCoarseScale.AddVisualizationComponentEngineeringPlasticStrain();
    myStructureCoarseScale.AddVisualizationComponentPrincipalEngineeringStress();
	myStructureCoarseScale.ElementTotalUpdateTmpStaticData();
    myStructureCoarseScale.ExportVtkDataFile("ConcurrentMultiscale.vtk");
#endif
*/
    // init some result data
    NuTo::FullMatrix<double> PlotData(1,9);

    // start analysis
    double maxDisp(1);
    double deltaDispFactor(0.0005);
    double maxDeltaDispFactor(0.05);
    double curDispFactor(0.00025);

    //update conre mat
	myStructureCoarseScale.NodeBuildGlobalDofs();

    //update tmpstatic data with zero displacements
	myStructureCoarseScale.ElementTotalUpdateTmpStaticData();

	//init some auxiliary variables
	NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixCSRVector2;
	NuTo::FullMatrix<double> dispForceVector;
	NuTo::FullMatrix<double> intForceVector;
	NuTo::FullMatrix<double> extForceVector;
	NuTo::FullMatrix<double> rhsVector;

	//allocate solver
	NuTo::SparseDirectSolverMUMPS mySolver;
	mySolver.SetShowTime(true);

    //calculate stiffness
	myStructureCoarseScale.BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
	NuTo::FullMatrix<double> fullStiffness(stiffnessMatrixCSRVector2);
	fullStiffness.Info(12,10);

	double curDisp(maxDisp*curDispFactor);
    myStructureCoarseScale.ConstraintSetRHS(ConstraintRHS,curDisp);

    //update conre mat
    myStructureCoarseScale.NodeBuildGlobalDofs();

    //update displacements of all nodes according to the new conre mat
    {
        NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
        NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
        myStructureCoarseScale.NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
        myStructureCoarseScale.NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
        myStructureCoarseScale.ElementTotalUpdateTmpStaticData();
    }

    // build global external load vector and RHS vector
	myStructureCoarseScale.BuildGlobalExternalLoadVector(extForceVector);
	myStructureCoarseScale.BuildGlobalGradientInternalPotentialVector(intForceVector);
    //intForceVector.Info(10,13);
    rhsVector = extForceVector + dispForceVector - intForceVector;

	//calculate absolute tolerance for matrix entries to be not considered as zero
	double maxValue, minValue, ToleranceZeroStiffness;
	stiffnessMatrixCSRVector2.Max(maxValue);
	stiffnessMatrixCSRVector2.Min(minValue);
	std::cout << "min and max " << minValue << " , " << maxValue << std::endl;

	ToleranceZeroStiffness = (1e-14) * (fabs(maxValue)>fabs(minValue) ?  fabs(maxValue) : fabs(minValue));
	myStructureCoarseScale.SetToleranceStiffnessEntries(ToleranceZeroStiffness);
	int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
	int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
	std::cout << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << std::endl;

	//repeat until max displacement is reached
	bool convergenceStatusLoadSteps(false);
    int loadstep(1);
    NuTo::FullMatrix<double> displacementsActiveDOFsLastConverged,displacementsDependentDOFsLastConverged;
    myStructureCoarseScale.NodeExtractDofValues(displacementsActiveDOFsLastConverged,displacementsDependentDOFsLastConverged);
	while (!convergenceStatusLoadSteps)
    {

        double normResidual(1);
        double maxResidual(1);
        int numNewtonIterations(0);
		double normRHS(1.);
		double alpha(1.);
		int convergenceStatus(0);
		//0 - not converged, continue Newton iteration
		//1 - converged
		//2 - stop iteration, decrease load step
		while(convergenceStatus==0)
		{
			numNewtonIterations++;

			if (numNewtonIterations>MAXNUMNEWTONITERATIONS)
			{
				if (PRINTRESULT)
		    	{
				    std::cout << "numNewtonIterations (" << numNewtonIterations << ") > MAXNUMNEWTONITERATIONS (" << MAXNUMNEWTONITERATIONS << ")" << std::endl;
		    	}
				convergenceStatus = 2; //decrease load step
				break;
			}

			normRHS = rhsVector.Norm();

			// solve
			NuTo::FullMatrix<double> deltaDisplacementsActiveDOFs;
			NuTo::FullMatrix<double> oldDisplacementsActiveDOFs;
			NuTo::FullMatrix<double> displacementsActiveDOFs;
			NuTo::FullMatrix<double> displacementsDependentDOFs;
			NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrixCSR(stiffnessMatrixCSRVector2);
			stiffnessMatrixCSR.SetOneBasedIndexing();
			mySolver.Solve(stiffnessMatrixCSR, rhsVector, deltaDisplacementsActiveDOFs);

			// write displacements to node
			myStructureCoarseScale.NodeExtractDofValues(oldDisplacementsActiveDOFs, displacementsDependentDOFs);

			//perform a linesearch
			alpha = 1.;
			do
			{
				//add new displacement state
				displacementsActiveDOFs = oldDisplacementsActiveDOFs + deltaDisplacementsActiveDOFs*alpha;
				myStructureCoarseScale.NodeMergeActiveDofValues(displacementsActiveDOFs);
				myStructureCoarseScale.ElementTotalUpdateTmpStaticData();

				// calculate residual
				myStructureCoarseScale.BuildGlobalGradientInternalPotentialVector(intForceVector);
				rhsVector = extForceVector - intForceVector;
				normResidual = rhsVector.Norm();
				std::cout << "alpha " << alpha << ", normResidual " << normResidual << ", normResidualInit "<< normRHS << ", normRHS*(1-0.5*alpha) " << normRHS*(1-0.5*alpha) << std::endl;
				alpha*=0.5;
			}
			while(alpha>1e-3 && normResidual>normRHS*(1-0.5*alpha) && normResidual>1e-5);
			if (normResidual>normRHS*(1-0.5*alpha) && normResidual>1e-5)
			{
			    convergenceStatus=2;
			    break;
			}

		    maxResidual = rhsVector.Max();

			std::cout << std::endl << "Newton iteration " << numNewtonIterations << ", final alpha " << 2*alpha << ", normResidual " << normResidual<< ", maxResidual " << maxResidual<<std::endl;
			//char cDummy[100]="";
			//std::cin.getline(cDummy, 100);

			//check convergence
			if (normResidual<1e-5 || maxResidual<1e-5)
			{
				if (PRINTRESULT)
		    	{
					std::cout << "Convergence after " << numNewtonIterations << " Newton iterations, curdispFactor " << curDispFactor << ", deltaDispFactor "<< deltaDispFactor << std::endl<< std::endl;
		    	}
				convergenceStatus=1;
				break;
			}

			//convergence status == 0 (continue Newton iteration)
			//build new stiffness matrix
			myStructureCoarseScale.BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
			int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
			int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
			std::cout << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << std::endl;
		}

		if (deltaDispFactor<1e-7)
		    throw NuTo::MechanicsException("Example ConcurrentMultiscale : No convergence, delta strain factor < 1e-7");

		if (convergenceStatus==1)
		{
			myStructureCoarseScale.ElementTotalUpdateStaticData();

            // visualize results
#ifdef ENABLE_VISUALIZE
            myStructureCoarseScale.ExportVtkDataFile("ConcurrentMultiscale.vtk");
#endif
 			//store result/plot data
			NuTo::FullMatrix<double> SinglePlotData(1,9);
			SinglePlotData(0,0) = curDisp;
            //Get Average Stress
            NuTo::FullMatrix<double> averageStress;
	        myStructureCoarseScale.ElementTotalGetAverageStress(lX*lY,averageStress);
            SinglePlotData(0,3) = averageStress(0,0);
            SinglePlotData(0,4) = averageStress(1,0);
            SinglePlotData(0,5) = averageStress(2,0);
            SinglePlotData(0,6) = averageStress(3,0);
            SinglePlotData(0,7) = averageStress(4,0);
            SinglePlotData(0,8) = averageStress(5,0);

            std::cout << "average stress " << std::endl << averageStress << std::endl;
            std::cout << "current Strain " << std::endl << curDisp << std::endl;

//			externalEnergy+=averageStress(0,0)*curStrain(0,0)+averageStress(1,0)*curStrain(1,0)+averageStress(2,0)*curStrain(2,0);

			PlotData.AppendRows(SinglePlotData);
		    PlotData.WriteToFile("ConcurrentMultiscaleLoadDisp.txt"," ","#load displacement curve, disp, stress, force, sxx in center element, syy in center element","  ");

		    myStructureCoarseScale.NodeExtractDofValues(displacementsActiveDOFsLastConverged,displacementsDependentDOFsLastConverged);
            if (curDispFactor==1)
                convergenceStatusLoadSteps=true;
            else
            {
                //eventually increase load step
                if (numNewtonIterations<MAXNUMNEWTONITERATIONS/3)
                {
                    deltaDispFactor*=1.5;
                }

                //increase displacement
                curDispFactor+=deltaDispFactor;
                if (curDispFactor>1)
                {
                    deltaDispFactor -= curDispFactor -1.;
                    curDispFactor=1;
                }

                curDisp = maxDisp*curDispFactor;

                //old stiffness matrix is used in first step of next load increment in order to prevent spurious problems at the boundary
                std::cout << "press enter to next load increment, delta disp factor " << deltaDispFactor << " max delta disp factor " <<  maxDeltaDispFactor << std::endl << std::endl;
                //char cDummy[100]="";
                //std::cin.getline(cDummy, 100);
            }
        }
		else
		{
            assert(convergenceStatus==2);
			//calculate stiffness of previous loadstep (used as initial stiffness in the next load step)
			//this is done within the loop in order to ensure, that for the first step the stiffness matrix of the previous step is used
			//otherwise, the additional boundary displacements will result in an artifical localization in elements at the boundary
			curDispFactor-=deltaDispFactor;
            curDisp = maxDisp*curDispFactor;

            myStructureCoarseScale.ConstraintSetRHS(ConstraintRHS,curDisp);

			// build global dof numbering
			myStructureCoarseScale.NodeBuildGlobalDofs();

            //set previous converged displacements
			myStructureCoarseScale.NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);
			myStructureCoarseScale.ElementTotalUpdateTmpStaticData();

			// calculate previous residual (should be almost zero)
			myStructureCoarseScale.BuildGlobalGradientInternalPotentialVector(intForceVector);

			//decrease load step
			deltaDispFactor*=0.5;
			curDispFactor+=deltaDispFactor;

			//check for minimum delta (this mostly indicates an error in the software
			if (deltaDispFactor<MIN_DELTA_STRAIN_FACTOR)
			{
			    throw NuTo::MechanicsException("Example ConcurrentMultiscale : No convergence, delta strain factor < 1e-7");
			}

			std::cout << "press enter to reduce load increment" << std::endl;
			//char cDummy[100]="";
			//std::cin.getline(cDummy, 100);;
		}

		if (!convergenceStatusLoadSteps)
		{
            //update new displacement of RHS
		    myStructureCoarseScale.ConstraintSetRHS(ConstraintRHS,curDisp);

            // build global dof numbering
            myStructureCoarseScale.NodeBuildGlobalDofs();

            //update stiffness in order to calculate new dispForceVector
            myStructureCoarseScale.BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
            int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
            int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
            std::cout << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << std::endl;

            //update displacements of all nodes according to the new conre mat
            NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
            NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
            myStructureCoarseScale.NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
            myStructureCoarseScale.NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
            myStructureCoarseScale.ElementTotalUpdateTmpStaticData();

            // calculate initial residual for next load step
            myStructureCoarseScale.BuildGlobalGradientInternalPotentialVector(intForceVector);

            //update rhs vector for next Newton iteration
            rhsVector = dispForceVector + extForceVector - intForceVector;

		}
    }

	if (PRINTRESULT)
	{
//        std::cout<< "numerical fracture energy "<< externalEnergy/(thickness*lY) << std::endl;
	}
}
catch (NuTo::Exception& e)
{
    std::cout << e.ErrorMessage() << std::endl;
}
return 0;
}
