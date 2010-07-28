#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/MechanicsException.h"

//just for test
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/NonlocalDamagePlasticity.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"
#include <eigen2/Eigen/Core>

#define MAXNUMNEWTONITERATIONS 30
#define PRINTRESULT true

int main()
{
try
{
	//
	double lX(100);
	double lY(100);

	//create structure
	NuTo::Structure myStructure(2);
	myStructure.SetShowTime(true);

	NuTo::FullMatrix<int> createdGroupIds;
	myStructure.ImportFromGmsh("ImportGmsh.msh","displacements", "ConstitutiveLawIpNonlocal", "StaticDataNonlocal",createdGroupIds);
	myStructure.Info();

	//create constitutive law nonlocal damage
	int myMatDamage = myStructure.ConstitutiveLawCreate("NonlocalDamagePlasticity");
	double YoungsModulusDamage(20000);
	myStructure.ConstitutiveLawSetYoungsModulus(myMatDamage,YoungsModulusDamage);
	myStructure.ConstitutiveLawSetPoissonsRatio(myMatDamage,0.2);
	myStructure.ConstitutiveLawSetNonlocalRadius(myMatDamage,10);
	double fct(2);
	myStructure.ConstitutiveLawSetTensileStrength(myMatDamage,fct);
	myStructure.ConstitutiveLawSetCompressiveStrength(myMatDamage,fct*10);
	myStructure.ConstitutiveLawSetBiaxialCompressiveStrength(myMatDamage,fct*12.5);
	myStructure.ConstitutiveLawSetFractureEnergy(myMatDamage,0.1);

	//create constitutive law linear elastic
	int myMatLinear = myStructure.ConstitutiveLawCreate("LinearElastic");
	double YoungsModulusLE(10000);
	myStructure.ConstitutiveLawSetYoungsModulus(myMatLinear,YoungsModulusLE);
	myStructure.ConstitutiveLawSetPoissonsRatio(myMatLinear,0.2);

	//create section
	double thickness(0.5);
	//int mySectionParticle = myStructure.SectionCreate("Plane_Strain");
	//myStructure.SectionSetThickness(mySectionParticle,thickness);

	int mySectionMatrix = myStructure.SectionCreate("Plane_Strain");
	myStructure.SectionSetThickness(mySectionMatrix,thickness);

	//assign constitutive law 
	//myStructure.ElementGroupSetSection(createdGroupIds(0,0),mySectionParticle);
	myStructure.ElementGroupSetSection(createdGroupIds(1,0),mySectionMatrix);
	myStructure.ElementGroupSetConstitutiveLaw(102,myMatDamage);
	//myStructure.ElementGroupSetConstitutiveLaw(101,myMatLinear);
	bool deleteNodes=true;
	myStructure.ElementGroupDelete(101,deleteNodes);

	//Build nonlocal elements
	myStructure.BuildNonlocalData(myMatDamage);

	//Create groups to apply the boundary conditions
	//left boundary
	int GrpNodes_LeftBoundary = myStructure.GroupCreate("Nodes");
	int direction = 0; //either 0,1,2
	double min(0.);
	double max(0.);
	myStructure.GroupAddNodeCoordinateRange(GrpNodes_LeftBoundary,direction,min,max);

	//lower left node
	int GrpNodes_LowerBoundary = myStructure.GroupCreate("Nodes");
	direction=1;
	min=0;
	max=0;
	myStructure.GroupAddNodeCoordinateRange(GrpNodes_LowerBoundary,direction,min,max);
	int GrpNodes_LowerLeftNode = myStructure.GroupIntersection(GrpNodes_LowerBoundary,GrpNodes_LeftBoundary);

	//right boundary
	int GrpNodes_RightBoundary = myStructure.GroupCreate("Nodes");
	direction = 0; //either 0,1,2
	min=lX;
	max=lX;
	myStructure.GroupAddNodeCoordinateRange(GrpNodes_RightBoundary,direction,min,max);


	//fix left support
	NuTo::FullMatrix<double> DirectionX(2,1);
	DirectionX.SetValue(0,0,1.0);
	DirectionX.SetValue(1,0,0.0);

	NuTo::FullMatrix<double>DirectionY(2,1);
	DirectionY.SetValue(0,0,0.0);
	DirectionY.SetValue(1,0,1.0);

	myStructure.ConstraintSetDisplacementNodeGroup(GrpNodes_LeftBoundary,DirectionX, 0);
	myStructure.ConstraintSetDisplacementNodeGroup(GrpNodes_LowerLeftNode,DirectionY, 0);

	// update the RHS of the constrain equation with myStructure.ConstraintSetRHS
	int ConstraintRHS = myStructure.ConstraintSetDisplacementNodeGroup(GrpNodes_RightBoundary,DirectionX, 0);

    // start analysis
    double maxDisp(20*fct/YoungsModulusDamage*lX);
    double deltaDisp(0.02*fct/YoungsModulusDamage*lX);
    double curDisp(0.4*fct/YoungsModulusDamage*lX);

#ifdef ENABLE_VISUALIZE
    myStructure.AddVisualizationComponentSection();
    myStructure.AddVisualizationComponentConstitutive();
    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentEngineeringStrain();
    myStructure.AddVisualizationComponentEngineeringStress();
    myStructure.AddVisualizationComponentDamage();
    myStructure.AddVisualizationComponentEngineeringPlasticStrain();
	myStructure.ElementTotalUpdateTmpStaticData();
    myStructure.ExportVtkDataFile("ImportGmsh.vtk");
#endif


    NuTo::FullMatrix<double> PlotData(1,6);
    double externalEnergy(0.);
    //increment for different load steps
	myStructure.ElementTotalUpdateTmpStaticData();
	NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixCSRVector2;
	NuTo::FullMatrix<double> dispForceVector;

	//allocate solver
	NuTo::SparseDirectSolverMUMPS mySolver;
	mySolver.SetShowTime(true);

    //calculate absolute tolerance for matrix entries to be not considered as zero
	myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
	double maxValue, minValue, ToleranceZeroStiffness;
	stiffnessMatrixCSRVector2.Max(maxValue);
	stiffnessMatrixCSRVector2.Min(minValue);
	std::cout << "min and max " << minValue << " , " << maxValue << std::endl;

	ToleranceZeroStiffness = (1e-14) * (fabs(maxValue)>fabs(minValue) ?  fabs(maxValue) : fabs(minValue));
	myStructure.SetToleranceStiffnessEntries(ToleranceZeroStiffness);
	int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
	int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
	std::cout << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << std::endl;

    while (curDisp<maxDisp)
    {
    	//update displacement of RHS
    	myStructure.ConstraintSetRHS(ConstraintRHS,curDisp);

    	if (PRINTRESULT)
    	    std::cout<<"curDisp " << curDisp << ", deltaDisp " << deltaDisp << ", maxDisp " <<  maxDisp << std::endl;
	    // build global dof numbering
	    myStructure.NodeBuildGlobalDofs();

        double normResidual(1);
        double maxResidual(1);
        int numNewtonIterations(0);
		double normRHS(1.);
		double alpha(1.);
		while(normRHS>1e-5 && maxResidual>1e-5 && numNewtonIterations<MAXNUMNEWTONITERATIONS)
		{
			numNewtonIterations++;
			// build global stiffness matrix and equivalent load vector which corresponds to prescribed boundary values
			myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
			//NuTo::FullMatrix<double> stiffnessMatrixFull(stiffnessMatrix);
			//std::cout<<"stiffnessMatrixFull" << std::endl;
			//stiffnessMatrixFull.Info();
			//std::cout<<"dispForceVector" << std::endl;
			//dispForceVector.Trans().Info();
			int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
			int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
			std::cout << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << std::endl;

			// build global external load vector
			NuTo::FullMatrix<double> extForceVector;
			myStructure.BuildGlobalExternalLoadVector(extForceVector);
			//std::cout<<"extForceVector" << std::endl;
			//extForceVector.Trans().Info();

			// build global internal load vector
			NuTo::FullMatrix<double> intForceVector;
			myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector);

			// calculate right hand side
			NuTo::FullMatrix<double> rhsVector = dispForceVector + extForceVector - intForceVector;
			//std::cout<<"rhsVector" << std::endl;
			//rhsVector.Trans().Info();

			// solve
			NuTo::FullMatrix<double> deltaDisplacementsActiveDOFs;
			NuTo::FullMatrix<double> oldDisplacementsActiveDOFs;
			NuTo::FullMatrix<double> displacementsActiveDOFs;
			NuTo::FullMatrix<double> displacementsDependentDOFs;
			NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrixCSR(stiffnessMatrixCSRVector2);
			stiffnessMatrixCSR.SetOneBasedIndexing();
			mySolver.Solve(stiffnessMatrixCSR, rhsVector, deltaDisplacementsActiveDOFs);
			//std::cout<<"deltaDisplacementsActiveDOFs" << std::endl;
			//deltaDisplacementsActiveDOFs.Trans().Info();
			//double normDeltaDisp = deltaDisplacementsActiveDOFs.Norm();
			//std::cout << "norm DeltaDisp: " << normDeltaDisp << std::endl;

			// write displacements to node
			myStructure.NodeExtractDofValues(oldDisplacementsActiveDOFs, displacementsDependentDOFs);
			normRHS = rhsVector.Norm();
			//perform a linesearch
			NuTo::FullMatrix<double> residualVector;
			alpha = 1.;
			do
			{
				displacementsActiveDOFs = oldDisplacementsActiveDOFs + deltaDisplacementsActiveDOFs*alpha;
				myStructure.NodeMergeActiveDofValues(displacementsActiveDOFs);
				myStructure.ElementTotalUpdateTmpStaticData();
				//std::cout<<"displacementsActiveDOFs" << std::endl;
				//displacementsActiveDOFs.Trans().Info();

				// calculate residual
				myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector);
				residualVector = extForceVector - intForceVector;
				normResidual = residualVector.Norm();
				alpha*=0.5;
				std::cout << "normResidual t " << normResidual << "normResidual 0 "<< normRHS << "normRHS*(1-0.5*alpha)" << normRHS*(1-0.5*alpha) << std::endl;
			}
			while(alpha>1e-3 && normResidual>normRHS*(1-0.5*alpha));
			maxResidual = residualVector.Max();

			std::cout << std::endl << "Newton iteration, final alpha=" << 2*alpha << ", normResidual" << normResidual<< ", maxResidual" << maxResidual<<std::endl;
			//char cDummy[100]="";
			//std::cin.getline(cDummy, 100);
		}
		if (numNewtonIterations<MAXNUMNEWTONITERATIONS || alpha>0.2)
		{

			//NuTo::FullMatrix<double> singleNodeDisp;
			//myStructure.NodeGetDisplacements(NumNodeX-1, singleNodeDisp);
			//std::cout<<"singleNodeDisp" << std::endl;
			//singleNodeDisp.Info();
			myStructure.ElementTotalUpdateStaticData();

			//store residual force
			NuTo::FullMatrix<double> SupportingForce;
			myStructure.NodeGroupInternalForce(GrpNodes_RightBoundary,SupportingForce);
			NuTo::FullMatrix<double> SinglePlotData(1,6);
			SinglePlotData(0,0) = curDisp;
			SinglePlotData(0,1) = SupportingForce(0,0)/(thickness*lY);
			SinglePlotData(0,2) = SupportingForce(0,0);
		    externalEnergy+=deltaDisp*SupportingForce(0,0);

			PlotData.AppendRows(SinglePlotData);
		    PlotData.WriteToFile("ImportGmshLoadDisp.txt"," ","#load displacement curve, disp, stress, force, sxx in center element, syy in center element","  ");

			// visualize results
#ifdef ENABLE_VISUALIZE
			myStructure.ExportVtkDataFile("ImportGmsh.vtk");
#endif
			if (numNewtonIterations<MAXNUMNEWTONITERATIONS/3)
				deltaDisp*=1.5;
	        curDisp+=deltaDisp;
			if (curDisp>maxDisp)
				curDisp=maxDisp;

			if (PRINTRESULT)
	    	{
				std::cout << "number of Newton iterations " << numNewtonIterations << ", delta disp "<< deltaDisp << ", maxResidual " << maxResidual << ", normResidual " << normResidual << std::endl;
	    	}
			//std::cout << "press enter to next load increment" << std::endl;
			//char cDummy[100]="";
			//std::cin.getline(cDummy, 100);;
		}
		else
		{
			deltaDisp*=0.5;
			curDisp-=deltaDisp;
			if (deltaDisp<1e-7)
				throw NuTo::MechanicsException("Example ImportGmsh : No convergence, delta disp < 1e-7");

		}
    }
	if (PRINTRESULT)
        std::cout<< "numerical fracture energy "<< externalEnergy/(thickness*lY) << std::endl;
}
catch (NuTo::Exception& e)
{
    std::cout << e.ErrorMessage() << std::endl;
}
return 0;
}
