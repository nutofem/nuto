#include "nuto/math/FullMatrix.h"
#include "nuto/math/FullVector.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

//just for test
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include <eigen3/Eigen/Core>

#define MAXNUMNEWTONITERATIONS 30
#define PRINTRESULT false

int main()
{
try
{
	//create structure
	NuTo::Structure myStructure(2);

	int NumElementX(5);
	int NumElementY(1);
	double lX(10);
	double lY(2);

	int NumNodeX(NumElementX+1);
	int NumNodeY(NumElementY+1);

	double deltaX=lX/NumElementX;
	double deltaY=lY/NumElementY;

	//create nodes
	NuTo::FullVector<double,Eigen::Dynamic> Coordinates(2);
    int LeftSupport = myStructure.GroupCreate("Nodes");
    int RightSupport = myStructure.GroupCreate("Nodes");
    double posX;
    double posY=0;
    int leftLowerNode;
	for (int theNodeY=0; theNodeY<NumNodeY;theNodeY++,posY+=deltaY)
	{
		posX=0.;
		for (int theNodeX=0 ; theNodeX<NumNodeX;theNodeX++, posX+=deltaX )
		{
			Coordinates(0) = posX;
			Coordinates(1) = posY;
			int theNode = myStructure.NodeCreate("displacements",Coordinates);
			if (posX==0)
			{
				if (posY==0)
					leftLowerNode = theNode;
				else
					myStructure.GroupAddNode(LeftSupport,theNode);
			}
			if (posX==lX)
				myStructure.GroupAddNode(RightSupport,theNode);
		}
	}

	//create elements
	NuTo::FullVector<int,Eigen::Dynamic> Incidence(4);
	for (int theElementY=0; theElementY<NumElementY;theElementY++)
	{
		for (int theElementX=0; theElementX<NumElementX;theElementX++)
		{
			Incidence(0) = theElementX+NumNodeX*theElementY;
			Incidence(1) = theElementX+1+NumNodeX*theElementY;
			Incidence(2) = theElementX+1+NumNodeX*(theElementY+1);
			Incidence(3) = theElementX+NumNodeX*(theElementY+1);
		    myStructure.ElementCreate("PLANE2D4N",Incidence,"ConstitutiveLawIpNonlocal","StaticDataNonlocal");
		}
	}

	//create constitutive law
	int myMatDamage = myStructure.ConstitutiveLawCreate("NonlocalDamagePlasticity");
	double YoungsModulus(20000);
	myStructure.ConstitutiveLawSetYoungsModulus(myMatDamage,YoungsModulus);
	myStructure.ConstitutiveLawSetPoissonsRatio(myMatDamage,0.2);
	myStructure.ConstitutiveLawSetNonlocalRadius(myMatDamage,deltaX*2);
	double fct(2);
	myStructure.ConstitutiveLawSetTensileStrength(myMatDamage,fct);
	myStructure.ConstitutiveLawSetCompressiveStrength(myMatDamage,fct*10);
	myStructure.ConstitutiveLawSetBiaxialCompressiveStrength(myMatDamage,fct*12.5);
	myStructure.ConstitutiveLawSetFractureEnergy(myMatDamage,0.1);

	//create section
	double thickness(0.5);
	int mySection = myStructure.SectionCreate("Plane_Strain");
	myStructure.SectionSetThickness(mySection,thickness);
	int mySectionweak = myStructure.SectionCreate("Plane_Strain");
	myStructure.SectionSetThickness(mySectionweak,0.99*thickness);

	//assign constitutive law 
	myStructure.ElementTotalSetSection(mySection);
	int weakElement((NumElementX-1)/2);
	//std::cout << "weak element " << weakElement << std::endl;
	myStructure.ElementSetSection(weakElement,mySectionweak);
	myStructure.ElementTotalSetConstitutiveLaw(myMatDamage);

	//Build nonlocal elements
	myStructure.BuildNonlocalData(myMatDamage);

	//fix left support
	NuTo::FullVector<double,Eigen::Dynamic> DirectionX(2);
	DirectionX.SetValue(0,1.0);
	DirectionX.SetValue(1,0.0);

	NuTo::FullVector<double,Eigen::Dynamic> DirectionY(2);
	DirectionY.SetValue(0,0.0);
	DirectionY.SetValue(1,1.0);

	myStructure.ConstraintLinearSetDisplacementNodeGroup(LeftSupport,DirectionX, 0);
	myStructure.ConstraintLinearSetDisplacementNode(leftLowerNode,DirectionX, 0);
	myStructure.ConstraintLinearSetDisplacementNode(leftLowerNode,DirectionY, 0);

	// update the RHS of the constrain equation with myStructure.ConstraintSetRHS
	int ConstraintRHS = myStructure.ConstraintLinearSetDisplacementNodeGroup(RightSupport,DirectionX, 0);

    // start analysis
    double maxDisp(200*fct/YoungsModulus*lX);
    double deltaMaxDisp(0.5*fct/YoungsModulus*lX);
    double curDisp(0.80*fct/YoungsModulus*lX);
	//char cDummy[100];

	#ifdef ENABLE_VISUALIZE
    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentEngineeringStrain();
    myStructure.AddVisualizationComponentEngineeringStress();
    myStructure.AddVisualizationComponentDamage();
    myStructure.AddVisualizationComponentEngineeringPlasticStrain();
#endif

    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> PlotData(1,6);
    double externalEnergy(0.);
    //increment for different load steps
	myStructure.ElementTotalUpdateTmpStaticData();
    while (curDisp<maxDisp)
    {
    	//update displacement of RHS
    	myStructure.ConstraintSetRHS(ConstraintRHS,curDisp);

    	if (PRINTRESULT)
    	    std::cout<<"curDisp " << curDisp << std::endl;
	    // build global dof numbering
	    myStructure.NodeBuildGlobalDofs();

        double normResidual(1);
        double maxResidual(1);
        int numNewtonIterations(0);
		while(maxResidual>1e-4 && numNewtonIterations<MAXNUMNEWTONITERATIONS)
		{
			numNewtonIterations++;
			// build global stiffness matrix and equivalent load vector which corresponds to prescribed boundary values
			NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrix;
			NuTo::FullVector<double,Eigen::Dynamic> dispForceVector;
			myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector);
			NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> stiffnessMatrixFull(stiffnessMatrix);
			//std::cout<<"stiffnessMatrixFull" << std::endl;
			//stiffnessMatrixFull.Info();
			//std::cout<<"dispForceVector" << std::endl;
			//dispForceVector.Trans().Info();
			stiffnessMatrix.RemoveZeroEntries(0,1e-14);

			// build global external load vector
			NuTo::FullVector<double,Eigen::Dynamic> extForceVector;
			myStructure.BuildGlobalExternalLoadVector(1,extForceVector);
			//std::cout<<"extForceVector" << std::endl;
			//extForceVector.Trans().Info();

			// build global internal load vector
			NuTo::FullVector<double,Eigen::Dynamic> intForceVector;
			myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector);

			// calculate right hand side
			NuTo::FullVector<double,Eigen::Dynamic> rhsVector = dispForceVector + extForceVector - intForceVector;
			//std::cout<<"rhsVector" << std::endl;
			//rhsVector.Trans().Info();

			// solve
			NuTo::SparseDirectSolverMUMPS mySolver;
			NuTo::FullVector<double,Eigen::Dynamic> deltaDisplacementsActiveDOFs;
			NuTo::FullVector<double,Eigen::Dynamic> oldDisplacementsActiveDOFs;
			NuTo::FullVector<double,Eigen::Dynamic> displacementsActiveDOFs;
			NuTo::FullVector<double,Eigen::Dynamic> displacementsDependentDOFs;
			stiffnessMatrix.SetOneBasedIndexing();
			mySolver.Solve(stiffnessMatrix, rhsVector, deltaDisplacementsActiveDOFs);
			//std::cout<<"deltaDisplacementsActiveDOFs" << std::endl;
			//deltaDisplacementsActiveDOFs.Trans().Info();
			//double normDeltaDisp = deltaDisplacementsActiveDOFs.Norm();
			//std::cout << "norm DeltaDisp: " << normDeltaDisp << std::endl;

			// write displacements to node
			myStructure.NodeExtractDofValues(oldDisplacementsActiveDOFs, displacementsDependentDOFs);

			//perform a linesearch
			NuTo::FullVector<double,Eigen::Dynamic> residualVector;
			double normRHS = rhsVector.Norm();
			double alpha(1);
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
			}
			while(alpha>1e-3 && normResidual>normRHS);
			maxResidual = residualVector.Max();

			//std::cout << "press enter to next step in Newton iteration, final alpha=" << 2*alpha << std::endl;
			//std::cin.getline(cDummy, 100);
			//std::cout << "Newton iteration, final alpha=" << 2*alpha << ", normResidual" << normResidual<< ", maxResidual" << maxResidual<<std::endl;
		}
		if (numNewtonIterations<MAXNUMNEWTONITERATIONS)
		{

			//NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> singleNodeDisp;
			//myStructure.NodeGetDisplacements(NumNodeX-1, singleNodeDisp);
			//std::cout<<"singleNodeDisp" << std::endl;
			//singleNodeDisp.Info();
			myStructure.ElementTotalUpdateStaticData();

			//calculate engineering plastic strain of myelement at all integration points
			NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> Stress;
			myStructure.ElementGetEngineeringStress(weakElement, Stress);
	    	if (PRINTRESULT)
	    	{
				std::cout << "stress in weak element " << weakElement << std::endl;
				Stress.Trans().Info();
				std::cout << "number of Newton iterations " << numNewtonIterations << ", delta disp "<< deltaMaxDisp << ", maxResidual " << maxResidual << ", normResidual " << normResidual << std::endl;
	    	}

			//store residual force
	    	NuTo::FullVector<double,Eigen::Dynamic> SupportingForce;
			myStructure.NodeGroupInternalForce(RightSupport,SupportingForce);
			NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> SinglePlotData(1,6);
			SinglePlotData(0,0) = curDisp;
			SinglePlotData(0,1) = SupportingForce(0,0)/(thickness*lY);
			SinglePlotData(0,2) = SupportingForce(0,0);
			SinglePlotData(0,3) = Stress(0,0);
			SinglePlotData(0,4) = Stress(1,0);
			SinglePlotData(0,5) = Stress(2,0);
		    externalEnergy+=deltaMaxDisp*SupportingForce(0,0);

			PlotData.AppendRows(SinglePlotData);

			// visualize results
#ifdef ENABLE_VISUALIZE
			myStructure.ExportVtkDataFile("NonlocalDamagePlasticityBeamView.vtk");
#endif
			//std::cout << "press enter to next load increment" << std::endl;
			//std::cin.getline(cDummy, 100);;
			if (numNewtonIterations<MAXNUMNEWTONITERATIONS/3)
				deltaMaxDisp*=1.5;
	        curDisp+=deltaMaxDisp;
			if (curDisp>maxDisp)
				curDisp=maxDisp;
		}
		else
		{
			deltaMaxDisp*=0.5;
			curDisp-=deltaMaxDisp;

		}
    }
	if (PRINTRESULT)
        std::cout<< "numerical fracture energy "<< externalEnergy/(thickness*lY) << std::endl;
    PlotData.WriteToFile("NonlocalDamagePlasticityBeamLoadDisp.txt"," ","#load displacement curve, disp, stress, force, sxx in center element, syy in center element","  ");
}
catch (NuTo::Exception& e)
{
    std::cout << e.ErrorMessage() << std::endl;
}
return 0;
}
