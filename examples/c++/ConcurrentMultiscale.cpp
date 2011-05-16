#include <stdlib.h>
#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/StructureMultiscale.h"
#include "nuto/base/Exception.h"

#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"

#define MAXNUMNEWTONITERATIONS 20
//#define MAXNORMRHS 100
#define PRINTRESULT true
#define MIN_DELTA_STRAIN_FACTOR 1e-7

//void InitStructure(NuTo::Structure& myStructureCoarseScale);
//void Solve(NuTo::Structure& myStructureCoarseScale);
//void SetLoadFactor(NuTo::Structure& myStructureCoarseScale, double rLoadFactor);
//void StoreResults(int rIteration, double rLoadfactor, NuTo::FullMatrix<double>& rResult, NuTo::Structure& myStructureCoarseScale);

class MyStructureClass : public NuTo::Structure
{
public:
    //! @brief ... default constructor
    MyStructureClass(int rDimension) : NuTo::Structure(rDimension)
    {
        //create the Structure
        mlX = 200;
        mlY = 200;

        //create structure
#ifdef SHOW_TIME
		SetShowTime(true);
#endif //SHOW_TIME
        //add visualization
        //myStructureCoarseScale.AddVisualizationComponentSection();
        //myStructureCoarseScale.AddVisualizationComponentConstitutive();
        AddVisualizationComponentDisplacements();
        AddVisualizationComponentEngineeringStrain();
        AddVisualizationComponentEngineeringStress();
        //myStructureCoarseScale.AddVisualizationComponentDamage();
        //myStructureCoarseScale.AddVisualizationComponentEngineeringPlasticStrain();
        //myStructureCoarseScale.AddVisualizationComponentPrincipalEngineeringStress();

        //2x2 nodes 1x1 element grid
        //create nodes
        NuTo::FullMatrix<double> Coordinates(2,1);
        Coordinates(0,0) = 0.0;
        Coordinates(1,0) = 0.0;
        mNode1 = NodeCreate("displacements",Coordinates);

        Coordinates(0,0) = mlX;
        Coordinates(1,0) = 0.0;
        mNode2 = NodeCreate("displacements",Coordinates);

        Coordinates(0,0) = mlX;
        Coordinates(1,0) = mlY;
        mNode3 = NodeCreate("displacements",Coordinates);

        Coordinates(0,0) = 0.0;
        Coordinates(1,0) = mlY;
        mNode4 = NodeCreate("displacements",Coordinates);

        //create elements
        NuTo::FullMatrix<int> Incidence(4,1);
        Incidence(0,0) = mNode1;
        Incidence(1,0) = mNode2;
        Incidence(2,0) = mNode3;
        Incidence(3,0) = mNode4;
        int myElement1 = ElementCreate("PLANE2D4N",Incidence,"ConstitutiveLawIp","StaticData");

        ElementSetIntegrationType(myElement1,"2D4NGauss1Ip","StaticData");

        //create constitutive law from fine scale model
        int myMatCoarseScale = ConstitutiveLawCreate("Multiscale");
        ConstitutiveLawSetYoungsModulus(myMatCoarseScale,20000);
        ConstitutiveLawSetPoissonsRatio(myMatCoarseScale,0.0);

        //create section
        double thickness(1);
        int mySectionCoarseScale = SectionCreate("Plane_Strain");
        SectionSetThickness(mySectionCoarseScale,thickness);

        ElementTotalSetSection(mySectionCoarseScale);
        ElementTotalSetConstitutiveLaw(myMatCoarseScale);

        //set fine scale model as ip data for all integration points
        //create the fine scale model and save it as a binary
        std::string nameOfBinaryFinescale("/home/unger3/develop/nuto_build/examples/c++/myStructureFineScale.bin");
        //std::string nameOfBinaryFinescale("myStructureFineScale.bin");
        CreateAndSaveFineScaleModel(nameOfBinaryFinescale);

        //myStructureCoarseScale.ElementTotalSetFineScaleModel("myStructureFineScale.bin");
        ElementTotalSetFineScaleModel(nameOfBinaryFinescale);

        resultMatrix.Resize(0,5);
    }

    //! @brief temporarily creates the fine scale model and saves it as a binary
    //! this is then restored for each integration point of the fine scale model
    void CreateAndSaveFineScaleModel(std::string rNameOfBinaryFinescale)
    {
    	try
    	{
    		//
    	    bool square(false);
    	    bool allNodesAsMultiscale(false);

    	    //create structure
    	    NuTo::StructureMultiscale myStructureFineScale(2);
#ifdef SHOW_TIME
   	    myStructureFineScale.SetShowTime(false);
#endif //SHOW_TIME
    	    int GroupNodesBoundaryDamage(0);
    	    int GroupNodesBoundaryHomogeneous(0);
    	    int GroupNodesMultiscaleDamage(0);
    	    int GroupNodesMultiscaleHomogeneous(0);
    	    int GroupElementsDamage(0);
    	    int GroupElementsHomogeneous(0);

			//create section
			double thickness(1);
			int mySectionParticle = myStructureFineScale.SectionCreate("Plane_Strain");
			myStructureFineScale.SectionSetThickness(mySectionParticle,thickness);

			int mySectionMatrix = myStructureFineScale.SectionCreate("Plane_Strain");
			myStructureFineScale.SectionSetThickness(mySectionMatrix,thickness);

			//create damaged (structureTypeCount=0) and (homogeneous(structureTypeCount=1) fine scale structure
    	    double shiftXDirection;
			for (int structureTypeCount=0; structureTypeCount<2; structureTypeCount++)
    	    {
				if (structureTypeCount==0)
					shiftXDirection = -60;
				else
					shiftXDirection = 60;

				NuTo::FullMatrix<int> createdGroupIds;

				if (structureTypeCount==0)
				{
					//damage model
					if (square)
					{
						system("gmsh -2 -order 1 /home/unger3/develop/nuto/examples/c++/ConcurrentMultiscaleFineScaleSquareDamage.geo");
						myStructureFineScale.ImportFromGmsh("/home/unger3/develop/nuto/examples/c++/ConcurrentMultiscaleFineScaleSquareDamage.msh","displacements", "ConstitutiveLawIpNonlocal", "StaticDataNonlocal",createdGroupIds);

				    }
					else
					{
						system("gmsh -2 -order 1 /home/unger3/develop/nuto/examples/c++/ConcurrentMultiscaleFineScaleRoundDamage.geo");
						myStructureFineScale.ImportFromGmsh("/home/unger3/develop/nuto/examples/c++/ConcurrentMultiscaleFineScaleRoundDamage.msh","displacements", "ConstitutiveLawIpNonlocal", "StaticDataNonlocal",createdGroupIds);
					}
				}
				else
				{
					//homogeneous model
					if (square)
					{
						system("gmsh -2 -order 1 /home/unger3/develop/nuto/examples/c++/ConcurrentMultiscaleFineScaleSquareHomogeneous.geo");
						myStructureFineScale.ImportFromGmsh("/home/unger3/develop/nuto/examples/c++/ConcurrentMultiscaleFineScaleSquareHomogeneous.msh","displacements", "ConstitutiveLawIpNonlocal", "StaticDataNonlocal",createdGroupIds);
					}
					else
					{
						system("gmsh -2 -order 1 /home/unger3/develop/nuto/examples/c++/ConcurrentMultiscaleFineScaleRoundHomogeneous.geo");
						myStructureFineScale.ImportFromGmsh("/home/unger3/develop/nuto/examples/c++/ConcurrentMultiscaleFineScaleRoundHomogeneous.msh","displacements", "ConstitutiveLawIpNonlocal", "StaticDataNonlocal",createdGroupIds);
					}
				}

				//create constitutive law nonlocal damage
				int myMatDamage = myStructureFineScale.ConstitutiveLawCreate("NonlocalDamagePlasticity");
				double YoungsModulusDamage(20000);
				myStructureFineScale.ConstitutiveLawSetYoungsModulus(myMatDamage,YoungsModulusDamage);
				myStructureFineScale.ConstitutiveLawSetPoissonsRatio(myMatDamage,0.0);
				double nonlocalRadius(8);
				myStructureFineScale.ConstitutiveLawSetNonlocalRadius(myMatDamage,nonlocalRadius);
				double fct(2);
				myStructureFineScale.ConstitutiveLawSetTensileStrength(myMatDamage,fct);
				myStructureFineScale.ConstitutiveLawSetCompressiveStrength(myMatDamage,fct*10);
				myStructureFineScale.ConstitutiveLawSetBiaxialCompressiveStrength(myMatDamage,fct*12.5);
				myStructureFineScale.ConstitutiveLawSetFractureEnergy(myMatDamage,.2);

				//create constitutive law linear elastic (finally not used, since the elements are deleted)
				int myMatLinear = myStructureFineScale.ConstitutiveLawCreate("LinearElastic");
				double YoungsModulusLE(20000);
				myStructureFineScale.ConstitutiveLawSetYoungsModulus(myMatLinear,YoungsModulusLE);
				myStructureFineScale.ConstitutiveLawSetPoissonsRatio(myMatLinear,0.0);

				int groupIdMatrix= createdGroupIds(0,0);
				//assign constitutive law
			    myStructureFineScale.ElementGroupSetSection(groupIdMatrix,mySectionMatrix);
			    //myStructureFineScale.ElementGroupSetConstitutiveLaw(groupIdMatrix,myMatLinear);
			    myStructureFineScale.ElementGroupSetConstitutiveLaw(groupIdMatrix,myMatDamage);

				//Build nonlocal elements
				myStructureFineScale.BuildNonlocalData(myMatDamage);

				int GrpNodes_Boundary;
				int GroupNodes_Multiscale;
				int GrpNodes_Model;
				if (square)
				{
					//Create groups to apply the boundary conditions
					double lX(100);
					double lY(100);
					//left boundary
					int GrpNodes_Left = myStructureFineScale.GroupCreate("Nodes");
					int direction = 0; //either 0,1,2
					double min(shiftXDirection+0.);
					double max(shiftXDirection+0.);
					myStructureFineScale.GroupAddNodeCoordinateRange(GrpNodes_Left,direction,min,max);

					//right boundary
					int GrpNodes_Right = myStructureFineScale.GroupCreate("Nodes");
					direction = 0; //either 0,1,2
					min=shiftXDirection+lX;
					max=shiftXDirection+lX;
					myStructureFineScale.GroupAddNodeCoordinateRange(GrpNodes_Right,direction,min,max);

					//top boundary
					int GrpNodes_Top = myStructureFineScale.GroupCreate("Nodes");
					direction=1;
					min=lY;
					max=lY;
					myStructureFineScale.GroupAddNodeCoordinateRange(GrpNodes_Top,direction,min,max);

					//bottom boundary
					int GrpNodes_Bottom = myStructureFineScale.GroupCreate("Nodes");
					direction=1;
					min=0;
					max=0;
					myStructureFineScale.GroupAddNodeCoordinateRange(GrpNodes_Bottom,direction,min,max);

					//join the groups
					int GrpNodes_BottomTop = myStructureFineScale.GroupUnion(GrpNodes_Bottom,GrpNodes_Top);
					int GrpNodes_LeftRight = myStructureFineScale.GroupUnion(GrpNodes_Left,GrpNodes_Right);
					int GrpNodes_Boundarytmp = myStructureFineScale.GroupUnion(GrpNodes_BottomTop,GrpNodes_LeftRight);

					//intersect with the nodes which are in the xrange from min to max
					GrpNodes_Model = myStructureFineScale.GroupCreate("Nodes");
					direction = 0; //either 0,1,2
					min=shiftXDirection+0.;
					max=shiftXDirection+lX;
					myStructureFineScale.GroupAddNodeCoordinateRange(GrpNodes_Model,direction,min,max);

					GrpNodes_Boundary = myStructureFineScale.GroupIntersection(GrpNodes_Boundarytmp,GrpNodes_Model);

/*					std::cout << "all nodes are boundary nodes!!!!!!!!!!" << std::endl;
					direction=1;
					min=-1;
					max=101;
					GrpNodes_Boundary = myStructureFineScale.GroupCreate("Nodes");
					myStructureFineScale.GroupAddNodeCoordinateRange(GrpNodes_Boundary,direction,min,max);
*/
				}
				else
				{
					GrpNodes_Boundary = myStructureFineScale.GroupCreate("Nodes");
					NuTo::FullMatrix<double> center(2,1);
					center(0,0) = shiftXDirection+0.;
					center(1,0) = 0.;
			        double rmin=49.999;
					//double rmin=0;
					double rmax=50.001;
					myStructureFineScale.GroupAddNodeRadiusRange(GrpNodes_Boundary,center,rmin,rmax);

					GrpNodes_Model = myStructureFineScale.GroupCreate("Nodes");
			        rmin=0.;
					rmax=50.001;
					myStructureFineScale.GroupAddNodeRadiusRange(GrpNodes_Model,center,rmin,rmax);
				}

				if (allNodesAsMultiscale)
				{
					GroupNodes_Multiscale = GrpNodes_Model;
				}
				else
				{
					GroupNodes_Multiscale = GrpNodes_Boundary;
				}

				if (structureTypeCount==0)
				{
					GroupNodesBoundaryDamage = GrpNodes_Boundary;
	    	        GroupElementsDamage = groupIdMatrix;
	    	        GroupNodesMultiscaleDamage = GroupNodes_Multiscale;
				}
				else
				{
					GroupNodesBoundaryHomogeneous = GrpNodes_Boundary;
    	            GroupElementsHomogeneous = groupIdMatrix;
	    	        GroupNodesMultiscaleHomogeneous = GroupNodes_Multiscale;
				}
    	    }
			myStructureFineScale.SetGroupBoundaryNodesElements(GroupNodesBoundaryDamage,GroupNodesBoundaryHomogeneous,
					GroupNodesMultiscaleDamage,GroupNodesMultiscaleHomogeneous,GroupElementsDamage,GroupElementsHomogeneous);
    		//update conre mat
    		myStructureFineScale.NodeBuildGlobalDofs();

    		myStructureFineScale.AddVisualizationComponentSection();
    		myStructureFineScale.AddVisualizationComponentConstitutive();
    		myStructureFineScale.AddVisualizationComponentDisplacements();
    		myStructureFineScale.AddVisualizationComponentEngineeringStrain();
    		myStructureFineScale.AddVisualizationComponentEngineeringStress();
    		myStructureFineScale.AddVisualizationComponentDamage();
    		myStructureFineScale.AddVisualizationComponentEngineeringPlasticStrain();
    		myStructureFineScale.AddVisualizationComponentPrincipalEngineeringStress();
    		myStructureFineScale.Info();

    	    //myStructureFineScale.Save("myStructureFineScale.xml","xml");
#ifdef ENABLE_SERIALIZATION
    		myStructureFineScale.Save(rNameOfBinaryFinescale,"binary");
#else
    		throw NuTo::Exception("ConcurrentMultiscale - example requires serialization");
#endif
    		//myStructureFineScale.Restore(rNameOfBinaryFinescale,"binary");
    	    myStructureFineScale.ElementTotalUpdateTmpStaticData();
    	    myStructureFineScale.ExportVtkDataFile(std::string("/home/unger3/develop/nuto_build/examples/c++/ConcurrentMultiscaleInitialFineScaleModel.vtk"));
    	}
    	catch (NuTo::Exception& e)
    	{
    	    std::cout << e.ErrorMessage() << std::endl;
    	}
    }

    //! @brief set the load factor (load or displacement control)
    //! @param load factor
    void SetLoadFactor(double rLoadFactor)
    {
        this->ConstraintSetRHS(mConstraint1X,0);
        this->ConstraintSetRHS(mConstraint1Y,0);
        this->ConstraintSetRHS(mConstraint2X,rLoadFactor*mStrain(0,0)*mlX);
        this->ConstraintSetRHS(mConstraint2Y,rLoadFactor*0.5*mStrain(2,0)*mlX);
        this->ConstraintSetRHS(mConstraint3X,rLoadFactor*(mStrain(0,0)*mlX+0.5*mStrain(2,0)*mlY));
        this->ConstraintSetRHS(mConstraint3Y,rLoadFactor*(mStrain(1,0)*mlY+0.5*mStrain(2,0)*mlX));
        this->ConstraintSetRHS(mConstraint4X,rLoadFactor*0.5*mStrain(2,0)*mlY);
        this->ConstraintSetRHS(mConstraint4Y,rLoadFactor*mStrain(1,0)*mlY);
    }

    void SetTotalStrain(NuTo::FullMatrix<double> rStrain)
    {
        mStrain = rStrain;
    }

    // set the initial constraints and the constraint number of the right boundary to be modified (used in the SetLoadFactorRoutine)
    // do the initial step using the linear elastic solution (no fine scale model just the calibration of the initial angles)
    void InitStructure()
    {
        //directionX
        NuTo::FullMatrix<double> DirectionX(2,1);
        DirectionX.SetValue(0,0,1.0);
        NuTo::FullMatrix<double> DirectionY(2,1);
        DirectionY.SetValue(1,0,1.0);

        mConstraint1X = this->ConstraintLinearSetDisplacementNode(mNode1,DirectionX,0);
        mConstraint1Y = this->ConstraintLinearSetDisplacementNode(mNode1,DirectionY,0);

        mConstraint2X = this->ConstraintLinearSetDisplacementNode(mNode2,DirectionX,0);
        mConstraint2Y = this->ConstraintLinearSetDisplacementNode(mNode2,DirectionY,0);

        mConstraint3X = this->ConstraintLinearSetDisplacementNode(mNode3,DirectionX,0);
        mConstraint3Y = this->ConstraintLinearSetDisplacementNode(mNode3,DirectionY,0);

        mConstraint4X = this->ConstraintLinearSetDisplacementNode(mNode4,DirectionX,0);
        mConstraint4Y = this->ConstraintLinearSetDisplacementNode(mNode4,DirectionY,0);

        this->SetLoadFactor(1.);
        //update conre mat
        this->NodeBuildGlobalDofs();

        //build stiffness matrix
        NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrix;
        NuTo::FullMatrix<double> dispForceVector;
        this->BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector);
        //NuTo::FullMatrix<double> stiffnessMatrixFull(stiffnessMatrix);
        //std::cout << "stiffness" << std::endl;
        //stiffnessMatrixFull.Info(12,3);

        // build global external load vector
        NuTo::FullMatrix<double> extForceVector;
        this->BuildGlobalExternalLoadVector(extForceVector);
        //std::cout<<"extForceVector" << std::endl;
        //extForceVector.Trans().Info();

        // build global internal load vector
        NuTo::FullMatrix<double> intForceVector;
        this->BuildGlobalGradientInternalPotentialVector(intForceVector);

        // calculate right hand side
        NuTo::FullMatrix<double> rhsVector = dispForceVector + extForceVector - intForceVector;
        //std::cout<<"rhsVector" << std::endl;
        //rhsVector.Trans().Info();

        // solve
        NuTo::SparseDirectSolverMUMPS mySolver;
        NuTo::FullMatrix<double> deltaDisplacementsActiveDOFs;
        NuTo::FullMatrix<double> oldDisplacementsActiveDOFs;
        NuTo::FullMatrix<double> displacementsActiveDOFs;
        NuTo::FullMatrix<double> displacementsDependentDOFs;
        stiffnessMatrix.SetOneBasedIndexing();
        if (stiffnessMatrix.GetNumRows()!=0)
            mySolver.Solve(stiffnessMatrix, rhsVector, deltaDisplacementsActiveDOFs);
        else
            deltaDisplacementsActiveDOFs.Resize(0,1);
        //std::cout<<"deltaDisplacementsActiveDOFs" << std::endl;
        //deltaDisplacementsActiveDOFs.Trans().Info();
        //double normDeltaDisp = deltaDisplacementsActiveDOFs.Norm();
        //std::cout << "norm DeltaDisp: " << normDeltaDisp << std::endl;

        // write displacements to node
        this->NodeExtractDofValues(oldDisplacementsActiveDOFs, displacementsDependentDOFs);

        //no linesearch necessary
        NuTo::FullMatrix<double> residualVector;
        displacementsActiveDOFs = oldDisplacementsActiveDOFs + deltaDisplacementsActiveDOFs;
        this->NodeMergeActiveDofValues(displacementsActiveDOFs);
        this->ElementTotalUpdateTmpStaticData();

        // calculate residual
        this->BuildGlobalGradientInternalPotentialVector(intForceVector);
        residualVector = extForceVector - intForceVector;
        double normResidual = residualVector.Norm();
        if (normResidual>1e-6)
        {
            throw NuTo::MechanicsException("[main::NewtonRaphsonAuxRoutines::InitStructure] norm of the residual should be zero in the first elastic iteration.");
        }
        this->ExportVtkDataFile(std::string("/home/unger3/develop/nuto_build/examples/c++/CoarseScaleInitialSolution.vtk"));

        //update structure - set total strain of the linear solution being the total strain of the fine scale solution
        this->ElementTotalUpdateStaticData();
        //initialize the nonlinear solution procedure, init the initial crackangle
        this->ElementTotalSetFineScaleParameter("UseNonlinearSolution",0);

    }

    //! @brief do a postprocessing step after each converged load step
    void PostProcessDataAfterConvergence(int rLoadStep, int rNumNewtonIterations, double rLoadFactor, double rDeltaLoadFactor)const
    {
        std::cout << " Macroscale Convergence after " << rNumNewtonIterations << " Newton iterations, curLoadFactor " << rLoadFactor << ", deltaLoadFactor "<< rDeltaLoadFactor << std::endl<< std::endl;
        std::stringstream ssLoadStep;
        ssLoadStep << rLoadStep;
        std::stringstream ssIteration;
        ssIteration << rNumNewtonIterations;
        ExportVtkDataFile(std::string("/home/unger3/develop/nuto_build/examples/c++/MacroscaleConcurrentConverged") + ssLoadStep.str()+"_" + ssIteration.str() + std::string(".vtk"));
    }

    //! @brief do a postprocessing step after each line search within the load step
    void PostProcessDataAfterLineSearch(int rLoadStep, int rNewtonIteration, double rLineSearchFactor, double rLoadFactor)const
    {
        std::cout << " Macroscale step, iteration " << rNewtonIteration <<
                     ", line search factor " << rLineSearchFactor <<
                     ", load factor " << rLoadFactor << std::endl;
        std::stringstream ssLoadStep;
        ssLoadStep << rLoadStep;
        std::stringstream ssIteration;
        ssIteration << rNewtonIteration;
        this->ExportVtkDataFile(std::string("/home/unger3/develop/nuto_build/examples/c++/MacroscaleConcurrentMultiscale") + ssLoadStep.str()+"_" + ssIteration.str() + std::string(".vtk"));
    }

protected:
    NuTo::FullMatrix<double> mStrain;
    int mConstraint1X;
    int mConstraint1Y;
    int mConstraint2X;
    int mConstraint2Y;
    int mConstraint3X;
    int mConstraint3Y;
    int mConstraint4X;
    int mConstraint4Y;
   double mlX;
    double mlY;
    int mNode1;
    int mNode2;
    int mNode3;
    int mNode4;
    NuTo::FullMatrix<double> resultMatrix;

};

int main()
{
try
{
    MyStructureClass myStructureCoarseScale(2);

    //crack transition zone (radius)
    myStructureCoarseScale.ElementTotalSetFineScaleParameter("CrackTransitionRadius",1);
    myStructureCoarseScale.ElementTotalSetFineScaleParameter("SquareCoarseScaleModel",0);

    //penalty stiffness for crack angle
    double PenaltyStiffnessCrackAngle(1);
    double PenaltyStiffnessScalingFactorCrackAngle(2.*M_PI);
    myStructureCoarseScale.ElementTotalSetFineScaleParameter("ConstraintPenaltyStiffnessCrackAngle", PenaltyStiffnessCrackAngle);
    myStructureCoarseScale.ElementTotalSetFineScaleParameter("PenaltyStiffnessScalingFactorCrackAngle", PenaltyStiffnessScalingFactorCrackAngle);

    //penalty stiffness for tangential crack opening
    //double PenaltyStiffnessTangentialCrackOpening(1);
    //double PenaltyStiffnessScalingFactorPenaltyStiffnessTangentialCrackOpening(0.1);
    //myStructureCoarseScale.ElementTotalSetFineScaleParameter("ConstraintPenaltyStiffnessTangentialCrackOpening", PenaltyStiffnessTangentialCrackOpening);
    //myStructureCoarseScale.ElementTotalSetFineScaleParameter("PenaltyStiffnessScalingFactorTangentialCrackOpening", PenaltyStiffnessScalingFactorPenaltyStiffnessTangentialCrackOpening);

    //augmented lagrangian for normal crack opening to be non negativ
    //double AugmentedLagrangeStiffnessCrackOpening(10000);
    //myStructureCoarseScale.ElementTotalSetFineScaleParameter("AugmentedLagrangeCrackOpening", AugmentedLagrangeStiffnessCrackOpening);

    //difference in both principal strains in order to solve the crack direction elastic from the total strain
    myStructureCoarseScale.ElementTotalSetFineScaleParameter("ToleranceElasticCrackAngleHigh",0.0000001);
    //difference in both principal strains in order to solve the crack direction elastic from the previous state
    myStructureCoarseScale.ElementTotalSetFineScaleParameter("ToleranceElasticCrackAngleLow", 0.00000005);

    // apply constraints set the strains and
    // do the initial step using the linear elastic solution (no fine scale model just the calibration of the initial angles)
    NuTo::FullMatrix<double> strain(3,1);
    // pure tension in diagonal direction for nu=0
    //strain(0,0) = 0.001;
    //strain(1,0) = 0.001;
    //strain(2,0) = 0.002;
    strain(0,0) = 0.001;
    strain(1,0) = 0.00;
    strain(2,0) = 0.00;
    //strain(0,0) = 0.00;
    //strain(1,0) = 0.000001;
    //strain(2,0) = 0.002;

    myStructureCoarseScale.SetTotalStrain(strain);
    myStructureCoarseScale.InitStructure();
    int numLoadStepMacro=50;
    for (int loadstepMacro=1;loadstepMacro <= numLoadStepMacro; loadstepMacro++)
    {
        NuTo::FullMatrix<double> curStrain(strain*((double)loadstepMacro/(double)numLoadStepMacro));
        std::cout << "new load step on the macroscale " << loadstepMacro << std::endl;
        curStrain.Trans().Info(12,4);
        myStructureCoarseScale.SetTotalStrain(curStrain);
        myStructureCoarseScale.NewtonRaphson();
       myStructureCoarseScale.ElementTotalUpdateStaticData();
       double energy = myStructureCoarseScale.ElementTotalGetTotalEnergy();
        std::cout << "total energy on the macroscale " << energy << std::endl;
        //system("paraview --state=/home/unger3/develop/nuto_build/examples/c++/test.pvsm");
    }

    myStructureCoarseScale.ExportVtkDataFile(std::string("/home/unger3/develop/nuto_build/examples/c++/MacroscaleConcurrentMultiscaleFinal.vtk"));

}
catch (NuTo::Exception& e)
{
    std::cout << "Error executing ConcurrentMultiscale "<< std::endl;
    std::cout << e.ErrorMessage() << std::endl;
    exit(1);
}
}

