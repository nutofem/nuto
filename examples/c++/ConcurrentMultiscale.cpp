#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
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
        mlX = 100;
        mlY = 100;

        //create structure
        SetShowTime(true);

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
        ConstitutiveLawSetYoungsModulus(myMatCoarseScale,1);
        ConstitutiveLawSetPoissonsRatio(myMatCoarseScale,0.2);

        //create section
        double thickness(1);
        int mySectionCoarseScale = SectionCreate("Plane_Strain");
        SectionSetThickness(mySectionCoarseScale,thickness);

        ElementTotalSetSection(mySectionCoarseScale);
        ElementTotalSetConstitutiveLaw(myMatCoarseScale);

        //set fine scale model as ip data for all integration points
        //myStructureCoarseScale.ElementTotalSetFineScaleModel("myStructureFineScale.bin");
        ElementTotalSetFineScaleModel("/home/unger3/develop/nuto_build/examples/c++/myStructureFineScale.bin");

        resultMatrix.Resize(0,5);
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
        this->ConstraintSetRHS(mConstraint4Y,rLoadFactor*0.5*mStrain(2,0)*mlX);
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


    //crack transition zone
    myStructureCoarseScale.ElementTotalSetFineScaleParameter("CrackTransitionZone",10);

    //penalty stiffness for crack angle
    double PenaltyStiffnessCrackAngle(1);
    double PenaltyStiffnessScalingFactorCrackAngle(2.*M_PI);
    myStructureCoarseScale.ElementTotalSetFineScaleParameter("ConstraintPenaltyStiffnessCrackAngle", PenaltyStiffnessCrackAngle);
    myStructureCoarseScale.ElementTotalSetFineScaleParameter("PenaltyStiffnessScalingFactorCrackAngle", PenaltyStiffnessScalingFactorCrackAngle);

    //penalty stiffness for tangential crack opening
    double PenaltyStiffnessTangentialCrackOpening(.1);
    double PenaltyStiffnessScalingFactorPenaltyStiffnessTangentialCrackOpening(0.1);
    myStructureCoarseScale.ElementTotalSetFineScaleParameter("ConstraintPenaltyStiffnessTangentialCrackOpening", PenaltyStiffnessTangentialCrackOpening);
    myStructureCoarseScale.ElementTotalSetFineScaleParameter("PenaltyStiffnessScalingFactorTangentialCrackOpening", PenaltyStiffnessScalingFactorPenaltyStiffnessTangentialCrackOpening);

    //augmented lagrangian for normal crack opening to be non negativ
    double AugmentedLagrangeStiffnessCrackOpening(1);
    //myStructureCoarseScale.ElementTotalSetFineScaleParameter("AugmentedLagrangeCrackOpening", AugmentedLagrangeStiffnessCrackOpening);

    //difference in both principal strains in order to solve the crack direction elastic from the total strain
    myStructureCoarseScale.ElementTotalSetFineScaleParameter("ToleranceElasticCrackAngleHigh",0.1);
    //difference in both principal strains in order to solve the crack direction elastic from the previous state
    myStructureCoarseScale.ElementTotalSetFineScaleParameter("ToleranceElasticCrackAngleLow",0.05);

    // apply constraints set the strains and
    // do the initial step using the linear elastic solution (no fine scale model just the calibration of the initial angles)
    NuTo::FullMatrix<double> strain(3,1);
    strain(0,0) = 0.1;
    strain(1,0) = 0.;
    strain(2,0) = 0.;
    myStructureCoarseScale.SetTotalStrain(strain);
    myStructureCoarseScale.InitStructure();

    strain(0,0) = 0.1;
    strain(1,0) = 0.;
    strain(2,0) = 0.0;
    myStructureCoarseScale.SetTotalStrain(strain);
    myStructureCoarseScale.NewtonRaphson();

    myStructureCoarseScale.ExportVtkDataFile(std::string("/home/unger3/develop/nuto_build/examples/c++/MacroscaleConcurrentMultiscaleFinal.vtk"));

}
catch (NuTo::Exception& e)
{
    std::cout << "Error executing ConcurrentMultiscale "<< std::endl;
    std::cout << e.ErrorMessage() << std::endl;
    exit(1);
}
}

