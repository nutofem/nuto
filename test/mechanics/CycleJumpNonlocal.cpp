#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/timeIntegration/JumpDirect.h"


#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>


constexpr unsigned int dimension = 2;
class Parameters
{
public:

	static constexpr bool EnableDisplacementControl = false;

    static const int mDimension = dimension;
    static constexpr double mRectangleSize = 250;

    static const bool mPerformLineSearch = true;		// original from Philip
    static const bool mAutomaticTimeStepping = true;

    static constexpr double mMatrixYoungsModulus = 4.0e4;   // concrete
    static constexpr double mMatrixPoissonsRatio = 0.2;
    static constexpr double mMatrixThickness = 10;
//    static constexpr double mMatrixNonlocalRadius = 2;     // original philip
    static constexpr double mMatrixNonlocalRadius = 18;
    static constexpr double mMatrixNonlocalRadiusParameter = 0.;
    static constexpr double mMatrixTensileStrength = 3;
//    static constexpr double mMatrixCompressiveStrength = 30;	// original
    static constexpr double mMatrixCompressiveStrength = 20;
    static constexpr double mMatrixFractureEnergy = 1;


//    static constexpr double mTimeStep = 1e-2;		// Philip
    static constexpr double mTimeStep = 1e-3;
//    static constexpr double mMinTimeStep = 1e-5;		// Philip
    static constexpr double mMinTimeStep = 1e-7;
    static constexpr double mMaxTimeStep = 1e-1;
//    static constexpr double mToleranceForce = 1e-8;
    static constexpr double mToleranceForce = 1e-7;				// for fine mesh
    static constexpr double mSimulationTime = 3.0;
//    static constexpr double mLoad = 10*0.02;
    static constexpr double mDisplacement = 10*0.05;
    static constexpr double mLoad = 3*1.5e4;			// original tests with 4*1.5e4


    static const boost::filesystem::path mOutputPath;
    static const boost::filesystem::path mMeshFilePathMatrix;
    static const boost::filesystem::path mMeshFilePathFibre;

    static const NuTo::FullVector<double, dimension> mDirectionX;
    static const NuTo::FullVector<double, dimension> mDirectionY;
};


const boost::filesystem::path Parameters::mOutputPath("./ResultsCycleJumpNonlocal");

const NuTo::FullVector<double, dimension> Parameters::mDirectionX = NuTo::FullVector<double, dimension>::UnitX();
const NuTo::FullVector<double, dimension> Parameters::mDirectionY = NuTo::FullVector<double, dimension>::UnitY();


int main(int argc, char* argv[])
{

    //**********************************************
    //          Structure
    //**********************************************

    NuTo::Structure myStructure(Parameters::mDimension);

    //**********************************************
    //         Integration Scheme
    //**********************************************

    NuTo::JumpDirect myIntegrationScheme(&myStructure);
    myIntegrationScheme.SetTimeStep(Parameters::mTimeStep);
    myIntegrationScheme.SetMinTimeStep(Parameters::mMinTimeStep);
    myIntegrationScheme.SetMaxTimeStep(Parameters::mMaxTimeStep);
    myIntegrationScheme.SetToleranceForce(Parameters::mToleranceForce);
    myIntegrationScheme.SetAutomaticTimeStepping(Parameters::mAutomaticTimeStepping);
    myIntegrationScheme.SetPerformLineSearch(Parameters::mPerformLineSearch);
    myIntegrationScheme.SetResultDirectory(Parameters::mOutputPath.string(), true);


    //**********************************************
    //          Section
    //**********************************************

    int mySection = myStructure.SectionCreate(NuTo::Section::PLANE_STRAIN);
    myStructure.SectionSetThickness(mySection, Parameters::mMatrixThickness);

    //**********************************************
    //          Material
    //**********************************************

    NuTo::FullVector<double, Eigen::Dynamic> myDamageLaw(1);
    myDamageLaw(0) = NuTo::Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING;


//    int matrixMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    int matrixMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS_FATIGUE);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mMatrixYoungsModulus);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mMatrixPoissonsRatio);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH, Parameters::mMatrixTensileStrength);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH, Parameters::mMatrixCompressiveStrength);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS, Parameters::mMatrixNonlocalRadius);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS_PARAMETER, Parameters::mMatrixNonlocalRadiusParameter);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY, Parameters::mMatrixFractureEnergy);
//    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::VISCOSITY_EXPONENT, 1.5);
    myStructure.ConstitutiveLawSetDamageLaw(matrixMaterial, myDamageLaw);


    //**********************************************
    //          Geometry
    //**********************************************

    // create nodes
    double size(250.);
    NuTo::FullVector<double,Eigen::Dynamic> nodeCoordinates(2);

    nodeCoordinates << 0., 0;
    myStructure.NodeCreate(0, nodeCoordinates);

    nodeCoordinates << size, 0.;
    myStructure.NodeCreate(1, nodeCoordinates);

    nodeCoordinates << size, size;
    myStructure.NodeCreate(2, nodeCoordinates);

    nodeCoordinates << 0., size;
    myStructure.NodeCreate(3, nodeCoordinates);

    // create interpolation type

    int myInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRIANGLE2D);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::NONLOCALEQSTRAIN, NuTo::Interpolation::EQUIDISTANT1);

    myStructure.InterpolationTypeSetIntegrationType(myInterpolationType, NuTo::IntegrationType::IntegrationType2D3NGauss3Ip, NuTo::IpData::STATICDATA);

    // create elements
    NuTo::FullVector<int,Eigen::Dynamic> elementIncidence(3);

    elementIncidence << 0, 1, 2;		// provide node numbers for the 1st element
    myStructure.ElementCreate(myInterpolationType, elementIncidence, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,NuTo::IpData::eIpDataType::STATICDATA);

    elementIncidence << 2, 3, 0;		// provide node numbers for the 1st element
    myStructure.ElementCreate(myInterpolationType, elementIncidence, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,NuTo::IpData::eIpDataType::STATICDATA);


    myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);
    myStructure.ElementTotalSetSection(mySection);
    myStructure.ElementTotalSetConstitutiveLaw(matrixMaterial);

    //**********************************************
    //          Boundary Conditions
    //**********************************************
	//bottom boundary
	int GrpNodes_BottomBoundary = myStructure.GroupCreate("Nodes");
	int directionGrp = 1; //either 0,1,2
	double min(0);
	double max(0);
	myStructure.GroupAddNodeCoordinateRange(GrpNodes_BottomBoundary,directionGrp,min,max);
	myStructure.ConstraintLinearSetDisplacementNodeGroup(GrpNodes_BottomBoundary, Parameters::mDirectionY, 0);

	//bottom left node
	int GrpNodes_BottomLeftNode = myStructure.GroupCreate("Nodes");
    NuTo::FullVector<double, 2> center;
    center[0] = 0;
    center[1] = 0;
	min = 0;
	max = 1e-6;
	myStructure.GroupAddNodeRadiusRange(GrpNodes_BottomLeftNode,center,min,max);
	myStructure.ConstraintLinearSetDisplacementNodeGroup(GrpNodes_BottomLeftNode, Parameters::mDirectionX, 0);

	//top boundary (loaded)
	int GrpNodes_TopBoundary = myStructure.GroupCreate("Nodes");
	min = Parameters::mRectangleSize;
	max = Parameters::mRectangleSize;
	myStructure.GroupAddNodeCoordinateRange(GrpNodes_TopBoundary,directionGrp,min,max);

    //**********************************************
    //          Loads
    //**********************************************
	const double LoadFraction(0.017);
	if (Parameters::EnableDisplacementControl) {
		// middle top boundary
		int PrescribedDisplacement = myStructure.ConstraintLinearSetDisplacementNodeGroup(GrpNodes_TopBoundary, Parameters::mDirectionY, 0);

	    NuTo::FullMatrix<double, 2, 2> dispRHS;
	    dispRHS(0, 0) = 0;
	    dispRHS(1, 0) = LoadFraction * Parameters::mSimulationTime;
	    dispRHS(0, 1) = 0;
	    dispRHS(1, 1) = LoadFraction * Parameters::mDisplacement;

	   	myIntegrationScheme.SetTimeDependentConstraint(PrescribedDisplacement, dispRHS);

    	myIntegrationScheme.AddResultGroupNodeForce("myforce", GrpNodes_TopBoundary);

	}	else	{
		// new block for load control
		//print set nodes
		std::cout << "#GrpNodes_TopBoundary = " << std::endl;
		NuTo::FullVector<int, Eigen::Dynamic> TopBoundaryNodesVector(myStructure.GroupGetMemberIds(GrpNodes_TopBoundary));
		// via constraint
		// loop through all the nodes from the top boundary except of the first node, which is loaded and not constrained
		// the constraint equations have the form: Uy(ConstrainedNodeNumber) - Uy(LoadedNodeNumber) = 0, where ConstrainedNodeNumber > 0.
		int LoadedNodeNumber(TopBoundaryNodesVector[0]);
		for (int var = 1; var < TopBoundaryNodesVector.size(); ++var) {
			int ConstrainedNodeNumber(TopBoundaryNodesVector[var]),
					ConstraintEquationNumber;
			std::cout << "Create a constraint equation for the node #" << ConstrainedNodeNumber << std::endl;
			// create constraint equation
			ConstraintEquationNumber = myStructure.ConstraintLinearEquationCreate(ConstrainedNodeNumber,"Y_DISPLACEMENT",1.,0.);
			// add the term: -1. * Uy(LoadedNodeNumber)
			myStructure.ConstraintLinearEquationAddTerm(ConstraintEquationNumber,LoadedNodeNumber,"Y_DISPLACEMENT",-1.);
		}

		// apply Load to the node #LoadedNodeNumber
		myStructure.LoadCreateNodeForce(0,LoadedNodeNumber,Parameters::mDirectionY,Parameters::mLoad);

		NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> loadingBC(2,2);
		loadingBC(0,0) = 0.;   				// t=0
        loadingBC(0,1) = 0.;   				// u(0)=0
        loadingBC(1,0) = LoadFraction;   	// t=LoadFraction
        loadingBC(1,1) = LoadFraction;   	// u(1)= LoadFraction
		myIntegrationScheme.SetTimeDependentLoadCase(0, loadingBC);
	}


//**********************************************
//          Visualisation
//**********************************************
#ifdef ENABLE_VISUALIZE
        // visualize results
        int visualizationGroup = myStructure.GroupCreate(NuTo::Groups::eGroupId::Elements);
        myStructure.GroupAddElementsTotal(visualizationGroup);

        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::DISPLACEMENTS);
        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRAIN);
        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRESS);
        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::DAMAGE);
        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::NONLOCAL_EQ_STRAIN);

        myStructure.ExportVtkDataFileElements("CycleJumpNonlocal.vtk");

#endif // ENABLE_VISUALIZE
//**********************************************
//          Solver
//**********************************************

    const int numElements(myStructure.GetNumElements());


    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();

    center[0] = 0;
    center[1] = Parameters::mRectangleSize;
    int grpNodes_output = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeRadiusRange(grpNodes_output, center, 0, 5e-1);

//    if (Parameters::EnableDisplacementControl) {
//    	myIntegrationScheme.AddResultGroupNodeForce("myforce", GrpNodes_TopBoundary);
//    }
    myIntegrationScheme.AddResultNodeDisplacements("mydisplacements3", 3);
    myIntegrationScheme.AddResultElementIpData("Out_ElementStress1",1,NuTo::IpData::eIpStaticDataType::ENGINEERING_STRESS);
    myIntegrationScheme.AddResultTime("time");

    //**********************************************
    //			Periodic Loading
    //**********************************************

    // set incrementation
	myIntegrationScheme.SetHarmonicIncrementation(20);


    // apply harmonic excitation
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> HarmonicConstraintFactor(1,3);
    const double AmplitudeFraction(0.1);
    if (Parameters::EnableDisplacementControl) {
    	HarmonicConstraintFactor(0,0) = AmplitudeFraction * LoadFraction * Parameters::mDisplacement;
    } else {
    	HarmonicConstraintFactor(0,0) = AmplitudeFraction * LoadFraction;
    }
    HarmonicConstraintFactor(0,1) = 100;					//	set frequency
    HarmonicConstraintFactor(0,2) = 2.9;


    myIntegrationScheme.SetHarmonicExcitation(HarmonicConstraintFactor);

    double ExtrapolationTolerance(0.04);					// bulo 0.004
	myIntegrationScheme.SetHarmonicExtrapolation(true,&ExtrapolationTolerance);
//	myIntegrationScheme.SetHarmonicExtrapolation(false);


    //**********************************************


    myIntegrationScheme.Solve(Parameters::mSimulationTime);


    std::cout << " ===> End CycleJumpNonlocal Test <=== " << std::endl;

    return 0;
}

