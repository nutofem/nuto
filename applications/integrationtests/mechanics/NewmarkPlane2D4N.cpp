#include "BoostUnitTest.h"
#include <sstream>
#include <iomanip>
#include <boost/filesystem.hpp>

#include "math/EigenCompanion.h"
#include "mechanics/MechanicsEnums.h"
#include "mechanics/groups/Group.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/mesh/MeshGenerator.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "visualize/VisualizeEnum.h"

void CheckResult(std::string resultDir, std::string file, Eigen::MatrixXd expected, double tolerance = 1.e-4)
{
    boost::filesystem::path resultFile = resultDir;
    resultFile /= std::string(file);
	Eigen::MatrixXd result = NuTo::EigenCompanion::ReadFromFile(resultFile.string());
    BoostUnitTest::CheckEigenMatrix(result, expected, tolerance);
}

BOOST_AUTO_TEST_CASE(NewmarkPlane2D4N)
{
    std::string resultDir = boost::filesystem::initial_path().string()+std::string("/ResultsNewmarkPlane2D4N");

    //delete result directory
    if (boost::filesystem::exists(resultDir))    // does p actually exist?
    {
        if (boost::filesystem::is_directory(resultDir))      // is p a directory?
        {
        	boost::filesystem::remove_all(resultDir);
        }
    }
    // create result directory
    boost::filesystem::create_directory(resultDir);

    double mL = 100;
    double mH = 10;
    double density = 1;       //N and mm^3

    NuTo::Structure myStructure(2);
    myStructure.SetNumTimeDerivatives(2);

	//create nodes
    auto meshInfo = NuTo::MeshGenerator::Grid(myStructure, {mL, mH}, {10, 10});
    myStructure.InterpolationTypeAdd(meshInfo.second, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.ElementTotalConvertToInterpolationType();

    //section
	double thickness(1);
    auto mySection = NuTo::SectionPlane::Create(thickness, false);
    myStructure.ElementTotalSetSection(mySection);

	//constitutive
    int myMatLattice = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    myStructure.ConstitutiveLawSetParameterDouble(myMatLattice,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,30000);
    myStructure.ConstitutiveLawSetParameterDouble(myMatLattice,NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,0.0);
    myStructure.ConstitutiveLawSetParameterDouble(myMatLattice,NuTo::Constitutive::eConstitutiveParameter::DENSITY,density);
    myStructure.ElementTotalSetConstitutiveLaw(myMatLattice);

	auto& grpNodes_Left = myStructure.GroupGetNodesAtCoordinate(NuTo::eDirection::X, 0);
	auto& grpNodes_Right = myStructure.GroupGetNodesAtCoordinate(NuTo::eDirection::X, mL);

	myStructure.CalculateMaximumIndependentSets();
    myStructure.Info();

    auto& nodeLeftBottom = myStructure.NodeGetAtCoordinate(Eigen::Vector2d::Zero());
	myStructure.LoadCreateNodeForce(0, myStructure.NodeGetId(&nodeLeftBottom),Eigen::Vector2d::UnitX(), 1);

	NuTo::NewmarkDirect myIntegrationScheme(&myStructure);
	myIntegrationScheme.SetDampingCoefficientMass(0.05);

	double simulationTime(10.);
    double period = 5;
    auto sinusoidalQuarterWave = [=](double time)
    {
        if (time < 1.25)
            return 0.1 * sin(time/period*2.*M_PI);
        return 0.1;
    };

    myStructure.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
            NuTo::Constraint::Component(grpNodes_Left, {NuTo::eDirection::X})); 
    
    myStructure.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
            NuTo::Constraint::Component(grpNodes_Right, {NuTo::eDirection::X}, sinusoidalQuarterWave)); 

    myStructure.SetShowTime(false);
    myStructure.SetNumProcessors(1);
   
    //add visualization
    int visualizationGroup = myStructure.GroupGetElementsTotal();
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    
    myIntegrationScheme.SetTimeStep(10);
    myIntegrationScheme.SetMaxTimeStep(10);
    myIntegrationScheme.SetMinTimeStep(0.001*myIntegrationScheme.GetMaxTimeStep());
	myIntegrationScheme.AddResultTime("Time");
	myIntegrationScheme.AddResultGroupNodeForce("Forces_GroupNodes_Left",myStructure.GroupGetId(&grpNodes_Left));
	myIntegrationScheme.AddResultGroupNodeForce("Forces_GroupNodes_Right",myStructure.GroupGetId(&grpNodes_Right));
    myIntegrationScheme.SetResultDirectory(resultDir, true);
    myIntegrationScheme.SetToleranceForce(1.e-10);
    myIntegrationScheme.Solve(simulationTime);

    Eigen::Matrix2d forceLeftRef;
    forceLeftRef.setZero();
    forceLeftRef(0,0) = -1; //disp on fixed node
    forceLeftRef(1,0) = -3.001682840791e+02;

    Eigen::Matrix2d forceRightRef;
    forceRightRef.setZero();
    forceRightRef(1,0) = 3.016648179801e+02;

    CheckResult(resultDir, "Forces_GroupNodes_Left.dat", forceLeftRef, 1.e-4);
    CheckResult(resultDir, "Forces_GroupNodes_Right.dat", forceRightRef, 1.e-4);
    CheckResult(resultDir, "Time.dat", Eigen::Vector2d(0, 10), 1.e-4);
}


