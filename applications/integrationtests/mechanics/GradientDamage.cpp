#include "BoostUnitTest.h"

#include <boost/filesystem.hpp>

#include "base/Timer.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/sections/SectionTruss.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/damageLaws/DamageLawExponential.h"
#include "mechanics/constitutive/laws/GradientDamageEngineeringStress.h"
#include "mechanics/elements/IpDataEnum.h"
#include "mechanics/groups/GroupBase.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/sections/SectionVariableTruss.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/timeIntegration/ImplEx.h"

#include "mechanics/mesh/MeshGenerator.h"

#include "mechanics/elements/ContinuumBoundaryElement.h"
#include "math/SparseMatrixCSRVector2General.h"
#include "visualize/VisualizeEnum.h"


int SetConstitutiveLaw(NuTo::Structure& rStructure)
{
    // create a damage law
    int lawId = rStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::DENSITY, 1.0);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 30000);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.2);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS, 3);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH, 4.);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH, 4. * 10);
    rStructure.ConstitutiveLawSetDamageLaw(lawId, NuTo::Constitutive::DamageLawExponential::Create(4./30000., 4. / 0.021));
    return lawId;
}

void ApplyDofValues(NuTo::Structure& rStructure)
{
    double boundaryDisplacement = 1.e-1;

    rStructure.CalculateMaximumIndependentSets();
    rStructure.NodeBuildGlobalDofs();

    int gAllNodes = rStructure.GroupGetNodesTotal();
    for (int nodeId : rStructure.GroupGetMemberIds(gAllNodes))
    {
        NuTo::NodeBase* node = rStructure.NodeGetNodePtr(nodeId);
        Eigen::VectorXd disps = node->Get(NuTo::Node::eDof::COORDINATES) / 100. * boundaryDisplacement;
        node->Set(NuTo::Node::eDof::DISPLACEMENTS, disps);

        if (node->GetNum(NuTo::Node::eDof::NONLOCALEQSTRAIN) > 0)
            node->Set(NuTo::Node::eDof::NONLOCALEQSTRAIN, disps(0, 0) / 10);
    }
}

void CheckStiffnesses(NuTo::Structure& rStructure)
{
    ApplyDofValues(rStructure);
    BOOST_CHECK(rStructure.ElementCheckHessian0(1.e-8, 1.e-4, true));
    BOOST_CHECK(rStructure.CheckHessian0(1.e-8, 1.e-4, true));
}

void Visualize(NuTo::Structure& rStructure, std::string rDir)
{
    int visualizationGroup = rStructure.GroupGetElementsTotal();

    rStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    rStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::NONLOCAL_EQ_STRAIN);
    rStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::LOCAL_EQ_STRAIN);
    rStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DAMAGE);


    std::string resultDir = "./ResultsGradientDamage";
    boost::filesystem::create_directory(resultDir);
    resultDir += "/" + rDir;
    boost::filesystem::create_directory(resultDir);

    rStructure.ExportVtkDataFileElements(resultDir+"/Elements.vtu", true);
}

void AddInterpolationType(NuTo::Structure& rS, int rITid)
{
    rS.InterpolationTypeAdd(rITid, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    rS.InterpolationTypeAdd(rITid, NuTo::Node::eDof::NONLOCALEQSTRAIN, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
}

template <int TDim>
int AddBoundaryElementsAtCoordinate(NuTo::Structure& rS, double rCoordinate)
{
    int gNodesBoundary = rS.GroupCreate(NuTo::eGroupId::Nodes);
    rS.GroupAddNodeCoordinateRange(gNodesBoundary, 0, rCoordinate - 1.e-6, rCoordinate + 1.e-6);

    int elemGroupBoundary = rS.GroupCreate(NuTo::eGroupId::Elements);
    rS.GroupAddElementsFromNodes(elemGroupBoundary, gNodesBoundary, false);

    int gBoundaryElements = rS.BoundaryElementsCreate(elemGroupBoundary, gNodesBoundary);
    for (int boundaryElementId : rS.GroupGetMemberIds(gBoundaryElements))
    {
        auto& boundaryElement = dynamic_cast<NuTo::ContinuumBoundaryElement<TDim>&>(*rS.ElementGetElementPtr(boundaryElementId));
        boundaryElement.SetAlpha(42.);
    }

    return rS.GroupGetNumMembers(gBoundaryElements);
}

template <int TDim>
void AddBoundaryElements(NuTo::Structure& s, double rLength, int rNumExpectedBoundaryElements)
{
    int numBoundaryElements = 0;
    numBoundaryElements += AddBoundaryElementsAtCoordinate<TDim>(s, 0);
    numBoundaryElements += AddBoundaryElementsAtCoordinate<TDim>(s, rLength);
    BOOST_CHECK_EQUAL(numBoundaryElements, rNumExpectedBoundaryElements);
}

void TestStructure1D(bool rUseRobinBoundaryElements)
{
    NuTo::Timer timer(__FUNCTION__);
    const int numElements = 4;
    const double length = 42;
    const double area = 10;

    NuTo::Structure s(1);
    s.SetShowTime(false);
    s.SetVerboseLevel(0);

    int interpolationType = NuTo::MeshGenerator::Grid(s, {length}, {numElements}).second;
    AddInterpolationType(s, interpolationType);

    // create sections
    const double xWeakSpot = 20;
    const double lWeakSpot = 10;
    const double alpha = 0.20;

    auto section = NuTo::SectionVariableTruss::Create(area, xWeakSpot, lWeakSpot, alpha);

    s.ElementTotalConvertToInterpolationType();
    s.ElementTotalSetConstitutiveLaw(SetConstitutiveLaw(s));
    s.ElementTotalSetSection(section);

    if (rUseRobinBoundaryElements)
    {
        AddBoundaryElements<1>(s, length, 2);
    }
    CheckStiffnesses(s);
    Visualize(s, "TRUSS1D");
}



void TestStructure2D(NuTo::Interpolation::eShapeType rShape, bool isPlaneStrain, bool rUseRobinBoundaryElements)
{
    NuTo::Timer timer(std::string(__FUNCTION__) + " " + NuTo::Interpolation::ShapeTypeToString(rShape));

    // define geometry
    double lX = 20.;
    double lY = 5.;
    double lZ = 2.;

    int numElementsX = 1;
    int numElementsY = 2;

    NuTo::Structure s(2);
    s.SetShowTime(false);
    s.SetVerboseLevel(0);

    int interpolationType = NuTo::MeshGenerator::Grid(s, {lX, lY}, {numElementsX, numElementsY}, rShape).second;
    AddInterpolationType(s, interpolationType);

    auto section = NuTo::SectionPlane::Create(lZ, isPlaneStrain);

    s.ElementTotalConvertToInterpolationType();
    s.ElementTotalSetConstitutiveLaw(SetConstitutiveLaw(s));
    s.ElementTotalSetSection(section);

    if (rUseRobinBoundaryElements)
    {
        AddBoundaryElements<2>(s, lX, 2*numElementsY);
    }
    CheckStiffnesses(s);
    Visualize(s, NuTo::Interpolation::ShapeTypeToString(rShape));
}


void TestStructure3D(NuTo::Interpolation::eShapeType rShape, bool rUseRobinBoundaryElements)
{
    NuTo::Timer timer(std::string(__FUNCTION__) + " " + NuTo::Interpolation::ShapeTypeToString(rShape));

    NuTo::Structure s(3);
    s.SetShowTime(false);
    s.SetVerboseLevel(0);

    double lX = 3, lY = 4, lZ = 5;

    int interpolationType = NuTo::MeshGenerator::Grid(s, {lX, lY, lZ}, {1,1,1}, rShape).second;
    AddInterpolationType(s, interpolationType);

    s.ElementTotalConvertToInterpolationType();
    s.ElementTotalSetConstitutiveLaw(SetConstitutiveLaw(s));

    if (rUseRobinBoundaryElements)
    {
        int numExpectedBoundaryElements = 2; // for brick
        if (rShape == NuTo::Interpolation::eShapeType::TETRAHEDRON3D)
            numExpectedBoundaryElements = 4;

        AddBoundaryElements<3>(s, lX, numExpectedBoundaryElements);
    }

    CheckStiffnesses(s);
    Visualize(s, NuTo::Interpolation::ShapeTypeToString(rShape));
}


void GroupRemoveNodesWithoutDisplacements(NuTo::Structure& rStructure, int rGroupNodeId)
{
    for (int nodeId : rStructure.GroupGetMemberIds(rGroupNodeId))
    {
        NuTo::NodeBase* node = rStructure.NodeGetNodePtr(nodeId);
        if (node->GetNum(NuTo::Node::eDof::DISPLACEMENTS) == 0)
        {
            NuTo::GroupBase* group = rStructure.GroupGetGroupPtr(rGroupNodeId);
            group->RemoveMember(nodeId);
        }
    }
}

void SetupNewmark(NuTo::NewmarkDirect& rTimeIntegration, int rBC, std::string rDir)
{
    double simulationTime = 1;
    double dispEnd = 0.01;
    int numLoadSteps = 3;

    Eigen::Matrix2d timeDepDisp;
    timeDepDisp << 0, 0, simulationTime, dispEnd;

    rTimeIntegration.AddTimeDependentConstraint(rBC, timeDepDisp);
    rTimeIntegration.SetTimeStep(simulationTime / numLoadSteps);
    rTimeIntegration.SetAutomaticTimeStepping(false);
    rTimeIntegration.SetToleranceForce(1e-8);
    rTimeIntegration.SetMaxNumIterations(12);
    rTimeIntegration.SetPerformLineSearch(true);

    std::string resultDir = "./ResultsGradientDamage";
    boost::filesystem::create_directory(resultDir);
    resultDir += "/" + rDir;
    boost::filesystem::create_directory(resultDir);
    rTimeIntegration.SetResultDirectory(resultDir, true);
}

void Check1D2D3D()
{
    NuTo::Timer timer(std::string(__FUNCTION__) + " Setup");
    double lx = 50;
    double ly = .5;
    double lz = .5;

    int numElements = 3;

    NuTo::Structure s1D(1);
    NuTo::Structure s2D(2);
    NuTo::Structure s3D(3);

    s1D.SetShowTime(false);
    s2D.SetShowTime(false);
    s3D.SetShowTime(false);

    int interpolationType1D = NuTo::MeshGenerator::Grid(s1D, {lx}, {numElements}).second;
    int interpolationType2D = NuTo::MeshGenerator::Grid(s2D, {lx, ly}, {numElements, 1}).second;
    int interpolationType3D = NuTo::MeshGenerator::Grid(s3D, {lx, ly, lz}, {numElements, 1, 1}).second;

    AddInterpolationType(s1D, interpolationType1D);
    AddInterpolationType(s2D, interpolationType2D);
    AddInterpolationType(s3D, interpolationType3D);

    s1D.InterpolationTypeSetIntegrationType(interpolationType1D, NuTo::eIntegrationType::IntegrationType1D2NGauss2Ip);
    s2D.InterpolationTypeSetIntegrationType(interpolationType2D, NuTo::eIntegrationType::IntegrationType2D4NGauss4Ip);
    s3D.InterpolationTypeSetIntegrationType(interpolationType3D, NuTo::eIntegrationType::IntegrationType3D8NGauss2x2x2Ip);

    s1D.ElementTotalConvertToInterpolationType();
    s2D.ElementTotalConvertToInterpolationType();
    s3D.ElementTotalConvertToInterpolationType();

    s1D.ElementTotalSetConstitutiveLaw(SetConstitutiveLaw(s1D));
    s2D.ElementTotalSetConstitutiveLaw(SetConstitutiveLaw(s2D));
    s3D.ElementTotalSetConstitutiveLaw(SetConstitutiveLaw(s3D));

    auto mySection1D = NuTo::SectionTruss::Create(lz*ly);
    auto mySection2D = NuTo::SectionPlane::Create(lz, false);

    s1D.ElementTotalSetSection(mySection1D);
    s2D.ElementTotalSetSection(mySection2D);

    int weakElementId = numElements / 2;
    double kappa = 0.001;
    {
        NuTo::ElementBase* element = s1D.ElementGetElementPtr(weakElementId);
        for (int i = 0; i < element->GetNumIntegrationPoints(); ++i)
            element->GetIPData().GetIPConstitutiveLaw(i).GetData<NuTo::GradientDamageEngineeringStress>().SetData(kappa);
    }
    {
        NuTo::ElementBase* element = s2D.ElementGetElementPtr(weakElementId);
        for (int i = 0; i < element->GetNumIntegrationPoints(); ++i)
            element->GetIPData().GetIPConstitutiveLaw(i).GetData<NuTo::GradientDamageEngineeringStress>().SetData(kappa);
    }
    {
        NuTo::ElementBase* element = s3D.ElementGetElementPtr(weakElementId);
        for (int i = 0; i < element->GetNumIntegrationPoints(); ++i)
            element->GetIPData().GetIPConstitutiveLaw(i).GetData<NuTo::GradientDamageEngineeringStress>().SetData(kappa);
    }


    int leftNodes1D = s1D.GroupCreate(NuTo::eGroupId::Nodes);
    int leftNodes2D = s2D.GroupCreate(NuTo::eGroupId::Nodes);
    int leftNodes3D = s3D.GroupCreate(NuTo::eGroupId::Nodes);

    int rightNodes1D = s1D.GroupCreate(NuTo::eGroupId::Nodes);
    int rightNodes2D = s2D.GroupCreate(NuTo::eGroupId::Nodes);
    int rightNodes3D = s3D.GroupCreate(NuTo::eGroupId::Nodes);

    s1D.GroupAddNodeCoordinateRange(leftNodes1D, 0, -1.e-4, 1.e-4);
    s2D.GroupAddNodeCoordinateRange(leftNodes2D, 0, -1.e-4, 1.e-4);
    s3D.GroupAddNodeCoordinateRange(leftNodes3D, 0, -1.e-4, 1.e-4);

    s1D.GroupAddNodeCoordinateRange(rightNodes1D, 0, lx - 1.e-4, lx + 1.e-4);
    s2D.GroupAddNodeCoordinateRange(rightNodes2D, 0, lx - 1.e-4, lx + 1.e-4);
    s3D.GroupAddNodeCoordinateRange(rightNodes3D, 0, lx - 1.e-4, lx + 1.e-4);

    GroupRemoveNodesWithoutDisplacements(s1D, leftNodes1D);
    GroupRemoveNodesWithoutDisplacements(s1D, rightNodes1D);
    GroupRemoveNodesWithoutDisplacements(s2D, leftNodes2D);
    GroupRemoveNodesWithoutDisplacements(s2D, rightNodes2D);
    GroupRemoveNodesWithoutDisplacements(s3D, leftNodes3D);
    GroupRemoveNodesWithoutDisplacements(s3D, rightNodes3D);

    s1D.ConstraintLinearSetDisplacementNodeGroup(leftNodes1D, Eigen::Matrix<double, 1, 1>::UnitX(), 0.0);
    s2D.ConstraintLinearSetDisplacementNodeGroup(leftNodes2D, Eigen::Matrix<double, 2, 1>::UnitX(), 0.0);
    s3D.ConstraintLinearSetDisplacementNodeGroup(leftNodes3D, Eigen::Matrix<double, 3, 1>::UnitX(), 0.0);

    int bc1D = s1D.ConstraintLinearSetDisplacementNodeGroup(rightNodes1D, Eigen::Matrix<double, 1, 1>::UnitX(), 0.0);
    int bc2D = s2D.ConstraintLinearSetDisplacementNodeGroup(rightNodes2D, Eigen::Matrix<double, 2, 1>::UnitX(), 0.0);
    int bc3D = s3D.ConstraintLinearSetDisplacementNodeGroup(rightNodes3D, Eigen::Matrix<double, 3, 1>::UnitX(), 0.0);

    s2D.ConstraintLinearSetDisplacementNode(0, Eigen::Vector2d::UnitY(), 0.);
    s3D.ConstraintLinearSetDisplacementNode(0, Eigen::Vector3d::UnitY(), 0.);

    int nFixRotation = s3D.NodeGetIdAtCoordinate(Eigen::Vector3d({0,0,lz}), 1.e-4);
    s3D.ConstraintLinearSetDisplacementNode(nFixRotation, Eigen::Vector3d::UnitY(), 0.);

    Visualize(s1D, "tmp");
    Visualize(s2D, "tmp");
    Visualize(s3D, "tmp");

    s1D.CalculateMaximumIndependentSets();
    s2D.CalculateMaximumIndependentSets();
    s3D.CalculateMaximumIndependentSets();


    NuTo::NewmarkDirect myIntegrationScheme1D(&s1D);
    NuTo::NewmarkDirect myIntegrationScheme2D(&s2D);
    NuTo::NewmarkDirect myIntegrationScheme3D(&s3D);

    SetupNewmark(myIntegrationScheme1D, bc1D, "Newmark1D");
    SetupNewmark(myIntegrationScheme2D, bc2D, "Newmark2D");
    SetupNewmark(myIntegrationScheme3D, bc3D, "Newmark3D");

    timer.Reset(std::string(__FUNCTION__) + " Solution via NewmarkDirect");

    myIntegrationScheme1D.Solve(1.);
    myIntegrationScheme2D.Solve(1.);
    myIntegrationScheme3D.Solve(1.);

    timer.Reset(std::string(__FUNCTION__) + " Checks");


    for (int i = 0; i < numElements; ++i)
    {
        double stress1D = s1D.ElementGetEngineeringStress(i)(0,0);
        double stress2D = s2D.ElementGetEngineeringStress(i)(0,0);
        double stress3D = s3D.ElementGetEngineeringStress(i)(0,0);

        double strain1D = s1D.ElementGetEngineeringStrain(i)(0,0);
        double strain2D = s2D.ElementGetEngineeringStrain(i)(0,0);
        double strain3D = s3D.ElementGetEngineeringStrain(i)(0,0);

        double damage1D = s1D.ElementGetDamage(i)(0,0);
        double damage2D = s2D.ElementGetDamage(i)(0,0);
        double damage3D = s3D.ElementGetDamage(i)(0,0);

        BOOST_CHECK_CLOSE_FRACTION(stress1D, stress2D, 1.e-3);
        BOOST_CHECK_CLOSE_FRACTION(stress1D, stress3D, 1.e-3);

        BOOST_CHECK_CLOSE_FRACTION(strain1D, strain2D, 1.e-3);
        BOOST_CHECK_CLOSE_FRACTION(strain1D, strain3D, 1.e-3);

        BOOST_CHECK_CLOSE_FRACTION(damage1D, damage2D, 1.e-3);
        BOOST_CHECK_CLOSE_FRACTION(damage1D, damage3D, 1.e-3);
    }
    timer.Reset(std::string(__FUNCTION__) + " Cleanup");
}

bool useRobinBoundaryElements = true;

BOOST_AUTO_TEST_CASE(GradientDamage1D)
{
    BOOST_CHECK_NO_THROW(TestStructure1D(useRobinBoundaryElements));
}

BOOST_AUTO_TEST_CASE(GradientDamage2D)
{
    BOOST_CHECK_NO_THROW(TestStructure2D(NuTo::Interpolation::eShapeType::QUAD2D,
                                         false,
                                         useRobinBoundaryElements));

    BOOST_CHECK_NO_THROW(TestStructure2D(NuTo::Interpolation::eShapeType::TRIANGLE2D,
                                         true,
                                         useRobinBoundaryElements));
}

BOOST_AUTO_TEST_CASE(GradientDamage3D)
{
    BOOST_CHECK_NO_THROW(TestStructure3D(NuTo::Interpolation::eShapeType::TETRAHEDRON3D,
                                         useRobinBoundaryElements));
}

BOOST_AUTO_TEST_CASE(GradientDamage1D2D3D)
{
    Check1D2D3D();
}
