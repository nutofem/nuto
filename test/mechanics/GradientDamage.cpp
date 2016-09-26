#include <boost/filesystem.hpp>

#include "nuto/math/FullMatrix.h"
#include "nuto/base/Timer.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/sections/SectionTruss.h"

#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/constitutive/staticData/Leaf.h"
#include "nuto/mechanics/constitutive/laws/GradientDamageEngineeringStress.h"
#include "nuto/mechanics/elements/ElementDataEnum.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/groups/GroupBase.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/sections/SectionEnum.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/mechanics/timeIntegration/ImplEx.h"

#include "nuto/mechanics/elements/ContinuumBoundaryElement.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"
#include "nuto/visualize/VisualizeEnum.h"

//#define PRINTRESULT
#include <eigen3/Eigen/Eigenvalues>

namespace NuToTest {
namespace GradientDamage {


bool CheckDamageLawsDerivatives(NuTo::GradientDamageEngineeringStress rConstitutiveLaw)
{
    double epsilon = 1.e-8;
    double E = rConstitutiveLaw.GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS);
    double e0 = rConstitutiveLaw.GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH) / E;
    double step = e0 / 5;
    for (int i = 1; i < 100; ++i)
    {
        double kappa = i * step + epsilon;
//        kappa = i*step;
        double sigma1 = (1 - rConstitutiveLaw.CalculateDamage(kappa)) * E * kappa;
        double sigma2 = (1 - rConstitutiveLaw.CalculateDamage(kappa + epsilon)) * E * (kappa + epsilon);

        double DsigmaDkappa = -rConstitutiveLaw.CalculateDerivativeDamage(kappa) * E * kappa + (1 - rConstitutiveLaw.CalculateDamage(kappa)) * E;
        double DsigmaDkappa_CDF = (sigma2 - sigma1) / epsilon;

        double differenceSigma = DsigmaDkappa - DsigmaDkappa_CDF;

#ifdef PRINTRESULT
        std::cout << "kappa:" << kappa << " | differenceSigma: " << differenceSigma << std::endl;
        std::cout << "Dsigma:" << DsigmaDkappa << " | Dsigma_CDF: " << DsigmaDkappa_CDF << std::endl;
        std::cout << "sigma1:" << sigma1 << " | sigma2: " << sigma2 << std::endl << std::endl;
#endif

        if (std::abs(differenceSigma) > 1.e-3)
            return false;
    }

    return true;
}

void CheckDamageLaws()
{
    NuTo::GradientDamageEngineeringStress myConstitutiveLaw;

    myConstitutiveLaw.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::DENSITY,1.0);
    myConstitutiveLaw.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,30000);
    myConstitutiveLaw.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,0.3);
    myConstitutiveLaw.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS,1.0);
    myConstitutiveLaw.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS_PARAMETER,10.0);
    myConstitutiveLaw.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH,4.);
    myConstitutiveLaw.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY,0.21);

    myConstitutiveLaw.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::DAMAGE_LAW, static_cast<double>(NuTo::Constitutive::eDamageLawType::ISOTROPIC_NO_SOFTENING));
    if (not CheckDamageLawsDerivatives(myConstitutiveLaw))
        throw NuTo::MechanicsException("DamageLaw::ISOTROPIC_NO_SOFTENING: wrong damage derivatives");

    myConstitutiveLaw.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::DAMAGE_LAW, static_cast<double>(NuTo::Constitutive::eDamageLawType::ISOTROPIC_LINEAR_SOFTENING));
    if (not CheckDamageLawsDerivatives(myConstitutiveLaw))
        throw NuTo::MechanicsException("DamageLaw::ISOTROPIC_LINEAR_SOFTENING: wrong damage derivatives");

    myConstitutiveLaw.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::DAMAGE_LAW, static_cast<double>(NuTo::Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING));
    if (not CheckDamageLawsDerivatives(myConstitutiveLaw))
        throw NuTo::MechanicsException("DamageLaw::ISOTROPIC_EXPONENTIAL_SOFTENING: wrong damage derivatives");

    myConstitutiveLaw.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::DAMAGE_LAW, static_cast<double>(NuTo::Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_RES_LOAD));
    if (not CheckDamageLawsDerivatives(myConstitutiveLaw))
        throw NuTo::MechanicsException("DamageLaw::ISOTROPIC_EXPONENTIAL_SOFTENING_RES_LOAD: wrong damage derivatives");


    myConstitutiveLaw.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::DAMAGE_LAW, static_cast<double>(NuTo::Constitutive::eDamageLawType::ISOTROPIC_CUBIC_HERMITE));
    if (not CheckDamageLawsDerivatives(myConstitutiveLaw))
        throw NuTo::MechanicsException("DamageLaw::ISOTROPIC_CUBIC_HERMITE: wrong damage derivatives");

}

int SetConstitutiveLaw(NuTo::Structure& rStructure)
{
    // create a damage law
    int lawId = rStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::DENSITY, 1.0);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 30000);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.2);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS, 3);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS_PARAMETER, 0.);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH, 4.);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH, 4. * 10);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY, 0.021);
    rStructure.ConstitutiveLawSetDamageLaw(lawId, NuTo::Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING);

    return lawId;
}

void ApplyDofValues(NuTo::Structure& rStructure)
{
    double boundaryDisplacement = 1.e-1;

    rStructure.CalculateMaximumIndependentSets();
    rStructure.NodeBuildGlobalDofs();

    int gAllNodes = rStructure.GroupGetNodesTotal();
    auto nodeIds = rStructure.GroupGetMemberIds(gAllNodes);
    for (int i = 0; i < nodeIds.GetNumRows(); ++i)
    {
        NuTo::NodeBase* node = rStructure.NodeGetNodePtr(nodeIds.GetValue(i));
        Eigen::VectorXd disps = node->Get(NuTo::Node::eDof::COORDINATES) / 100. * boundaryDisplacement;
        node->Set(NuTo::Node::eDof::DISPLACEMENTS, disps);

        if (node->GetNum(NuTo::Node::eDof::NONLOCALEQSTRAIN) > 0)
            node->Set(NuTo::Node::eDof::NONLOCALEQSTRAIN, disps(0, 0) / 10);
    }
}

void CheckStiffnesses(NuTo::Structure& rStructure)
{
    ApplyDofValues(rStructure);

//    rStructure.Info();

    bool elementStiffnessCorrect = rStructure.ElementCheckHessian0(1.e-8, 1.e-4, true);
    if (not elementStiffnessCorrect)
        throw NuTo::Exception(__PRETTY_FUNCTION__, "element stiffness matrices incorrect!");

    bool globalStiffnessCorrect = rStructure.CheckHessian0(1.e-8, 1.e-4, true);
    if (not globalStiffnessCorrect)
        throw NuTo::Exception(__PRETTY_FUNCTION__,  "global stiffness matrix incorrect!");
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

int AddInterpolationType(NuTo::Structure& rStructure, NuTo::Interpolation::eShapeType rShape)
{
    int it = rStructure.InterpolationTypeCreate(rShape);
    rStructure.InterpolationTypeAdd(it, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    rStructure.InterpolationTypeAdd(it, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    rStructure.InterpolationTypeAdd(it, NuTo::Node::eDof::NONLOCALEQSTRAIN, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    return it;
}

int AddBoundaryElementsAtCoordinate(NuTo::Structure &s, double rCoordinate)
{
    int gNodesBoundary = s.GroupCreate("NODES");
    s.GroupAddNodeCoordinateRange(gNodesBoundary, 0, rCoordinate - 1.e-6, rCoordinate + 1.e-6);

    int elemGroupBoundary = s.GroupCreate("ELEMENTS");
    s.GroupAddElementsFromNodes(elemGroupBoundary, gNodesBoundary, false);

    int gBoundaryElements = s.BoundaryElementsCreate(elemGroupBoundary, gNodesBoundary);
    return s.GroupGetNumMembers(gBoundaryElements);
}

void AddBoundaryElements(NuTo::Structure& s, double rLength, int rNumExpectedBoundaryElements)
{
    int numBoundaryElements = 0;
    numBoundaryElements += AddBoundaryElementsAtCoordinate(s, 0);
    numBoundaryElements += AddBoundaryElementsAtCoordinate(s, rLength);

    if (numBoundaryElements != rNumExpectedBoundaryElements)
    {
        throw NuTo::MechanicsException("Wrong number of boundary elements. \nCreated: " + std::to_string(numBoundaryElements) + "\nExpected: " + std::to_string(rNumExpectedBoundaryElements));
    }
}

void TestStructure1D(bool rUseRobinBoundaryElements)
{
    NuTo::Timer timer(__FUNCTION__);
    const int numElements = 2;
    const double length = 100;
    const double area = 10;

    NuTo::Structure s(1);
    s.SetShowTime(false);
    s.SetVerboseLevel(0);

    // create nodes
    int numNodes = numElements + 1;
    double lengthElement = length / numElements;

    NuTo::FullVector<double, Eigen::Dynamic> nodeCoordinates(1);
    for (int iNode = 0; iNode < numNodes; ++iNode)
    {
        nodeCoordinates(0) = iNode * lengthElement;
        s.NodeCreate(iNode, nodeCoordinates);
    }

    int interpolationType = AddInterpolationType(s, NuTo::Interpolation::eShapeType::TRUSS1D);

    // create elements
    NuTo::FullVector<int, Eigen::Dynamic> nodes(2);
    for (int iElement = 0; iElement < numElements; iElement++)
    {
        nodes(0) = iElement;
        nodes(1) = iElement + 1;
        s.ElementCreate(interpolationType, nodes, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
    }


    // create sections
    int mySection = s.SectionCreate("Truss");
    s.SectionSetArea(mySection, area);

    // set function for area reduction
    NuTo::SectionTruss* secTruss = s.SectionGetSectionPtr(mySection)->AsSectionTruss();
    const double xWeakSpot = 25;
    const double lWeakSpot = 5;
    const double alpha = 0.10;
    const double exponent = 4;
    double areaParameters[4];
    areaParameters[0] = xWeakSpot;
    areaParameters[1] = lWeakSpot;
    areaParameters[2] = alpha;
    areaParameters[3] = exponent;
    secTruss->SetAreaParameters(areaParameters);

    s.ElementTotalConvertToInterpolationType();
    s.ElementTotalSetConstitutiveLaw(SetConstitutiveLaw(s));
    s.ElementTotalSetSection(mySection);

    if (rUseRobinBoundaryElements)
    {
        AddBoundaryElements(s, length, 2);
    }
    CheckStiffnesses(s);
    Visualize(s, "TRUSS1D");
}



void TestStructure2D(NuTo::Interpolation::eShapeType rShape, NuTo::eSectionType rSection, bool rUseRobinBoundaryElements)
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

    //create nodes
    int numNodesX = numElementsX + 1;
    int numNodesY = numElementsY + 1;
    double deltaX = lX / (numElementsX);
    double deltaY = lY / (numElementsY);

    int nodeNum = 0;
    for (int countY = 0; countY < numNodesY; countY++)
    {
        for (int countX = 0; countX < numNodesX; countX++)
        {
            NuTo::FullVector<double, Eigen::Dynamic> coordinates(2);
            coordinates(0) = countX * deltaX;
            coordinates(1) = countY * deltaY;
            s.NodeCreate(nodeNum, coordinates);
            nodeNum++;
        }
    }

    int myInterpolationType = AddInterpolationType(s, rShape);

    //create elements
    for (int countY = 0; countY < numElementsY; countY++)
    {
        for (int countX = 0; countX < numElementsX; countX++)
        {
            if (rShape == NuTo::Interpolation::eShapeType::QUAD2D)
            {
                NuTo::FullVector<int, Eigen::Dynamic> nodes(4);
                nodes(0) = countX + countY * numNodesX;
                nodes(1) = countX + 1 + countY * numNodesX;
                nodes(2) = countX + 1 + (countY + 1) * numNodesX;
                nodes(3) = countX + (countY + 1) * numNodesX;
                s.ElementCreate(myInterpolationType, nodes, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
            }
            if (rShape == NuTo::Interpolation::eShapeType::TRIANGLE2D)
            {
                NuTo::FullVector<int, Eigen::Dynamic> nodes(3);
                nodes(0) = countX + countY * numNodesX;
                nodes(1) = countX + 1 + countY * numNodesX;
                nodes(2) = countX + 1 + (countY + 1) * numNodesX;
                s.ElementCreate(myInterpolationType, nodes, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);

                nodes(0) = countX + countY * numNodesX;
                nodes(1) = countX + 1 + (countY + 1) * numNodesX;
                nodes(2) = countX + (countY + 1) * numNodesX;
                s.ElementCreate(myInterpolationType, nodes, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
            }
        }
    }

    int myConstitutiveLaw = SetConstitutiveLaw(s);
    int mySection = s.SectionCreate(rSection);
    s.SectionSetThickness(mySection, lZ);

    s.ElementTotalConvertToInterpolationType();
    s.ElementTotalSetConstitutiveLaw(myConstitutiveLaw);
    s.ElementTotalSetSection(mySection);

    if (rUseRobinBoundaryElements)
    {
        AddBoundaryElements(s, lX, 2*numElementsY);
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
    NuTo::FullVector<int, Eigen::Dynamic> nodeIds(8);
    nodeIds[0] = s.NodeCreate(NuTo::FullVector<double, 3> ({ 0, 0, 0}));
    nodeIds[1] = s.NodeCreate(NuTo::FullVector<double, 3> ({lX, 0, 0}));
    nodeIds[2] = s.NodeCreate(NuTo::FullVector<double, 3> ({lX,lY, 0}));
    nodeIds[3] = s.NodeCreate(NuTo::FullVector<double, 3> ({ 0,lY, 0}));
    nodeIds[4] = s.NodeCreate(NuTo::FullVector<double, 3> ({ 0, 0,lZ}));
    nodeIds[5] = s.NodeCreate(NuTo::FullVector<double, 3> ({lX, 0,lZ}));
    nodeIds[6] = s.NodeCreate(NuTo::FullVector<double, 3> ({lX,lY,lZ}));
    nodeIds[7] = s.NodeCreate(NuTo::FullVector<double, 3> ({ 0,lY,lZ}));

    int myInterpolationType = AddInterpolationType(s, rShape);

    if (rShape == NuTo::Interpolation::eShapeType::BRICK3D)
        s.ElementCreate(myInterpolationType, nodeIds, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
    if (rShape == NuTo::Interpolation::eShapeType::TETRAHEDRON3D)
    {
        NuTo::FullVector<int,Eigen::Dynamic> nodesTet0(4);
        NuTo::FullVector<int,Eigen::Dynamic> nodesTet1(4);
        NuTo::FullVector<int,Eigen::Dynamic> nodesTet2(4);
        NuTo::FullVector<int,Eigen::Dynamic> nodesTet3(4);
        NuTo::FullVector<int,Eigen::Dynamic> nodesTet4(4);
        NuTo::FullVector<int,Eigen::Dynamic> nodesTet5(4);

        nodesTet0 << nodeIds(0), nodeIds(1), nodeIds(3), nodeIds(7);
        nodesTet1 << nodeIds(0), nodeIds(1), nodeIds(7), nodeIds(4);
        nodesTet2 << nodeIds(5), nodeIds(4), nodeIds(7), nodeIds(1);
        nodesTet3 << nodeIds(6), nodeIds(5), nodeIds(7), nodeIds(1);
        nodesTet4 << nodeIds(2), nodeIds(7), nodeIds(1), nodeIds(6);
        nodesTet5 << nodeIds(2), nodeIds(3), nodeIds(1), nodeIds(7);

        s.ElementCreate(myInterpolationType, nodesTet0, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
        s.ElementCreate(myInterpolationType, nodesTet1, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
        s.ElementCreate(myInterpolationType, nodesTet2, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
        s.ElementCreate(myInterpolationType, nodesTet3, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
        s.ElementCreate(myInterpolationType, nodesTet4, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
        s.ElementCreate(myInterpolationType, nodesTet5, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
    }
    int mySection = s.SectionCreate("VOLUME");

    s.ElementTotalConvertToInterpolationType();
    s.ElementTotalSetConstitutiveLaw(SetConstitutiveLaw(s));
    s.ElementTotalSetSection(mySection);


    if (rUseRobinBoundaryElements)
    {
        int numExpectedBoundaryElements = 2; // for brick
        if (rShape == NuTo::Interpolation::eShapeType::TETRAHEDRON3D)
            numExpectedBoundaryElements = 4;

        AddBoundaryElements(s, lX, numExpectedBoundaryElements);
    }

    CheckStiffnesses(s);
    Visualize(s, NuTo::Interpolation::ShapeTypeToString(rShape));
}


void GroupRemoveNodesWithoutDisplacements(NuTo::Structure& rStructure, int rGroupNodeId)
{
    auto ids = rStructure.GroupGetMemberIds(rGroupNodeId);
    for (int i = 0; i < ids.GetNumRows(); ++i)
    {
        int nodeId = ids.GetValue(i);
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

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> timeDepDisp(2, 2);
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
    int numNodesX = numElements + 1;
    double lengthElementX = lx / numElements;


    NuTo::Structure s1D(1);
    NuTo::Structure s2D(2);
    NuTo::Structure s3D(3);

    s1D.SetShowTime(false);
    s2D.SetShowTime(false);
    s3D.SetShowTime(false);

    int interpolationType1D = AddInterpolationType(s1D, NuTo::Interpolation::eShapeType::TRUSS1D);
    s1D.InterpolationTypeSetIntegrationType(interpolationType1D, NuTo::eIntegrationType::IntegrationType1D2NGauss2Ip, NuTo::IpData::eIpDataType::STATICDATA);

    int interpolationType2D = AddInterpolationType(s2D, NuTo::Interpolation::eShapeType::QUAD2D);
    s2D.InterpolationTypeSetIntegrationType(interpolationType2D, NuTo::eIntegrationType::IntegrationType2D4NGauss4Ip, NuTo::IpData::eIpDataType::STATICDATA);

    int interpolationType3D = AddInterpolationType(s3D, NuTo::Interpolation::eShapeType::BRICK3D);
    s3D.InterpolationTypeSetIntegrationType(interpolationType3D, NuTo::eIntegrationType::IntegrationType3D8NGauss2x2x2Ip, NuTo::IpData::eIpDataType::STATICDATA);


    NuTo::FullVector<double, Eigen::Dynamic> nodeCoordinates(1);
    for (int iNode = 0; iNode < numNodesX; ++iNode)
    {
        nodeCoordinates(0) = iNode * lengthElementX; // two nodes per element
        s1D.NodeCreate(iNode, nodeCoordinates);
    }
    NuTo::FullVector<int, Eigen::Dynamic> elementIncidence(2);
    for (int iElement = 0; iElement < numElements; iElement++)
    {
        elementIncidence(0) = iElement;
        elementIncidence(1) = iElement + 1;
        s1D.ElementCreate(interpolationType1D, elementIncidence, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
    }

    int nodeNum = 0;
    for (int countY = 0; countY < 2; countY++)
        for (int countX = 0; countX < numNodesX; countX++)
        {
            NuTo::FullVector<double, Eigen::Dynamic> coordinates(2);
            coordinates(0) = countX * lengthElementX;
            coordinates(1) = countY * ly;
            s2D.NodeCreate(nodeNum, coordinates);
            nodeNum++;
        }

    //create elements
    for (int countX = 0; countX < numElements; countX++)
    {
        NuTo::FullVector<int, Eigen::Dynamic> nodes(4);
        nodes(0) = countX;
        nodes(1) = countX + 1;
        nodes(2) = countX + 1 + numNodesX;
        nodes(3) = countX + numNodesX;
        s2D.ElementCreate(interpolationType2D, nodes, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
    }

    nodeNum = 0;
    for (int iZ=0; iZ<2; iZ++)
        for (int iY=0; iY<2; iY++)
            for (int iX=0; iX<numNodesX; iX++)
            {
                NuTo::FullVector<double,Eigen::Dynamic> coordinates(3);
                coordinates(0) = iX*lengthElementX;
                coordinates(1) = iY*ly;
                coordinates(2) = iZ*lz;
                s3D.NodeCreate(nodeNum, coordinates);
                nodeNum++;
            }

    for (int iX=0; iX<numElements; iX++)
    {
        NuTo::FullVector<int,Eigen::Dynamic> nodes(8);
        nodes(0) = iX;
        nodes(1) = iX+1;
        nodes(2) = iX+1 +  2 * numNodesX;
        nodes(3) = iX   +  2 * numNodesX;
        nodes(4) = iX                    + numNodesX;
        nodes(5) = iX+1                  + numNodesX;
        nodes(6) = iX+1 +  2 * numNodesX + numNodesX;
        nodes(7) = iX   +  2 * numNodesX + numNodesX;

        s3D.ElementCreate(interpolationType3D, nodes, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
    }
    s1D.ElementTotalConvertToInterpolationType();
    s2D.ElementTotalConvertToInterpolationType();
    s3D.ElementTotalConvertToInterpolationType();

    s1D.ElementTotalSetConstitutiveLaw(SetConstitutiveLaw(s1D));
    s2D.ElementTotalSetConstitutiveLaw(SetConstitutiveLaw(s2D));
    s3D.ElementTotalSetConstitutiveLaw(SetConstitutiveLaw(s3D));

    int mySection1D  = s1D.SectionCreate("Truss");
    int mySection2D  = s2D.SectionCreate("Plane_Stress");
    int mySection3D  = s3D.SectionCreate("Volume");

    s1D.SectionSetArea(mySection1D, lz*ly);
    s2D.SectionSetThickness(mySection2D, lz);

    s1D.ElementTotalSetSection(mySection1D);
    s2D.ElementTotalSetSection(mySection2D);
    s3D.ElementTotalSetSection(mySection3D);

    int weakElementId = numElements / 2;
    double kappa = 0.001;
    {
        NuTo::ElementBase* element = s1D.ElementGetElementPtr(weakElementId);
        for (int i = 0; i < element->GetNumIntegrationPoints(); ++i)
            dynamic_cast<NuTo::Constitutive::StaticData::Leaf<double>*>(element->GetConstitutiveStaticData(i))->SetData(kappa);
    }
    {
        NuTo::ElementBase* element = s2D.ElementGetElementPtr(weakElementId);
        for (int i = 0; i < element->GetNumIntegrationPoints(); ++i)
            dynamic_cast<NuTo::Constitutive::StaticData::Leaf<double>*>(element->GetConstitutiveStaticData(i))->SetData(kappa);
    }
    {
        NuTo::ElementBase* element = s3D.ElementGetElementPtr(weakElementId);
        for (int i = 0; i < element->GetNumIntegrationPoints(); ++i)
            dynamic_cast<NuTo::Constitutive::StaticData::Leaf<double>*>(element->GetConstitutiveStaticData(i))->SetData(kappa);
    }


    int leftNodes1D = s1D.GroupCreate("NODES");
    int leftNodes2D = s2D.GroupCreate("NODES");
    int leftNodes3D = s3D.GroupCreate("NODES");

    int rightNodes1D = s1D.GroupCreate("NODES");
    int rightNodes2D = s2D.GroupCreate("NODES");
    int rightNodes3D = s3D.GroupCreate("NODES");

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

    s1D.ConstraintLinearSetDisplacementNodeGroup(leftNodes1D, NuTo::FullVector<double, 1>::UnitX(), 0.0);
    s2D.ConstraintLinearSetDisplacementNodeGroup(leftNodes2D, NuTo::FullVector<double, 2>::UnitX(), 0.0);
    s3D.ConstraintLinearSetDisplacementNodeGroup(leftNodes3D, NuTo::FullVector<double, 3>::UnitX(), 0.0);

    int bc1D = s1D.ConstraintLinearSetDisplacementNodeGroup(rightNodes1D, NuTo::FullVector<double, 1>::UnitX(), 0.0);
    int bc2D = s2D.ConstraintLinearSetDisplacementNodeGroup(rightNodes2D, NuTo::FullVector<double, 2>::UnitX(), 0.0);
    int bc3D = s3D.ConstraintLinearSetDisplacementNodeGroup(rightNodes3D, NuTo::FullVector<double, 3>::UnitX(), 0.0);

    s2D.ConstraintLinearSetDisplacementNode(0, NuTo::FullVector<double, 2>::UnitY(), 0.);
    s3D.ConstraintLinearSetDisplacementNode(0, NuTo::FullVector<double, 3>::UnitY(), 0.);

    int nFixRotation = s3D.NodeGetIdAtCoordinate(NuTo::FullVector<double, 3>({0,0,lz}), 1.e-4);
    s3D.ConstraintLinearSetDisplacementNode(nFixRotation, NuTo::FullVector<double, 3>::UnitY(), 0.);

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

        // assert that everything is almost equal
        auto notEqual = [](double a, double b) { return std::abs(a - b) > 1.e-3; };

        if (notEqual(stress1D, stress2D) or notEqual(stress1D, stress3D) )
            throw NuTo::MechanicsException("Stresses not equal. Better check the plots in " + myIntegrationScheme1D.GetResultDirectory());

        if (notEqual(strain1D, strain2D) or notEqual(strain1D, strain3D) )
            throw NuTo::MechanicsException("Strains not equal. Better check the plots in " + myIntegrationScheme1D.GetResultDirectory());

        if (notEqual(damage1D, damage2D) or notEqual(damage1D, damage3D) )
            throw NuTo::MechanicsException("Damage not equal. Better check the plots in " + myIntegrationScheme1D.GetResultDirectory());
    }

    return;
    timer.Reset(std::string(__FUNCTION__) + " Cleanup");
}

}  // namespace NuToTest

void CheckBoundaryElements()
{

}

}  // namespace GradientDamage

int main()
{


    try
    {
        NuTo::Timer Timer("GradientDamage", true);

        NuToTest::GradientDamage::CheckDamageLaws();

        bool useRobinBoundaryElements = false;

        NuToTest::GradientDamage::TestStructure1D(useRobinBoundaryElements);
        NuToTest::GradientDamage::TestStructure2D(NuTo::Interpolation::eShapeType::QUAD2D, NuTo::eSectionType::PLANE_STRESS, useRobinBoundaryElements);
        NuToTest::GradientDamage::TestStructure2D(NuTo::Interpolation::eShapeType::TRIANGLE2D, NuTo::eSectionType::PLANE_STRAIN, useRobinBoundaryElements);
        NuToTest::GradientDamage::TestStructure3D(NuTo::Interpolation::eShapeType::BRICK3D, useRobinBoundaryElements);
//      NuToTest::GradientDamage::TestStructure3D(NuTo::Interpolation::eShapeType::TETRAHEDRON3D, useRobinBoundaryElements);

        NuToTest::GradientDamage::Check1D2D3D();


    }
    catch (NuTo::MechanicsException& e)
    {
        std::cout << e.ErrorMessage();
        return EXIT_FAILURE;
    }
    catch (...)
    {
        std::cout << "Something else went wrong." << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << std::endl;
    std::cout << "#####################################" << std::endl;
    std::cout << "##  Congratulations! Everything    ##" << std::endl;
    std::cout << "##   went better than expected!    ##" << std::endl;
    std::cout << "#####################################" << std::endl;

    return EXIT_SUCCESS;

}
