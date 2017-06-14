
#include <mpi.h>


#include "mechanics/feti/StructureFeti.h"

#include <ctime>
#include <chrono>
#include <iomanip>
#include "mechanics/feti/NewmarkFeti.h"

#include "mechanics/nodes/NodeBase.h"

#include "boost/filesystem.hpp"

#include "mechanics/DirectionEnum.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/groups/Group.h"
#include "visualize/VisualizeEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/elements/IpDataEnum.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"

using std::cout;
using std::endl;
using NuTo::Constitutive::ePhaseFieldEnergyDecomposition;
using NuTo::Constitutive::eConstitutiveType;
using NuTo::Constitutive::eConstitutiveParameter;
using NuTo::Node::eDof;
using NuTo::Interpolation::eTypeOrder;
using NuTo::Interpolation::eShapeType;
using NuTo::eGroupId;
using NuTo::eVisualizeWhat;
using NuTo::eDirection;
using NuTo::Constraint::Component;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Vector2d;
using Eigen::Matrix2d;
using EigenSolver = Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>;
using FetiPreconditioner = NuTo::NewmarkFeti<EigenSolver>::eFetiPreconditioner;
using FetiIterativeSolver = NuTo::NewmarkFeti<EigenSolver>::eIterativeSolver;
using FetiScaling = NuTo::NewmarkFeti<EigenSolver>::eFetiScaling;
using namespace NuTo;

// geometry
constexpr int dimension = 2;
constexpr double thickness = 1.0;

// material
constexpr double youngsModulus = 4.0e4;
constexpr double poissonsRatio = 0.2;

// integration
constexpr double timeStep = 1.0;
constexpr double simulationTime = 1.0;
constexpr double loadFactor = 13.37;

void AssignSection(NuTo::Structure& structure);
void AssignMaterial(NuTo::Structure& structure);

std::map<int, VectorXd> GetNodeAndDisplacementsMap(NuTo::Structure& structure);

std::map<int, VectorXd> ComputeReferenceSolution(int rank);

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    NuTo::StructureFeti structure(dimension);
    structure.GetLogger().OpenFile("FetiTestOutputRank_" + std::to_string(rank));
    structure.GetLogger().SetQuiet(true);

    std::string meshFile = "feti_beam_coarse_2_subdomains_24_ele.mesh" + std::to_string(rank);

    const int interpolationTypeId = structure.InterpolationTypeCreate(eShapeType::QUAD2D);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES, eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT1);

    structure.ImportMeshJson(meshFile, interpolationTypeId);

    AssignMaterial(structure);
    AssignSection(structure);


    structure.GetLogger() << "*********************************** \n"
                          << "**      node groups              ** \n"
                          << "*********************************** \n\n";

    Vector2d coordinate;

    coordinate[0] = 0;
    coordinate[1] = 0;
    const auto& groupNodesAtBottomLeft = structure.GroupGetNodeRadiusRange(coordinate, 0., 1.e-6);

    coordinate[0] = 60;
    coordinate[1] = 0;
    const auto& groupNodeAtBottomRight = structure.GroupGetNodeRadiusRange(coordinate, 0., 1.e-6);

    const auto groupNodesAtBoundaries = Group<NodeBase>::Unite(groupNodeAtBottomRight, groupNodesAtBottomLeft);

    coordinate[0] = 20;
    coordinate[1] = 10;
    const auto& groupNodesLoad = structure.GroupGetNodeRadiusRange(coordinate, 0, 1.e-6);

    structure.GetLogger() << "*********************************** \n"
                          << "**      virtual constraints      ** \n"
                          << "*********************************** \n\n";

    std::vector<int> nodeIdsBoundaries = groupNodesAtBoundaries.GetMemberIds();

    std::vector<int> nodeIdsLoads = groupNodesLoad.GetMemberIds();

    structure.ApplyVirtualConstraints(nodeIdsBoundaries, nodeIdsLoads);

    structure.GetLogger() << "*********************************** \n"
                          << "**      real constraints         ** \n"
                          << "*********************************** \n\n";

    structure.NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

    structure.ApplyConstraintsTotalFeti(groupNodesAtBoundaries);

    structure.GetLogger() << "*********************************** \n"
                          << "**      load                     ** \n"
                          << "*********************************** \n\n";

    // prescribe displacement of groupNodesLoad in Y direction
    std::map<int, double> dofIdAndPrescribedDisplacementMap;
    for (auto const& nodeId : nodeIdsLoads)
    {
        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);
        dofIdAndPrescribedDisplacementMap.emplace(dofIds[1], 1.);
    }

    structure.ApplyPrescribedDisplacements(dofIdAndPrescribedDisplacementMap);

    int loadId = structure.LoadCreateNodeGroupForce(&groupNodesLoad, Vector2d::UnitY(), 0.);


    structure.GetLogger() << "*********************************** \n"
                          << "**      visualization            ** \n"
                          << "*********************************** \n\n";

    int groupAllElements = 9999;
    structure.GroupCreate(groupAllElements, eGroupId::Elements);
    structure.GroupAddElementsTotal(groupAllElements);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DISPLACEMENTS);

    structure.GetLogger() << "*********************************** \n"
                          << "**      integration scheme       ** \n"
                          << "*********************************** \n\n";

    NuTo::NewmarkFeti<EigenSolver> myIntegrationScheme(&structure);
    boost::filesystem::path resultPath(boost::filesystem::initial_path().string() + "/FetiSolutionRank_" +
                                       std::to_string(structure.mRank));

    myIntegrationScheme.SetTimeStep(timeStep);
    myIntegrationScheme.SetResultDirectory(resultPath.string(), true);
    myIntegrationScheme.SetToleranceIterativeSolver(1.e-7);

    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    myIntegrationScheme.SetTimeDependentLoadCase(loadId, dispRHS);

    structure.GetLogger() << "*********************************** \n"
                          << "**      solve                    ** \n"
                          << "*********************************** \n\n";

    myIntegrationScheme.Solve(simulationTime);

    std::map<int, VectorXd> nodeIdAndDisplacements = GetNodeAndDisplacementsMap(structure);
    std::map<int, VectorXd> nodeIdAndDisplacementsReference = ComputeReferenceSolution(rank);


    for (const auto& nodeIdAndDisplacementsPair : nodeIdAndDisplacements)
    {
        VectorXd vec00 = nodeIdAndDisplacementsPair.second;
        VectorXd vec01 = nodeIdAndDisplacementsReference[nodeIdAndDisplacementsPair.first - 1];

        if (not(vec00 - vec01).isMuchSmallerThan(1.e-4, 1.))
        {
            std::cout << "nodeId: \n" << nodeIdAndDisplacementsPair.first << std::endl;
            std::cout << std::setprecision(6) << std::endl;
            std::cout << "vec00 \n" << vec00 << std::endl;

            std::cout << "vec01 \n" << vec01 << std::endl;

            std::cout << "Solutions are not equal" << std::endl;

            return EXIT_FAILURE;
        }
    }


    MPI_Finalize();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AssignSection(NuTo::Structure& structure)
{
    structure.GetLogger() << "*********************************** \n"
                          << "**      section                  ** \n"
                          << "*********************************** \n\n";

    auto section00 = NuTo::SectionPlane::Create(thickness, false);
    structure.ElementTotalSetSection(section00);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AssignMaterial(NuTo::Structure& structure)
{
    structure.GetLogger() << "*********************************** \n"
                          << "**      material                 ** \n"
                          << "*********************************** \n\n";

    int material00 = structure.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);

    structure.ElementTotalSetConstitutiveLaw(material00);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::map<int, VectorXd> GetNodeAndDisplacementsMap(NuTo::Structure& structure)
{
    std::map<int, VectorXd> nodeIdAndDisplacements;
    for (const auto& nodeIdPtrPair : structure.NodeGetNodeMap())
    {
        VectorXd displacements;
        structure.NodeGetDisplacements(nodeIdPtrPair.first, displacements);
        nodeIdAndDisplacements.emplace(nodeIdPtrPair.first, displacements);
    }

    return nodeIdAndDisplacements;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::map<int, VectorXd> ComputeReferenceSolution(int rank)
{

    NuTo::Structure structure(dimension);
    structure.SetShowTime(false);
    structure.GetLogger().OpenFile("FetiTestOutputReferenceSolution_rank" + std::to_string(rank));

    std::string meshFile = "feti_beam_coarse_2_subdomains_24_ele_compare.msh";
    auto eleGroupAndInterpolationTypeList = structure.ImportFromGmsh(meshFile);

    const int interpolationTypeId = eleGroupAndInterpolationTypeList[0].second;
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeSetIntegrationType(interpolationTypeId, eIntegrationType::IntegrationType2D4NGauss4Ip);

    structure.ElementTotalConvertToInterpolationType();

    AssignSection(structure);
    AssignMaterial(structure);

    structure.NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

    structure.GetLogger() << "*********************************** \n"
                          << "**      real constraints         ** \n"
                          << "*********************************** \n\n";

    Vector2d coordinate;
    coordinate[0] = 0;
    coordinate[1] = 0;

    auto& groupNodesLeftBoundary = structure.GroupGetNodeRadiusRange(coordinate, 0, 1.e-6);
    structure.Constraints().Add(eDof::DISPLACEMENTS,
                                Constraint::Component(groupNodesLeftBoundary, {eDirection::X, eDirection::Y}));

    coordinate[0] = 60;
    coordinate[1] = 0;

    auto& groupNodesRightBoundary = structure.GroupGetNodeRadiusRange(coordinate, 0, 1.e-6);
    structure.Constraints().Add(eDof::DISPLACEMENTS,
                                Constraint::Component(groupNodesRightBoundary, {eDirection::X, eDirection::Y}));

    structure.GetLogger() << "*********************************** \n"
                          << "**      load                     ** \n"
                          << "*********************************** \n\n";

    coordinate[0] = 20;
    coordinate[1] = 10;

    auto& loadNodeGroup = structure.GroupGetNodeRadiusRange(coordinate, 0, 1.e-6);
    structure.Constraints().Add(
            eDof::DISPLACEMENTS,
            Constraint::Component(loadNodeGroup, {eDirection::Y}, Constraint::RhsRamp(simulationTime, loadFactor)));

    structure.GetLogger() << "*********************************** \n"
                          << "**      intergration scheme      ** \n"
                          << "*********************************** \n\n";

    NuTo::NewmarkDirect newmarkDirect(&structure);
    boost::filesystem::path resultPath(boost::filesystem::initial_path().string() + "/FetiReferenceSolution_rank" +
                                       std::to_string(rank));

    std::cout << "resultPath" << resultPath.string() << std::endl;
    newmarkDirect.SetTimeStep(timeStep);
    newmarkDirect.SetResultDirectory(resultPath.string(), true);

    structure.GetLogger() << "*********************************** \n"
                          << "**      solve                    ** \n"
                          << "*********************************** \n\n";

    newmarkDirect.Solve(simulationTime);

    structure.GetLogger() << "*********************************** \n"
                          << "**      post process             ** \n"
                          << "*********************************** \n\n";


    return GetNodeAndDisplacementsMap(structure);
}
