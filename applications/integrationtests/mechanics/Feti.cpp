
#include <mpi.h>
#include <ctime>
#include <chrono>
#include <iomanip>

#include "boost/filesystem.hpp"

#include "mechanics/feti/StructureFeti.h"
#include "mechanics/feti/NewmarkFeti.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/DirectionEnum.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/groups/Group.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/elements/IpDataEnum.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"

#include "visualize/VisualizeEnum.h"

using std::cout;
using std::endl;
using namespace NuTo;
using namespace Constitutive;
using namespace Interpolation;
using Node::eDof;
using Constraint::Component;
using Eigen::VectorXd;
using Eigen::Vector2d;
using Eigen::Matrix2d;
using EigenSolver = Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>;
using FetiPreconditioner = NewmarkFeti<EigenSolver>::eFetiPreconditioner;
using FetiIterativeSolver = NewmarkFeti<EigenSolver>::eIterativeSolver;
using FetiScaling = NewmarkFeti<EigenSolver>::eFetiScaling;

// geometry
constexpr int dimension = 2;
constexpr double thickness = 1.0;
const Vector2d coordinateAtBottomLeft(0., 0.);
const Vector2d coordinateAtBottomRight(60., 0.);
const Vector2d coordinateAtLoad(20., 10.);

// material
constexpr double youngsModulus = 4.0e4;
constexpr double poissonsRatio = 0.2;

// integration
constexpr double timeStep = 1.0;
constexpr double simulationTime = 1.0;
constexpr double loadFactor = 13.37;

void AssignSection(NuTo::Structure& structure);
void AssignMaterial(NuTo::Structure& structure);

std::map<int, VectorXd> GetNodeToDisplacementsMap(NuTo::Structure& structure);

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

    const auto& groupNodesAtBottomLeft = structure.GroupGetNodeRadiusRange(coordinateAtBottomLeft);
    const auto& groupNodeAtBottomRight = structure.GroupGetNodeRadiusRange(coordinateAtBottomRight);

    const auto groupNodesAtBoundaries = Group<NodeBase>::Unite(groupNodeAtBottomRight, groupNodesAtBottomLeft);

    const auto& groupNodesLoad = structure.GroupGetNodeRadiusRange(coordinateAtLoad, 0, 1.e-6);

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

    structure.AddVisualizationComponent(structure.GroupGetElementsTotal(), eVisualizeWhat::DISPLACEMENTS);

    structure.GetLogger() << "*********************************** \n"
                          << "**      integration scheme       ** \n"
                          << "*********************************** \n\n";

    NuTo::NewmarkFeti<EigenSolver> newmark(&structure);
    boost::filesystem::path resultPath(boost::filesystem::initial_path().string() + "/FetiSolutionRank_" +
                                       std::to_string(structure.mRank));

    newmark.SetTimeStep(timeStep);
    newmark.SetResultDirectory(resultPath.string(), true);
    newmark.SetToleranceIterativeSolver(1.e-7);

    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    newmark.SetTimeDependentLoadCase(loadId, dispRHS);

    structure.GetLogger() << "*********************************** \n"
                          << "**      solve                    ** \n"
                          << "*********************************** \n\n";

    newmark.Solve(simulationTime);

    std::map<int, VectorXd> nodeIdAndDisplacements = GetNodeToDisplacementsMap(structure);
    std::map<int, VectorXd> nodeIdAndDisplacementsReference = ComputeReferenceSolution(rank);
    
    structure.GetLogger() << "*********************************** \n"
                          << "**      compare results          ** \n"
                          << "*********************************** \n\n";

    for (const auto& nodeIdAndDisplacementsPair : nodeIdAndDisplacements)
    {
        VectorXd vec00 = nodeIdAndDisplacementsPair.second;
        VectorXd vec01 = nodeIdAndDisplacementsReference[nodeIdAndDisplacementsPair.first - 1];

        if (not(vec00 - vec01).isMuchSmallerThan(1.e-4, 1.))
        {
            cout << "nodeId: \n" << nodeIdAndDisplacementsPair.first << endl;
            cout << std::setprecision(6) << endl;
            cout << "vec00 \n" << vec00 << endl;

            cout << "vec01 \n" << vec01 << endl;

            cout << "Solutions are not equal" << endl;

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

    auto section = NuTo::SectionPlane::Create(thickness, false);
    structure.ElementTotalSetSection(section);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AssignMaterial(NuTo::Structure& structure)
{
    structure.GetLogger() << "*********************************** \n"
                          << "**      material                 ** \n"
                          << "*********************************** \n\n";

    const int materialId = structure.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);

    structure.ElementTotalSetConstitutiveLaw(materialId);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::map<int, VectorXd> GetNodeToDisplacementsMap(NuTo::Structure& structure)
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
    structure.GetLogger().SetQuiet(true);

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

    auto& groupNodesLeftBoundary = structure.GroupGetNodeRadiusRange(coordinateAtBottomLeft);
    structure.Constraints().Add(eDof::DISPLACEMENTS, Component(groupNodesLeftBoundary, {eDirection::X, eDirection::Y}));

    auto& groupNodesRightBoundary = structure.GroupGetNodeRadiusRange(coordinateAtBottomRight);
    structure.Constraints().Add(eDof::DISPLACEMENTS,
                                Component(groupNodesRightBoundary, {eDirection::X, eDirection::Y}));

    structure.GetLogger() << "*********************************** \n"
                          << "**      load                     ** \n"
                          << "*********************************** \n\n";

    auto& loadNodeGroup = structure.GroupGetNodeRadiusRange(coordinateAtLoad, 0, 1.e-6);
    structure.Constraints().Add(eDof::DISPLACEMENTS, Component(loadNodeGroup, {eDirection::Y},
                                                               Constraint::RhsRamp(simulationTime, loadFactor)));

    structure.GetLogger() << "*********************************** \n"
                          << "**      intergration scheme      ** \n"
                          << "*********************************** \n\n";

    boost::filesystem::path resultPath(boost::filesystem::initial_path().string() + "/FetiReferenceSolutionRank_" +
                                       std::to_string(rank));

    NuTo::NewmarkDirect newmarkDirect(&structure);
    newmarkDirect.SetResultDirectory(resultPath.string(), true);
    newmarkDirect.SetTimeStep(timeStep);

    structure.GetLogger() << "*********************************** \n"
                          << "**      solve                    ** \n"
                          << "*********************************** \n\n";

    newmarkDirect.Solve(simulationTime);

    structure.GetLogger() << "*********************************** \n"
                          << "**      post process             ** \n"
                          << "*********************************** \n\n";

    return GetNodeToDisplacementsMap(structure);
}
