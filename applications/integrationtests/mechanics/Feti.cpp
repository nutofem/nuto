
#include <mpi.h>


#include "mechanics/feti/StructureFeti.h"

#include <ctime>
#include <chrono>
#include "mechanics/feti/NewmarkFeti.h"

#include "mechanics/nodes/NodeBase.h"

#include "boost/filesystem.hpp"

#include "mechanics/groups/GroupEnum.h"
#include "visualize/VisualizeEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/elements/IpDataEnum.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"

using NuTo::Constitutive::ePhaseFieldEnergyDecomposition;
using NuTo::Constitutive::eConstitutiveType;
using NuTo::Constitutive::eConstitutiveParameter;
using NuTo::Node::eDof;
using NuTo::Interpolation::eTypeOrder;
using NuTo::Interpolation::eShapeType;
using NuTo::eGroupId;
using NuTo::eVisualizeWhat;

using Eigen::VectorXd;
using Eigen::MatrixXd;
using EigenSolver = Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>>;

// geometry
constexpr   int         dimension                   = 2;
constexpr   double      thickness                   = 1.0;

// material
constexpr   double      youngsModulus               = 4.0e4;
constexpr   double      poissonsRatio               = 0.2;

// integration
constexpr   bool        automaticTimeStepping       = false;
constexpr   double      timeStep                    = 1.0;
constexpr   double      toleranceDisp              = 1e-6;
constexpr   double      simulationTime              = 1.0;
constexpr   double      loadFactor                  = 13.37;


const Eigen::Vector2d directionX    = Eigen::Vector2d::UnitX();
const Eigen::Vector2d directionY    = Eigen::Vector2d::UnitY();

void AssignSection(NuTo::Structure& structure);
void AssignMaterial(NuTo::Structure& structure);
void AddVisualization(NuTo::Structure& structure);
std::map<int, VectorXd> GetNodeAndDisplacementsMap(NuTo::Structure& structure);

std::map<int, VectorXd> ComputeReferenceSolution();

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    NuTo::StructureFeti structure(dimension);
    structure.SetNumTimeDerivatives(0);
    structure.SetVerboseLevel(10);
    structure.SetShowTime(false);
    structure.GetLogger().OpenFile("FetiTestOutputRank_" + std::to_string(rank));
    structure.GetLogger().SetQuiet(true);

    std::string meshFile = "feti_beam_coarse_2_subdomains_24_ele.mesh" + std::to_string(rank);
    structure.GetLogger() << meshFile << "\n";

    const int interpolationTypeId = structure.InterpolationTypeCreate(eShapeType::QUAD2D);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES,      eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS,    eTypeOrder::EQUIDISTANT1);

    structure.ImportMeshJson(meshFile,interpolationTypeId);


    AssignMaterial(structure);
    AssignSection(structure);

    structure.GetLogger() << "**********************************************" << "\n";
    structure.GetLogger() << "**  virtual constraints                     **" << "\n";
    structure.GetLogger() << "**********************************************" << "\n\n";

    Eigen::VectorXd nodeCoords(2);


    if (structure.mRank == 1)
    {
        int groupNodesFakeConstraints00 = structure.GroupCreate(eGroupId::Nodes);
        int groupNodesFakeConstraints01 = structure.GroupCreate(eGroupId::Nodes);

        nodeCoords[0] = 10;
        nodeCoords[1] = 0;
        structure.GroupAddNodeRadiusRange(groupNodesFakeConstraints00, nodeCoords, 0, 1.e-6);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesFakeConstraints00, directionX, 0);


        nodeCoords[0] = 10;
        nodeCoords[1] = 10;
        structure.GroupAddNodeRadiusRange(groupNodesFakeConstraints01, nodeCoords, 0, 1.e-6);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesFakeConstraints01, directionX, 0);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesFakeConstraints01, directionY, 0);

        structure.GetLogger() << "Number of nodes that are constraint in 1st group: \t"
                              << structure.GroupGetNumMembers(groupNodesFakeConstraints00) << "\n";

        structure.GetLogger() << "Number of nodes that are constraint in 2nd group: \t"
                              << structure.GroupGetNumMembers(groupNodesFakeConstraints01) << "\n";


    }

    if (structure.mRank == 0)
    {
        int groupNodesFakeConstraints00 = structure.GroupCreate(eGroupId::Nodes);
        int groupNodesFakeConstraints01 = structure.GroupCreate(eGroupId::Nodes);

        nodeCoords[0] = 40;
        nodeCoords[1] = 0;
        structure.GroupAddNodeRadiusRange(groupNodesFakeConstraints00, nodeCoords, 0, 1.e-6);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesFakeConstraints00, directionX, 0);


        nodeCoords[0] = 50;
        nodeCoords[1] = 10;
        structure.GroupAddNodeRadiusRange(groupNodesFakeConstraints01, nodeCoords, 0, 1.e-6);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesFakeConstraints01, directionX, 0);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesFakeConstraints01, directionY, 0);

        structure.GetLogger() << "Number of nodes that are constraint in 1st group: \t"
                              << structure.GroupGetNumMembers(groupNodesFakeConstraints00) << "\n";

        structure.GetLogger() << "Number of nodes that are constraint in 2nd group: \t"
                              << structure.GroupGetNumMembers(groupNodesFakeConstraints01) << "\n";


    }


    structure.NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

    structure.GetLogger() << "**********************************************" << "\n";
    structure.GetLogger() << "**  real constraints                        **" << "\n";
    structure.GetLogger() << "**********************************************" << "\n\n";


    if (rank == 1)
    {
        nodeCoords[0] = 0;
        nodeCoords[1] = 0;

        int groupNodesLeftBoundary = structure.GroupCreate(eGroupId::Nodes);

        structure.GroupAddNodeRadiusRange(groupNodesLeftBoundary, nodeCoords, 0, 1.e-6);
        structure.ApplyConstraintsTotalFeti(groupNodesLeftBoundary);
    }


    if (rank == 0)
    {

        nodeCoords[0] = 60;
        nodeCoords[1] = 0;

        int groupNodesRightBoundary = structure.GroupCreate(eGroupId::Nodes);

        structure.GroupAddNodeRadiusRange(groupNodesRightBoundary, nodeCoords, 0, 1.e-6);
        structure.ApplyConstraintsTotalFeti(groupNodesRightBoundary);


    }

    structure.GetLogger() << "**********************************************" << "\n";
    structure.GetLogger() << "**  load                                    **" << "\n";
    structure.GetLogger() << "**********************************************" << "\n\n";

    int loadNodeGroup = structure.GroupCreate(eGroupId::Nodes);

    nodeCoords[0] = 20;
    nodeCoords[1] = 10;
    structure.GroupAddNodeRadiusRange(loadNodeGroup, nodeCoords, 0, 1.e-6);


    // prescribe displacement of loadNodeGroup in Y direction
    std::map<int, double> dofIdAndPrescribedDisplacementMap;
    std::vector<int> nodeIds = structure.GroupGetMemberIds(loadNodeGroup);
    for (auto const& nodeId : nodeIds)
    {
        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);
        dofIdAndPrescribedDisplacementMap.emplace(dofIds[1], 1.);
    }

    structure.ApplyPrescribedDisplacements(dofIdAndPrescribedDisplacementMap);

    int loadId = structure.LoadCreateNodeGroupForce(0, loadNodeGroup, directionY, 0);


    structure.GetLogger() << "***********************************" << "\n";
    structure.GetLogger() << "**      Visualization            **" << "\n";
    structure.GetLogger() << "***********************************" << "\n\n";

    int groupAllElements = 9999;
    structure.GroupCreate(groupAllElements, eGroupId::Elements);
    structure.GroupAddElementsTotal(groupAllElements);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DISPLACEMENTS);

    structure.GetLogger() << "**********************************************" << "\n";
    structure.GetLogger() << "**  integration sheme                       **" << "\n";
    structure.GetLogger() << "**********************************************" << "\n\n";


    NuTo::NewmarkFeti<EigenSolver> myIntegrationScheme(&structure);
    boost::filesystem::path resultPath(boost::filesystem::initial_path().string() + "/FetiSolutionRank_" + std::to_string(structure.mRank));

    myIntegrationScheme.SetTimeStep                 ( timeStep                  );
    myIntegrationScheme.SetAutomaticTimeStepping    ( automaticTimeStepping     );
    myIntegrationScheme.SetResultDirectory          ( resultPath.string(), true );
    myIntegrationScheme.SetToleranceResidual        ( eDof::DISPLACEMENTS,      toleranceDisp );

    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    myIntegrationScheme.SetTimeDependentLoadCase(loadId, dispRHS);

    structure.GetLogger() << "***********************************" << "\n";
    structure.GetLogger() << "**      Solve                    **" << "\n";
    structure.GetLogger() << "***********************************" << "\n\n";

    myIntegrationScheme.Solve(simulationTime);

    std::map<int, VectorXd> nodeIdAndDisplacements          = GetNodeAndDisplacementsMap(structure);
    std::map<int, VectorXd> nodeIdAndDisplacementsReference = ComputeReferenceSolution();


    for (const auto& nodeIdAndDisplacementsPair : nodeIdAndDisplacements)
    {
        VectorXd vec00 = nodeIdAndDisplacementsPair.second;
        VectorXd vec01 = nodeIdAndDisplacementsReference[nodeIdAndDisplacementsPair.first - 1];

        if (not (vec00 - vec01).isMuchSmallerThan(1.e-4,1.e-1))
        {
            std::cout << "nodeIdAndDisplacementsPair.first \n" << nodeIdAndDisplacementsPair.first << std::endl;
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
    structure.GetLogger() << "***********************************" << "\n";
    structure.GetLogger() << "**      Section                  **" << "\n";
    structure.GetLogger() << "***********************************" << "\n\n";

    auto section00 = NuTo::SectionPlane::Create(thickness, false);

    structure.ElementTotalSetSection(section00);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AssignMaterial(NuTo::Structure& structure)
{
    structure.GetLogger() << "***********************************" << "\n";
    structure.GetLogger() << "**      Material                 **" << "\n";
    structure.GetLogger() << "***********************************" << "\n\n";

    int material00 = structure.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);

    structure.ElementTotalSetConstitutiveLaw(material00);

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AddVisualization(NuTo::Structure& structure)
{
    structure.GetLogger() << "***********************************" << "\n";
    structure.GetLogger() << "**      Visualization            **" << "\n";
    structure.GetLogger() << "***********************************" << "\n\n";

    int groupAllElements = 9999;
    structure.GroupCreate(groupAllElements, eGroupId::Elements);
    structure.GroupAddElementsTotal(groupAllElements);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DISPLACEMENTS);
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

std::map<int, VectorXd> ComputeReferenceSolution()
{

    NuTo::Structure structure(dimension);
    structure.SetShowTime(false);
    structure.GetLogger().OpenFile("FetiTestOutputReferenceSolution");

    std::string meshFile = "feti_beam_coarse_2_subdomains_24_ele_compare.msh";
    auto eleGroupAndInterpolationTypeList = structure.ImportFromGmsh(meshFile);

    const int interpolationTypeId = eleGroupAndInterpolationTypeList[0].second;
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS,    eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeSetIntegrationType(interpolationTypeId,NuTo::eIntegrationType::IntegrationType2D4NGauss4Ip);

    structure.ElementTotalConvertToInterpolationType();

    AssignSection(structure);
    AssignMaterial(structure);

    structure.NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

    structure.GetLogger() << "**********************************************" << "\n";
    structure.GetLogger() << "**  real constraints                        **" << "\n";
    structure.GetLogger() << "**********************************************" << "\n\n";

    Eigen::VectorXd nodeCoords(2);
    nodeCoords[0] = 0;
    nodeCoords[1] = 0;

    int groupNodesLeftBoundary = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeRadiusRange(groupNodesLeftBoundary, nodeCoords, 0, 1.e-6);

    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, directionX, 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, directionY, 0.0);

    nodeCoords[0] = 60;
    nodeCoords[1] = 0;

    int groupNodesRightBoundary = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeRadiusRange(groupNodesRightBoundary, nodeCoords, 0, 1.e-6);

    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary, directionX, 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary, directionY, 0.0);

    structure.GetLogger() << "**********************************************" << "\n";
    structure.GetLogger() << "**  load                                    **" << "\n";
    structure.GetLogger() << "**********************************************" << "\n\n";

    int loadNodeGroup = structure.GroupCreate(eGroupId::Nodes);
    nodeCoords[0] = 20;
    nodeCoords[1] = 10;
    structure.GroupAddNodeRadiusRange(loadNodeGroup, nodeCoords, 0, 1.e-6);
    int loadId = structure.ConstraintLinearSetDisplacementNodeGroup(loadNodeGroup, directionY, 1);

//    AddVisualization(structure);

    structure.GetLogger() << "**********************************************" << "\n";
    structure.GetLogger() << "**  integration sheme                       **" << "\n";
    structure.GetLogger() << "**********************************************" << "\n\n";

    NuTo::NewmarkDirect myIntegrationScheme(&structure);
    boost::filesystem::path resultPath(boost::filesystem::initial_path().string() + "/FetiReferenceSolution");

    std::cout << "resultPath" << resultPath.string() << std::endl;
    myIntegrationScheme.SetTimeStep                 ( timeStep                  );
    myIntegrationScheme.SetAutomaticTimeStepping    ( automaticTimeStepping     );
    myIntegrationScheme.SetResultDirectory          ( resultPath.string(), true );
    myIntegrationScheme.SetToleranceResidual        ( eDof::DISPLACEMENTS, toleranceDisp );

    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    myIntegrationScheme.AddTimeDependentConstraint(loadId, dispRHS);

    structure.GetLogger() << "***********************************" << "\n";
    structure.GetLogger() << "**      Solve                    **" << "\n";
    structure.GetLogger() << "***********************************" << "\n\n";

    myIntegrationScheme.Solve(simulationTime);

    structure.GetLogger() << "***********************************" << "\n";
    structure.GetLogger() << "**      Post process             **" << "\n";
    structure.GetLogger() << "***********************************" << "\n\n";

    return GetNodeAndDisplacementsMap(structure);
}
