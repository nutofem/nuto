
#include <mpi.h>
#include <boost/mpi.hpp>
#include "mechanics/constitutive/damageLaws/DamageLawExponential.h"
#include "mechanics/feti/StructureFeti.h"
#include <ctime>
#include <chrono>
#include "mechanics/feti/NewmarkFeti.h"
#include "mechanics/MechanicsEnums.h"
#include "visualize/VisualizeEnum.h"
#include "mechanics/nodes/NodeBase.h"

#include "boost/filesystem.hpp"
#include "mechanics/sections/SectionPlane.h"

#include "mechanics/mesh/MeshGenerator.h"


using NuTo::Constitutive::eConstitutiveType;
using NuTo::Constitutive::eConstitutiveParameter;
using NuTo::Node::eDof;
using NuTo::Interpolation::eTypeOrder;
using NuTo::Interpolation::eShapeType;
using NuTo::eGroupId;
using NuTo::eVisualizeWhat;
using NuTo::eDirection;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Vector2d;
using Eigen::Matrix2d;
using EigenSolver = Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>;


constexpr int dim = 2;
constexpr double thickness = 1.0;

// material
constexpr double nonlocalRadius = 4; // mm
constexpr double youngsModulus = 4.0e4;
constexpr double poissonsRatio = 0.2;
constexpr double tensileStrength = 3;
constexpr double compressiveStrength = 30;
constexpr double fractureEnergy = 0.01;
constexpr double alpha = 0.3;

// integration
constexpr bool performLineSearch = true;
constexpr bool automaticTimeStepping = true;
constexpr double timeStep = 1e-3;
constexpr double minTimeStep = 1e-5;
constexpr double maxTimeStep = 1e-1;

constexpr double toleranceDisp = 1e-8;
constexpr double toleranceNlEqStrain = 1e-8;
constexpr double tolerance = 1e-5;

constexpr double simulationTime = 1.0;
constexpr double loadFactor = -0.5;
constexpr double maxIterations = 10;

const Eigen::Vector2d directionX = Eigen::Vector2d::UnitX();
const Eigen::Vector2d directionY = Eigen::Vector2d::UnitY();

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
//
//    int rank = 0;
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//    NuTo::StructureFeti structure(dim);
//    structure.SetNumTimeDerivatives(0);
//    structure.SetVerboseLevel(5);
//    structure.SetShowTime(true);
//    structure.GetLogger().OpenFile("FetiGradientDamage_output_rank_" + std::to_string(rank));
//    structure.GetLogger().SetQuiet(true);
//
//    std::vector<double> meshDimensions;
//    meshDimensions.push_back(60.);
//    meshDimensions.push_back(10.);
//
//    std::vector<int> numElements;
//    numElements.push_back(20);
//    numElements.push_back(10);
//
//    auto importContainer = structure.CreateRectangularMesh2D(meshDimensions, numElements);
//
//    const int interpolationTypeId = importContainer.second;
//    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT1);
//    structure.InterpolationTypeAdd(interpolationTypeId, eDof::NONLOCALEQSTRAIN, eTypeOrder::EQUIDISTANT1);
//    structure.ElementTotalConvertToInterpolationType();
//
//
//    structure.GetLogger() << "*********************************** \n"
//                          << "**      material                 ** \n"
//                          << "*********************************** \n\n";
//
//    const int materialId = structure.ConstitutiveLawCreate(eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
//    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
//    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);
//    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::TENSILE_STRENGTH, tensileStrength);
//    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::COMPRESSIVE_STRENGTH,
//                                                compressiveStrength);
//    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::NONLOCAL_RADIUS, nonlocalRadius);
//    structure.ConstitutiveLawSetDamageLaw(
//            materialId, NuTo::Constitutive::DamageLawExponential::Create(tensileStrength / youngsModulus,
//                                                                         tensileStrength / fractureEnergy, alpha));
//
//    structure.ElementTotalSetConstitutiveLaw(materialId);
//
//
//    structure.GetLogger() << "*********************************** \n"
//                          << "**      section                  ** \n"
//                          << "*********************************** \n\n";
//
//    auto section = NuTo::SectionPlane::Create(thickness, true);
//    structure.ElementTotalSetSection(section);
//
//
//    structure.GetLogger() << "*********************************** \n"
//                          << "**      node groups              ** \n"
//                          << "*********************************** \n\n";
//
//    Eigen::VectorXd coordinates(dim);
//
//    int groupNodesLeftBoundary = structure.GroupCreate(eGroupId::Nodes);
//    structure.GroupAddNodeCoordinateRange(groupNodesLeftBoundary, 0, -1.e-6, +1.e-6);
//
//    int loadNodeGroup = structure.GroupCreate(eGroupId::Nodes);
//    coordinates[0] = 60;
//    coordinates[1] = 0;
//    structure.GroupAddNodeRadiusRange(loadNodeGroup, coordinates, 0, 1.e-6);
//
//
//    structure.GetLogger() << "*********************************** \n"
//                          << "**      virtual constraints      ** \n"
//                          << "*********************************** \n\n";
//
//    std::vector<int> nodeIdsBoundaries = structure.GroupGetMemberIds(groupNodesLeftBoundary);
//    std::vector<int> nodeIdsLoads = structure.GroupGetMemberIds(loadNodeGroup);
//
//    //    std::cout << nodeIdsBoundaries << "\n";
//
//    structure.ApplyVirtualConstraints(nodeIdsBoundaries, nodeIdsLoads);
//
//    structure.GetLogger() << "*********************************** \n"
//                          << "**      real constraints         ** \n"
//                          << "*********************************** \n\n";
//
//    structure.NodeBuildGlobalDofs(__PRETTY_FUNCTION__);
//    structure.ApplyConstraintsTotalFeti(groupNodesLeftBoundary);
//
//    structure.GetLogger() << "*********************************** \n"
//                          << "**      load                     ** \n"
//                          << "*********************************** \n\n";
//
//
//    // prescribe displacement of loadNodeGroup in Y direction
//    std::map<int, double> dofIdAndPrescribedDisplacementMap;
//    std::vector<int> nodeIds = structure.GroupGetMemberIds(loadNodeGroup);
//    for (auto const& nodeId : nodeIds)
//    {
//        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);
//        dofIdAndPrescribedDisplacementMap.emplace(dofIds[1], 1.);
//    }
//
//    structure.ApplyPrescribedDisplacements(dofIdAndPrescribedDisplacementMap);
//
//    int loadId = structure.LoadCreateNodeGroupForce(0, loadNodeGroup, directionY, 0);
//
//    structure.GetLogger() << "*********************************** \n"
//                          << "**      visualization            ** \n"
//                          << "*********************************** \n\n";
//
//    int groupAllElements = 9999;
//    structure.GroupCreate(groupAllElements, eGroupId::Elements);
//    structure.GroupAddElementsTotal(groupAllElements);
//    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DISPLACEMENTS);
//    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::NONLOCAL_EQ_STRAIN);
//    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::ENGINEERING_STRAIN);
//    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::ENGINEERING_STRESS);
//
//    structure.GetLogger() << "*********************************** \n"
//                          << "**      integration scheme       ** \n"
//                          << "*********************************** \n\n";
//
//
//    NuTo::NewmarkFeti<EigenSolver> myIntegrationScheme(&structure);
//
//    boost::filesystem::path resultPath(boost::filesystem::path(getenv("HOME")).string() +
//                                       std::string("/results/feti/") + std::to_string(structure.mRank));
//
//
//    myIntegrationScheme.SetTimeStep(timeStep);
//    myIntegrationScheme.SetMaxNumIterations(maxIterations);
//    myIntegrationScheme.SetMinTimeStep(minTimeStep);
//    myIntegrationScheme.SetMaxTimeStep(maxTimeStep);
//    myIntegrationScheme.SetAutomaticTimeStepping(automaticTimeStepping);
//    myIntegrationScheme.SetResultDirectory(resultPath.string(), true);
//    myIntegrationScheme.SetPerformLineSearch(performLineSearch);
//    myIntegrationScheme.SetToleranceResidual(eDof::DISPLACEMENTS, toleranceDisp);
//    myIntegrationScheme.SetToleranceIterativeSolver(1.e-4);
//    myIntegrationScheme.SetIterativeSolver(NuTo::NewmarkFeti<EigenSolver>::eIterativeSolver::BiconjugateGradientStabilized);
//    myIntegrationScheme.SetFetiPreconditioner(NuTo::NewmarkFeti<EigenSolver>::eFetiPreconditioner::Lumped);
//
//    Eigen::Matrix2d dispRHS;
//    dispRHS(0, 0) = 0;
//    dispRHS(1, 0) = simulationTime;
//    dispRHS(0, 1) = 0;
//    dispRHS(1, 1) = loadFactor;
//
//    myIntegrationScheme.SetTimeDependentLoadCase(loadId, dispRHS);
//
//    structure.GetLogger() << "*********************************** \n"
//                          << "**      solve                    ** \n"
//                          << "*********************************** \n\n";
//
//
//    myIntegrationScheme.Solve(simulationTime);
//
    MPI_Finalize();
}
