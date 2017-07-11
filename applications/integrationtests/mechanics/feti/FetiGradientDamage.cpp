
#include <mpi.h>
#include <boost/mpi.hpp>
#include "mechanics/constitutive/damageLaws/DamageLawExponential.h"
#include "mechanics/feti/StructureFeti.h"
#include <ctime>
#include <chrono>
#include "mechanics/feti/NewmarkFeti.h"
#include "visualize/VisualizeEnum.h"
#include "boost/filesystem.hpp"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/groups/Group.h"


using std::cout;
using std::endl;
using namespace NuTo;
using namespace Constitutive;
using namespace Interpolation;
using Node::eDof;
using Eigen::VectorXd;
using Eigen::Vector2d;
using Eigen::Matrix2d;
using EigenSolver = Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>;

// geometry
constexpr int dim = 2;
constexpr double thickness = 1.0;
constexpr double lengthX = 60.;
constexpr double lengthY = 10.;
const Vector2d coordinateAtBottomLeft(0., 0.);
const Vector2d coordinateAtBottomRight(lengthX, 0.);
const Vector2d coordinateAtLoad(0.5 * lengthX, lengthY);

// material
constexpr double nonlocalRadius = 4; // mm
constexpr double youngsModulus = 4.0e4;
constexpr double poissonsRatio = 0.2;
constexpr double tensileStrength = 3;
constexpr double compressiveStrength = 30;
constexpr double fractureEnergy = 0.01;
constexpr double alpha = 0.99;

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

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    NuTo::StructureFeti structure(dim);
    structure.SetNumTimeDerivatives(0);
    structure.SetVerboseLevel(5);
    structure.SetShowTime(true);
    structure.GetLogger().OpenFile("FetiGradientDamageLogFile_" + std::to_string(rank));
    structure.GetLogger().SetQuiet(true);

    std::vector<double> meshDimensions;
    meshDimensions.push_back(lengthX);
    meshDimensions.push_back(lengthY);

    std::vector<int> numElements;
    numElements.push_back(20);
    numElements.push_back(10);

    auto importContainer = structure.CreateRectangularMesh2D(meshDimensions, numElements);

    const int interpolationTypeId = importContainer.second;
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::NONLOCALEQSTRAIN, eTypeOrder::EQUIDISTANT1);
    structure.ElementTotalConvertToInterpolationType();


    structure.GetLogger() << "*********************************** \n"
                          << "**      material                 ** \n"
                          << "*********************************** \n\n";

    const int materialId = structure.ConstitutiveLawCreate(eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::TENSILE_STRENGTH, tensileStrength);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::COMPRESSIVE_STRENGTH,
                                                compressiveStrength);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::NONLOCAL_RADIUS, nonlocalRadius);
    structure.ConstitutiveLawSetDamageLaw(
            materialId, NuTo::Constitutive::DamageLawExponential::Create(tensileStrength / youngsModulus,
                                                                         tensileStrength / fractureEnergy, alpha));

    structure.ElementTotalSetConstitutiveLaw(materialId);

    structure.GetLogger() << "*********************************** \n"
                          << "**      section                  ** \n"
                          << "*********************************** \n\n";

    auto section = NuTo::SectionPlane::Create(thickness, true);
    structure.ElementTotalSetSection(section);

    structure.GetLogger() << "*********************************** \n"
                          << "**      node groups              ** \n"
                          << "*********************************** \n\n";

    const auto& groupNodesAtBottomLeft = structure.GroupGetNodeRadiusRange(coordinateAtBottomLeft);
    const auto& groupNodeAtBottomRight = structure.GroupGetNodeRadiusRange(coordinateAtBottomRight);

    const auto groupNodesAtBoundaries = Group<NodeBase>::Unite(groupNodeAtBottomRight, groupNodesAtBottomLeft);

    const auto& groupNodeLoad = structure.GroupGetNodeRadiusRange(coordinateAtLoad);

    structure.GetLogger() << "*********************************** \n"
                          << "**      virtual constraints      ** \n"
                          << "*********************************** \n\n";

    std::vector<int> nodeIdsBoundaries = groupNodesAtBoundaries.GetMemberIds();
    std::vector<int> nodeIdsLoads = groupNodeLoad.GetMemberIds();


    structure.ApplyVirtualConstraints(nodeIdsBoundaries, nodeIdsLoads);

    structure.GetLogger() << "*********************************** \n"
                          << "**      real constraints         ** \n"
                          << "*********************************** \n\n";

    structure.NodeBuildGlobalDofs(__PRETTY_FUNCTION__);
    structure.ApplyConstraintsTotalFeti(groupNodesAtBoundaries);

    structure.GetLogger() << "*********************************** \n"
                          << "**      load                     ** \n"
                          << "*********************************** \n\n";

    std::map<int, double> dofIdAndPrescribedDisplacementMap;
    for (auto const& nodeId : nodeIdsLoads)
    {
        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);
        dofIdAndPrescribedDisplacementMap.emplace(dofIds[1], 1.);
    }

    structure.ApplyPrescribedDisplacements(dofIdAndPrescribedDisplacementMap);

    int loadId = structure.LoadCreateNodeGroupForce(&groupNodeLoad, Vector2d::UnitY(), 0.);

    structure.GetLogger() << "*********************************** \n"
                          << "**      visualization            ** \n"
                          << "*********************************** \n\n";

    int groupAllElements = 9999;
    structure.GroupCreate(groupAllElements, eGroupId::Elements);
    structure.GroupAddElementsTotal(groupAllElements);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DAMAGE);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::NONLOCAL_EQ_STRAIN);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::ENGINEERING_STRESS);

    structure.GetLogger() << "*********************************** \n"
                          << "**      integration scheme       ** \n"
                          << "*********************************** \n\n";

    NuTo::NewmarkFeti<EigenSolver> newmarkFeti(&structure);

    boost::filesystem::path resultPath(boost::filesystem::initial_path().string() + "/FetiGradientDamageResultDir_" +
                                       std::to_string(structure.mRank));

    newmarkFeti.SetTimeStep(timeStep);
    newmarkFeti.SetMaxNumIterations(maxIterations);
    newmarkFeti.SetMinTimeStep(minTimeStep);
    newmarkFeti.SetMaxTimeStep(maxTimeStep);
    newmarkFeti.SetAutomaticTimeStepping(automaticTimeStepping);
    newmarkFeti.SetResultDirectory(resultPath.string(), true);
    newmarkFeti.SetPerformLineSearch(performLineSearch);
    newmarkFeti.SetToleranceResidual(eDof::DISPLACEMENTS, toleranceDisp);
    newmarkFeti.SetToleranceResidual(eDof::NONLOCALEQSTRAIN, toleranceNlEqStrain);
    newmarkFeti.SetToleranceIterativeSolver(1.e-4);
    newmarkFeti.SetIterativeSolver(NuTo::NewmarkFeti<EigenSolver>::eIterativeSolver::ProjectedGmres);

    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    newmarkFeti.SetTimeDependentLoadCase(loadId, dispRHS);

    structure.GetLogger() << "*********************************** \n"
                          << "**      solve                    ** \n"
                          << "*********************************** \n\n";

    newmarkFeti.Solve(simulationTime);

    MPI_Finalize();
}
