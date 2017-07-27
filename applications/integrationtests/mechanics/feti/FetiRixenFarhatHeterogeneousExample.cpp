
#include <mpi.h>
#include <boost/mpi.hpp>
#include "boost/filesystem.hpp"

#include <ctime>
#include <chrono>

#include "mechanics/feti/NewmarkFeti.h"
#include "mechanics/feti/StructureFeti.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/groups/Group.h"

#include "visualize/VisualizeEnum.h"

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
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Vector2d;
using Eigen::Matrix2d;
using EigenSolver = Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>;
using FetiIterativeSolver = NuTo::NewmarkFeti<EigenSolver>::eIterativeSolver;
using FetiScaling = NuTo::NewmarkFeti<EigenSolver>::eFetiScaling;

using EigenSolver = Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>;

// geometry
constexpr int dim = 2;
constexpr double lengthX = 20.;
constexpr double thickness = 1.0;

// material
double youngsModulus = 4094;
constexpr double poissonsRatio = 0.3;

// integration
constexpr bool performLineSearch = false;
constexpr bool automaticTimeStepping = false;
constexpr double timeStep = 1e-0;
constexpr double minTimeStep = 1e-1;
constexpr double maxTimeStep = 1e-0;
constexpr double toleranceDisp = 1e-5;
constexpr double simulationTime = 1.0;
constexpr double loadFactor = 1.0;
constexpr double maxIterations = 1;
constexpr double tol = 1.e-6;

int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;
    const int rank = world.rank();

    NuTo::StructureFeti structure(dim);
    structure.SetNumTimeDerivatives(0);
    structure.SetVerboseLevel(5);
    structure.SetShowTime(false);
    structure.GetLogger().OpenFile("FetiRixenFarhatHeterogeneousExampleLogFile_" + std::to_string(rank));
    structure.GetLogger().SetQuiet(true);

    std::string meshFile = "FetiRixenFarhatHeterogeneousExample.mesh" + std::to_string(rank);

    const int interpolationTypeId = structure.InterpolationTypeCreate(eShapeType::QUAD2D);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES, eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT1);

    structure.ImportMeshJson(meshFile, interpolationTypeId);

    if (rank == 8 or rank == 9 or rank == 10 or rank == 11)
    {
        youngsModulus = 1;
    }


    const int materialId = structure.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);
    structure.ElementTotalSetConstitutiveLaw(materialId);


    auto section = NuTo::SectionPlane::Create(thickness, true);
    structure.ElementTotalSetSection(section);

    structure.ExportVtkDataFileElements("elements_" + std::to_string(rank) + ".vtk");

    structure.GetLogger() << "*********************************** \n"
                          << "**      create node groups       ** \n"
                          << "*********************************** \n\n";

    const auto& groupNodesLeftBoundary = structure.GroupGetNodesAtCoordinate(eDirection::X, 0.);
    const auto& groupNodesLoad = structure.GroupGetNodesAtCoordinate(eDirection::X, lengthX);


    structure.GetLogger() << "*********************************** \n"
                          << "**      virtual constraints      ** \n"
                          << "*********************************** \n\n";


    std::vector<int> nodeIdsBoundaries = groupNodesLeftBoundary.GetMemberIds();
    std::vector<int> nodeIdsLoads = groupNodesLoad.GetMemberIds();


    structure.ApplyVirtualConstraints(nodeIdsBoundaries, nodeIdsLoads);

    structure.GetLogger() << "*********************************** \n"
                          << "**      real constraints         ** \n"
                          << "*********************************** \n\n";


    structure.NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

    structure.NodeInfo(10);

    std::vector<int> boundaryDofIds;
    for (const int nodeId : nodeIdsBoundaries)
    {
        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);
        boundaryDofIds.push_back(dofIds[0]);
        boundaryDofIds.push_back(dofIds[1]);
    }

    structure.ApplyConstraintsTotalFeti(boundaryDofIds);


    structure.GetLogger() << "*********************************** \n"
                          << "**      load                     ** \n"
                          << "*********************************** \n\n";


    std::map<int, double> dofIdAndPrescribedDisplacementMap;

    for (auto const& nodeId : nodeIdsLoads)
    {
        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);
        dofIdAndPrescribedDisplacementMap.emplace(dofIds[0], 1.);
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
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::ENGINEERING_STRESS);


    structure.GetLogger() << "*********************************** \n"
                          << "**      integration scheme       ** \n"
                          << "*********************************** \n\n";

    NuTo::NewmarkFeti<EigenSolver> newmarkFeti(&structure);
    boost::filesystem::path resultPath(boost::filesystem::initial_path().string() +
                                       "/FetiRixenFarhatHeterogeneousExampleResultDir_" + std::to_string(rank));

    newmarkFeti.SetTimeStep(timeStep);

    newmarkFeti.SetMinTimeStep(minTimeStep);
    newmarkFeti.SetMaxTimeStep(maxTimeStep);
    newmarkFeti.SetMaxNumIterations(maxIterations);
    newmarkFeti.SetAutomaticTimeStepping(automaticTimeStepping);
    newmarkFeti.PostProcessing().SetResultDirectory(resultPath.string(), true);
    newmarkFeti.SetPerformLineSearch(performLineSearch);
    newmarkFeti.SetToleranceResidual(eDof::DISPLACEMENTS, toleranceDisp);
    newmarkFeti.SetIterativeSolver(NuTo::NewmarkFeti<EigenSolver>::eIterativeSolver::ConjugateGradient);
    newmarkFeti.SetMaxNumberOfFetiIterations(1000);
    newmarkFeti.SetToleranceIterativeSolver(1.e-6);

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
    structure.GetLogger() << "Total number of Dofs: \t" << structure.GetNumTotalDofs() << "\n\n";
}
