
#include <mpi.h>
#include <boost/mpi.hpp>
#include "boost/filesystem.hpp"

#include <chrono>
#include <mechanics/feti/FetiLumpedPreconditioner.h>
#include <mechanics/feti/FetiDirichletPreconditioner.h>

#include "mechanics/feti/NewmarkFeti.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/groups/Group.h"

#include "visualize/VisualizeEnum.h"

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

// geometry
constexpr int dim = 2;
constexpr double thickness = 1.0;
constexpr double lengthX = 40.;

// material
constexpr double youngsModulus = 2.1e5;
constexpr double poissonsRatio = 0.0;

// integration
constexpr double timeStep = 1.;
constexpr double simulationTime = 1.;
constexpr double loadFactor = 10.0;

// auxiliary
const Vector2d directionX = Vector2d::UnitX();

int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;

    const int rank = world.rank();

    NuTo::StructureFeti structure(dim);
    structure.SetNumTimeDerivatives(0);
    structure.GetLogger().OpenFile("output" + std::to_string(rank));
    structure.GetLogger().SetQuiet(true);

    std::string meshFile = "FetiMultiplicityScaling.mesh" + std::to_string(rank);

    const int interpolationTypeId = structure.InterpolationTypeCreate(eShapeType::QUAD2D);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES, eTypeOrder::EQUIDISTANT1);

    structure.ImportMeshJson(meshFile, interpolationTypeId);

    structure.GetLogger() << "*********************************** \n"
                          << "**      material                 ** \n"
                          << "*********************************** \n\n";

    const int materialId = structure.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);
    structure.ElementTotalSetConstitutiveLaw(materialId);

    structure.GetLogger() << "*********************************** \n"
                          << "**      section                  ** \n"
                          << "*********************************** \n\n";

    auto section = NuTo::SectionPlane::Create(thickness, true);
    structure.ElementTotalSetSection(section);

    structure.GetLogger() << "*********************************** \n"
                          << "**      node groups              ** \n"
                          << "*********************************** \n\n";

    const auto& groupNodesLeftBoundary = structure.GroupGetNodesAtCoordinate(eDirection::X, 0.);
    const auto& groupNodesLoad = structure.GroupGetNodesAtCoordinate(eDirection::X, lengthX);

    structure.GetLogger() << "*********************************** \n"
                          << "**      virtual constraints      ** \n"
                          << "*********************************** \n\n";

    structure.ApplyVirtualConstraints(groupNodesLeftBoundary.GetMemberIds(), groupNodesLoad.GetMemberIds());

    structure.GetLogger() << "*********************************** \n"
                          << "**      real constraints         ** \n"
                          << "*********************************** \n\n";

    structure.NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

    structure.ApplyConstraintsTotalFeti(groupNodesLeftBoundary);

    structure.GetLogger() << "*********************************** \n"
                          << "**      load                     ** \n"
                          << "*********************************** \n\n";

    // prescribe displacement of groupNodesLoad in Y direction
    std::map<int, double> dofIdAndPrescribedDisplacementMap;
    for (auto const& nodeId : groupNodesLoad.GetMemberIds())
    {
        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);
        dofIdAndPrescribedDisplacementMap.emplace(dofIds[0], 1.);
    }

    structure.ApplyPrescribedDisplacements(dofIdAndPrescribedDisplacementMap);

    int loadId = structure.LoadCreateNodeGroupForce(&groupNodesLoad, directionX, 0.);


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


    NuTo::NewmarkFeti<EigenSolver> newmarkFeti(&structure);

    boost::filesystem::path resultPath(boost::filesystem::initial_path().string() + "/feti_" + std::to_string(rank));

    newmarkFeti.SetTimeStep(timeStep);
    newmarkFeti.SetResultDirectory(resultPath.string(), true);
    newmarkFeti.SetToleranceIterativeSolver(1.e-8);
    newmarkFeti.SetMaxNumberOfFetiIterations(100);
    newmarkFeti.SetFetiPreconditioner(std::make_unique<NuTo::FetiLumpedPreconditioner>());
    newmarkFeti.SetIterativeSolver(FetiIterativeSolver::ConjugateGradient);
    newmarkFeti.SetFetiScaling(FetiScaling::None);


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

    std::ifstream file(newmarkFeti.GetResultDirectory() + "/FetiSolverInfo.txt");
    std::cout << newmarkFeti.GetResultDirectory() << std::endl;

    int numIterations = 0;
    file >> numIterations >> numIterations;
    file.close();

    std::cout << "#iterations: \t" << numIterations << std::endl;
}
