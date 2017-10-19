
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
#include "base/Exception.h"
#include "typedefs.h"

using EigenSolver = Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>;
using FetiIterativeSolver = NuTo::NewmarkFeti<EigenSolver>::eIterativeSolver;
using FetiScaling = NuTo::NewmarkFeti<EigenSolver>::eFetiScaling;

// geometry
constexpr int dim = 2;
constexpr double thickness = 1.0;
constexpr double lengthX = 40.;

// material
double youngsModulus = 1000;
constexpr double poissonsRatio = 0.0;

// integration
constexpr double timeStep = 1.;
constexpr double simulationTime = 1.;
constexpr double loadFactor = 10.0;

// auxiliary
const Vector2d directionX = Vector2d::UnitX();

void InitializeStructure(NuTo::StructureFeti& structure);
void InitializeNewmarkFeti(NuTo::NewmarkFeti<EigenSolver>& newmarkFeti);
int ReadNumIterationsFromFile(const std::string& file);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc, argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int numIterationsLumpedNoScaling = 0;
    int numIterationsLumpedSuperlumpedScaling = 0;
    int numIterationsLumpedMultiplicityScaling = 0;

    if (rank == 0)
        std::cout << "Conjugate gradient, lumped preconditioner, no scaling." << std::endl;
    {
        NuTo::StructureFeti structure(dim);
        InitializeStructure(structure);

        NuTo::NewmarkFeti<EigenSolver> newmarkFeti(&structure);
        InitializeNewmarkFeti(newmarkFeti);
        newmarkFeti.SetFetiPreconditioner(std::make_unique<NuTo::FetiLumpedPreconditioner>());
        newmarkFeti.SetFetiScaling(FetiScaling::None);

        newmarkFeti.Solve(simulationTime);

        numIterationsLumpedNoScaling =
                ReadNumIterationsFromFile(newmarkFeti.PostProcessing().GetResultDirectory() + "/FetiSolverInfo.txt");
    }

    if (rank == 0)
        std::cout << "Conjugate gradient, lumped preconditioner, multiplicity scaling." << std::endl;
    {
        NuTo::StructureFeti structure(dim);
        InitializeStructure(structure);

        NuTo::NewmarkFeti<EigenSolver> newmarkFeti(&structure);
        InitializeNewmarkFeti(newmarkFeti);
        newmarkFeti.SetFetiPreconditioner(std::make_unique<NuTo::FetiLumpedPreconditioner>());
        newmarkFeti.SetFetiScaling(FetiScaling::Multiplicity);

        newmarkFeti.Solve(simulationTime);

        numIterationsLumpedMultiplicityScaling =
                ReadNumIterationsFromFile(newmarkFeti.PostProcessing().GetResultDirectory() + "/FetiSolverInfo.txt");

    }

    if (rank == 0)
        std::cout << "Conjugate gradient, lumped preconditioner, k scaling." << std::endl;
    {
        NuTo::StructureFeti structure(dim);
        InitializeStructure(structure);

        NuTo::NewmarkFeti<EigenSolver> newmarkFeti(&structure);
        InitializeNewmarkFeti(newmarkFeti);
        newmarkFeti.SetFetiPreconditioner(std::make_unique<NuTo::FetiLumpedPreconditioner>());
        newmarkFeti.SetFetiScaling(FetiScaling::Superlumped);

        newmarkFeti.Solve(simulationTime);

        numIterationsLumpedSuperlumpedScaling =
                ReadNumIterationsFromFile(newmarkFeti.PostProcessing().GetResultDirectory() + "/FetiSolverInfo.txt");

    }

    if (numIterationsLumpedSuperlumpedScaling >= numIterationsLumpedNoScaling)
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Scaling should improve convergence. iterations with scaling: " +
                                                           std::to_string(numIterationsLumpedSuperlumpedScaling) +
                                                           " iterations withour scaling: " +
                                                           std::to_string(numIterationsLumpedNoScaling));

    int numIterationsDirichletNoScaling = 0;
    int numIterationsDirichletSuperlumpedScaling = 0;
    int numIterationsDirichletMultiplicityScaling = 0;

    if (rank == 0)
        std::cout << "Conjugate gradient, Dirichlet preconditioner, no scaling." << std::endl;
    {
        NuTo::StructureFeti structure(dim);
        InitializeStructure(structure);

        NuTo::NewmarkFeti<EigenSolver> newmarkFeti(&structure);
        InitializeNewmarkFeti(newmarkFeti);
        newmarkFeti.SetFetiPreconditioner(std::make_unique<NuTo::FetiDirichletPreconditioner>());
        newmarkFeti.SetFetiScaling(FetiScaling::None);

        newmarkFeti.Solve(simulationTime);

        numIterationsDirichletNoScaling =
                ReadNumIterationsFromFile(newmarkFeti.PostProcessing().GetResultDirectory() + "/FetiSolverInfo.txt");
    }


    if (rank == 0)
        std::cout << "Conjugate gradient, Dirichlet preconditioner, multiplicity scaling." << std::endl;
    {
        NuTo::StructureFeti structure(dim);
        InitializeStructure(structure);

        NuTo::NewmarkFeti<EigenSolver> newmarkFeti(&structure);
        InitializeNewmarkFeti(newmarkFeti);
        newmarkFeti.SetFetiPreconditioner(std::make_unique<NuTo::FetiLumpedPreconditioner>());
        newmarkFeti.SetFetiScaling(FetiScaling::Multiplicity);

        newmarkFeti.Solve(simulationTime);

        numIterationsDirichletMultiplicityScaling =
                ReadNumIterationsFromFile(newmarkFeti.PostProcessing().GetResultDirectory() + "/FetiSolverInfo.txt");

    }

    if (rank == 0)
        std::cout << "Conjugate gradient, Dirichlet preconditioner, k scaling." << std::endl;
    {
        NuTo::StructureFeti structure(dim);
        InitializeStructure(structure);

        NuTo::NewmarkFeti<EigenSolver> newmarkFeti(&structure);
        InitializeNewmarkFeti(newmarkFeti);
        newmarkFeti.SetFetiPreconditioner(std::make_unique<NuTo::FetiDirichletPreconditioner>());
        newmarkFeti.SetFetiScaling(FetiScaling::Superlumped);

        newmarkFeti.Solve(simulationTime);

        numIterationsDirichletSuperlumpedScaling =
                ReadNumIterationsFromFile(newmarkFeti.PostProcessing().GetResultDirectory() + "/FetiSolverInfo.txt");

    }

    if (numIterationsDirichletSuperlumpedScaling >= numIterationsDirichletNoScaling)
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Scaling should improve convergence. iterations with scaling: " +
                                                           std::to_string(numIterationsDirichletSuperlumpedScaling) +
                                                           " iterations withour scaling: " +
                                                           std::to_string(numIterationsDirichletNoScaling));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int ReadNumIterationsFromFile(const std::string& fileName)
{
    int numIterations = 0;
    std::ifstream file(fileName);
    std::string line;
    std::getline(file, line);
    file >> numIterations;
    file.close();
    return numIterations;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void InitializeNewmarkFeti(NuTo::NewmarkFeti<EigenSolver>& newmarkFeti)
{
    boost::mpi::communicator world;
    const int rank = world.rank();

    Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    boostFs::path resultPath(boostFs::initial_path().string() + "/feti_" + std::to_string(rank));

    newmarkFeti.SetTimeStep(timeStep);
    newmarkFeti.PostProcessing().SetResultDirectory(resultPath.string(), true);
    newmarkFeti.SetToleranceIterativeSolver(1.e-8);
    newmarkFeti.SetMaxNumberOfFetiIterations(1000);
    newmarkFeti.SetIterativeSolver(FetiIterativeSolver::ConjugateGradient);
    newmarkFeti.SetFetiScaling(FetiScaling::None);
    newmarkFeti.SetTimeDependentLoadCase(0, dispRHS);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void InitializeStructure(NuTo::StructureFeti& structure)
{
    boost::mpi::communicator world;

    const int rank = world.rank();
    structure.SetNumTimeDerivatives(0);
    structure.GetLogger().OpenFile("output" + std::to_string(rank));
    structure.GetLogger().SetQuiet(true);

    std::string meshFile = "meshes/FetiMultiplicityScaling.mesh" + std::to_string(rank);

    const int interpolationTypeId = structure.InterpolationTypeCreate(eShapeType::QUAD2D);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES, eTypeOrder::EQUIDISTANT1);

    structure.ImportMeshJson(meshFile, interpolationTypeId);

    structure.GetLogger() << "*********************************** \n"
                          << "**      material                 ** \n"
                          << "*********************************** \n\n";

    if (structure.mRank == 0)
        youngsModulus = 10.;

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

    structure.LoadCreateNodeGroupForce(&groupNodesLoad, directionX, 0.);


    structure.GetLogger() << "*********************************** \n"
                          << "**      visualization            ** \n"
                          << "*********************************** \n\n";

    int groupAllElements = 9999;
    structure.GroupCreate(groupAllElements, eGroupId::Elements);
    structure.GroupAddElementsTotal(groupAllElements);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DISPLACEMENTS);
}
