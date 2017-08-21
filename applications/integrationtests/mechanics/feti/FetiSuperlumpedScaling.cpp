
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

    int numIterationsLumpedNoScaling = 0;
    int numIterationsLumpedSuperlumpedScaling = 0;

    // Conjugate gradient, lumped preconditioner, no scaling
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

    // Conjugate gradient, lumped preconditioner, Multiplicity scaling
    {
        NuTo::StructureFeti structure(dim);
        InitializeStructure(structure);

        if (structure.mRank == 0)
            std::cout << "Multiplicity-scaling *************************** \n";

        NuTo::NewmarkFeti<EigenSolver> newmarkFeti(&structure);
        InitializeNewmarkFeti(newmarkFeti);
        newmarkFeti.SetFetiPreconditioner(std::make_unique<NuTo::FetiLumpedPreconditioner>());
        newmarkFeti.SetFetiScaling(FetiScaling::Multiplicity);

        newmarkFeti.Solve(simulationTime);

        numIterationsLumpedSuperlumpedScaling =
                ReadNumIterationsFromFile(newmarkFeti.PostProcessing().GetResultDirectory() + "/FetiSolverInfo.txt");

        if (structure.mRank == 0)
            std::cout << "Multiplicity-scaling *************************** \n";
    }

    // Conjugate gradient, lumped preconditioner, superlumped scaling
    {
        NuTo::StructureFeti structure(dim);
        InitializeStructure(structure);

        if (structure.mRank == 0)
            std::cout << "K-scaling *************************** \n";

        NuTo::NewmarkFeti<EigenSolver> newmarkFeti(&structure);
        InitializeNewmarkFeti(newmarkFeti);
        newmarkFeti.SetFetiPreconditioner(std::make_unique<NuTo::FetiLumpedPreconditioner>());
        newmarkFeti.SetFetiScaling(FetiScaling::Superlumped);

        newmarkFeti.Solve(simulationTime);

        numIterationsLumpedSuperlumpedScaling =
                ReadNumIterationsFromFile(newmarkFeti.PostProcessing().GetResultDirectory() + "/FetiSolverInfo.txt");

        if (structure.mRank == 0)
            std::cout << "K-scaling *************************** \n";
    }

    assert(numIterationsLumpedSuperlumpedScaling < numIterationsLumpedNoScaling and "Scaling should improve convergence");

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

    std::string meshFile = "FetiMultiplicityScaling.mesh" + std::to_string(rank);

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
