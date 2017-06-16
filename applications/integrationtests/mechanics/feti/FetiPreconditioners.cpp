//
// Created by phuschke on 6/16/17.
//



#include <mpi.h>
#include <boost/mpi.hpp>
#include "mechanics/constitutive/damageLaws/DamageLawExponential.h"
#include "mechanics/feti/StructureFeti.h"
#include <ctime>
#include <chrono>
#include "mechanics/feti/NewmarkFeti.h"
#include "visualize/VisualizeEnum.h"

#include "mechanics/nodes/NodeBase.h"
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
using FetiPreconditioner = NewmarkFeti<EigenSolver>::eFetiPreconditioner;
using FetiIterativeSolver = NewmarkFeti<EigenSolver>::eIterativeSolver;

// geometry
constexpr int dim = 2;
constexpr double thickness = 1.0;
constexpr double lengthX = 60.;
constexpr double lengthY = 10.;
const Vector2d coordinateAtBottomLeft(0., 0.);
const Vector2d coordinateAtBottomRight(lengthX, 0.);
const Vector2d coordinateAtLoad(0.5*lengthX, lengthY);

// material
constexpr double youngsModulus = 123456;
constexpr double poissonsRatio = 0.2;

// integration
constexpr double timeStep = 1.;
constexpr double simulationTime = 1.0;
constexpr double loadFactor = 13.37;

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    NuTo::StructureFeti structure(dim);
    structure.SetNumTimeDerivatives(0);
    structure.SetVerboseLevel(5);
    structure.SetShowTime(true);
    structure.GetLogger().OpenFile("FetiPreconditionersLogFile_" + std::to_string(rank));
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
    structure.ElementTotalConvertToInterpolationType();


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

    structure.AddVisualizationComponent(structure.GroupGetElementsTotal(), eVisualizeWhat::DISPLACEMENTS);

    structure.GetLogger() << "*********************************** \n"
                          << "**      integration scheme       ** \n"
                          << "*********************************** \n\n";

    NuTo::NewmarkFeti<EigenSolver> newmarkFeti(&structure);

    boost::filesystem::path resultPath(boost::filesystem::initial_path().string() + "/FetiPreconditionersResultDir_" +
                                       std::to_string(structure.mRank));

    newmarkFeti.SetTimeStep(timeStep);
    MPI_Barrier(MPI_COMM_WORLD);
    newmarkFeti.SetResultDirectory(resultPath.string(), true);
    MPI_Barrier(MPI_COMM_WORLD);
    newmarkFeti.SetToleranceIterativeSolver(1.e-8);
    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;
    newmarkFeti.SetTimeDependentLoadCase(loadId, dispRHS);

    newmarkFeti.SetIterativeSolver(FetiIterativeSolver::ConjugateGradient);
    newmarkFeti.SetFetiPreconditioner(FetiPreconditioner::None);
    newmarkFeti.Solve(simulationTime);

    MPI_Barrier(MPI_COMM_WORLD);
    newmarkFeti.SetResultDirectory(resultPath.string(), true);
    MPI_Barrier(MPI_COMM_WORLD);

    newmarkFeti.SetIterativeSolver(FetiIterativeSolver::ConjugateGradient);
    newmarkFeti.SetFetiPreconditioner(FetiPreconditioner::Lumped);
    newmarkFeti.Solve(2*simulationTime);

    MPI_Barrier(MPI_COMM_WORLD);
    newmarkFeti.SetResultDirectory(resultPath.string(), true);
    MPI_Barrier(MPI_COMM_WORLD);

    newmarkFeti.SetIterativeSolver(FetiIterativeSolver::ConjugateGradient);
    newmarkFeti.SetFetiPreconditioner(FetiPreconditioner::Dirichlet);
    newmarkFeti.Solve(3*simulationTime);

    MPI_Finalize();
}
