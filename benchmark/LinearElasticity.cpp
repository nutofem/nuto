#include "Benchmark.h"
#include "mechanics/tools/MeshGenerator.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/MechanicsEnums.h"
#include "math/FullMatrix.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "visualize/VisualizeEnum.h"

using namespace NuTo;

void SolveAMediumSizedProblem()
{
    Structure structure(3);
    structure.SetNumTimeDerivatives(2);
    int section = structure.SectionCreate(eSectionType::VOLUME);
    int constitutiveLaw = structure.ConstitutiveLawCreate(Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(constitutiveLaw, Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 1.0);
    int interpolationType = structure.InterpolationTypeCreate(Interpolation::eShapeType::BRICK3D);
    structure.InterpolationTypeAdd(interpolationType, Node::eDof::COORDINATES, Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationType, Node::eDof::DISPLACEMENTS, Interpolation::eTypeOrder::EQUIDISTANT1);

    std::array<int, 3> numElements {10, 10, 100}; // i.e. 10k Elements
    std::array<double, 3> length {1.0, 1.0, 10.0};
    MeshGenerator::MeshCuboid(structure, section, constitutiveLaw, interpolationType, numElements, length);

    // why, oh why...
    structure.ElementTotalConvertToInterpolationType();

    int bottomNodes = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(bottomNodes, 2, 0.0, 0.0);

    int topNodes = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(topNodes, 2, 10.0, 10.0);
    std::cout << structure.GroupGetNumMembers(topNodes) << std::endl;

    int topElements = structure.GroupCreate(eGroupId::Elements);
    structure.GroupAddElementsFromNodes(topElements, topNodes, false);
    std::cout << structure.GroupGetNumMembers(topElements) << std::endl;

    structure.ConstraintLinearSetDisplacementNodeGroup(bottomNodes, Eigen::Vector3d::UnitX(), 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(bottomNodes, Eigen::Vector3d::UnitY(), 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(bottomNodes, Eigen::Vector3d::UnitZ(), 0.0);

    structure.SetNumLoadCases(1);
    structure.LoadSurfacePressureCreate3D(0, topElements, topNodes, 10.0);

    int visualizationGroup = structure.GroupCreate(NuTo::eGroupId::Elements);
    structure.GroupAddElementsTotal(visualizationGroup);

    structure.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::ENGINEERING_STRESS);
    structure.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);

    NewmarkDirect newmark(&structure);
    Eigen::Matrix<double, 2, 2> loadFactor;
    loadFactor << 0, 0, 100.0, 1.0;
    newmark.SetTimeDependentLoadCase(0, loadFactor);
    newmark.SetResultDirectory("LinearElasticityResults", true);
    newmark.SetTimeStep(10.0);
    newmark.Solve(100.0);
}

BENCHMARK(LinearElasticity, CompleteRun, runner)
{
    while(runner.KeepRunningIterations(1))
    {
        SolveAMediumSizedProblem();
    }
}
