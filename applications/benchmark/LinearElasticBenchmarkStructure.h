#pragma once
#include "mechanics/mesh/MeshGenerator.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/MechanicsEnums.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/timeIntegration/postProcessing/PostProcessor.h"
#include "visualize/VisualizeEnum.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/groups/Group.h"

namespace NuTo
{
namespace Benchmark
{

class LinearElasticBenchmarkStructure
{
public:
    LinearElasticBenchmarkStructure(std::vector<int> numElements, int numProc = 1)
        : mS(3)
    {
        mS.SetNumProcessors(numProc);
        mS.SetShowTime(false);
        mS.SetNumTimeDerivatives(2);
        SetupMesh(numElements);
        SetupLaw();
    }

    void SetupVisualization()
    {
        int visualizationGroup = mS.GroupGetElementsTotal();
        mS.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::DISPLACEMENTS);
        mS.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::ENGINEERING_STRAIN);
        mS.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::ENGINEERING_STRESS);
        mS.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
    }

    void SolveWithNewmark(int numSolves, std::string resultDir)
    {
        SetupBCs();
        NewmarkDirect newmark(&mS);
        Eigen::Matrix<double, 2, 2> loadFactor;
        loadFactor << 0, 0, 100.0, 1.0;
        newmark.SetTimeDependentLoadCase(0, loadFactor);
        newmark.PostProcessing().SetResultDirectory(resultDir, true);
        double timeStep = 3.1415;
        newmark.SetTimeStep(timeStep);
        newmark.Solve(timeStep * numSolves);
    }

    Structure& GetStructure()
    {
        return mS;
    }

    void SetupBCs()
    {
        auto& bottomNodes = mS.GroupGetNodesAtCoordinate(NuTo::eDirection::Z, 0);
        mS.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
                             NuTo::Constraint::Component(
                                     bottomNodes, {NuTo::eDirection::X, NuTo::eDirection::Y, NuTo::eDirection::Z}));

        auto& topNodes = mS.GroupGetNodesAtCoordinate(NuTo::eDirection::Z, lz);
        int topNodesId = mS.GroupGetId(&topNodes);
        int topElements = mS.GroupCreate(eGroupId::Elements);
        mS.GroupAddElementsFromNodes(topElements, topNodesId, false);
        mS.LoadSurfacePressureCreate3D(topElements, topNodesId, 10.0);
    }

private:
    void SetupMesh(std::vector<int> numElements)
    {
        auto meshInfo = MeshGenerator::Grid(mS, {lx, ly, lz}, numElements);
        mS.InterpolationTypeAdd(meshInfo.second, Node::eDof::DISPLACEMENTS, Interpolation::eTypeOrder::EQUIDISTANT2);
        mS.ElementTotalConvertToInterpolationType();
    }

    void SetupLaw()
    {
        int lawId = mS.ConstitutiveLawCreate(Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
        mS.ConstitutiveLawSetParameterDouble(lawId, Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 1.0);
        mS.ElementTotalSetConstitutiveLaw(lawId);
    }

    static constexpr double lx = 25;
    static constexpr double ly = 35;
    static constexpr double lz = 45;

    Structure mS;
};


} // Benchmark
} // NuTo
