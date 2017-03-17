#pragma once
#include "mechanics/mesh/MeshGenerator.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/MechanicsEnums.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "visualize/VisualizeEnum.h"

namespace NuTo
{
namespace Benchmark
{

class LinearElasticBenchmarkStructure
{
public:
    LinearElasticBenchmarkStructure(std::vector<int> rNumElements, int rNumProc = 1) : mS(3)
    {
        mS.SetNumProcessors(rNumProc);
        mS.SetShowTime(false);
        mS.SetNumTimeDerivatives(2);
        SetupMesh(rNumElements);
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

    void SolveWithNewmark(int rNumSolves, std::string rResultDir)
    {
        SetupBCs();
        NewmarkDirect newmark(&mS);
        Eigen::Matrix<double, 2, 2> loadFactor;
        loadFactor << 0, 0, 100.0, 1.0;
        newmark.SetTimeDependentLoadCase(0, loadFactor);
        newmark.SetResultDirectory(rResultDir, true);
        double timeStep = 3.1415;
        newmark.SetTimeStep(timeStep);
        newmark.Solve(timeStep * rNumSolves);
    }

    Structure& GetStructure()
    {
        return mS;
    }

    void SetupBCs()
    {
        int bottomNodes = mS.GroupCreate(eGroupId::Nodes);
        mS.GroupAddNodeCoordinateRange(bottomNodes, 2, 0.0, 0.0);

        int topNodes = mS.GroupCreate(eGroupId::Nodes);
        mS.GroupAddNodeCoordinateRange(topNodes, 2, lz, lz);

        int topElements = mS.GroupCreate(eGroupId::Elements);
        mS.GroupAddElementsFromNodes(topElements, topNodes, false);

        mS.ConstraintLinearSetDisplacementNodeGroup(bottomNodes, Eigen::Vector3d::UnitX(), 0.0);
        mS.ConstraintLinearSetDisplacementNodeGroup(bottomNodes, Eigen::Vector3d::UnitY(), 0.0);
        mS.ConstraintLinearSetDisplacementNodeGroup(bottomNodes, Eigen::Vector3d::UnitZ(), 0.0);

        mS.SetNumLoadCases(1);
        mS.LoadSurfacePressureCreate3D(0, topElements, topNodes, 10.0);
    }

private:

    void SetupMesh(std::vector<int> rNumElements)
    {
        auto meshInfo = MeshGenerator::Grid(mS, {lx, ly, lz}, rNumElements);
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
