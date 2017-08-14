#include "base/Exception.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constraints/Constraints.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/DirectionEnum.h"
#include "mechanics/groups/Group.h"
#include "mechanics/mesh/MeshGenerator.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/sections/Section.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/sections/SectionTruss.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/timeIntegration/postProcessing/PostProcessor.h"
#include "visualize/VisualizeEnum.h"

#include <boost/filesystem.hpp>

#define TIMESTEP 2000.0
#define SIMULATIONTIME 10000.0
#define TOLERANCE 1.e-2
#define MAXITERATION 20
#define EXTERNALFORCE -1.e7;


using namespace NuTo;

template <int TDim>
void TestCreepModel()
{
    constexpr int numElementsDirection = 3;
    Structure S(TDim);
    NewmarkDirect NM(&S);

    // Mesh
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    std::pair<int, int> meshData;
    switch (TDim)
    {
    case 1:
        meshData = MeshGenerator::Grid(S, {1.0}, {numElementsDirection});
        break;
    case 2:
        meshData = MeshGenerator::Grid(S, {1.0, 1.0}, {numElementsDirection, numElementsDirection});
        break;
    case 3:
        meshData = MeshGenerator::Grid(S, {1.0, 1.0, 1.0},
                                       {numElementsDirection, numElementsDirection, numElementsDirection});
        break;
    default:
        throw Exception(__PRETTY_FUNCTION__, "Invalid dimension");
    }
    int elementGroupID = meshData.first;
    int interpolationTypeID = meshData.second;


    // Section
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch (TDim)
    {
    case 1:
        S.ElementTotalSetSection(SectionTruss::Create(1.0));
        break;
    case 2:
        S.ElementTotalSetSection(SectionPlane::Create(1.0, true));
        break;
    default:
        break;
    }


    // Constitutive law
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    int lawID = S.ConstitutiveLawCreate(Constitutive::eConstitutiveType::CREEP);

    Eigen::VectorXd kelvinChainStiffness;
    Eigen::VectorXd kelvinChainRetardationTime;
    S.ConstitutiveLawSetParameterFullVectorDouble(lawID, Constitutive::eConstitutiveParameter::KELVIN_CHAIN_STIFFNESS,
                                                  (Eigen::VectorXd(1) << 5.0e9).finished());
    S.ConstitutiveLawSetParameterFullVectorDouble(lawID,
                                                  Constitutive::eConstitutiveParameter::KELVIN_CHAIN_RETARDATIONTIME,
                                                  (Eigen::VectorXd(1) << 20000.0).finished());
    S.ConstitutiveLawSetParameterDouble(lawID, Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.0);

    //    int lawID = S.ConstitutiveLawCreate(Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    //    S.ConstitutiveLawSetParameterDouble(lawID, Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 60 * 10e9);
    //    S.ConstitutiveLawSetParameterDouble(lawID, Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.0);

    S.ElementGroupSetConstitutiveLaw(elementGroupID, lawID);
    S.InterpolationTypeAdd(interpolationTypeID, Node::eDof::DISPLACEMENTS, Interpolation::eTypeOrder::EQUIDISTANT1);
    S.ElementTotalConvertToInterpolationType();


    // Constraints
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    auto& leftNodesGroup = S.GroupGetNodeCoordinateRange(eDirection::X, 0.0, 0.0);
    assert(leftNodesGroup.GetNumMembers() > 0);

    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(leftNodesGroup, {eDirection::X}));

    // Additional 2D/3D constraints
    if (TDim > 1)
    {

        auto& nodeOrigin = S.NodeGetAtCoordinate(Eigen::VectorXd::Zero(TDim));
        S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(nodeOrigin, {eDirection::Y}));

        if (TDim > 2)
        {
            S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(nodeOrigin, {eDirection::Z}));

            auto& additionalNode = S.NodeGetAtCoordinate((Eigen::VectorXd(TDim) << 0.0, 1.0, 0.0).finished());
            S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(additionalNode, {eDirection::Z}));
        }
    }


    // Loads
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Eigen::MatrixXd timeDependentLoad(3, 2);
    timeDependentLoad(0, 0) = 0;
    timeDependentLoad(1, 0) = TIMESTEP;
    timeDependentLoad(2, 0) = SIMULATIONTIME;
    timeDependentLoad(0, 1) = 0;
    timeDependentLoad(1, 1) = EXTERNALFORCE;
    timeDependentLoad(2, 1) = EXTERNALFORCE;

    auto& rightNodesGroup = S.GroupGetNodeCoordinateRange(eDirection::X, 1.0, 1.0);
    assert(rightNodesGroup.GetNumMembers() > 0);

    Eigen::VectorXd direction = Eigen::VectorXd::Zero(TDim);
    direction[0] = 1;
    int load = S.LoadCreateNodeGroupForce(S.GroupGetId(&rightNodesGroup), direction, 1);
    NM.SetTimeDependentLoadCase(load, timeDependentLoad);

    //    if (TDim == 2)
    //    {
    //        int rightEdgeNodeGroup = S.GroupCreate(eGroupId::Nodes);
    //        S.GroupAddNodeFunction(rightEdgeNodeGroup, [](NodeBase* node) -> bool {
    //            if (node->Get(Node::eDof::COORDINATES)[0] == 1.0 &&
    //                (node->Get(Node::eDof::COORDINATES)[1] == 0.0 || node->Get(Node::eDof::COORDINATES)[1] == 1.0))
    //                return true;
    //            else
    //                return false;
    //        });
    //        assert(S.GroupGetNumMembers(rightEdgeNodeGroup) == 2);
    //        int loadEdge = S.LoadCreateNodeGroupForce(rightEdgeNodeGroup, direction, 1);
    //        NM.SetTimeDependentLoadCase(loadEdge, timeDependentLoad);
    //    }


    // Visualization
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    int visualizeGroup = S.GroupCreate(eGroupId::Elements);
    S.GroupAddElementsTotal(visualizeGroup);

    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::DISPLACEMENTS);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::ENGINEERING_STRAIN);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::ENGINEERING_STRESS);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);


    // Solve
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    std::string resultDir = "CreepModelResults";
    boost::filesystem::create_directory(resultDir);

    resultDir.append("/");
    resultDir.append(std::to_string(TDim));
    resultDir.append("D");

    NM.SetAutomaticTimeStepping(false);
    NM.SetTimeStep(TIMESTEP);
    NM.SetPerformLineSearch(false);
    NM.SetToleranceResidual(Node::eDof::DISPLACEMENTS, TOLERANCE);
    NM.SetMaxNumIterations(MAXITERATION);
    NM.PostProcessing().SetResultDirectory(resultDir, true);
    NM.Solve(SIMULATIONTIME);
}


int main(int argc, char* argv[])
{
    // TODO: Add direction option to test Y and Z as well
    // TODO: Test Poisson Ratio 0 and some other value
    TestCreepModel<1>();
    TestCreepModel<2>();
    TestCreepModel<3>();
    return 0;
}
