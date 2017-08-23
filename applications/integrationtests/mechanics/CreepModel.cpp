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
#define SIMULATIONTIME 100000.0
#define TOLERANCE 1.e-2
#define MAXITERATION 20
#define EXTERNALFORCE -1.e9;
#define TOTALYOUNGSMODULUS 2.e9;

using namespace NuTo;
using namespace NuTo::Constraint;

template <int TDim>
void TestCreepModel(std::string testName, const std::array<eDirection, TDim> directions, double YoungsModulus,
                    Eigen::VectorXd kelvinChainStiffness, Eigen::VectorXd kelvinChainRetardationTimes,
                    double poissonRatio)
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
        S.ElementTotalSetSection(SectionPlane::Create(1.0, false));
        break;
    default:
        break;
    }


    // Constitutive law
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //    int lawID = S.ConstitutiveLawCreate(Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    int lawID = S.ConstitutiveLawCreate(Constitutive::eConstitutiveType::CREEP);


    S.ConstitutiveLawSetParameterDouble(lawID, Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, YoungsModulus);
    S.ConstitutiveLawSetParameterDouble(lawID, Constitutive::eConstitutiveParameter::POISSONS_RATIO, poissonRatio);

    S.ConstitutiveLawSetParameterFullVectorDouble(lawID, Constitutive::eConstitutiveParameter::KELVIN_CHAIN_STIFFNESS,
                                                  kelvinChainStiffness);
    S.ConstitutiveLawSetParameterFullVectorDouble(
            lawID, Constitutive::eConstitutiveParameter::KELVIN_CHAIN_RETARDATIONTIME, kelvinChainRetardationTimes);

    S.ElementGroupSetConstitutiveLaw(elementGroupID, lawID);
    S.InterpolationTypeAdd(interpolationTypeID, Node::eDof::DISPLACEMENTS, Interpolation::eTypeOrder::EQUIDISTANT1);
    S.ElementTotalConvertToInterpolationType();


    // Constraints
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    auto virtualNodePtr = S.NodeGetNodePtr(S.NodeCreate(Eigen::VectorXd::Ones(TDim) * -1, {Node::eDof::DISPLACEMENTS}));

    auto& leftNodesGroup = S.GroupGetNodeCoordinateRange(directions[0], 0.0, 0.0);
    auto& rightNodesGroup = S.GroupGetNodeCoordinateRange(directions[0], 1.0, 1.0);
    assert(rightNodesGroup.GetNumMembers() > 0);
    assert(leftNodesGroup.GetNumMembers() > 0);

    S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(leftNodesGroup, {directions[0]}));
    for (auto& itNode : rightNodesGroup)
    {
        S.Constraints().Add(Node::eDof::DISPLACEMENTS,
                            Equation({Term(*virtualNodePtr, ToComponentIndex(directions[0]), 1),
                                      Term(*S.NodeGetNodePtr(itNode.first), ToComponentIndex(directions[0]), -1)}));
    }

    // Additional 2D/3D constraints
    if (TDim > 1)
    {
        auto& nodeOrigin = S.NodeGetAtCoordinate(Eigen::VectorXd::Zero(TDim));
        S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(nodeOrigin, {directions[1]}));
        S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(*virtualNodePtr, {directions[1]}));

        // Additional 3D constraints
        if (TDim > 2)
        {
            S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(nodeOrigin, {directions[2]}));
            S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(*virtualNodePtr, {directions[2]}));

            Eigen::VectorXd additionalNodeCoordinates = Eigen::VectorXd::Zero(TDim);
            additionalNodeCoordinates[ToComponentIndex(directions[1])] = 1.0;
            auto& additionalNode = S.NodeGetAtCoordinate(additionalNodeCoordinates);
            S.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(additionalNode, {directions[2]}));
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


    Eigen::VectorXd direction = Eigen::VectorXd::Zero(TDim);
    direction[ToComponentIndex(directions[0])] = 1;
    int load = S.LoadCreateNodeForce(virtualNodePtr, direction, 1);
    NM.SetTimeDependentLoadCase(load, timeDependentLoad);


    // Visualization
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    int visualizeGroup = S.GroupCreate(eGroupId::Elements);
    S.GroupAddElementsTotal(visualizeGroup);

    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::DISPLACEMENTS);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::ENGINEERING_STRAIN);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::ENGINEERING_STRESS);
    S.AddVisualizationComponent(visualizeGroup, eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);


    // Set result directory
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    std::string resultDir = "CreepModelResults";
    boost::filesystem::create_directory(resultDir);

    resultDir.append("/");
    resultDir.append(std::to_string(TDim));
    resultDir.append("D");
    boost::filesystem::create_directory(resultDir);

    resultDir.append("/");
    resultDir.append(testName);
    resultDir.append("_direction=");
    switch (directions[0])
    {
    case eDirection::X:
        resultDir.append("X");
        break;
    case eDirection::Y:
        resultDir.append("Y");
        break;
    case eDirection::Z:
        resultDir.append("Z");
        break;
    default:
        resultDir.append("UNKNOWN");
        break;
    }
    resultDir.append("_nu=");
    resultDir.append(std::to_string(poissonRatio));

    // Solve
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    NM.SetAutomaticTimeStepping(false);
    NM.SetTimeStep(TIMESTEP);
    NM.SetPerformLineSearch(false);
    NM.SetToleranceResidual(Node::eDof::DISPLACEMENTS, TOLERANCE);
    NM.SetMaxNumIterations(MAXITERATION);
    NM.PostProcessing().SetResultDirectory(resultDir, true);
    NM.Solve(SIMULATIONTIME);
}


void PerformTestSeries(std::string testName, double YoungsModulus, Eigen::VectorXd kelvinChainStiffness,
                       Eigen::VectorXd kelvinChainRetardationTimes, double poissonRatio)
{
    TestCreepModel<1>(testName, {eDirection::X}, YoungsModulus, kelvinChainStiffness, kelvinChainRetardationTimes,
                      poissonRatio);
    TestCreepModel<2>(testName, {eDirection::X, eDirection::Y}, YoungsModulus, kelvinChainStiffness,
                      kelvinChainRetardationTimes, poissonRatio);
    TestCreepModel<2>(testName, {eDirection::Y, eDirection::X}, YoungsModulus, kelvinChainStiffness,
                      kelvinChainRetardationTimes, poissonRatio);
    TestCreepModel<3>(testName, {eDirection::X, eDirection::Y, eDirection::Z}, YoungsModulus, kelvinChainStiffness,
                      kelvinChainRetardationTimes, poissonRatio);
    TestCreepModel<3>(testName, {eDirection::Y, eDirection::Z, eDirection::X}, YoungsModulus, kelvinChainStiffness,
                      kelvinChainRetardationTimes, poissonRatio);
    TestCreepModel<3>(testName, {eDirection::Z, eDirection::X, eDirection::Y}, YoungsModulus, kelvinChainStiffness,
                      kelvinChainRetardationTimes, poissonRatio);
}

int main(int argc, char* argv[])
{
    // Poisson Ratio = 0.0
    PerformTestSeries("TwoChainElementsWithSpring", 4.e9, (Eigen::VectorXd(2) << 20.e9, 5.e9).finished(),
                      (Eigen::VectorXd(2) << 5000., 10000.).finished(), 0.0);

    // Poisson Ratio = 0.2
    PerformTestSeries("TwoChainElementsWithSpring", 4.e9, (Eigen::VectorXd(2) << 20.e9, 5.e9).finished(),
                      (Eigen::VectorXd(2) << 5000., 10000.).finished(), 0.2);

    return 0;
}
