/*
 * MisesPlasticity.cpp
 *
 *  Created on: 28 Sep 2015
 *      Author: ttitsche
 */
#include <boost/filesystem.hpp>
#include <string>
#include "mechanics/MechanicsEnums.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/groups/Group.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "visualize/VisualizeEnum.h"
#include "mechanics/mesh/MeshGenerator.h"


void SetConstitutiveLaw(NuTo::Structure& s)
{
    using namespace NuTo::Constitutive;
    int myMat = s.ConstitutiveLawCreate(eConstitutiveType::MISES_PLASTICITY_ENGINEERING_STRESS);

    s.ConstitutiveLawSetParameterDouble(myMat, eConstitutiveParameter::YOUNGS_MODULUS, 100);
    s.ConstitutiveLawSetParameterDouble(myMat, eConstitutiveParameter::POISSONS_RATIO, 0.0);
    s.ConstitutiveLawSetParameterDouble(myMat, eConstitutiveParameter::INITIAL_YIELD_STRENGTH, 100);
    s.ConstitutiveLawSetParameterDouble(myMat, eConstitutiveParameter::INITIAL_HARDENING_MODULUS, 10);
    s.ConstitutiveLawInfo(10);
    s.ElementTotalSetConstitutiveLaw(myMat);

    int visualizationGroup = s.GroupGetElementsTotal();
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_PLASTIC_STRAIN);
}


void Mises2D(const std::string& rDir)
{
    boost::filesystem::create_directory(rDir + "2D/");

    std::string gnuplotFileName = rDir + "LoadDisp2D.plt";
    std::ofstream gnuplotFile;
    gnuplotFile.open(gnuplotFileName);
    gnuplotFile << "#!/usr/bin/gnuplot \n";
    gnuplotFile << "data = '< paste " + rDir + "2D/Displ.dat " + rDir + "2D/Force.dat' \n";
    gnuplotFile << "plot data using 1:3 with lp\n";
    gnuplotFile << "pause 2; refresh; reread;\n";
    gnuplotFile.close();

    NuTo::Structure s(2);
    s.SetVerboseLevel(10);

    int interpol = NuTo::MeshGenerator::Grid(s, {1., 1.}, {1, 1}).second;
    s.InterpolationTypeAdd(interpol, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    s.ElementTotalConvertToInterpolationType();

    SetConstitutiveLaw(s);

    auto section = NuTo::SectionPlane::Create(3.1415, true);
    s.ElementTotalSetSection(section);

    // boundary conditions
    auto& nodeOrigin = s.NodeGetAtCoordinate(Eigen::Vector2d::Zero());
    auto& nodesLeft = s.GroupGetNodesAtCoordinate(NuTo::eDirection::X, 0);
    auto& nodesRight = s.GroupGetNodesAtCoordinate(NuTo::eDirection::X, 1);

    double simulationTime = 1;
    double deltaD = 5;

    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
                        NuTo::Constraint::Component(nodeOrigin, {NuTo::eDirection::Y}));
    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Constraint::Component(nodesLeft, {NuTo::eDirection::X}));
    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
                        NuTo::Constraint::Component(nodesRight, {NuTo::eDirection::X},
                                                    NuTo::Constraint::RhsRamp(simulationTime, deltaD)));

    s.NodeBuildGlobalDofs();
    s.CalculateMaximumIndependentSets();

    NuTo::NewmarkDirect myIntegrationScheme(&s);

    myIntegrationScheme.SetResultDirectory(rDir + "2D/", true);
    myIntegrationScheme.AddResultGroupNodeForce("Force", s.GroupGetId(&nodesRight));
    myIntegrationScheme.AddResultNodeDisplacements("Displ", nodesRight.GetMemberIds()[0]);
    myIntegrationScheme.SetTimeStep(0.02);
    myIntegrationScheme.SetToleranceForce(1e-6);
    myIntegrationScheme.SetAutomaticTimeStepping(true);
    myIntegrationScheme.SetMaxTimeStep(0.02);
    myIntegrationScheme.SetVerboseLevel(0);
    myIntegrationScheme.SetShowTime(true);
    myIntegrationScheme.SetPerformLineSearch(true);

    s.SetVerboseLevel(0);
    myIntegrationScheme.Solve(simulationTime);
}

void Mises3D(const std::string& rDir)
{
    boost::filesystem::create_directory(rDir + "3D/");

    std::string gnuplotFileName = rDir + "LoadDisp3D.plt";
    std::ofstream gnuplotFile;
    gnuplotFile.open(gnuplotFileName);
    gnuplotFile << "#!/usr/bin/gnuplot \n";
    gnuplotFile << "data = '< paste " + rDir + "3D/Displ.dat " + rDir + "3D/Force.dat' \n";
    gnuplotFile << "plot data using 1:4 with lp\n";
    gnuplotFile << "pause 2; refresh; reread;\n";
    gnuplotFile.close();

    NuTo::Structure s(3);
    s.SetVerboseLevel(10);

    int interpol = NuTo::MeshGenerator::Grid(s, {1., 1., 1.}, {1, 1, 1}, NuTo::Interpolation::eShapeType::TETRAHEDRON3D)
                           .second;
    s.InterpolationTypeAdd(interpol, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    s.ElementTotalConvertToInterpolationType();

    SetConstitutiveLaw(s);

    // boundary conditions
    auto& nodeOrigin = s.NodeGetAtCoordinate(Eigen::Vector3d::Zero());
    auto& nodeFix = s.NodeGetAtCoordinate(Eigen::Vector3d(0, 1, 0));
    auto& nodesLeft = s.GroupGetNodesAtCoordinate(NuTo::eDirection::X, 0);
    auto& nodesRight = s.GroupGetNodesAtCoordinate(NuTo::eDirection::X, 1);

    double simulationTime = 1;
    double deltaD = 5;

    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
                        NuTo::Constraint::Component(nodeOrigin, {NuTo::eDirection::Y, NuTo::eDirection::Z}));
    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Constraint::Component(nodeFix, {NuTo::eDirection::Z}));
    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Constraint::Component(nodesLeft, {NuTo::eDirection::X}));
    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
                        NuTo::Constraint::Component(nodesRight, {NuTo::eDirection::X},
                                                    NuTo::Constraint::RhsRamp(simulationTime, deltaD)));


    s.NodeBuildGlobalDofs();
    s.CalculateMaximumIndependentSets();

    NuTo::NewmarkDirect myIntegrationScheme(&s);

    myIntegrationScheme.AddResultGroupNodeForce("Force", s.GroupGetId(&nodesLeft));
    myIntegrationScheme.AddResultNodeDisplacements("Displ", nodesLeft.GetMemberIds()[0]);
    myIntegrationScheme.SetTimeStep(0.02);
    myIntegrationScheme.SetToleranceForce(1e-6);
    myIntegrationScheme.SetAutomaticTimeStepping(true);
    myIntegrationScheme.SetMaxTimeStep(0.02);
    myIntegrationScheme.SetVerboseLevel(0);
    myIntegrationScheme.SetShowTime(true);
    myIntegrationScheme.SetPerformLineSearch(true);

    bool deleteDirectory = true;
    myIntegrationScheme.SetResultDirectory(rDir + "3D/", deleteDirectory);

    s.SetVerboseLevel(0);
    myIntegrationScheme.Solve(simulationTime);
}

int main(int, char* argv[])
{
    std::string resultDir = std::string(argv[0]) + "Out/";
    boost::filesystem::create_directory(resultDir);

    Mises2D(resultDir);
    Mises3D(resultDir);
}
