/*
 * MisesPlasticity.cpp
 *
 *  Created on: 28 Sep 2015
 *      Author: ttitsche
 */


#include <boost/filesystem.hpp>
#include <string>
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/elements/IpDataEnum.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "visualize/VisualizeEnum.h"
#include "mechanics/mesh/MeshGenerator.h"


void SetConstitutiveLaw(NuTo::Structure& rStructure)
{
    int myMat = rStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::MISES_PLASTICITY_ENGINEERING_STRESS);

    rStructure.ConstitutiveLawSetParameterDouble(myMat,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,100);
    rStructure.ConstitutiveLawSetParameterDouble(myMat,NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,0.0);
    rStructure.ConstitutiveLawSetParameterDouble(myMat,NuTo::Constitutive::eConstitutiveParameter::INITIAL_YIELD_STRENGTH,100);
//    rStructure.ConstitutiveLawAddYieldStrength(myMat,0.25,150);
//    rStructure.ConstitutiveLawAddYieldStrength(myMat,0.3,150);

    rStructure.ConstitutiveLawSetParameterDouble(myMat,NuTo::Constitutive::eConstitutiveParameter::INITIAL_HARDENING_MODULUS,10);
    rStructure.ConstitutiveLawInfo(10);

    rStructure.ElementTotalSetConstitutiveLaw(myMat);

    int visualizationGroup = rStructure.GroupCreate(NuTo::eGroupId::Elements);
    rStructure.GroupAddElementsTotal(visualizationGroup);

    rStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    rStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    rStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    rStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_PLASTIC_STRAIN);


}


void Mises2D(const std::string& rDir)
{

    boost::filesystem::create_directory(rDir+"2D/");


    std::string gnuplotFileName = rDir+"LoadDisp2D.plt";
    std::ofstream gnuplotFile;
    gnuplotFile.open(gnuplotFileName);
    gnuplotFile << "#!/usr/bin/gnuplot \n";
    gnuplotFile << "data = '< paste " + rDir + "2D/Displ.dat "+ rDir +"2D/Force.dat' \n";
    gnuplotFile << "plot data using 1:3 with lp\n";
    gnuplotFile << "pause 2; refresh; reread;\n";
    gnuplotFile.close();

    NuTo::Structure myStructure(2);
    myStructure.SetVerboseLevel(10);

    int interpol = NuTo::MeshGenerator::Grid(myStructure, {1., 1.}, {1, 1}).second;
    myStructure.InterpolationTypeAdd(interpol, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myStructure.ElementTotalConvertToInterpolationType();

    SetConstitutiveLaw(myStructure);

    auto section = NuTo::SectionPlane::Create(3.1415, true);
    myStructure.ElementTotalSetSection(section);

    // boundary conditions
    int nOrigin = myStructure.NodeGetIdAtCoordinate(Eigen::Vector2d::Zero(), 1.e-6);

    int gNodesXminus = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    int gNodesXplus  = myStructure.GroupCreate(NuTo::eGroupId::Nodes);

    myStructure.GroupAddNodeCoordinateRange(gNodesXminus, 0, 0, 0);
    myStructure.GroupAddNodeCoordinateRange(gNodesXplus,  0, 1, 1);

    myStructure.ConstraintLinearSetDisplacementNodeGroup(gNodesXminus, Eigen::Vector2d::UnitX(), 0);
    myStructure.ConstraintLinearSetDisplacementNode(nOrigin, Eigen::Vector2d::UnitY(), 0);

    int bc = myStructure.ConstraintLinearSetDisplacementNodeGroup(gNodesXplus, Eigen::Vector2d::UnitX(), 0);

    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);


    myIntegrationScheme.SetResultDirectory(rDir+"2D/", true);
    myIntegrationScheme.AddResultGroupNodeForce("Force", gNodesXplus);
    myIntegrationScheme.AddResultNodeDisplacements("Displ", myStructure.GroupGetMemberIds(gNodesXplus)[0]);

    double simulationTime = 1;
    double deltaD = 5;


    Eigen::Matrix2d dispRHS;
    dispRHS << 0, 0, simulationTime, deltaD;


    myIntegrationScheme.AddTimeDependentConstraint(bc, dispRHS);
    myIntegrationScheme.SetTimeStep(0.02);

    myIntegrationScheme.SetToleranceForce(1e-6);
    myIntegrationScheme.SetAutomaticTimeStepping(true);
    myIntegrationScheme.SetMaxTimeStep(0.02);
    myIntegrationScheme.SetVerboseLevel(0);
    myIntegrationScheme.SetShowTime(true);
    myIntegrationScheme.SetPerformLineSearch(true);


    myStructure.SetVerboseLevel(0);

    myIntegrationScheme.Solve(simulationTime);

}




void Mises3D(const std::string& rDir)
{

    boost::filesystem::create_directory(rDir+"3D/");


    std::string gnuplotFileName = rDir+"LoadDisp3D.plt";
    std::ofstream gnuplotFile;
    gnuplotFile.open(gnuplotFileName);
    gnuplotFile << "#!/usr/bin/gnuplot \n";
    gnuplotFile << "data = '< paste " + rDir + "3D/Displ.dat "+ rDir +"3D/Force.dat' \n";
    gnuplotFile << "plot data using 1:4 with lp\n";
    gnuplotFile << "pause 2; refresh; reread;\n";
    gnuplotFile.close();

    NuTo::Structure myStructure(3);
    myStructure.SetVerboseLevel(10);

    int interpol = NuTo::MeshGenerator::Grid(myStructure, {1., 1., 1.}, {1, 1, 1},
                                             NuTo::Interpolation::eShapeType::TETRAHEDRON3D).second;
    myStructure.InterpolationTypeAdd(interpol, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myStructure.ElementTotalConvertToInterpolationType();

    SetConstitutiveLaw(myStructure);

    // boundary conditions
    int nOrigin = myStructure.NodeGetIdAtCoordinate(Eigen::Vector3d::Zero(), 1.e-6);
    int nFix = myStructure.NodeGetIdAtCoordinate(Eigen::Vector3d(0,1,0), 1.e-6);

    int gNodesXminus = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    int gNodesXplus  = myStructure.GroupCreate(NuTo::eGroupId::Nodes);

    myStructure.GroupAddNodeCoordinateRange(gNodesXminus, 0, 0, 0);
    myStructure.GroupAddNodeCoordinateRange(gNodesXplus,  0, 1, 1);

    myStructure.ConstraintLinearSetDisplacementNodeGroup(gNodesXminus, Eigen::Vector3d::UnitX(), 0);

    myStructure.ConstraintLinearSetDisplacementNode(nOrigin, Eigen::Vector3d::UnitY(), 0);
    myStructure.ConstraintLinearSetDisplacementNode(nOrigin, Eigen::Vector3d::UnitZ(), 0);
    myStructure.ConstraintLinearSetDisplacementNode(nFix, Eigen::Vector3d::UnitZ(), 0);

    int bc = myStructure.ConstraintLinearSetDisplacementNodeGroup(gNodesXplus, Eigen::Vector3d::UnitX(), 0);

    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);

    myIntegrationScheme.AddResultGroupNodeForce("Force", gNodesXplus);
    myIntegrationScheme.AddResultNodeDisplacements("Displ", myStructure.GroupGetMemberIds(gNodesXplus)[0]);

    double simulationTime = 1;
    double deltaD = 5;


    Eigen::Matrix2d dispRHS;
    dispRHS << 0, 0, simulationTime, deltaD;


    myIntegrationScheme.AddTimeDependentConstraint(bc, dispRHS);
    myIntegrationScheme.SetTimeStep(0.02);

    myIntegrationScheme.SetToleranceForce(1e-6);
    myIntegrationScheme.SetAutomaticTimeStepping(true);
    myIntegrationScheme.SetMaxTimeStep(0.02);
    myIntegrationScheme.SetVerboseLevel(0);
    myIntegrationScheme.SetShowTime(true);
    myIntegrationScheme.SetPerformLineSearch(true);

    bool deleteDirectory = true;
    myIntegrationScheme.SetResultDirectory(rDir+"3D/", deleteDirectory);

    myStructure.SetVerboseLevel(0);

    myIntegrationScheme.Solve(simulationTime);

}


int main(int argc, char* argv[])
{
    std::string resultDir = std::string(argv[0]) + "Out/";
    boost::filesystem::create_directory(resultDir);

    Mises2D(resultDir);
    Mises3D(resultDir);


    std::cout << std::endl;
    std::cout << "#####################################" << std::endl;
    std::cout << "##  Congratulations! Everything    ##" << std::endl;
    std::cout << "##   went better than expected!    ##" << std::endl;
    std::cout << "#####################################" << std::endl;

    return 0;
}

