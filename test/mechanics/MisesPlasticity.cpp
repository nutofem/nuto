/*
 * MisesPlasticity.cpp
 *
 *  Created on: 28 Sep 2015
 *      Author: ttitsche
 */


#include <boost/filesystem.hpp>
#include <string>
#include "math/FullMatrix.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/elements/IpDataEnum.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/sections/SectionEnum.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "visualize/VisualizeEnum.h"


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

    std::vector<int> nodeIds(4);
    nodeIds[0] = myStructure.NodeCreate(NuTo::FullVector<double, 2>({1,1}));
    nodeIds[1] = myStructure.NodeCreate(NuTo::FullVector<double, 2>({0,1}));
    nodeIds[2] = myStructure.NodeCreate(NuTo::FullVector<double, 2>({0,0}));
    nodeIds[3] = myStructure.NodeCreate(NuTo::FullVector<double, 2>({1,0}));

    int interpol = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
    myStructure.InterpolationTypeAdd(interpol, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(interpol, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);


    myStructure.ElementCreate(interpol, nodeIds);
    myStructure.ElementTotalConvertToInterpolationType();

    SetConstitutiveLaw(myStructure);

    int section = myStructure.SectionCreate(NuTo::eSectionType::PLANE_STRAIN);
    myStructure.SectionSetThickness(section, 3.1415);
    myStructure.ElementTotalSetSection(section);

    // boundary conditions
    int nOrigin = myStructure.NodeGetIdAtCoordinate(NuTo::FullVector<double, 2>({0,0}), 1.e-6);

    int gNodesXminus = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    int gNodesXplus  = myStructure.GroupCreate(NuTo::eGroupId::Nodes);

    myStructure.GroupAddNodeCoordinateRange(gNodesXminus, 0, 0, 0);
    myStructure.GroupAddNodeCoordinateRange(gNodesXplus,  0, 1, 1);

    NuTo::FullVector<double, Eigen::Dynamic> dirX = NuTo::FullVector<double, 2>::UnitX();
    NuTo::FullVector<double, Eigen::Dynamic> dirY = NuTo::FullVector<double, 2>::UnitY();

    myStructure.ConstraintLinearSetDisplacementNodeGroup(gNodesXminus, dirX, 0);
    myStructure.ConstraintLinearSetDisplacementNode(nOrigin, dirY, 0);

    int bc = myStructure.ConstraintLinearSetDisplacementNodeGroup(gNodesXplus, dirX, 0);

    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);


    myIntegrationScheme.SetResultDirectory(rDir+"2D/", true);
    myIntegrationScheme.AddResultGroupNodeForce("Force", gNodesXplus);
    myIntegrationScheme.AddResultNodeDisplacements("Displ", myStructure.GroupGetMemberIds(gNodesXplus).GetValue(0));

    double simulationTime = 1;
    double deltaD = 5;


    NuTo::FullMatrix<double, 2, 2> dispRHS;
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

    NuTo::FullVector<int, Eigen::Dynamic> nodeIds(8);


    nodeIds[0] = myStructure.NodeCreate(NuTo::FullVector<double, 3>({1,1,1}));
    nodeIds[1] = myStructure.NodeCreate(NuTo::FullVector<double, 3>({0,1,1}));
    nodeIds[2] = myStructure.NodeCreate(NuTo::FullVector<double, 3>({0,0,1}));
    nodeIds[3] = myStructure.NodeCreate(NuTo::FullVector<double, 3>({1,0,1}));
    nodeIds[4] = myStructure.NodeCreate(NuTo::FullVector<double, 3>({1,1,0}));
    nodeIds[5] = myStructure.NodeCreate(NuTo::FullVector<double, 3>({0,1,0}));
    nodeIds[6] = myStructure.NodeCreate(NuTo::FullVector<double, 3>({0,0,0}));
    nodeIds[7] = myStructure.NodeCreate(NuTo::FullVector<double, 3>({1,0,0}));


    std::vector<int> nodesTet0 = {nodeIds[0], nodeIds[1], nodeIds[3], nodeIds[7]};
    std::vector<int> nodesTet1 = {nodeIds[0], nodeIds[1], nodeIds[7], nodeIds[4]};
    std::vector<int> nodesTet2 = {nodeIds[5], nodeIds[4], nodeIds[7], nodeIds[1]};
    std::vector<int> nodesTet3 = {nodeIds[6], nodeIds[5], nodeIds[7], nodeIds[1]};
    std::vector<int> nodesTet4 = {nodeIds[2], nodeIds[7], nodeIds[1], nodeIds[6]};
    std::vector<int> nodesTet5 = {nodeIds[2], nodeIds[3], nodeIds[1], nodeIds[7]};




    int interpol = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TETRAHEDRON3D);
    myStructure.InterpolationTypeAdd(interpol, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(interpol, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);


    myStructure.ElementCreate(interpol, nodesTet0);
    myStructure.ElementCreate(interpol, nodesTet1);
    myStructure.ElementCreate(interpol, nodesTet2);
    myStructure.ElementCreate(interpol, nodesTet3);
    myStructure.ElementCreate(interpol, nodesTet4);
    myStructure.ElementCreate(interpol, nodesTet5);



    myStructure.ElementTotalConvertToInterpolationType();

    SetConstitutiveLaw(myStructure);

    myStructure.ElementTotalSetSection(myStructure.SectionCreate(NuTo::eSectionType::VOLUME));

    // boundary conditions
    int nOrigin = myStructure.NodeGetIdAtCoordinate(NuTo::FullVector<double, 3>({0,0,0}), 1.e-6);

    int gNodesXminus = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    int gNodesXplus  = myStructure.GroupCreate(NuTo::eGroupId::Nodes);

    myStructure.GroupAddNodeCoordinateRange(gNodesXminus, 0, 0, 0);
    myStructure.GroupAddNodeCoordinateRange(gNodesXplus,  0, 1, 1);

    NuTo::FullVector<double, Eigen::Dynamic> dirX = NuTo::FullVector<double, 3>::UnitX();
    NuTo::FullVector<double, Eigen::Dynamic> dirY = NuTo::FullVector<double, 3>::UnitY();
    NuTo::FullVector<double, Eigen::Dynamic> dirZ = NuTo::FullVector<double, 3>::UnitZ();

    myStructure.ConstraintLinearSetDisplacementNodeGroup(gNodesXminus, dirX, 0);
    myStructure.ConstraintLinearSetDisplacementNode(nOrigin, dirY, 0);
    myStructure.ConstraintLinearSetDisplacementNode(nOrigin, dirZ, 0);

    int bc = myStructure.ConstraintLinearSetDisplacementNodeGroup(gNodesXplus, dirX, 0);

    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);

    myIntegrationScheme.AddResultGroupNodeForce("Force", gNodesXplus);
    myIntegrationScheme.AddResultNodeDisplacements("Displ", myStructure.GroupGetMemberIds(gNodesXplus).GetValue(0));

    double simulationTime = 1;
    double deltaD = 5;


    NuTo::FullMatrix<double, 2, 2> dispRHS;
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

