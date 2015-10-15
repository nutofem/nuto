/*
 * MisesPlasticity.cpp
 *
 *  Created on: 28 Sep 2015
 *      Author: ttitsche
 */


#include <boost/filesystem.hpp>
#include <string>
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"


void SetConstitutiveLaw(NuTo::Structure& rStructure)
{
    int myMat = rStructure.ConstitutiveLawCreate(NuTo::Constitutive::MISES_PLASTICITY_ENGINEERING_STRESS);

    rStructure.ConstitutiveLawSetParameterDouble(myMat,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,100);
    rStructure.ConstitutiveLawSetParameterDouble(myMat,NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,0.0);
    rStructure.ConstitutiveLawSetParameterDouble(myMat,NuTo::Constitutive::eConstitutiveParameter::INITIAL_YIELD_STRENGTH,100);
//    rStructure.ConstitutiveLawAddYieldStrength(myMat,0.25,150);
//    rStructure.ConstitutiveLawAddYieldStrength(myMat,0.3,150);

    rStructure.ConstitutiveLawSetParameterDouble(myMat,NuTo::Constitutive::eConstitutiveParameter::INITIAL_HARDENING_MODULUS,10);
    rStructure.ConstitutiveLawInfo(10);

    rStructure.ElementTotalSetConstitutiveLaw(myMat);

    rStructure.AddVisualizationComponentEngineeringStrain();
    rStructure.AddVisualizationComponentDisplacements();
    rStructure.AddVisualizationComponentEngineeringStress();
    rStructure.AddVisualizationComponentEngineeringPlasticStrain();

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

    NuTo::FullVector<int, Eigen::Dynamic> nodeIds(4);
    nodeIds[0] = myStructure.NodeCreate(NuTo::FullVector<double, 2>({1,1}));
    nodeIds[1] = myStructure.NodeCreate(NuTo::FullVector<double, 2>({0,1}));
    nodeIds[2] = myStructure.NodeCreate(NuTo::FullVector<double, 2>({0,0}));
    nodeIds[3] = myStructure.NodeCreate(NuTo::FullVector<double, 2>({1,0}));

    int interpol = myStructure.InterpolationTypeCreate(NuTo::Interpolation::QUAD2D);
    myStructure.InterpolationTypeAdd(interpol, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(interpol, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);


    myStructure.ElementCreate(interpol, nodeIds, NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::STATICDATA);
    myStructure.ElementTotalConvertToInterpolationType();

    SetConstitutiveLaw(myStructure);

    int section = myStructure.SectionCreate(NuTo::Section::PLANE_STRAIN);
    myStructure.SectionSetThickness(section, 3.1415);
    myStructure.ElementTotalSetSection(section);

    // boundary conditions
    int nOrigin = myStructure.NodeGetIdAtCoordinate(NuTo::FullVector<double, 2>({0,0}), 1.e-6);

    int gNodesXminus = myStructure.GroupCreate(NuTo::Groups::Nodes);
    int gNodesXplus  = myStructure.GroupCreate(NuTo::Groups::Nodes);

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


    myIntegrationScheme.SetTimeDependentConstraint(bc, dispRHS);
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

    int interpol = myStructure.InterpolationTypeCreate(NuTo::Interpolation::BRICK3D);
    myStructure.InterpolationTypeAdd(interpol, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(interpol, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);


    myStructure.ElementCreate(interpol, nodeIds, NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::STATICDATA);
    myStructure.ElementTotalConvertToInterpolationType();

    SetConstitutiveLaw(myStructure);

    myStructure.ElementTotalSetSection(myStructure.SectionCreate(NuTo::Section::VOLUME));

    // boundary conditions
    int nOrigin = myStructure.NodeGetIdAtCoordinate(NuTo::FullVector<double, 3>({0,0,0}), 1.e-6);

    int gNodesXminus = myStructure.GroupCreate(NuTo::Groups::Nodes);
    int gNodesXplus  = myStructure.GroupCreate(NuTo::Groups::Nodes);

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


    myIntegrationScheme.SetTimeDependentConstraint(bc, dispRHS);
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


    std::cout << std::endl;
    std::cout << "#####################################" << std::endl;
    std::cout << "##  Congratulations! Everything    ##" << std::endl;
    std::cout << "##   went better than expected!    ##" << std::endl;
    std::cout << "#####################################" << std::endl;

    return 0;
}

