/*
 * GeometryConcrete2D.cpp
 *
 *  Created on: 3 Sep 2015
 *      Author: ttitsche
 */

#include <boost/filesystem.hpp>
#include <fstream>
#include "nuto/math/FullMatrix.h"

#include "nuto/geometryConcrete/GeometryConcrete.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"

void CreateMesoscaleGeometryMesh(std::string rGmshFile, double rLX, double rLY)
{
    // ***************************************************************************
    // ***********   define the geometry      ************************************

    NuTo::GeometryConcrete geometry;
    geometry.SetSeed(1337);
    geometry.SetSpecimenBox(0, rLX, 0, rLY, 0, rLY);
    geometry.SetGradingCurve(NuTo::GeometryConcrete::B16, 3);
    geometry.SetParticleVolumeFraction(0.4);
    geometry.SetAbsoluteGrowthRate(0.1);

    geometry.MaximizeParticleDistance(0.75);

    geometry.ExportGmshGeo2D(rGmshFile, 0.75, rLY/2.);

    std::cout << "Meshing..." << std::endl;
    system(("gmsh -2 -order 2 " + rGmshFile + ".geo -o " + rGmshFile + ".msh -v 2").c_str());

}




int main(int argc, char* argv[])
{
    boost::filesystem::path outputPath = std::string(argv[0]) + "Out/";
    boost::filesystem::path resultPath = outputPath.string() + "Results/";
    boost::filesystem::create_directory(outputPath);
    boost::filesystem::create_directory(resultPath);

    std::string gmshFile = outputPath.string() + "geometry";

    std::cout << "Gmsh File:  " << gmshFile << ".msh" << std::endl;
    std::cout << "Result dir: " << resultPath.string() << std::endl;

    NuTo::Structure myStructure(2);
    myStructure.SetNumTimeDerivatives(0);
    myStructure.SetNumProcessors(4);
    myStructure.LoggerOpenFile(outputPath.string() + "Log.dat");
    myStructure.LoggerSetQuiet(false);
    myStructure.SetVerboseLevel(10);
    myStructure.SetShowTime(false);


    // **********  Create and import a mesoscale geometry via gmsh  **************
    double lX = 32;
    double lY = 16;

    CreateMesoscaleGeometryMesh(gmshFile, lX, lY);

    auto groupIndices = myStructure.ImportFromGmsh(gmshFile+".msh", NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::STATICDATA);
    assert (groupIndices.GetNumRows() == 2); // two physical groups
    assert (groupIndices.GetNumColumns() == 2); // 1st col: group, 2nd col: interpolation type

    int gMatrix = groupIndices.GetValue(0,0);
    int gAggreg = groupIndices.GetValue(1,0);

    int interpolationMatrix = groupIndices.GetValue(0,1);
    int interpolationAggreg = groupIndices.GetValue(1,1);

    myStructure.InterpolationTypeAdd(interpolationMatrix, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);
    myStructure.InterpolationTypeAdd(interpolationAggreg, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);

    myStructure.ElementTotalConvertToInterpolationType();

    // **********  Define and set sections and constitutive laws    **************
    double thickness = 17;

    int mySection = myStructure.SectionCreate(NuTo::Section::PLANE_STRAIN);
    myStructure.SectionSetThickness(mySection, thickness);

    int myConstitutiveLawAggreg = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::LINEAR_ELASTIC_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetParameterDouble(myConstitutiveLawAggreg, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 1.);
    myStructure.ConstitutiveLawSetParameterDouble(myConstitutiveLawAggreg, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, .2);



    int myConstitutiveLawMatrix = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::MISES_PLASTICITY_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetParameterDouble(myConstitutiveLawMatrix, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 75000.);
    myStructure.ConstitutiveLawSetParameterDouble(myConstitutiveLawMatrix, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, .3);
    myStructure.ConstitutiveLawSetParameterDouble(myConstitutiveLawMatrix,NuTo::Constitutive::eConstitutiveParameter::INITIAL_YIELD_STRENGTH,75);
    myStructure.ConstitutiveLawSetParameterDouble(myConstitutiveLawMatrix,NuTo::Constitutive::eConstitutiveParameter::INITIAL_HARDENING_MODULUS,5000);


    myStructure.ElementTotalSetSection(mySection);
    myStructure.ElementGroupSetConstitutiveLaw(gMatrix, myConstitutiveLawMatrix);
    myStructure.ElementGroupSetConstitutiveLaw(gAggreg, myConstitutiveLawAggreg);

    // **********  Set boundary conditions  **************************************
    double deltaD = 0.01;

    int gNodesWest = myStructure.GroupCreate(NuTo::Groups::Nodes);
    int gNodesEast = myStructure.GroupCreate(NuTo::Groups::Nodes);

    int iNodeOrigin = myStructure.NodeGetIdAtCoordinate(NuTo::FullVector<double, 2>({0.,0.}), 1.e-6);
    myStructure.GroupAddNodeCoordinateRange(gNodesWest, 0, 0.,0.);
    myStructure.GroupAddNodeCoordinateRange(gNodesEast, 0, lX,lX);


    myStructure.ConstraintLinearSetDisplacementNodeGroup(gNodesWest, NuTo::FullVector<double,2>::UnitX(), 0.);
    myStructure.ConstraintLinearSetDisplacementNode(iNodeOrigin, NuTo::FullVector<double,2>::UnitY(), 0.);

    int bc = myStructure.ConstraintLinearSetDisplacementNodeGroup(gNodesEast, NuTo::FullVector<double,2>::UnitX(), 0);

    //**********************************************
    //          Visualisation
    //**********************************************
    myStructure.AddVisualizationComponent(gAggreg, NuTo::VisualizeBase::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(gAggreg, NuTo::VisualizeBase::ENGINEERING_STRAIN);
    myStructure.AddVisualizationComponent(gAggreg, NuTo::VisualizeBase::ENGINEERING_STRESS);
    myStructure.AddVisualizationComponent(gAggreg, NuTo::VisualizeBase::PRINCIPAL_ENGINEERING_STRESS);

    myStructure.AddVisualizationComponent(gMatrix, NuTo::VisualizeBase::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(gMatrix, NuTo::VisualizeBase::ENGINEERING_STRAIN);
    myStructure.AddVisualizationComponent(gMatrix, NuTo::VisualizeBase::ENGINEERING_STRESS);
    myStructure.AddVisualizationComponent(gMatrix, NuTo::VisualizeBase::PRINCIPAL_ENGINEERING_STRESS);
    myStructure.AddVisualizationComponent(gMatrix, NuTo::VisualizeBase::ENGINEERING_PLASTIC_STRAIN);

    //**********************************************
    //          Solver
    //**********************************************

    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();
    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);


    double simulationTime = 1;

    NuTo::FullMatrix<double, 2, 2> dispRHS;
    dispRHS << 0, 0, simulationTime, deltaD;

    myIntegrationScheme.AddTimeDependentConstraint(bc, dispRHS);
    myIntegrationScheme.SetTimeStep(.1*simulationTime);
    myIntegrationScheme.SetToleranceForce(1e-6);
    myIntegrationScheme.SetAutomaticTimeStepping(true);
    myIntegrationScheme.SetVerboseLevel(0);
    myIntegrationScheme.SetShowTime(true);

    bool deleteDirectory = true;
    myIntegrationScheme.SetResultDirectory(resultPath.string(), deleteDirectory);

    try
    {
        myIntegrationScheme.Solve(simulationTime);
    } catch (NuTo::Exception& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
        std::cout << "\n\n\n Errors occured! \n\n\n" << std::endl;
    }

    std::cout << "I'm done. Thanks for using NuTo!" << std::endl;

    return 0;
}

