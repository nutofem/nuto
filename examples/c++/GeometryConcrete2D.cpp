#include <boost/filesystem.hpp>
#include <fstream>

#include "geometryConcrete/GeometryConcrete.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/MechanicsEnums.h"
#include "visualize/VisualizeEnum.h"

void CreateMesoscaleGeometryMesh(std::string rGmshFile, double rLX, double rLY)
{
    // define the geometry
    NuTo::GeometryConcrete geometry;
    geometry.SetSeed(1337);
    geometry.SetSpecimenBox(0, rLX, 0, rLY, 0, rLY);
    geometry.SetGradingCurve(NuTo::GeometryConcrete::B16, 3);
    geometry.SetParticleVolumeFraction(0.4);
    geometry.SetAbsoluteGrowthRate(0.1);

    geometry.MaximizeParticleDistance(0.75);

    geometry.ExportGmshGeo2D(rGmshFile, 0.75, rLY / 2., 0.75);

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

    // Create and import a mesoscale geometry via gmsh
    double lX = 32;
    double lY = 16;

    CreateMesoscaleGeometryMesh(gmshFile, lX, lY);

    auto groupIndices = myStructure.ImportFromGmsh(gmshFile + ".msh");
    assert(groupIndices.GetNumRows() == 2); // two physical groups
    assert(groupIndices.GetNumColumns() == 2); // 1st col: group, 2nd col: interpolation type

    int gMatrix = groupIndices(0, 0);
    int gAggreg = groupIndices(1, 0);

    int interpolationMatrix = groupIndices(0, 1);
    int interpolationAggreg = groupIndices(1, 1);

    myStructure.InterpolationTypeAdd(
            interpolationMatrix, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myStructure.InterpolationTypeAdd(
            interpolationAggreg, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);

    myStructure.ElementTotalConvertToInterpolationType();

    // Define and set sections and constitutive laws
    double thickness = 17;

    int mySection = myStructure.SectionCreate(NuTo::eSectionType::PLANE_STRAIN);
    myStructure.SectionSetThickness(mySection, thickness);

    int myConstitutiveLawAggreg =
            myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetParameterDouble(
            myConstitutiveLawAggreg, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 1.);
    myStructure.ConstitutiveLawSetParameterDouble(
            myConstitutiveLawAggreg, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, .2);

    int myConstitutiveLawMatrix = myStructure.ConstitutiveLawCreate(
            NuTo::Constitutive::eConstitutiveType::MISES_PLASTICITY_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetParameterDouble(
            myConstitutiveLawMatrix, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 75000.);
    myStructure.ConstitutiveLawSetParameterDouble(
            myConstitutiveLawMatrix, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, .3);
    myStructure.ConstitutiveLawSetParameterDouble(
            myConstitutiveLawMatrix, NuTo::Constitutive::eConstitutiveParameter::INITIAL_YIELD_STRENGTH, 75);
    myStructure.ConstitutiveLawSetParameterDouble(
            myConstitutiveLawMatrix, NuTo::Constitutive::eConstitutiveParameter::INITIAL_HARDENING_MODULUS, 5000);

    myStructure.ElementTotalSetSection(mySection);
    myStructure.ElementGroupSetConstitutiveLaw(gMatrix, myConstitutiveLawMatrix);
    myStructure.ElementGroupSetConstitutiveLaw(gAggreg, myConstitutiveLawAggreg);

    // Set boundary conditions
    double deltaD = 0.01;

    int gNodesWest = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    int gNodesEast = myStructure.GroupCreate(NuTo::eGroupId::Nodes);

    int iNodeOrigin = myStructure.NodeGetIdAtCoordinate(Eigen::Vector2d({0., 0.}), 1.e-6);
    myStructure.GroupAddNodeCoordinateRange(gNodesWest, 0, 0., 0.);
    myStructure.GroupAddNodeCoordinateRange(gNodesEast, 0, lX, lX);

    myStructure.ConstraintLinearSetDisplacementNodeGroup(gNodesWest, Eigen::Vector2d::UnitX(), 0.);
    myStructure.ConstraintLinearSetDisplacementNode(iNodeOrigin, Eigen::Vector2d::UnitY(), 0.);

    int bc = myStructure.ConstraintLinearSetDisplacementNodeGroup(gNodesEast, Eigen::Vector2d::UnitX(), 0);

    // Visualisation
    myStructure.AddVisualizationComponent(gAggreg, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(gAggreg, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    myStructure.AddVisualizationComponent(gAggreg, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    myStructure.AddVisualizationComponent(gAggreg, NuTo::eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);

    myStructure.AddVisualizationComponent(gMatrix, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(gMatrix, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    myStructure.AddVisualizationComponent(gMatrix, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    myStructure.AddVisualizationComponent(gMatrix, NuTo::eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
    myStructure.AddVisualizationComponent(gMatrix, NuTo::eVisualizeWhat::ENGINEERING_PLASTIC_STRAIN);

    // Solver
    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();
    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);

    double simulationTime = 1;

    Eigen::Matrix2d dispRHS;
    dispRHS << 0, 0, simulationTime, deltaD;

    myIntegrationScheme.AddTimeDependentConstraint(bc, dispRHS);
    myIntegrationScheme.SetTimeStep(.1 * simulationTime);
    myIntegrationScheme.SetToleranceForce(1e-6);
    myIntegrationScheme.SetAutomaticTimeStepping(true);
    myIntegrationScheme.SetVerboseLevel(0);
    myIntegrationScheme.SetShowTime(true);

    bool deleteDirectory = true;
    myIntegrationScheme.SetResultDirectory(resultPath.string(), deleteDirectory);

    try
    {
        myIntegrationScheme.Solve(simulationTime);
    }
    catch (NuTo::Exception& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
        std::cout << "\n\n\n Errors occured! \n\n\n" << std::endl;
    }

    std::cout << "I'm done. Thank you for using NuTo!" << std::endl;

    return 0;
}
