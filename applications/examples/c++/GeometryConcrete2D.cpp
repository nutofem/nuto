#include <boost/filesystem.hpp>
#include <fstream>

#include "geometryConcrete/GeometryConcrete.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/groups/Group.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/MechanicsEnums.h"
#include "mechanics/sections/SectionPlane.h"
#include "visualize/VisualizeEnum.h"

using namespace NuTo;

void CreateMesoscaleGeometryMesh(std::string rGmshFile, double rLX, double rLY)
{
    // define the geometry
    GeometryConcrete geometry;
    geometry.SetSeed(1337);
    geometry.SetSpecimenBox(0, rLX, 0, rLY, 0, rLY);
    geometry.SetGradingCurve(GeometryConcrete::B16, 3);
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

    Structure structure(2);
    structure.SetNumTimeDerivatives(0);
    structure.SetNumProcessors(4);
    structure.LoggerOpenFile(outputPath.string() + "Log.dat");
    structure.LoggerSetQuiet(false);
    structure.SetVerboseLevel(10);
    structure.SetShowTime(false);

    // Create and import a mesoscale geometry via gmsh
    double lX = 32;
    double lY = 16;

    CreateMesoscaleGeometryMesh(gmshFile, lX, lY);

    auto groupIndices = structure.ImportFromGmsh(gmshFile + ".msh");
    assert(groupIndices.size() == 2); // two physical groups

    int gMatrix = groupIndices[0].first;
    int gAggreg = groupIndices[1].first;

    int interpolationMatrix = groupIndices[0].second;
    int interpolationAggreg = groupIndices[1].second;

    structure.InterpolationTypeAdd(interpolationMatrix, Node::eDof::DISPLACEMENTS,
                                   Interpolation::eTypeOrder::EQUIDISTANT2);
    structure.InterpolationTypeAdd(interpolationAggreg, Node::eDof::DISPLACEMENTS,
                                   Interpolation::eTypeOrder::EQUIDISTANT2);

    structure.ElementTotalConvertToInterpolationType();

    // Define and set sections and constitutive laws
    double thickness = 17;

    auto mySection = SectionPlane::Create(thickness, true);

    int myConstitutiveLawAggreg =
            structure.ConstitutiveLawCreate(Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(myConstitutiveLawAggreg,
                                                Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 1.);
    structure.ConstitutiveLawSetParameterDouble(myConstitutiveLawAggreg,
                                                Constitutive::eConstitutiveParameter::POISSONS_RATIO, .2);

    int myConstitutiveLawMatrix =
            structure.ConstitutiveLawCreate(Constitutive::eConstitutiveType::MISES_PLASTICITY_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(myConstitutiveLawMatrix,
                                                Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 75000.);
    structure.ConstitutiveLawSetParameterDouble(myConstitutiveLawMatrix,
                                                Constitutive::eConstitutiveParameter::POISSONS_RATIO, .3);
    structure.ConstitutiveLawSetParameterDouble(myConstitutiveLawMatrix,
                                                Constitutive::eConstitutiveParameter::INITIAL_YIELD_STRENGTH, 75);
    structure.ConstitutiveLawSetParameterDouble(myConstitutiveLawMatrix,
                                                Constitutive::eConstitutiveParameter::INITIAL_HARDENING_MODULUS, 5000);

    structure.ElementTotalSetSection(mySection);
    structure.ElementGroupSetConstitutiveLaw(gMatrix, myConstitutiveLawMatrix);
    structure.ElementGroupSetConstitutiveLaw(gAggreg, myConstitutiveLawAggreg);

    // Set boundary conditions
    double deltaD = 0.01;

    auto& origin = structure.NodeGetAtCoordinate(Eigen::Vector2d({0., 0.}), 1.e-6);
    auto nodesWest = structure.GroupGetNodeCoordinateRange(eDirection::X, 0., 0.);
    auto nodesEast = structure.GroupGetNodeCoordinateRange(eDirection::X, lX, lX);

    structure.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(nodesWest, {eDirection::X}));
    structure.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Component(origin, {eDirection::Y}));

    double simulationTime = 1;
    structure.Constraints().Add(
            Node::eDof::DISPLACEMENTS,
            Constraint::Component(nodesEast, {eDirection::X}, Constraint::RhsRamp(simulationTime, deltaD)));

    // Visualisation
    structure.AddVisualizationComponent(gAggreg, eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(gAggreg, eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(gAggreg, eVisualizeWhat::ENGINEERING_STRESS);
    structure.AddVisualizationComponent(gAggreg, eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);

    structure.AddVisualizationComponent(gMatrix, eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(gMatrix, eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(gMatrix, eVisualizeWhat::ENGINEERING_STRESS);
    structure.AddVisualizationComponent(gMatrix, eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
    structure.AddVisualizationComponent(gMatrix, eVisualizeWhat::ENGINEERING_PLASTIC_STRAIN);

    // Solver
    structure.NodeBuildGlobalDofs();
    structure.CalculateMaximumIndependentSets();
    NewmarkDirect integrationScheme(&structure);
    integrationScheme.SetTimeStep(.1 * simulationTime);
    integrationScheme.SetToleranceForce(1e-6);
    integrationScheme.SetAutomaticTimeStepping(true);
    integrationScheme.SetVerboseLevel(0);
    integrationScheme.SetShowTime(true);

    bool deleteDirectory = true;
    integrationScheme.PostProcessing().SetResultDirectory(resultPath.string(), deleteDirectory);

    integrationScheme.Solve(simulationTime);

    std::cout << "I'm done. Thank you for using NuTo!" << std::endl;

    return 0;
}
