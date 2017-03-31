//! \file Coupled2DMesoTransient.cpp
//! \brief Simulation of a plate under thermal loading. *Beware:* To run the
//! example, you need to create the needed mesh with gmsh by running 
//! `gmsh -2 -order 2 -o HeatedPlate.msh HeatedPlate.geo`.
//! Thermal BCs: 
//!     \f[T|_{x = 0} = 0\f]
//!     \f[T|_{x = 100\mathrm{mm}} = f(t) \text{ where $f$ is ISO 834 fire curve }\f]
//!     \f[q|_{y=0, y = 20\mathrm{mm}} = 0\f]
//!
#include <iostream>
#include <cmath>
#include "math/SparseMatrixCSRGeneral.h"
#include "math/SparseDirectSolverMUMPS.h"
#include "math/LinearInterpolation.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "mechanics/constitutive/laws/AdditiveOutput.h"
#include "mechanics/constitutive/laws/ThermalStrains.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/MechanicsEnums.h"
#include "visualize/VisualizeEnum.h"

using namespace NuTo;

struct Properties
{
    double youngsModulus;
    double capacity;
    double conductivity;
    double density;
    double expansionCoeff;
};

double iso_temperature_curve(double seconds)
{
    return 345.0*std::log10(8.0*seconds/60.0 + 1.0);
}

std::array<double, 2> SandstoneExpansion(double temperature)
{
    static std::vector<std::array<double, 2>> values = {{{0.0, 0.0},
                                                   {200.0, 0.25e-2},
                                                   {400.0, 0.6e-2},
                                                   {600.0, 1.3e-2},
                                                   {800.0, 1.5e-2},
                                                  {1600.0, 1.5e-2}}};
    static auto interpolation = Math::LinearInterpolation(values);
    return {interpolation(temperature), interpolation.derivative(temperature)};
}

std::array<double,2> CruzGillenCement(double temperature)
{
    static std::vector<std::array<double, 2>> values = {{{0.0, 0.0},
                                                   {100.0, 0.2e-2},
                                                   {200.0, 0.2e-2},
                                                   {300.0, -0.2e-2},
                                                   {400.0, -0.6e-2},
                                                   {500.0, -1.1e-2},
                                                   {600.0, -1.5e-2},
                                                   {700.0, -1.7e-2},
                                                   {800.0, -1.8e-2},
                                                  {1600.0, -1.8e-2}}};
    static auto interpolation = Math::LinearInterpolation(values);
    return {interpolation(temperature), interpolation.derivative(temperature)};
}

void SetConstitutiveLaws(Structure &structure, int group, Properties properties,
        std::function<std::array<double, 2>(double)> ExpansionFunction)
{
    int additive_input_id = structure.ConstitutiveLawCreate(
            Constitutive::eConstitutiveType::ADDITIVE_INPUT_EXPLICIT);
    int additive_output_id = structure.ConstitutiveLawCreate(
            Constitutive::eConstitutiveType::ADDITIVE_OUTPUT);

    int lin_elastic_id = structure.ConstitutiveLawCreate(
            Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(lin_elastic_id,
            Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, properties.youngsModulus);
    structure.ConstitutiveLawSetParameterDouble(lin_elastic_id,
            Constitutive::eConstitutiveParameter::POISSONS_RATIO, .2);

    int heat_conduction_id = structure.ConstitutiveLawCreate(
            Constitutive::eConstitutiveType::HEAT_CONDUCTION);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id,
            Constitutive::eConstitutiveParameter::HEAT_CAPACITY, properties.capacity);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id,
            Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY, properties.conductivity);
    structure.ConstitutiveLawSetParameterDouble(heat_conduction_id,
            Constitutive::eConstitutiveParameter::DENSITY, properties.density);

    int thermal_strains_id = structure.ConstitutiveLawCreate(
            Constitutive::eConstitutiveType::THERMAL_STRAINS);
    structure.ConstitutiveLawSetParameterDouble(thermal_strains_id,
            Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT, properties.expansionCoeff);

    AdditiveInputExplicit* additive_input =
        static_cast<AdditiveInputExplicit*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_input_id));
    AdditiveOutput* additive_output =
        static_cast<AdditiveOutput*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_output_id));
    ConstitutiveBase* lin_elastic = structure.ConstitutiveLawGetConstitutiveLawPtr(lin_elastic_id);
    ConstitutiveBase* thermal_strains = structure.ConstitutiveLawGetConstitutiveLawPtr(thermal_strains_id);
    ConstitutiveBase* heat_conduction = structure.ConstitutiveLawGetConstitutiveLawPtr(heat_conduction_id);

    thermal_strains->SetParameterFunction(ExpansionFunction);

    additive_input->AddConstitutiveLaw(*lin_elastic);
    additive_input->AddConstitutiveLaw(*thermal_strains, Constitutive::eInput::ENGINEERING_STRAIN);

    additive_output->AddConstitutiveLaw(*additive_input);
    additive_output->AddConstitutiveLaw(*heat_conduction);

    structure.ElementGroupSetConstitutiveLaw(group, additive_output_id);
}

void SetInterpolation(Structure& structure, int group)
{
    structure.InterpolationTypeAdd(group, Node::eDof::DISPLACEMENTS,
            Interpolation::eTypeOrder::EQUIDISTANT2);
    structure.InterpolationTypeAdd(group, Node::eDof::TEMPERATURE,
            Interpolation::eTypeOrder::EQUIDISTANT2);
    structure.InterpolationTypeSetIntegrationType(group, eIntegrationType::IntegrationType2D3NGauss4Ip);
}

void SetVisualization(Structure& structure)
{
    int visualizationGroup = structure.GroupCreate(eGroupId::Elements);
    structure.GroupAddElementsTotal(visualizationGroup);

    structure.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::TEMPERATURE);
    structure.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::ENGINEERING_STRESS);
    structure.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::THERMAL_STRAIN);
    structure.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
}

int main()
{
    // create one-dimensional structure
    Structure structure(2);
    structure.SetNumTimeDerivatives(1);

    // import mesh (needs to be generated by gmsh first!)
    auto groupIndices = structure.ImportFromGmsh("./HeatedPlate.msh");

    // create section
    double thickness = 20.0;
    auto section = SectionPlane::Create(thickness, true);
    structure.ElementTotalSetSection(section);

    auto matrix_group = groupIndices[0].first;
    auto aggregate_group = groupIndices[1].first;

    // set constitutive laws
    //Properties concrete = {41071.0, 893e-6, 1.89, 2899.0, 17.3e-6};
    //SetConstitutiveLaws(structure, matrix_group, concrete, CruzGillenCement);
    
    Properties matrix_properties = {25e3, 1000e-6, 1.1, 3120.0, 20.0e-6};
    SetConstitutiveLaws(structure, matrix_group, matrix_properties, CruzGillenCement);

    Properties aggregate_properties = {70e3, 700e-6, 3.3, 2500.0, 12.5e-6};
    SetConstitutiveLaws(structure, aggregate_group, aggregate_properties, SandstoneExpansion);

    // set interpolation types
    auto interpolationMatrix = groupIndices[0].second;
    auto interpolationAggreg = groupIndices[1].second;

    SetInterpolation(structure, interpolationMatrix);
    SetInterpolation(structure, interpolationAggreg);

    structure.ElementTotalConvertToInterpolationType();

    SetVisualization(structure);

    // set boundary conditions and loads
    auto nodesWest = structure.GroupCreate("Nodes");
    auto nodesEast = structure.GroupCreate("Nodes");
    auto nodesSouth = structure.GroupCreate(eGroupId::Nodes);
    auto nodesNorth = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(nodesWest, 0, 0.0, 0.0);
    structure.GroupAddNodeCoordinateRange(nodesEast, 0, 32.0, 32.0);
    structure.GroupAddNodeCoordinateRange(nodesSouth, 1, 0.0, 0.0);
    structure.GroupAddNodeCoordinateRange(nodesNorth, 1, 16.0, 16.0);

    // displacement BC
    structure.ConstraintLinearSetDisplacementNodeGroup(nodesWest, Eigen::Vector2d::UnitX(), 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(nodesSouth, Eigen::Vector2d::UnitY(), 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(nodesNorth, Eigen::Vector2d::UnitY(), 0.0);

    // temperature BC
    structure.ConstraintLinearSetTemperatureNodeGroup(nodesWest, 0.0);
    auto east_bc = structure.ConstraintLinearSetTemperatureNodeGroup(nodesEast, 0.0);

    auto test = SandstoneExpansion(300.0);
    std::cout << test[0] << " " << test[1] << std::endl;

    // solve system
    NewmarkDirect newmark(&structure);
    double simulationTime = 3600.0;
    newmark.AddTimeDependentConstraintFunction(east_bc, iso_temperature_curve);
    newmark.SetPerformLineSearch(false);
    newmark.SetTimeStep(45);
    newmark.SetMaxTimeStep(0.2*simulationTime);
    newmark.SetToleranceResidual(Node::eDof::TEMPERATURE, 1e-4);
    newmark.SetToleranceResidual(Node::eDof::DISPLACEMENTS, 1e-3);
    newmark.SetAutomaticTimeStepping(true);

    bool deleteDirectory = true;
    newmark.SetResultDirectory("HeatedPlateResults", deleteDirectory);
    newmark.Solve(simulationTime);

    return 0;
}
