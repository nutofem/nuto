#include "visualize/ComponentName.h"
#include "base/Exception.h"

std::map<NuTo::eVisualizeWhat, std::string> GetMap()
{
    std::map<NuTo::eVisualizeWhat, std::string> compMap = {
            {NuTo::eVisualizeWhat::ACCELERATION, "Accelerations"},
            {NuTo::eVisualizeWhat::ANGULAR_ACCELERATION, "AngularAccelerations"},
            {NuTo::eVisualizeWhat::ANGULAR_VELOCITY, "AngularVelocities"},
            {NuTo::eVisualizeWhat::BOND_STRESS, "BondStress"},
            {NuTo::eVisualizeWhat::DAMAGE, "Damage"},
            {NuTo::eVisualizeWhat::CRACK_PHASE_FIELD, "CrackPhaseField"},
            {NuTo::eVisualizeWhat::DISPLACEMENTS, "Displacements"},
            {NuTo::eVisualizeWhat::ENGINEERING_PLASTIC_STRAIN, "EngineeringPlasticStrain"},
            {NuTo::eVisualizeWhat::ENGINEERING_STRAIN, "EngineeringStrain"},
            {NuTo::eVisualizeWhat::ENGINEERING_STRESS, "EngineeringStress"},
            {NuTo::eVisualizeWhat::HEAT_FLUX, "HeatFlux"},
            {NuTo::eVisualizeWhat::LATTICE_STRAIN, "LatticeStrain"},
            {NuTo::eVisualizeWhat::LATTICE_STRESS, "LatticeStress"},
            {NuTo::eVisualizeWhat::LATTICE_PLASTIC_STRAIN, "LatticePlasticStrain"},
            {NuTo::eVisualizeWhat::LOCAL_EQ_STRAIN, "LocalEqStrain"},
            {NuTo::eVisualizeWhat::NONLOCAL_EQ_STRAIN, "NonlocalEqStrain"},
            {NuTo::eVisualizeWhat::PARTICLE_RADIUS, "ParticleRadius"},
            {NuTo::eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS, "PrincipalEngineeringStress"},
            {NuTo::eVisualizeWhat::RELATIVE_HUMIDITY, "RelativeHumidity"},
            {NuTo::eVisualizeWhat::ROTATION, "Rotations"},
            {NuTo::eVisualizeWhat::SHRINKAGE_STRAIN, "ShrinkageStrains"},
            {NuTo::eVisualizeWhat::SLIP, "Slip"},
            {NuTo::eVisualizeWhat::TEMPERATURE, "Temperature"},
            {NuTo::eVisualizeWhat::THERMAL_STRAIN, "ThermalStrain"},
            {NuTo::eVisualizeWhat::TOTAL_INELASTIC_EQ_STRAIN, "TotalInelasticEqStrain"},
            {NuTo::eVisualizeWhat::VELOCITY, "Velocities"},
            {NuTo::eVisualizeWhat::WATER_VOLUME_FRACTION, "WaterVolumeFraction"},
            {NuTo::eVisualizeWhat::ELECTRIC_POTENTIAL, "ElectricPotential"},
            {NuTo::eVisualizeWhat::ELECTRIC_FIELD, "ElectricField"},
            {NuTo::eVisualizeWhat::ELECTRIC_DISPLACEMENT, "ElectricDisplacement"}};
    return compMap;
}

std::string NuTo::GetComponentName(NuTo::eVisualizeWhat component)
{
    try
    {
        return GetMap().at(component);
    }
    catch (const std::out_of_range& e)
    {
        throw Exception(__PRETTY_FUNCTION__, "Enum undefined or not implemented.");
    }
}

NuTo::eVisualizeWhat NuTo::GetComponentEnum(std::string component)
{
    const auto compMap = GetMap();
    for (auto entry : compMap)
        if (entry.second == component)
            return entry.first;

    throw Exception(__PRETTY_FUNCTION__, component + " has no enum equivalent or is not implemented.");
}
