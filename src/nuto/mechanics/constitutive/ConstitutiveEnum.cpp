#include "ConstitutiveEnum.h"

#include <algorithm>
#include <boost/assign/list_of.hpp>
#include "nuto/mechanics/MechanicsException.h"

const std::map<NuTo::Constitutive::eConstitutiveType, std::string> NuTo::Constitutive::GetConstitutiveTypeMap()
{
    const std::map<eConstitutiveType, std::string> map =
    {{eConstitutiveType::ADDITIVE_INPUT_EXPLICIT,                                   "ADDITIVE_INPUT_EXPLICIT"},
        {eConstitutiveType::ADDITIVE_OUTPUT,                      "ADDITIVE_OUTPUT"},
        {eConstitutiveType::DAMAGE_VISCO_PLASTICITY_ENGINEERING_STRESS,             "DAMAGE_VISCO_PLASTICITY_ENGINEERING_STRESS"},
        {eConstitutiveType::DAMAGE_VISCO_PLASTICITY_HARDENING_ENGINEERING_STRESS,   "DAMAGE_VISCO_PLASTICITY_HARDENING_ENGINEERING_STRESS"},
        {eConstitutiveType::FIBRE_MATRIX_BOND_STRESS_SLIP,                          "FIBRE_MATRIX_BOND_STRESS_SLIP"},
        {eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS,                     "GRADIENT_DAMAGE_ENGINEERING_STRESS"},
        {eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS_FATIGUE,             "GRADIENT_DAMAGE_ENGINEERING_STRESS_FATIGUE"},
        {eConstitutiveType::GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS,          "GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS"},
        {eConstitutiveType::HEAT_CONDUCTION,                                        "HEAT_CONDUCTION"},
        {eConstitutiveType::LATTICE_CONCRETE,                                       "LATTICE_CONCRETE"},
        {eConstitutiveType::LINEAR_DAMPING_ENGINEERING_STRESS,                      "LINEAR_DAMPING_ENGINEERING_STRESS"},
        {eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS,                      "LINEAR_ELASTIC_ENGINEERING_STRESS"},
        {eConstitutiveType::LINEAR_SPRING,                                          "LINEAR_SPRING"},
        {eConstitutiveType::MISES_PLASTICITY_ENGINEERING_STRESS,                    "MISES_PLASTICITY_ENGINEERING_STRESS"},
        {eConstitutiveType::MOISTURE_TRANSPORT,                                     "MOISTURE_TRANSPORT"},
        {eConstitutiveType::MULTISCALE,                                             "MULTISCALE"},
        {eConstitutiveType::NONLOCAL_DAMAGE_PLASTICITY_ENGINEERING_STRESS,          "NONLOCAL_DAMAGE_PLASTICITY_ENGINEERING_STRESS"},
        {eConstitutiveType::PHASE_FIELD,                                            "PHASE_FIELD"},
        {eConstitutiveType::SHRINKAGE_CAPILLARY_STRESS_BASED,                       "SHRINKAGE_CAPILLARY_STRESS_BASED"},
        {eConstitutiveType::STRAIN_GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS,   "STRAIN_GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS"},
        {eConstitutiveType::THERMAL_STRAINS,                                        "THERMAL_STRAINS"}};
    return map;
}

std::string NuTo::Constitutive::ConstitutiveTypeToString(NuTo::Constitutive::eConstitutiveType rOutput)
{
    try
    {
        return GetConstitutiveTypeMap().find(rOutput)->second;
    }
    catch (const std::out_of_range& e)
    {
        throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Enum undefined or not implemented.");
    }
}

NuTo::Constitutive::eConstitutiveType NuTo::Constitutive::ConstitutiveTypeToEnum(std::string rOutput)
{
    std::transform(rOutput.begin(), rOutput.end(),rOutput.begin(), ::toupper);

    for(auto entry : GetConstitutiveTypeMap())
        if (entry.second == rOutput)
            return entry.first;

    throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Enum undefined or not implemented.");
}

const std::map<NuTo::Constitutive::eDamageLawType, std::string> NuTo::Constitutive::GetDamageLawMap()
{
    const std::map<eDamageLawType, std::string> map =
       {{eDamageLawType::ISOTROPIC_NO_SOFTENING,                    "ISOTROPIC_NO_SOFTENING"},
        {eDamageLawType::ISOTROPIC_LINEAR_SOFTENING,                "ISOTROPIC_LINEAR_SOFTENING"},
        {eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING,           "ISOTROPIC_EXPONENTIAL_SOFTENING"},
        {eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_RES_LOAD,  "ISOTROPIC_EXPONENTIAL_SOFTENING_RES_LOAD"},
        {eDamageLawType::ISOTROPIC_CUBIC_HERMITE,                   "ISOTROPIC_CUBIC_HERMITE"}};
    return map;
}

std::string NuTo::Constitutive::DamageLawToString(NuTo::Constitutive::eDamageLawType rDamageLaw)
{
    try
    {
        return GetDamageLawMap().at(rDamageLaw);
    }
    catch (const std::out_of_range& e)
    {
        throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Enum undefined or not implemented.");
    }
}

NuTo::Constitutive::eDamageLawType NuTo::Constitutive::DamageLawToEnum(std::string rDamageLaw)
{
    std::transform(rDamageLaw.begin(), rDamageLaw.end(),rDamageLaw.begin(), ::toupper);

    for(auto entry : GetDamageLawMap())
        if (entry.second == rDamageLaw)
            return entry.first;
    throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Enum undefined or not implemented.");
}

const std::map<NuTo::Constitutive::eConstitutiveParameter, std::string> NuTo::Constitutive::GetConstitutiveParameterMap()
{
    const std::map<eConstitutiveParameter, std::string> map =
       {{eConstitutiveParameter::BIAXIAL_COMPRESSIVE_STRENGTH,              "BIAXIAL_COMPRESSIVE_STRENGTH"},
        {eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_RH,         "BOUNDARY_DIFFUSION_COEFFICIENT_RH"},
        {eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_WV,         "BOUNDARY_DIFFUSION_COEFFICIENT_WV"},
        {eConstitutiveParameter::COMPRESSIVE_STRENGTH,                      "COMPRESSIVE_STRENGTH"},
        {eConstitutiveParameter::DAMPING_COEFFICIENT,                       "DAMPING_COEFFICIENT"},
        {eConstitutiveParameter::DAMAGE_DISTRIBUTION,                       "DAMAGE_DISTRIBUTION"},
        {eConstitutiveParameter::DAMAGE_LAW,                                "DAMAGE_LAW"},
        {eConstitutiveParameter::DENSITY,                                   "DENSITY"},
        {eConstitutiveParameter::DENSITY_WATER,                             "DENSITY_WATER"},
        {eConstitutiveParameter::DIFFUSION_COEFFICIENT_RH,                  "DIFFUSION_COEFFICIENT_RH"},
        {eConstitutiveParameter::DIFFUSION_COEFFICIENT_WV,                  "DIFFUSION_COEFFICIENT_WV"},
        {eConstitutiveParameter::DIFFUSION_EXPONENT_RH,                     "DIFFUSION_EXPONENT_RH"},
        {eConstitutiveParameter::DIFFUSION_EXPONENT_WV,                     "DIFFUSION_EXPONENT_WV"},
        {eConstitutiveParameter::ENABLE_MODIFIED_TANGENTIAL_STIFFNESS,      "ENABLE_MODIFIED_TANGENTIAL_STIFFNESS"},
        {eConstitutiveParameter::ENABLE_SORPTION_HYSTERESIS,                "ENABLE_SORPTION_HYSTERESIS"},
        {eConstitutiveParameter::FATIGUE_EXTRAPOLATION,                     "FATIGUE_EXTRAPOLATION"},
        {eConstitutiveParameter::FRACTURE_ENERGY,                           "FRACTURE_ENERGY"},
        {eConstitutiveParameter::GRADIENT_CORRECTION_ADSORPTION_DESORPTION, "GRADIENT_CORRECTION_ADSORPTION_DESORPTION"},
        {eConstitutiveParameter::GRADIENT_CORRECTION_DESORPTION_ADSORPTION, "GRADIENT_CORRECTION_DESORPTION_ADSORPTION"},
        {eConstitutiveParameter::HARDENING_EXPONENT,                        "HARDENING_EXPONENT"},
        {eConstitutiveParameter::HARDENING_VALUE,                           "HARDENING_VALUE"},
        {eConstitutiveParameter::HEAT_CAPACITY,                             "HEAT_CAPACITY"},
        {eConstitutiveParameter::INITIAL_HARDENING_MODULUS,                 "INITIAL_HARDENING_MODULUS"},
        {eConstitutiveParameter::INITIAL_YIELD_STRENGTH,                    "INITIAL_YIELD_STRENGTH"},
        {eConstitutiveParameter::MACROSCOPIC_BULK_MODULUS,                  "MACROSCOPIC_BULK_MODULUS"},
        {eConstitutiveParameter::MASS_EXCHANGE_RATE,                        "MASS_EXCHANGE_RATE"},
        {eConstitutiveParameter::NONLOCAL_RADIUS,                           "NONLOCAL_RADIUS"},
        {eConstitutiveParameter::NONLOCAL_RADIUS_PARAMETER,                 "NONLOCAL_RADIUS_PARAMETER"},
        {eConstitutiveParameter::POISSONS_RATIO,                            "POISSONS_RATIO"},
        {eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_ADSORPTION,        "POLYNOMIAL_COEFFICIENTS_ADSORPTION"},
        {eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION,        "POLYNOMIAL_COEFFICIENTS_DESORPTION"},
        {eConstitutiveParameter::PORE_VOLUME_FRACTION,                      "PORE_VOLUME_FRACTION"},
        {eConstitutiveParameter::DENSITY_SATURATED_WATER_VAPOR,             "DENSITY_SATURATED_WATER_VAPOR"},
        {eConstitutiveParameter::SOLID_PHASE_BULK_MODULUS,                  "SOLID_PHASE_BULK_MODULUS"},
        {eConstitutiveParameter::SPRING_STIFFNESS,                          "SPRING_STIFFNESS"},
        {eConstitutiveParameter::TEMPERATURE,                               "TEMPERATURE"},
        {eConstitutiveParameter::TENSILE_STRENGTH,                          "TENSILE_STRENGTH"},
        {eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT,             "THERMAL_EXPANSION_COEFFICIENT"},
        {eConstitutiveParameter::THERMAL_CONDUCTIVITY,                      "THERMAL_CONDUCTIVITY"},
        {eConstitutiveParameter::VISCOPLASTIC_YIELD_SURFACE_OFFSET,         "VISCOPLASTIC_YIELD_SURFACE_OFFSET"},
        {eConstitutiveParameter::VISCOSITY,                                 "VISCOSITY"},
        {eConstitutiveParameter::VISCOSITY_EXPONENT,                        "VISCOSITY_EXPONENT"},
        {eConstitutiveParameter::YOUNGS_MODULUS,                            "YOUNGS_MODULUS"}};

    return map;
}

std::string NuTo::Constitutive::ConstitutiveParameterToString(NuTo::Constitutive::eConstitutiveParameter rParameter)
{
    try
    {
        return GetConstitutiveParameterMap().find(rParameter)->second;
    }
    catch (const std::out_of_range& e)
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Enum undefined or not implemented.");
    }
}

NuTo::Constitutive::eConstitutiveParameter NuTo::Constitutive::ConstitutiveParameterToEnum(std::string rParameter)
{
    std::transform(rParameter.begin(), rParameter.end(), rParameter.begin(), ::toupper);
    auto map = GetConstitutiveParameterMap();
    for(auto entry : map)
        if (entry.second == rParameter)
            return entry.first;

    throw NuTo::MechanicsException(__PRETTY_FUNCTION__, rParameter + " has no enum equivalent or is not implemented.");
}

std::string NuTo::Constitutive::InputToString(const NuTo::Constitutive::eInput &e)
{
    const std::map<eInput, std::string > lut =
            boost::assign::map_list_of(eInput::ENGINEERING_STRAIN, "ENGINEERING_STRAIN")
            (eInput::TEMPERATURE,"TEMPERATURE")
            (eInput::TEMPERATURE_GRADIENT,"TEMPERATURE_GRADIENT")
            (eInput::TEMPERATURE_CHANGE,"TEMPERATURE_CHANGE")
            (eInput::NONLOCAL_EQ_PLASTIC_STRAIN,"NONLOCAL_EQ_PLASTIC_STRAIN")
            (eInput::NONLOCAL_EQ_STRAIN,"NONLOCAL_EQ_STRAIN")
            (eInput::NONLOCAL_TOTAL_STRAIN_1D,"NONLOCAL_TOTAL_STRAIN_1D")
            (eInput::ENGINEERING_STRESS_1D,"ENGINEERING_STRESS_1D")
            (eInput::DEFORMATION_GRADIENT_REAL_1D,"DEFORMATION_GRADIENT_REAL_1D")
            (eInput::NONLOCAL_TOTAL_STRAIN_REAL_1D,"NONLOCAL_TOTAL_STRAIN_REAL_1D")
            (eInput::NONLOCAL_TOTAL_STRAIN_VIRT_1D,"NONLOCAL_TOTAL_STRAIN_VIRT_1D")
            (eInput::RELATIVE_HUMIDITY,"RELATIVE_HUMIDITY")
            (eInput::RELATIVE_HUMIDITY_BOUNDARY,"RELATIVE_HUMIDITY_BOUNDARY")
            (eInput::RELATIVE_HUMIDITY_DT1,"RELATIVE_HUMIDITY_DT1")
            (eInput::RELATIVE_HUMIDITY_GRADIENT,"RELATIVE_HUMIDITY_GRADIENT")
            (eInput::WATER_VOLUME_FRACTION,"WATER_VOLUME_FRACTION")
            (eInput::WATER_VOLUME_FRACTION_BOUNDARY,"WATER_VOLUME_FRACTION_BOUNDARY")
            (eInput::WATER_VOLUME_FRACTION_DT1,"WATER_VOLUME_FRACTION_DT1")
            (eInput::WATER_VOLUME_FRACTION_GRADIENT,"WATER_VOLUME_FRACTION_GRADIENT")
            (eInput::CALCULATE_STATIC_DATA, "CALCULATE_STATIC_DATA")
            (eInput::TIME_STEP,"TIME_STEP");
    std::map< eInput, std::string >::const_iterator it = lut.find( e );
    if ( lut.end() != it )
        return it->second;

    return std::string("undefined");
}

std::string NuTo::Constitutive::OutputToString(const NuTo::Constitutive::eOutput e)
{
    const std::map< eOutput, std::string > lut =
            boost::assign::map_list_of(eOutput::ENGINEERING_STRAIN, "ENGINEERING_STRAIN" )
            (eOutput::D_ENGINEERING_STRAIN_D_RELATIVE_HUMIDITY,"D_ENGINEERING_STRAIN_D_RELATIVE_HUMIDITY")
            (eOutput::D_ENGINEERING_STRAIN_D_WATER_VOLUME_FRACTION,"D_ENGINEERING_STRAIN_D_WATER_VOLUME_FRACTION")
            (eOutput::ENGINEERING_STRAIN_VISUALIZE, "ENGINEERING_STRAIN_VISUALIZE" )
            (eOutput::SHRINKAGE_STRAIN_VISUALIZE, "SHRINKAGE_STRAIN_VISUALIZE" )
            (eOutput::THERMAL_STRAIN, "THERMAL_STRAIN" )
            (eOutput::ENGINEERING_PLASTIC_STRAIN_VISUALIZE, "ENGINEERING_PLASTIC_STRAIN_VISUALIZE" )
            (eOutput::ENGINEERING_VISCOPLASTIC_STRAIN_3D, "ENGINEERING_VISCOPLASTIC_STRAIN_3D" )
            (eOutput::ENGINEERING_TOTAL_INELASTIC_STRAIN_3D, "ENGINEERING_TOTAL_INELASTIC_STRAIN_3D" )
            (eOutput::ENGINEERING_STRESS,"ENGINEERING_STRESS")
            (eOutput::ENGINEERING_STRESS_VISUALIZE,"ENGINEERING_STRESS_VISUALIZE")
            (eOutput::ENGINEERING_STRESS_ELASTIC_1D,"ENGINEERING_STRESS_ELASTIC_1D")
            (eOutput::ENGINEERING_STRESS_ELASTIC_2D,"ENGINEERING_STRESS_ELASTIC_2D")
            (eOutput::ENGINEERING_STRESS_ELASTIC_3D,"ENGINEERING_STRESS_ELASTIC_3D")
            (eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN,"D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN")
            (eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_DT1, "D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_DT1")
            (eOutput::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY,"D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY")
            (eOutput::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION,"D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION")
            (eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_ELASTIC_1D,"D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_ELASTIC_1D")
            (eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_ELASTIC_2D,"D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_ELASTIC_2D")
            (eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_ELASTIC_3D,"D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_ELASTIC_3D")
            (eOutput::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_PLASTIC_STRAIN_1D,"D_ENGINEERING_STRESS_D_NONLOCAL_EQ_PLASTIC_STRAIN_1D")
            (eOutput::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN,"D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN")
            (eOutput::D_ENGINEERING_STRESS_D_NONLOCAL_TOTAL_STRAIN_1D,"D_ENGINEERING_STRESS_D_NONLOCAL_TOTAL_STRAIN_1D")
            (eOutput::D_ENGINEERING_STRESS_D_THERMAL_STRAIN,"D_ENGINEERING_STRESS_D_THERMAL_STRAIN")
            (eOutput::D_ENGINEERING_STRESS_D_TEMPERATURE,"D_ENGINEERING_STRESS_D_TEMPERATURE")
            (eOutput::D_ENGINEERING_STRESS_D_TEMPERATURE_1D,"D_ENGINEERING_STRESS_D_TEMPERATURE_1D")
            (eOutput::D_ENGINEERING_STRESS_D_TEMPERATURE_2D,"D_ENGINEERING_STRESS_D_TEMPERATURE_2D")
            (eOutput::D_ENGINEERING_STRESS_D_TEMPERATURE_3D,"D_ENGINEERING_STRESS_D_TEMPERATURE_3D")
            (eOutput::D_STRAIN_D_TEMPERATURE, "D_STRAIN_D_TEMPERATURE")
            (eOutput::HEAT_FLUX,"HEAT_FLUX")
            (eOutput::HEAT_CHANGE,"HEAT_CHANGE")
            (eOutput::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT,"D_HEAT_FLUX_D_TEMPERATURE_GRADIENT")
            (eOutput::D_HEAT_D_TEMPERATURE,"D_HEAT_D_TEMPERATURE")
            (eOutput::DAMAGE,"DAMAGE")
            (eOutput::EXTRAPOLATION_ERROR,"EXTRAPOLATION_ERROR")
            (eOutput::UPDATE_STATIC_DATA,"UPDATE_STATIC_DATA")
            (eOutput::UPDATE_TMP_STATIC_DATA,"UPDATE_TMP_STATIC_DATA")
            (eOutput::LOCAL_EQ_PLASTIC_STRAIN,"LOCAL_EQ_PLASTIC_STRAIN")
            (eOutput::LOCAL_EQ_TOTAL_INELASTIC_STRAIN,"LOCAL_EQ_TOTAL_INELASTIC_STRAIN")
            (eOutput::LOCAL_EQ_STRAIN,"LOCAL_EQ_STRAIN")
            (eOutput::D_LOCAL_EQ_PLASTIC_STRAIN_D_STRAIN_1D,"D_LOCAL_EQ_PLASTIC_STRAIN_D_STRAIN_1D")
            (eOutput::D_LOCAL_EQ_STRAIN_D_STRAIN,"D_LOCAL_EQ_STRAIN_D_STRAIN")
            (eOutput::D_LOCAL_EQ_STRAIN_XI_D_STRAIN,"D_LOCAL_EQ_STRAIN_XI_D_STRAIN")
            (eOutput::ENGINEERING_STRESS_REAL_1D,"ENGINEERING_STRESS_REAL_1D")
            (eOutput::ENGINEERING_STRAIN_VIRT_1D,"ENGINEERING_STRAIN_VIRT_1D")
            (eOutput::D_ENGINEERING_STRESS_REAL_D_ENGINEERING_STRAIN_REAL_1D,"D_ENGINEERING_STRESS_REAL_D_ENGINEERING_STRAIN_REAL_1D")
            (eOutput::D_ENGINEERING_STRESS_REAL_D_NONLOCAL_TOTAL_STRAIN_REAL_1D,"D_ENGINEERING_STRESS_REAL_D_NONLOCAL_TOTAL_STRAIN_REAL_1D")
            (eOutput::D_ENGINEERING_STRAIN_VIRT_D_STRESS_REAL_1D,"D_ENGINEERING_STRAIN_VIRT_D_STRESS_REAL_1D")
            (eOutput::D_ENGINEERING_STRAIN_VIRT_D_NONLOCAL_TOTAL_STRAIN_VIRT_1D,"D_ENGINEERING_STRAIN_VIRT_D_NONLOCAL_TOTAL_STRAIN_VIRT_1D")
            (eOutput::D_INTERNAL_GRADIENT_RH_D_RH_BB_H0, "D_INTERNAL_GRADIENT_RH_D_RH_BB_H0")
            (eOutput::D_INTERNAL_GRADIENT_RH_D_RH_NN_H0, "D_INTERNAL_GRADIENT_RH_D_RH_NN_H0")
            (eOutput::D_INTERNAL_GRADIENT_RH_D_WV_BN_H0, "D_INTERNAL_GRADIENT_RH_D_WV_BN_H0")
            (eOutput::D_INTERNAL_GRADIENT_RH_D_WV_NN_H0, "D_INTERNAL_GRADIENT_RH_D_WV_NN_H0")
            (eOutput::D_INTERNAL_GRADIENT_WV_D_RH_NN_H0, "D_INTERNAL_GRADIENT_WV_D_RH_NN_H0")
            (eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BB_H0, "D_INTERNAL_GRADIENT_WV_D_WV_BB_H0")
            (eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BN_H0, "D_INTERNAL_GRADIENT_WV_D_WV_BN_H0")
            (eOutput::D_INTERNAL_GRADIENT_WV_D_WV_NN_H0, "D_INTERNAL_GRADIENT_WV_D_WV_NN_H0")
            (eOutput::D_INTERNAL_GRADIENT_RH_D_RH_NN_H1, "D_INTERNAL_GRADIENT_RH_D_RH_NN_H1")
            (eOutput::D_INTERNAL_GRADIENT_RH_D_WV_NN_H1, "D_INTERNAL_GRADIENT_RH_D_WV_NN_H1")
            (eOutput::D_INTERNAL_GRADIENT_WV_D_WV_NN_H1, "D_INTERNAL_GRADIENT_WV_D_WV_NN_H1")
            (eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_B,     "INTERNAL_GRADIENT_RELATIVE_HUMIDITY_B")
            (eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_N,     "INTERNAL_GRADIENT_RELATIVE_HUMIDITY_N")
            (eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_B, "INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_B")
            (eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_N, "INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_N")
            (eOutput::D_INTERNAL_GRADIENT_RH_D_RH_BOUNDARY_NN_H0, "D_INTERNAL_GRADIENT_RH_D_RH_BOUNDARY_NN_H0")
            (eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BOUNDARY_NN_H0, "D_INTERNAL_GRADIENT_WV_D_WV_BOUNDARY_NN_H0")
            (eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_BOUNDARY_N,     "INTERNAL_GRADIENT_RELATIVE_HUMIDITY_BOUNDARY_N")
            (eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_BOUNDARY_N, "INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_BOUNDARY_N")
            (eOutput::NONLOCAL_PARAMETER_XI, "NONLOCAL_PARAMETER_XI");
    std::map< eOutput, std::string >::const_iterator it = lut.find( e );
    if ( lut.end() != it )
        return it->second;

    return std::string("undefined");
}
