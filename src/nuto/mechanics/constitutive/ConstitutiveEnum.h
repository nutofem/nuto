#pragma once

#include <map>
#include <algorithm>
#include <boost/assign/list_of.hpp>
#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{

namespace Constitutive
{
enum eConstitutiveType
{
    ADDITIVE_INPUT_EXPLICIT,                                //!< container for multiple constitutive laws linked by a shared input
    ADDITIVE_INPUT_IMPLICIT,                                //!< container for multiple constitutive laws linked by a shared input
    CONSTITUTIVE_LAWS_ADDITIVE_OUTPUT,                      //!< container for multiple constitutive laws linked by addition of their outputs
    DAMAGE_VISCO_PLASTICITY_ENGINEERING_STRESS,             //!< viscoplastic damage model
    DAMAGE_VISCO_PLASTICITY_HARDENING_ENGINEERING_STRESS,   //!< viscoplastic damage model with hardening
    FIBRE_MATRIX_BOND_STRESS_SLIP,                          //!< material model for the matrix-fibre interface
    GRADIENT_DAMAGE_ENGINEERING_STRESS,                     //!< gradient damage model
    GRADIENT_DAMAGE_ENGINEERING_STRESS_FATIGUE,             //!< gradient damage model for fatigued
    GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS,          //!< gradient damage plasticity model
    HEAT_CONDUCTION,                                        //!< Heat conduction
    LATTICE_CONCRETE,                                       //!< material law for lattice model
    LINEAR_ELASTIC_ENGINEERING_STRESS,                      //!< linear-elastic behavior
    LINEAR_SPRING,                                          //!< linear spring model
    MISES_PLASTICITY_ENGINEERING_STRESS,                    //!< mises plasticity with isotropic and kinematic hardening
    MOISTURE_TRANSPORT,                                     //!< moisture transport model
    MULTISCALE,                                             //!< multiscale model, where the average stress is calculated from a full fine scale model
    NONLOCAL_DAMAGE_PLASTICITY_ENGINEERING_STRESS,          //!< nonlocal damage model with plasticity in the effective stress space
    PHASE_FIELD,                                            //!< phase field model
    SHRINKAGE_CAPILLARY_STRAIN_BASED,                       //!< strain based drying shrinkage - capillary term
    SHRINKAGE_CAPILLARY_STRESS_BASED,                       //!< stress based drying shrinkage - capillary term
    STRAIN_GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS,   //!< strain gradient damage plasticity model (damage and plasticity are function of nonlocal total strain)
    THERMAL_STRAINS                                         //!< strain induced by temperature change
};

static inline std::map<eConstitutiveType, std::string> GetConstitutiveTypeMap()
{
    std::map<eConstitutiveType, std::string> map;
    map[ADDITIVE_INPUT_EXPLICIT]                                = "ADDITIVE_INPUT_EXPLICIT";
    map[CONSTITUTIVE_LAWS_ADDITIVE_OUTPUT]                      = "CONSTITUTIVE_LAWS_ADDITIVE_OUTPUT";
    map[DAMAGE_VISCO_PLASTICITY_ENGINEERING_STRESS]             = "DAMAGE_VISCO_PLASTICITY_ENGINEERING_STRESS";
    map[DAMAGE_VISCO_PLASTICITY_HARDENING_ENGINEERING_STRESS]   = "DAMAGE_VISCO_PLASTICITY_HARDENING_ENGINEERING_STRESS";
    map[FIBRE_MATRIX_BOND_STRESS_SLIP]                          = "FIBRE_MATRIX_BOND_STRESS_SLIP";
    map[GRADIENT_DAMAGE_ENGINEERING_STRESS]                     = "GRADIENT_DAMAGE_ENGINEERING_STRESS";
    map[GRADIENT_DAMAGE_ENGINEERING_STRESS_FATIGUE]             = "GRADIENT_DAMAGE_ENGINEERING_STRESS_FATIGUE";
    map[GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS]          = "GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS";
    map[HEAT_CONDUCTION]                                        = "HEAT_CONDUCTION";
    map[LATTICE_CONCRETE]                                       = "LATTICE_CONCRETE";
    map[LINEAR_ELASTIC_ENGINEERING_STRESS]                      = "LINEAR_ELASTIC_ENGINEERING_STRESS";
    map[LINEAR_SPRING]                                          = "LINEAR_SPRING";
    map[MISES_PLASTICITY_ENGINEERING_STRESS]                    = "MISES_PLASTICITY_ENGINEERING_STRESS";
    map[MOISTURE_TRANSPORT]                                     = "MOISTURE_TRANSPORT";
    map[MULTISCALE]                                             = "MULTISCALE";
    map[NONLOCAL_DAMAGE_PLASTICITY_ENGINEERING_STRESS]          = "NONLOCAL_DAMAGE_PLASTICITY_ENGINEERING_STRESS";
    map[PHASE_FIELD]                                            = "PHASE_FIELD";
    map[SHRINKAGE_CAPILLARY_STRESS_BASED]                       = "SHRINKAGE_CAPILLARY_STRESS_BASED";
    map[STRAIN_GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS]   = "STRAIN_GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS";
    map[THERMAL_STRAINS]                                        = "THERMAL_STRAINS";
    return map;
}

static inline std::string ConstitutiveTypeToString(eConstitutiveType rOutput)
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

static inline eConstitutiveType ConstitutiveTypeToEnum(std::string rOutput)
{
    std::transform(rOutput.begin(), rOutput.end(),rOutput.begin(), ::toupper);

    for(auto entry : GetConstitutiveTypeMap())
        if (entry.second == rOutput)
            return entry.first;

    throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Enum undefined or not implemented.");
}


enum class eConstitutiveStaticDataType
{
    MOISTURE_TRANSPORT,
    MULTI_PHYSICS
};



enum eNonlocalDamageYieldSurface
{
	RANKINE_ROUNDED,            //!< Rankine yield surface (rounded in tension
	ONLY_DP,                    //!< only Drucker-Prager yield surface
	COMBINED_ROUNDED,           //!< Drucker-Prager and rounded Rankine
	COMBINED_SHARP,             //!< Drucker-Prager and Rankine (no rounding)
	TWO_DP                      //!< two Drucker Prager cones
};

enum eSolutionPhaseType
{
	HOMOGENIZED_LINEAR_ELASTIC, //!< linear (homogenized tangent stiffness is used to calculate stiffness and stress)
	NONLINEAR_NO_CRACK,         //!< nonlinear with disp boundary conditions, but no crack can localize
	NONLINEAR_CRACKED           //!< nonlinear with crack enrichment
};

enum eDamageLawType
{
    ISOTROPIC_NO_SOFTENING,                         //!< constant post peak behaviour
    ISOTROPIC_LINEAR_SOFTENING,                     //!< linear
    ISOTROPIC_EXPONENTIAL_SOFTENING,                //!< exponential
    ISOTROPIC_EXPONENTIAL_SOFTENING_RES_LOAD,       //!< exponential with residual loading capacity
    ISOTROPIC_CUBIC_HERMITE                         //!< cubic hermite h00
};

enum class ePhaseFieldEnergyDecomposition
{
    ISOTROPIC,                                              //!< isotropic degradation
    ANISOTROPIC_SPECTRAL_DECOMPOSITION,                     //!< spectral decompostion of the elastic strain tensor proposed by Miehe et al.
};

static inline std::map<eDamageLawType, std::string> GetDamageLawMap()
{
    std::map<eDamageLawType, std::string> map;
    map[ISOTROPIC_NO_SOFTENING]                     = "ISOTROPIC_NO_SOFTENING";
    map[ISOTROPIC_LINEAR_SOFTENING]                 = "ISOTROPIC_LINEAR_SOFTENING";
    map[ISOTROPIC_EXPONENTIAL_SOFTENING]            = "ISOTROPIC_EXPONENTIAL_SOFTENING";
    map[ISOTROPIC_EXPONENTIAL_SOFTENING_RES_LOAD]   = "ISOTROPIC_EXPONENTIAL_SOFTENING_RES_LOAD";
    map[ISOTROPIC_CUBIC_HERMITE]                    = "ISOTROPIC_CUBIC_HERMITE";
    return map;
}
static inline std::string DamageLawToString(eDamageLawType rDamageLaw)
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

static inline eDamageLawType DamageLawToEnum(std::string rDamageLaw)
{
    std::transform(rDamageLaw.begin(), rDamageLaw.end(),rDamageLaw.begin(), ::toupper);

    for(auto entry : GetDamageLawMap())
        if (entry.second == rDamageLaw)
            return entry.first;
    throw NuTo::MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Enum undefined or not implemented.");
}

enum class eConstitutiveParameter
{
    ALPHA,                                      //!<
    ARTIFICIAL_VISCOSITY,                       //!<
    BIAXIAL_COMPRESSIVE_STRENGTH,               //!<
    BOUNDARY_DIFFUSION_COEFFICIENT_RH,          //!<
    BOUNDARY_DIFFUSION_COEFFICIENT_WV,          //!<
    COMPRESSIVE_STRENGTH,                       //!<
    DAMAGE_DISTRIBUTION,                        //!<
    DAMAGE_LAW,                                 //!<
    DENSITY,                                    //!<
    DENSITY_WATER,                              //!<
    DIFFUSION_COEFFICIENT_RH,                   //!<
    DIFFUSION_COEFFICIENT_WV,                   //!<
    DIFFUSION_EXPONENT_RH,                      //!<
    DIFFUSION_EXPONENT_WV,                      //!<
    ENABLE_MODIFIED_TANGENTIAL_STIFFNESS,       //!<
    ENABLE_SORPTION_HYSTERESIS,                 //!<
    FATIGUE_EXTRAPOLATION,                      //!<
    FRACTURE_ENERGY,                            //!<
    GRADIENT_CORRECTION_ADSORPTION_DESORPTION,  //!<
    GRADIENT_CORRECTION_DESORPTION_ADSORPTION,  //!<
    HARDENING_EXPONENT,                         //!<
    HARDENING_VALUE,                            //!<
    HEAT_CAPACITY,                              //!< specific heat capacity \f$c_T\f$
    INITIAL_HARDENING_MODULUS,                  //!<
    INITIAL_YIELD_STRENGTH,                     //!<
    LENGTH_SCALE_PARAMETER,                     //!<
    MACROSCOPIC_BULK_MODULUS,                   //!<
    MASS_EXCHANGE_RATE,                         //!<
    MAX_BOND_STRESS,                            //!<
    NONLOCAL_RADIUS,                            //!<
    NONLOCAL_RADIUS_PARAMETER,                  //!<
    NORMAL_STIFFNESS,                           //!<
    POISSONS_RATIO,                             //!<
    POLYNOMIAL_COEFFICIENTS_ADSORPTION,         //!<
    POLYNOMIAL_COEFFICIENTS_DESORPTION,         //!<
    PORE_VOLUME_FRACTION,                       //!<
    RESIDUAL_BOND_STRESS,                       //!<
    DENSITY_SATURATED_WATER_VAPOR,              //!<
    SLIP_AT_MAX_BOND_STRESS,                    //!<
    SLIP_AT_RESIDUAL_BOND_STRESS,               //!<
    SOLID_PHASE_BULK_MODULUS,                   //!<
    SPRING_STIFFNESS,                           //!<
    SPRING_DIRECTION,                           //!<
    TENSILE_STRENGTH,                           //!<
    TEMPERATURE,                                //!<
    THERMAL_EXPANSION_COEFFICIENT,              //!<
    THERMAL_CONDUCTIVITY,                       //!< \f$k \text{ in } \mathbf{q} = - k \nabla T \f$
    VISCOPLASTIC_YIELD_SURFACE_OFFSET,          //!<
    VISCOSITY,                                  //!<
    VISCOSITY_EXPONENT,                         //!<
    YOUNGS_MODULUS                              //!<

};

static inline std::map<eConstitutiveParameter, std::string> GetConstitutiveParameterMap()
{
    std::map<eConstitutiveParameter, std::string> map;

    map[eConstitutiveParameter::BIAXIAL_COMPRESSIVE_STRENGTH]               = "BIAXIAL_COMPRESSIVE_STRENGTH";
    map[eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_RH]          = "BOUNDARY_DIFFUSION_COEFFICIENT_RH";
    map[eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_WV]          = "BOUNDARY_DIFFUSION_COEFFICIENT_WV";
    map[eConstitutiveParameter::COMPRESSIVE_STRENGTH]                       = "COMPRESSIVE_STRENGTH";
    map[eConstitutiveParameter::DAMAGE_DISTRIBUTION]                        = "DAMAGE_DISTRIBUTION";
    map[eConstitutiveParameter::DAMAGE_LAW]                                 = "DAMAGE_LAW";
    map[eConstitutiveParameter::DENSITY]                                    = "DENSITY";
    map[eConstitutiveParameter::DENSITY_WATER]                              = "DENSITY_WATER";
    map[eConstitutiveParameter::DIFFUSION_COEFFICIENT_RH]                   = "DIFFUSION_COEFFICIENT_RH";
    map[eConstitutiveParameter::DIFFUSION_COEFFICIENT_WV]                   = "DIFFUSION_COEFFICIENT_WV";
    map[eConstitutiveParameter::DIFFUSION_EXPONENT_RH]                      = "DIFFUSION_EXPONENT_RH";
    map[eConstitutiveParameter::DIFFUSION_EXPONENT_WV]                      = "DIFFUSION_EXPONENT_WV";
    map[eConstitutiveParameter::ENABLE_MODIFIED_TANGENTIAL_STIFFNESS]       = "ENABLE_MODIFIED_TANGENTIAL_STIFFNESS";
    map[eConstitutiveParameter::ENABLE_SORPTION_HYSTERESIS]                 = "ENABLE_SORPTION_HYSTERESIS";
    map[eConstitutiveParameter::FATIGUE_EXTRAPOLATION]                      = "FATIGUE_EXTRAPOLATION";
    map[eConstitutiveParameter::FRACTURE_ENERGY]                            = "FRACTURE_ENERGY";
    map[eConstitutiveParameter::GRADIENT_CORRECTION_ADSORPTION_DESORPTION]  = "GRADIENT_CORRECTION_ADSORPTION_DESORPTION";
    map[eConstitutiveParameter::GRADIENT_CORRECTION_DESORPTION_ADSORPTION]  = "GRADIENT_CORRECTION_DESORPTION_ADSORPTION";
    map[eConstitutiveParameter::HARDENING_EXPONENT]                         = "HARDENING_EXPONENT";
    map[eConstitutiveParameter::HARDENING_VALUE]                            = "HARDENING_VALUE";
    map[eConstitutiveParameter::HEAT_CAPACITY]                              = "HEAT_CAPACITY";
    map[eConstitutiveParameter::INITIAL_HARDENING_MODULUS]                  = "INITIAL_HARDENING_MODULUS";
    map[eConstitutiveParameter::INITIAL_YIELD_STRENGTH]                     = "INITIAL_YIELD_STRENGTH";
    map[eConstitutiveParameter::MACROSCOPIC_BULK_MODULUS]                   = "MACROSCOPIC_BULK_MODULUS";
    map[eConstitutiveParameter::MASS_EXCHANGE_RATE]                         = "MASS_EXCHANGE_RATE";
    map[eConstitutiveParameter::NONLOCAL_RADIUS]                            = "NONLOCAL_RADIUS";
    map[eConstitutiveParameter::NONLOCAL_RADIUS_PARAMETER]                  = "NONLOCAL_RADIUS_PARAMETER";
    map[eConstitutiveParameter::POISSONS_RATIO]                             = "POISSONS_RATIO";
    map[eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_ADSORPTION]         = "POLYNOMIAL_COEFFICIENTS_ADSORPTION";
    map[eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION]         = "POLYNOMIAL_COEFFICIENTS_DESORPTION";
    map[eConstitutiveParameter::PORE_VOLUME_FRACTION]                       = "PORE_VOLUME_FRACTION";
    map[eConstitutiveParameter::DENSITY_SATURATED_WATER_VAPOR]              = "DENSITY_SATURATED_WATER_VAPOR";
    map[eConstitutiveParameter::SOLID_PHASE_BULK_MODULUS]                   = "SOLID_PHASE_BULK_MODULUS";
    map[eConstitutiveParameter::SPRING_STIFFNESS]                           = "SPRING_STIFFNESS";
    map[eConstitutiveParameter::TEMPERATURE]                                = "TEMPERATURE";
    map[eConstitutiveParameter::TENSILE_STRENGTH]                           = "TENSILE_STRENGTH";
    map[eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT]              = "THERMAL_EXPANSION_COEFFICIENT";
    map[eConstitutiveParameter::THERMAL_CONDUCTIVITY]                       = "THERMAL_CONDUCTIVITY";
    map[eConstitutiveParameter::VISCOPLASTIC_YIELD_SURFACE_OFFSET]          = "VISCOPLASTIC_YIELD_SURFACE_OFFSET";
    map[eConstitutiveParameter::VISCOSITY]                                  = "VISCOSITY";
    map[eConstitutiveParameter::VISCOSITY_EXPONENT]                         = "VISCOSITY_EXPONENT";
    map[eConstitutiveParameter::YOUNGS_MODULUS]                             = "YOUNGS_MODULUS";

    return map;
}

static inline std::string ConstitutiveParameterToString(eConstitutiveParameter rParameter)
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

static inline eConstitutiveParameter ConstitutiveParameterToEnum(std::string rParameter)
{
    std::transform(rParameter.begin(), rParameter.end(), rParameter.begin(), ::toupper);
    auto map = GetConstitutiveParameterMap();
    for(auto entry : map)
        if (entry.second == rParameter)
            return entry.first;

    throw NuTo::MechanicsException(__PRETTY_FUNCTION__, rParameter + " has no enum equivalent or is not implemented.");
}


namespace Input
{
enum eInput
{
    NONE,                               //!< for constitutive law additive input explicit (AddConstitutiveLaw-function)
    ENGINEERING_STRAIN,                 //!<
    TEMPERATURE,                        //!< Temperature \f$T\f$
    TEMPERATURE_GRADIENT,               //!< Temperature gradient \f$\nabla T\f$
    TEMPERATURE_CHANGE,                 //!< First time derivative \f$\frac{\partial T}{\partial t}\f$
    NONLOCAL_EQ_PLASTIC_STRAIN,         //!<
    NONLOCAL_EQ_STRAIN,                 //!<
    NONLOCAL_TOTAL_STRAIN_1D,           //!<
    ENGINEERING_STRESS_1D,              //!< usually the stress is an output, that 's why the additional input term is required
    DEFORMATION_GRADIENT_REAL_1D,       //!<
    NONLOCAL_TOTAL_STRAIN_REAL_1D,      //!<
    NONLOCAL_TOTAL_STRAIN_VIRT_1D,      //!<
    RELATIVE_HUMIDITY,                  //!<
    RELATIVE_HUMIDITY_BOUNDARY,         //!<
    RELATIVE_HUMIDITY_DT1,              //!< first time derivative
    RELATIVE_HUMIDITY_GRADIENT,         //!<
    WATER_VOLUME_FRACTION,              //!<
    WATER_VOLUME_FRACTION_BOUNDARY,     //!<
    WATER_VOLUME_FRACTION_DT1,          //!< first time derivative
    WATER_VOLUME_FRACTION_GRADIENT,     //!<
    INTERFACE_SLIP,                     //!<
    CRACK_PHASE_FIELD,                  //!<
    ELASTIC_ENERGY_DENSITY,             //!<
    CALCULATE_STATIC_DATA,
    PLANE_STATE,                        //!< for telling the law whether there is plane stress or plane strain
    TIME_STEP
};
}

static inline std::string InputToString ( const Input::eInput& e )
{
    const std::map< Input::eInput, std::string > lut =
    boost::assign::map_list_of(Input::ENGINEERING_STRAIN, "ENGINEERING_STRAIN")
                              (Input::TEMPERATURE,"TEMPERATURE")
                              (Input::TEMPERATURE_GRADIENT,"TEMPERATURE_GRADIENT")
                              (Input::TEMPERATURE_CHANGE,"TEMPERATURE_CHANGE")
                              (Input::NONLOCAL_EQ_PLASTIC_STRAIN,"NONLOCAL_EQ_PLASTIC_STRAIN")
                              (Input::NONLOCAL_EQ_STRAIN,"NONLOCAL_EQ_STRAIN")
                              (Input::NONLOCAL_TOTAL_STRAIN_1D,"NONLOCAL_TOTAL_STRAIN_1D")
                              (Input::ENGINEERING_STRESS_1D,"ENGINEERING_STRESS_1D")
                              (Input::DEFORMATION_GRADIENT_REAL_1D,"DEFORMATION_GRADIENT_REAL_1D")
                              (Input::NONLOCAL_TOTAL_STRAIN_REAL_1D,"NONLOCAL_TOTAL_STRAIN_REAL_1D")
                              (Input::NONLOCAL_TOTAL_STRAIN_VIRT_1D,"NONLOCAL_TOTAL_STRAIN_VIRT_1D")
                              (Input::RELATIVE_HUMIDITY,"RELATIVE_HUMIDITY")
                              (Input::RELATIVE_HUMIDITY_BOUNDARY,"RELATIVE_HUMIDITY_BOUNDARY")
                              (Input::RELATIVE_HUMIDITY_DT1,"RELATIVE_HUMIDITY_DT1")
                              (Input::RELATIVE_HUMIDITY_GRADIENT,"RELATIVE_HUMIDITY_GRADIENT")
                              (Input::WATER_VOLUME_FRACTION,"WATER_VOLUME_FRACTION")
                              (Input::WATER_VOLUME_FRACTION_BOUNDARY,"WATER_VOLUME_FRACTION_BOUNDARY")
                              (Input::WATER_VOLUME_FRACTION_DT1,"WATER_VOLUME_FRACTION_DT1")
                              (Input::WATER_VOLUME_FRACTION_GRADIENT,"WATER_VOLUME_FRACTION_GRADIENT")
                              (Input::CALCULATE_STATIC_DATA, "CALCULATE_STATIC_DATA")
                              (Input::TIME_STEP,"TIME_STEP");
 std::map< Input::eInput, std::string >::const_iterator it = lut.find( e );
  if ( lut.end() != it )
    return it->second;

  return std::string("undefined");
}

namespace Output
{
enum eOutput
{
    ENGINEERING_STRESS,
    ENGINEERING_STRESS_VISUALIZE,
    D_ENGINEERING_STRESS_D_PHASE_FIELD,
    ENGINEERING_STRAIN,
    D_ENGINEERING_STRAIN_D_RELATIVE_HUMIDITY,
    D_ENGINEERING_STRAIN_D_WATER_VOLUME_FRACTION,
    ENGINEERING_STRAIN_VISUALIZE,
    SHRINKAGE_STRAIN_VISUALIZE,
    THERMAL_STRAIN, //!< \f$\varepsilon_{th} = \alpha \Delta T I\f$
    D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN,
    D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY,
    D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION,
    ENGINEERING_PLASTIC_STRAIN_VISUALIZE,

    D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN,
    D_LOCAL_EQ_STRAIN_D_STRAIN,
    D_LOCAL_EQ_STRAIN_XI_D_STRAIN,

	ENGINEERING_VISCOPLASTIC_STRAIN_3D,  	//!<
        ENGINEERING_TOTAL_INELASTIC_STRAIN_3D,  //!<
	ENGINEERING_STRESS_ELASTIC_1D,	 //!<
	ENGINEERING_STRESS_ELASTIC_2D,	 //!<
	ENGINEERING_STRESS_ELASTIC_3D,	 //!<
	D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_ELASTIC_1D,
	D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_ELASTIC_2D,
	D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_ELASTIC_3D,
	D_ENGINEERING_STRESS_D_THERMAL_STRAIN,
    D_ENGINEERING_STRESS_D_NONLOCAL_EQ_PLASTIC_STRAIN_1D,
    D_ENGINEERING_STRESS_D_NONLOCAL_TOTAL_STRAIN_1D,
	D_ENGINEERING_STRESS_D_TEMPERATURE,
    D_STRAIN_D_TEMPERATURE,
	D_ENGINEERING_STRESS_D_TEMPERATURE_1D,
	D_ENGINEERING_STRESS_D_TEMPERATURE_2D,
	D_ENGINEERING_STRESS_D_TEMPERATURE_3D,
    HEAT_FLUX,
    HEAT_CHANGE, //!< First time derivative of heat, i.e. 
                 //!< \f$\frac{\partial Q}{\partial t} = \rho c_T \frac{\partial T}{\partial t}\f$
    D_HEAT_FLUX_D_TEMPERATURE_GRADIENT, //!< conductivity matrix
    D_HEAT_D_TEMPERATURE, //!< heat capacity
	DAMAGE,
	EXTRAPOLATION_ERROR,
	UPDATE_STATIC_DATA,
	UPDATE_TMP_STATIC_DATA,
    LOCAL_EQ_PLASTIC_STRAIN,
	LOCAL_EQ_TOTAL_INELASTIC_STRAIN,
    LOCAL_EQ_STRAIN,
    D_LOCAL_EQ_PLASTIC_STRAIN_D_STRAIN_1D,
    ENGINEERING_STRESS_REAL_1D,
	ENGINEERING_STRAIN_VIRT_1D,
	D_ENGINEERING_STRESS_REAL_D_ENGINEERING_STRAIN_REAL_1D,
	D_ENGINEERING_STRESS_REAL_D_NONLOCAL_TOTAL_STRAIN_REAL_1D,
	D_ENGINEERING_STRAIN_VIRT_D_STRESS_REAL_1D,
    D_ENGINEERING_STRAIN_VIRT_D_NONLOCAL_TOTAL_STRAIN_VIRT_1D,
    D_INTERNAL_GRADIENT_RH_D_RH_BB_H0,
    D_INTERNAL_GRADIENT_RH_D_RH_NN_H0,
    D_INTERNAL_GRADIENT_RH_D_WV_BN_H0,
    D_INTERNAL_GRADIENT_RH_D_WV_NN_H0,
    D_INTERNAL_GRADIENT_WV_D_RH_NN_H0,
    D_INTERNAL_GRADIENT_WV_D_WV_BB_H0,
    D_INTERNAL_GRADIENT_WV_D_WV_BN_H0,
    D_INTERNAL_GRADIENT_WV_D_WV_NN_H0,
    D_INTERNAL_GRADIENT_RH_D_RH_NN_H1,
    D_INTERNAL_GRADIENT_RH_D_WV_NN_H1,
    D_INTERNAL_GRADIENT_WV_D_WV_NN_H1,
    INTERNAL_GRADIENT_RELATIVE_HUMIDITY_B,
    INTERNAL_GRADIENT_RELATIVE_HUMIDITY_N,
    INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_B,
    INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_N,
    D_INTERNAL_GRADIENT_RH_D_RH_BOUNDARY_NN_H0,
    D_INTERNAL_GRADIENT_WV_D_WV_BOUNDARY_NN_H0,
    INTERNAL_GRADIENT_RELATIVE_HUMIDITY_BOUNDARY_N,
    INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_BOUNDARY_N,
	FATIGUE_SAVE_STATIC_DATA,
	FATIGUE_RESTORE_STATIC_DATA,
	FATIGUE_EXTRAPOLATE_STATIC_DATA,
	NONLOCAL_PARAMETER_XI,
    INTERFACE_CONSTITUTIVE_MATRIX,
    BOND_STRESS,
    SLIP,
    ELASTIC_ENERGY_DAMAGED_PART,
    D_ELASTIC_ENERGY_DAMAGED_PART_D_ENGINEERING_STRAIN,
};
}

static inline std::string OutputToString( const Output::eOutput& e )
{
	const std::map< Output::eOutput, std::string > lut =
    boost::assign::map_list_of(Output::ENGINEERING_STRAIN, "ENGINEERING_STRAIN" )
                              (Output::D_ENGINEERING_STRAIN_D_RELATIVE_HUMIDITY,"D_ENGINEERING_STRAIN_D_RELATIVE_HUMIDITY")
                              (Output::D_ENGINEERING_STRAIN_D_WATER_VOLUME_FRACTION,"D_ENGINEERING_STRAIN_D_WATER_VOLUME_FRACTION")
                              (Output::ENGINEERING_STRAIN_VISUALIZE, "ENGINEERING_STRAIN_VISUALIZE" )
                              (Output::SHRINKAGE_STRAIN_VISUALIZE, "SHRINKAGE_STRAIN_VISUALIZE" )
                              (Output::THERMAL_STRAIN, "THERMAL_STRAIN" )
                              (Output::ENGINEERING_PLASTIC_STRAIN_VISUALIZE, "ENGINEERING_PLASTIC_STRAIN_VISUALIZE" )
                              (Output::ENGINEERING_VISCOPLASTIC_STRAIN_3D, "ENGINEERING_VISCOPLASTIC_STRAIN_3D" )
                              (Output::ENGINEERING_TOTAL_INELASTIC_STRAIN_3D, "ENGINEERING_TOTAL_INELASTIC_STRAIN_3D" )
                              (Output::ENGINEERING_STRESS,"ENGINEERING_STRESS")
                              (Output::ENGINEERING_STRESS_VISUALIZE,"ENGINEERING_STRESS_VISUALIZE")
                              (Output::ENGINEERING_STRESS_ELASTIC_1D,"ENGINEERING_STRESS_ELASTIC_1D")
                              (Output::ENGINEERING_STRESS_ELASTIC_2D,"ENGINEERING_STRESS_ELASTIC_2D")
                              (Output::ENGINEERING_STRESS_ELASTIC_3D,"ENGINEERING_STRESS_ELASTIC_3D")
                              (Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN,"D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN")
                              (Output::D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY,"D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY")
                              (Output::D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION,"D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION")
                              (Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_ELASTIC_1D,"D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_ELASTIC_1D")
                              (Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_ELASTIC_2D,"D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_ELASTIC_2D")
                              (Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_ELASTIC_3D,"D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_ELASTIC_3D")
                              (Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_PLASTIC_STRAIN_1D,"D_ENGINEERING_STRESS_D_NONLOCAL_EQ_PLASTIC_STRAIN_1D")
                              (Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN,"D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN")
                              (Output::D_ENGINEERING_STRESS_D_NONLOCAL_TOTAL_STRAIN_1D,"D_ENGINEERING_STRESS_D_NONLOCAL_TOTAL_STRAIN_1D")
                              (Output::D_ENGINEERING_STRESS_D_THERMAL_STRAIN,"D_ENGINEERING_STRESS_D_THERMAL_STRAIN")
                              (Output::D_ENGINEERING_STRESS_D_TEMPERATURE,"D_ENGINEERING_STRESS_D_TEMPERATURE")
                              (Output::D_ENGINEERING_STRESS_D_TEMPERATURE_1D,"D_ENGINEERING_STRESS_D_TEMPERATURE_1D")
                              (Output::D_ENGINEERING_STRESS_D_TEMPERATURE_2D,"D_ENGINEERING_STRESS_D_TEMPERATURE_2D")
                              (Output::D_ENGINEERING_STRESS_D_TEMPERATURE_3D,"D_ENGINEERING_STRESS_D_TEMPERATURE_3D")
                              (Output::D_STRAIN_D_TEMPERATURE, "D_STRAIN_D_TEMPERATURE")
                              (Output::HEAT_FLUX,"HEAT_FLUX")
                              (Output::HEAT_CHANGE,"HEAT_CHANGE")
                              (Output::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT,"D_HEAT_FLUX_D_TEMPERATURE_GRADIENT")
                              (Output::D_HEAT_D_TEMPERATURE,"D_HEAT_D_TEMPERATURE")
                              (Output::DAMAGE,"DAMAGE")
                              (Output::EXTRAPOLATION_ERROR,"EXTRAPOLATION_ERROR")
                              (Output::UPDATE_STATIC_DATA,"UPDATE_STATIC_DATA")
                              (Output::UPDATE_TMP_STATIC_DATA,"UPDATE_TMP_STATIC_DATA")
                              (Output::LOCAL_EQ_PLASTIC_STRAIN,"LOCAL_EQ_PLASTIC_STRAIN")
                              (Output::LOCAL_EQ_TOTAL_INELASTIC_STRAIN,"LOCAL_EQ_TOTAL_INELASTIC_STRAIN")
                              (Output::LOCAL_EQ_STRAIN,"LOCAL_EQ_STRAIN")
                              (Output::D_LOCAL_EQ_PLASTIC_STRAIN_D_STRAIN_1D,"D_LOCAL_EQ_PLASTIC_STRAIN_D_STRAIN_1D")
                              (Output::D_LOCAL_EQ_STRAIN_D_STRAIN,"D_LOCAL_EQ_STRAIN_D_STRAIN")
                              (Output::D_LOCAL_EQ_STRAIN_XI_D_STRAIN,"D_LOCAL_EQ_STRAIN_XI_D_STRAIN")
                              (Output::ENGINEERING_STRESS_REAL_1D,"ENGINEERING_STRESS_REAL_1D")
                              (Output::ENGINEERING_STRAIN_VIRT_1D,"ENGINEERING_STRAIN_VIRT_1D")
                              (Output::D_ENGINEERING_STRESS_REAL_D_ENGINEERING_STRAIN_REAL_1D,"D_ENGINEERING_STRESS_REAL_D_ENGINEERING_STRAIN_REAL_1D")
                              (Output::D_ENGINEERING_STRESS_REAL_D_NONLOCAL_TOTAL_STRAIN_REAL_1D,"D_ENGINEERING_STRESS_REAL_D_NONLOCAL_TOTAL_STRAIN_REAL_1D")
                              (Output::D_ENGINEERING_STRAIN_VIRT_D_STRESS_REAL_1D,"D_ENGINEERING_STRAIN_VIRT_D_STRESS_REAL_1D")
                              (Output::D_ENGINEERING_STRAIN_VIRT_D_NONLOCAL_TOTAL_STRAIN_VIRT_1D,"D_ENGINEERING_STRAIN_VIRT_D_NONLOCAL_TOTAL_STRAIN_VIRT_1D")
                              (Output::D_INTERNAL_GRADIENT_RH_D_RH_BB_H0, "D_INTERNAL_GRADIENT_RH_D_RH_BB_H0")
                              (Output::D_INTERNAL_GRADIENT_RH_D_RH_NN_H0, "D_INTERNAL_GRADIENT_RH_D_RH_NN_H0")
                              (Output::D_INTERNAL_GRADIENT_RH_D_WV_BN_H0, "D_INTERNAL_GRADIENT_RH_D_WV_BN_H0")
                              (Output::D_INTERNAL_GRADIENT_RH_D_WV_NN_H0, "D_INTERNAL_GRADIENT_RH_D_WV_NN_H0")
                              (Output::D_INTERNAL_GRADIENT_WV_D_RH_NN_H0, "D_INTERNAL_GRADIENT_WV_D_RH_NN_H0")
                              (Output::D_INTERNAL_GRADIENT_WV_D_WV_BB_H0, "D_INTERNAL_GRADIENT_WV_D_WV_BB_H0")
                              (Output::D_INTERNAL_GRADIENT_WV_D_WV_BN_H0, "D_INTERNAL_GRADIENT_WV_D_WV_BN_H0")
                              (Output::D_INTERNAL_GRADIENT_WV_D_WV_NN_H0, "D_INTERNAL_GRADIENT_WV_D_WV_NN_H0")
                              (Output::D_INTERNAL_GRADIENT_RH_D_RH_NN_H1, "D_INTERNAL_GRADIENT_RH_D_RH_NN_H1")
                              (Output::D_INTERNAL_GRADIENT_RH_D_WV_NN_H1, "D_INTERNAL_GRADIENT_RH_D_WV_NN_H1")
                              (Output::D_INTERNAL_GRADIENT_WV_D_WV_NN_H1, "D_INTERNAL_GRADIENT_WV_D_WV_NN_H1")
                              (Output::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_B,     "INTERNAL_GRADIENT_RELATIVE_HUMIDITY_B")
                              (Output::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_N,     "INTERNAL_GRADIENT_RELATIVE_HUMIDITY_N")
                              (Output::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_B, "INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_B")
                              (Output::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_N, "INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_N")
                              (Output::D_INTERNAL_GRADIENT_RH_D_RH_BOUNDARY_NN_H0, "D_INTERNAL_GRADIENT_RH_D_RH_BOUNDARY_NN_H0")
                              (Output::D_INTERNAL_GRADIENT_WV_D_WV_BOUNDARY_NN_H0, "D_INTERNAL_GRADIENT_WV_D_WV_BOUNDARY_NN_H0")
                              (Output::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_BOUNDARY_N,     "INTERNAL_GRADIENT_RELATIVE_HUMIDITY_BOUNDARY_N")
                              (Output::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_BOUNDARY_N, "INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_BOUNDARY_N")
                              (Output::NONLOCAL_PARAMETER_XI, "NONLOCAL_PARAMETER_XI");

  std::map< Output::eOutput, std::string >::const_iterator it = lut.find( e );
  if ( lut.end() != it )
    return it->second;

  return std::string("undefined");
}
}//Constitutive


}//NuTo
