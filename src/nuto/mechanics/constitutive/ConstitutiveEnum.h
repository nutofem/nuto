// $Id$
#ifndef CONSTITUTIVEENUM_H_
#define CONSTITUTIVEENUM_H_

#include <map>
#include <boost/assign/list_of.hpp>
#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{
namespace Constitutive
{
enum eConstitutiveType
{
    LINEAR_ELASTIC,                                //!< linear-elastic behavior
    LINEAR_ELASTIC_ENGINEERING_STRESS,             //!< linear-elastic behavior
    MISES_PLASTICITY_ENGINEERING_STRESS,           //!< mises plasticity with isotropic and kinematic hardening
    NONLOCAL_DAMAGE_PLASTICITY_ENGINEERING_STRESS, //!< nonlocal damage model with plasticity in the effective stress space
    MULTISCALE,                                    //!< multiscale model, where the average stress is calculated from a full fine scale model
    LATTICE_CONCRETE,                              //!< material law for lattice model
    LINEAR_HEAT_FLUX,                              //!< material law for lattice model
    GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS, //!< gradient damage plasticity model
    GRADIENT_DAMAGE_ENGINEERING_STRESS,            //!< gradient damage model
    STRAIN_GRADIENT_DAMAGE_PLASTICITY_ENGINEERING_STRESS, //!< strain gradient damage plasticity model (damage and plasticity are function of nonlocal total strain)
    DAMAGE_VISCO_PLASTICITY_ENGINEERING_STRESS,    //!< viscoplastic damage model
    DAMAGE_VISCO_PLASTICITY_HARDENING_ENGINEERING_STRESS, //!< viscoplastic damage model with hardening
    MOISTURE_TRANSPORT,                             //!< moisture transport model
    LINEAR_SPRING,                                   //!< linear spring model
    MULTI_PHYSICS,                                   //!< multi physics
    INTERFACE_GOODMAN,                              //!< interface model proposed by Goodman et al.
    DRYING_SHRINKAGE                                //!< drying shrinkage
};

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
    ISOTROPIC_EXPONENTIAL_SOFTENING_SMOOTH,         //!< exponential with polynomial smoothing near peak
    ISOTROPIC_CUBIC_HERMITE                         //!< cubic hermite h00
};

enum class eConstitutiveParameter
{
    ALPHA,                                      //!<
    BIAXIAL_COMPRESSIVE_STRENGTH,               //!<
    BOUNDARY_TRANSPORT_CONSTANT_GAS_PHASE,      //!<
    BOUNDARY_TRANSPORT_CONSTANT_WATER_PHASE,    //!<
    COMPRESSIVE_STRENGTH,                       //!<
    DAMAGE_DISTRIBUTION,                        //!<
    DAMAGE_LAW,                                 //!<
    DENSITY,                                    //!<
    DENSITY_WATER_PHASE,                        //!<
    DIFFUSION_CONSTANT_GAS_PHASE,               //!<
    DIFFUSION_CONSTANT_WATER_PHASE,             //!<
    DIFFUSION_EXPONENT_GAS_PHASE,               //!<
    DIFFUSION_EXPONENT_WATER_PHASE,             //!<
    ENABLE_MODIFIED_TANGENTIAL_STIFFNESS,       //!<
    ENABLE_SORPTION_HYSTERESIS,                 //!<
    FATIGUE_EXTRAPOLATION,                      //!<
    FRACTURE_ENERGY,                            //!<
    GRADIENT_CORRECTION_ADSORPTION_DESORPTION,  //!<
    GRADIENT_CORRECTION_DESORPTION_ADSORPTION,  //!<
    HARDENING_EXPONENT,                         //!<
    HARDENING_VALUE,                            //!<
    INITIAL_HARDENING_MODULUS,                  //!<
    INITIAL_YIELD_STRENGTH,                     //!<
    MASS_EXCHANGE_RATE,                         //!<
    MAX_BOND_STRESS,                            //!<
    NONLOCAL_RADIUS,                            //!<
    NONLOCAL_RADIUS_PARAMETER,                  //!<
    NORMAL_STIFFNESS,                           //!<
    POISSONS_RATIO,                             //!<
    POLYNOMIAL_COEFFICIENTS_ADSORPTION,         //!<
    POLYNOMIAL_COEFFICIENTS_DESORPTION,         //!<
    POROSITY,                                   //!<
    RESIDUAL_BOND_STRESS,                       //!<
    SATURATION_DENSITY_GAS_PHASE,               //!<
    SLIP_AT_MAX_BOND_STRESS,                    //!<
    SLIP_AT_RESIDUAL_BOND_STRESS,               //!<
    SPRING_STIFFNESS,                           //!<
    SPRING_DIRECTION,                           //!<
    TENSILE_STRENGTH,                           //!<
    THERMAL_EXPANSION_COEFFICIENT,              //!<
    VISCOPLASTIC_YIELD_SURFACE_OFFSET,          //!<
    VISCOSITY,                                  //!<
    VISCOSITY_EXPONENT,                         //!<
    YOUNGS_MODULUS                              //!<

};

static inline eConstitutiveParameter GetConstitutiveVariableFromString(const std::string& rVariableName)
{
    std::map<std::string,eConstitutiveParameter> MapStringToConstitutiveVariable;

    MapStringToConstitutiveVariable["BIAXIALCOMPRESSIVESTRENGTH"]                   = eConstitutiveParameter::BIAXIAL_COMPRESSIVE_STRENGTH;
    MapStringToConstitutiveVariable["BOUNDARY_TRANSPORT_CONSTANT_GAS_PHASE"]        = eConstitutiveParameter::BOUNDARY_TRANSPORT_CONSTANT_GAS_PHASE;
    MapStringToConstitutiveVariable["BOUNDARY_TRANSPORT_CONSTANT_WATER_PHASE"]      = eConstitutiveParameter::BOUNDARY_TRANSPORT_CONSTANT_WATER_PHASE;
    MapStringToConstitutiveVariable["COMPRESSIVESTRENGTH"]                          = eConstitutiveParameter::COMPRESSIVE_STRENGTH;
    MapStringToConstitutiveVariable["DAMAGEDISTRIBUTION"]                           = eConstitutiveParameter::DAMAGE_DISTRIBUTION;
    MapStringToConstitutiveVariable["DAMAGELAW"]                                    = eConstitutiveParameter::DAMAGE_LAW;
    MapStringToConstitutiveVariable["DENSITY"]                                      = eConstitutiveParameter::DENSITY;
    MapStringToConstitutiveVariable["DENSITY_WATER_PHASE"]                          = eConstitutiveParameter::DENSITY_WATER_PHASE;
    MapStringToConstitutiveVariable["DIFFUSION_CONSTANT_GAS_PHASE"]                 = eConstitutiveParameter::DIFFUSION_CONSTANT_GAS_PHASE;
    MapStringToConstitutiveVariable["DIFFUSION_CONSTANT_WATER_PHASE"]               = eConstitutiveParameter::DIFFUSION_CONSTANT_WATER_PHASE;
    MapStringToConstitutiveVariable["DIFFUSION_EXPONENT_GAS_PHASE"]                 = eConstitutiveParameter::DIFFUSION_EXPONENT_GAS_PHASE;
    MapStringToConstitutiveVariable["DIFFUSION_EXPONENT_WATER_PHASE"]               = eConstitutiveParameter::DIFFUSION_EXPONENT_WATER_PHASE;
    MapStringToConstitutiveVariable["ENABLE_MODIFIED_TANGENTIAL_STIFFNESS"]         = eConstitutiveParameter::ENABLE_MODIFIED_TANGENTIAL_STIFFNESS;
    MapStringToConstitutiveVariable["ENABLE_SORPTION_HYSTERESIS"]                   = eConstitutiveParameter::ENABLE_SORPTION_HYSTERESIS;
    MapStringToConstitutiveVariable["FATIGUEEXTRAPOLATION"]                         = eConstitutiveParameter::FATIGUE_EXTRAPOLATION;
    MapStringToConstitutiveVariable["FRACTUREENERGY"]                               = eConstitutiveParameter::FRACTURE_ENERGY;
    MapStringToConstitutiveVariable["GRADIENT_CORRECTION_ADSORPTION_DESORPTION"]    = eConstitutiveParameter::GRADIENT_CORRECTION_ADSORPTION_DESORPTION;
    MapStringToConstitutiveVariable["GRADIENT_CORRECTION_DESORPTION_ADSORPTION"]    = eConstitutiveParameter::GRADIENT_CORRECTION_DESORPTION_ADSORPTION;
    MapStringToConstitutiveVariable["HARDENINGEXPONENT"]                            = eConstitutiveParameter::HARDENING_EXPONENT;
    MapStringToConstitutiveVariable["HARDENINGVALUE"]                               = eConstitutiveParameter::HARDENING_VALUE;
    MapStringToConstitutiveVariable["INITIALHARDENINGMODULUS"]                      = eConstitutiveParameter::INITIAL_HARDENING_MODULUS;
    MapStringToConstitutiveVariable["INITIALYIELDSTRENGTH"]                         = eConstitutiveParameter::INITIAL_YIELD_STRENGTH;
    MapStringToConstitutiveVariable["MASS_EXCHANGE_RATE"]                           = eConstitutiveParameter::MASS_EXCHANGE_RATE;
    MapStringToConstitutiveVariable["NONLOCALRADIUS"]                               = eConstitutiveParameter::NONLOCAL_RADIUS;
    MapStringToConstitutiveVariable["NONLOCALRADIUSPARAMETER"]                      = eConstitutiveParameter::NONLOCAL_RADIUS_PARAMETER;
    MapStringToConstitutiveVariable["POISSONSRATIO"]                                = eConstitutiveParameter::POISSONS_RATIO;
    MapStringToConstitutiveVariable["POLYNOMIAL_COEFFICIENTS_ADSORPTION"]           = eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_ADSORPTION;
    MapStringToConstitutiveVariable["POLYNOMIAL_COEFFICIENTS_DESORPTION"]           = eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION;
    MapStringToConstitutiveVariable["POROSITY"]                                     = eConstitutiveParameter::POROSITY;
    MapStringToConstitutiveVariable["SATURATION_DENSITY_GAS_PHASE"]                 = eConstitutiveParameter::SATURATION_DENSITY_GAS_PHASE;
    MapStringToConstitutiveVariable["SPRING_STIFFNESS"]                             = eConstitutiveParameter::SPRING_STIFFNESS;
    MapStringToConstitutiveVariable["TENSILESTRENGTH"]                              = eConstitutiveParameter::TENSILE_STRENGTH;
    MapStringToConstitutiveVariable["THERMALEXPANSIONCOEFFICIENT"]                  = eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT;
    MapStringToConstitutiveVariable["VISCOPLASTICYIELDSURFACEOFFSET"]               = eConstitutiveParameter::VISCOPLASTIC_YIELD_SURFACE_OFFSET;
    MapStringToConstitutiveVariable["VISCOSITY"]                                    = eConstitutiveParameter::VISCOSITY;
    MapStringToConstitutiveVariable["VISCOSITYEXPONENT"]                            = eConstitutiveParameter::VISCOSITY_EXPONENT;
    MapStringToConstitutiveVariable["YOUNGSMODULUS"]                                = eConstitutiveParameter::YOUNGS_MODULUS;
    // find element in map

    auto itResult = MapStringToConstitutiveVariable.find(rVariableName);
    if (itResult!=MapStringToConstitutiveVariable.end())
    {
        return itResult->second;
    }
    else
    {
        throw NuTo::MechanicsException("[NuTo::Constitutive::GetConstitutiveVariableFromString] There is no entry for constitutive variable identifier " +
                                       rVariableName + " in file ConstitutiveEnum.h.");
    }
}

namespace Input
{
enum eInput
{
    DEFORMATION_GRADIENT_1D,            //!<
    DEFORMATION_GRADIENT_2D,            //!<
    DEFORMATION_GRADIENT_3D,            //!<
    ENGINEERING_STRAIN_1D,              //!<
    ENGINEERING_STRAIN_2D,              //!<
    ENGINEERING_STRAIN_3D,              //!<
    TEMPERATURE,                        //!<
    TEMPERATURE_GRADIENT_1D,            //!<
    TEMPERATURE_GRADIENT_2D,            //!<
    TEMPERATURE_GRADIENT_3D,            //!<
    NONLOCAL_EQ_PLASTIC_STRAIN,         //!<
    NONLOCAL_EQ_STRAIN,                 //!<
    NONLOCAL_TOTAL_STRAIN_1D,           //!<
    ENGINEERING_STRESS_1D,              //!< usually the stress is an output, that 's why the additional input term is required
    DEFORMATION_GRADIENT_REAL_1D,       //!<
    NONLOCAL_TOTAL_STRAIN_REAL_1D,      //!<
    NONLOCAL_TOTAL_STRAIN_VIRT_1D,      //!<
    RELATIVE_HUMIDITY,                  //!<
    RELATIVE_HUMIDITY_BOUNDARY,         //!<
    RELATIVE_HUMIDITY_D1,               //!< first time derivative
    RELATIVE_HUMIDITY_GRADIENT,         //!<
    WATER_VOLUME_FRACTION,              //!<
    WATER_VOLUME_FRACTION_BOUNDARY,     //!<
    WATER_VOLUME_FRACTION_D1,           //!< first time derivative
    WATER_VOLUME_FRACTION_GRADIENT,     //!<
    INTERFACE_SLIP,                     //!<
};
}

static inline std::string InputToString ( const Input::eInput& e )
{
	const std::map< Input::eInput, std::string > lut =
    boost::assign::map_list_of(Input::DEFORMATION_GRADIENT_1D, "DEFORMATION_GRADIENT_1D")
                              (Input::DEFORMATION_GRADIENT_2D, "DEFORMATION_GRADIENT_2D")
                              (Input::DEFORMATION_GRADIENT_3D,"DEFORMATION_GRADIENT_3D")
                              (Input::ENGINEERING_STRAIN_1D, "ENGINEERING_STRAIN_1D")
                              (Input::ENGINEERING_STRAIN_2D, "ENGINEERING_STRAIN_2D")
                              (Input::ENGINEERING_STRAIN_3D,"ENGINEERING_STRAIN_3D")
                              (Input::TEMPERATURE,"TEMPERATURE")
                              (Input::TEMPERATURE_GRADIENT_1D,"TEMPERATURE_GRADIENT_1D")
                              (Input::TEMPERATURE_GRADIENT_2D,"TEMPERATURE_GRADIENT_2D")
                              (Input::TEMPERATURE_GRADIENT_3D,"TEMPERATURE_GRADIENT_3D")
                              (Input::NONLOCAL_EQ_PLASTIC_STRAIN,"NONLOCAL_EQ_PLASTIC_STRAIN")
                              (Input::NONLOCAL_EQ_STRAIN,"NONLOCAL_EQ_STRAIN")
                              (Input::NONLOCAL_TOTAL_STRAIN_1D,"NONLOCAL_TOTAL_STRAIN_1D")
                              (Input::ENGINEERING_STRESS_1D,"ENGINEERING_STRESS_1D")
                              (Input::DEFORMATION_GRADIENT_REAL_1D,"DEFORMATION_GRADIENT_REAL_1D")
                              (Input::NONLOCAL_TOTAL_STRAIN_REAL_1D,"NONLOCAL_TOTAL_STRAIN_REAL_1D")
                              (Input::NONLOCAL_TOTAL_STRAIN_VIRT_1D,"NONLOCAL_TOTAL_STRAIN_VIRT_1D");
 std::map< Input::eInput, std::string >::const_iterator it = lut.find( e );
  if ( lut.end() != it )
    return it->second;

  return std::string("undefined");
}

namespace Output
{
enum eOutput
{
	ENGINEERING_STRAIN_1D,           //!<
	ENGINEERING_STRAIN_2D,           //!<
	ENGINEERING_STRAIN_3D,           //!<
	ENGINEERING_PLASTIC_STRAIN_3D,   //!<
	ENGINEERING_VISCOPLASTIC_STRAIN_3D,  	//!<
	ENGINEERING_TOTAL_INELASTIC_STRAIN_3D,  //!<
	ENGINEERING_STRESS_1D,           //!<
	ENGINEERING_STRESS_2D,           //!<
	ENGINEERING_STRESS_3D,           //!<
        ENGINEERING_STRESS_2D_PORE_PRESSURE, //!<
	ENGINEERING_STRESS_ELASTIC_3D,	 //!<
	D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_1D,
	D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_2D,
	D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_3D,
	D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_ELASTIC_3D,
    D_ENGINEERING_STRESS_D_NONLOCAL_EQ_PLASTIC_STRAIN_1D,
    D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_1D,
    D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_2D,
    D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_3D,
    D_ENGINEERING_STRESS_D_NONLOCAL_TOTAL_STRAIN_1D,
	D_ENGINEERING_STRESS_D_TEMPERATURE_1D,
	D_ENGINEERING_STRESS_D_TEMPERATURE_2D,
	D_ENGINEERING_STRESS_D_TEMPERATURE_3D,
    D_ENGINEERING_STRESS_D_RELATIVE_HUMIDITY_2D,
    D_ENGINEERING_STRESS_D_WATER_VOLUME_FRACTION_2D,
	D_HEAT_FLUX_D_TEMPERATURE_RATE_1D,
	D_HEAT_FLUX_D_TEMPERATURE_RATE_2D,
	D_HEAT_FLUX_D_TEMPERATURE_RATE_3D,
	HEAT_FLUX_1D,
	HEAT_FLUX_2D,
	HEAT_FLUX_3D,
	D_HEAT_FLUX_D_TEMPERATURE_GRADIENT_1D,
	D_HEAT_FLUX_D_TEMPERATURE_GRADIENT_2D,
	D_HEAT_FLUX_D_TEMPERATURE_GRADIENT_3D,
	DAMAGE,
	UPDATE_STATIC_DATA,
	UPDATE_TMP_STATIC_DATA,
    LOCAL_EQ_PLASTIC_STRAIN,
    LOCAL_EQ_STRAIN,
    D_LOCAL_EQ_PLASTIC_STRAIN_D_STRAIN_1D,
    D_LOCAL_EQ_STRAIN_D_STRAIN_1D,
    D_LOCAL_EQ_STRAIN_D_STRAIN_2D,
    D_LOCAL_EQ_STRAIN_D_STRAIN_3D,
    D_LOCAL_EQ_STRAIN_XI_D_STRAIN_1D,
    D_LOCAL_EQ_STRAIN_XI_D_STRAIN_2D,
    D_LOCAL_EQ_STRAIN_XI_D_STRAIN_3D,
    ENGINEERING_STRESS_REAL_1D,
	ENGINEERING_STRAIN_VIRT_1D,
	D_ENGINEERING_STRESS_REAL_D_ENGINEERING_STRAIN_REAL_1D,
	D_ENGINEERING_STRESS_REAL_D_NONLOCAL_TOTAL_STRAIN_REAL_1D,
	D_ENGINEERING_STRAIN_VIRT_D_STRESS_REAL_1D,
    D_ENGINEERING_STRAIN_VIRT_D_NONLOCAL_TOTAL_STRAIN_VIRT_1D,
    BOUNDARY_SURFACE_RELATIVE_HUMIDIY_TRANSPORT_COEFFICIENT,
    BOUNDARY_SURFACE_WATER_VOLUME_FRACTION_TRANSPORT_COEFFICIENT,
    BOUNDARY_SURFACE_VAPOR_PHASE_RESIDUAL,
    BOUNDARY_SURFACE_WATER_PHASE_RESIDUAL,
    RESIDUAL_WATER_PHASE_B,
    RESIDUAL_WATER_PHASE_N,
    RESIDUAL_VAPOR_PHASE_B,
    RESIDUAL_VAPOR_PHASE_N,
    D_RESIDUAL_RH_D_RH_H0_BB,
    D_RESIDUAL_RH_D_RH_H0_NN,
    D_RESIDUAL_RH_D_WV_H0_BN,
    D_RESIDUAL_RH_D_WV_H0_NN,
    D_RESIDUAL_WV_D_RH_H0_NN,
    D_RESIDUAL_WV_D_WV_H0_BB,
    D_RESIDUAL_WV_D_WV_H0_BN,
    D_RESIDUAL_WV_D_WV_H0_NN,
    D_RESIDUAL_RH_D_RH_H1_NN,
    D_RESIDUAL_RH_D_WV_H1_NN,
    D_RESIDUAL_WV_D_WV_H1_NN,
	FATIGUE_SAVE_STATIC_DATA,
	FATIGUE_RESTORE_STATIC_DATA,
	FATIGUE_EXTRAPOLATE_STATIC_DATA,
    RESIDUAL_NORM_FACTOR_DISPLACEMENTS,
    RESIDUAL_NORM_FACTOR_RELATIVE_HUMIDITY,
    RESIDUAL_NORM_FACTOR_WATER_VOLUME_FRACTION,
	NONLOCAL_PARAMETER_XI,
    INTERFACE_CONSTITUTIVE_MATRIX,
    INTERFACE_STRESSES,
};
}

static inline std::string OutputToString( const Output::eOutput& e )
{
	const std::map< Output::eOutput, std::string > lut =
    boost::assign::map_list_of(Output::ENGINEERING_STRAIN_1D, "ENGINEERING_STRAIN_1D" )
                              (Output::ENGINEERING_STRAIN_2D, "ENGINEERING_STRAIN_2D" )
                              (Output::ENGINEERING_STRAIN_3D, "ENGINEERING_STRAIN_3D" )
                              (Output::ENGINEERING_PLASTIC_STRAIN_3D, "ENGINEERING_PLASTIC_STRAIN_3D" )
                              (Output::ENGINEERING_VISCOPLASTIC_STRAIN_3D, "ENGINEERING_VISCOPLASTIC_STRAIN_3D" )
                              (Output::ENGINEERING_TOTAL_INELASTIC_STRAIN_3D, "ENGINEERING_TOTAL_INELASTIC_STRAIN_3D" )
                              (Output::ENGINEERING_STRESS_1D,"ENGINEERING_STRESS_1D")
                              (Output::ENGINEERING_STRESS_2D,"ENGINEERING_STRESS_2D")
                              (Output::ENGINEERING_STRESS_3D,"ENGINEERING_STRESS_3D")
                              (Output::ENGINEERING_STRESS_2D_PORE_PRESSURE,"ENGINEERING_STRESS_2D_PORE_PRESSURE")
                              (Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_1D,"D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_1D")
                              (Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_2D,"D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_2D")
                              (Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_3D,"D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_3D")
                              (Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_PLASTIC_STRAIN_1D,"D_ENGINEERING_STRESS_D_NONLOCAL_EQ_PLASTIC_STRAIN_1D")
                              (Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_1D,"D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_1D")
                              (Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_2D,"D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_2D")
                              (Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_3D,"D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_3D")
                              (Output::D_ENGINEERING_STRESS_D_NONLOCAL_TOTAL_STRAIN_1D,"D_ENGINEERING_STRESS_D_NONLOCAL_TOTAL_STRAIN_1D")
                              (Output::D_ENGINEERING_STRESS_D_TEMPERATURE_1D,"D_ENGINEERING_STRESS_D_TEMPERATURE_1D")
                              (Output::D_ENGINEERING_STRESS_D_TEMPERATURE_2D,"D_ENGINEERING_STRESS_D_TEMPERATURE_2D")
                              (Output::D_ENGINEERING_STRESS_D_TEMPERATURE_3D,"D_ENGINEERING_STRESS_D_TEMPERATURE_3D")
                              (Output::D_HEAT_FLUX_D_TEMPERATURE_RATE_1D,"D_HEAT_FLUX_D_TEMPERATURE_RATE_1D")
                              (Output::D_HEAT_FLUX_D_TEMPERATURE_RATE_2D,"D_HEAT_FLUX_D_TEMPERATURE_RATE_2D")
                              (Output::D_HEAT_FLUX_D_TEMPERATURE_RATE_3D,"D_HEAT_FLUX_D_TEMPERATURE_RATE_3D")
                              (Output::HEAT_FLUX_1D,"HEAT_FLUX_1D")
                              (Output::HEAT_FLUX_2D,"HEAT_FLUX_2D")
                              (Output::HEAT_FLUX_3D,"HEAT_FLUX_3D")
                              (Output::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT_1D,"D_HEAT_FLUX_D_TEMPERATURE_GRADIENT_1D")
                              (Output::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT_2D,"D_HEAT_FLUX_D_TEMPERATURE_GRADIENT_2D")
                              (Output::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT_3D,"D_HEAT_FLUX_D_TEMPERATURE_GRADIENT_3D")
                              (Output::DAMAGE,"DAMAGE")
                              (Output::UPDATE_STATIC_DATA,"UPDATE_STATIC_DATA")
                              (Output::UPDATE_TMP_STATIC_DATA,"UPDATE_TMP_STATIC_DATA")
                              (Output::LOCAL_EQ_PLASTIC_STRAIN,"LOCAL_EQ_PLASTIC_STRAIN")
                              (Output::LOCAL_EQ_STRAIN,"LOCAL_EQ_STRAIN")
                              (Output::D_LOCAL_EQ_PLASTIC_STRAIN_D_STRAIN_1D,"D_LOCAL_EQ_PLASTIC_STRAIN_D_STRAIN_1D")
                              (Output::D_LOCAL_EQ_STRAIN_D_STRAIN_1D,"D_LOCAL_EQ_STRAIN_D_STRAIN_1D")
                              (Output::D_LOCAL_EQ_STRAIN_D_STRAIN_2D,"D_LOCAL_EQ_STRAIN_D_STRAIN_2D")
                              (Output::D_LOCAL_EQ_STRAIN_D_STRAIN_3D,"D_LOCAL_EQ_STRAIN_D_STRAIN_3D")
                              (Output::D_LOCAL_EQ_STRAIN_XI_D_STRAIN_1D,"D_LOCAL_EQ_STRAIN_XI_D_STRAIN_1D")
                              (Output::D_LOCAL_EQ_STRAIN_XI_D_STRAIN_2D,"D_LOCAL_EQ_STRAIN_XI_D_STRAIN_2D")
                              (Output::D_LOCAL_EQ_STRAIN_XI_D_STRAIN_3D,"D_LOCAL_EQ_STRAIN_XI_D_STRAIN_3D")
                              (Output::ENGINEERING_STRESS_REAL_1D,"ENGINEERING_STRESS_REAL_1D")
                              (Output::ENGINEERING_STRAIN_VIRT_1D,"ENGINEERING_STRAIN_VIRT_1D")
                              (Output::D_ENGINEERING_STRESS_REAL_D_ENGINEERING_STRAIN_REAL_1D,"D_ENGINEERING_STRESS_REAL_D_ENGINEERING_STRAIN_REAL_1D")
                              (Output::D_ENGINEERING_STRESS_REAL_D_NONLOCAL_TOTAL_STRAIN_REAL_1D,"D_ENGINEERING_STRESS_REAL_D_NONLOCAL_TOTAL_STRAIN_REAL_1D")
                              (Output::D_ENGINEERING_STRAIN_VIRT_D_STRESS_REAL_1D,"D_ENGINEERING_STRAIN_VIRT_D_STRESS_REAL_1D")
                              (Output::D_ENGINEERING_STRAIN_VIRT_D_NONLOCAL_TOTAL_STRAIN_VIRT_1D,"D_ENGINEERING_STRAIN_VIRT_D_NONLOCAL_TOTAL_STRAIN_VIRT_1D")
                              (Output::RESIDUAL_WATER_PHASE_B, "RESIDUAL_WATER_PHASE_B")
                              (Output::RESIDUAL_WATER_PHASE_N, "RESIDUAL_WATER_PHASE_N")
                              (Output::RESIDUAL_VAPOR_PHASE_B, "RESIDUAL_VAPOR_PHASE_B")
                              (Output::RESIDUAL_VAPOR_PHASE_N, "RESIDUAL_VAPOR_PHASE_N")
                              (Output::D_RESIDUAL_RH_D_RH_H0_BB, "D_RESIDUAL_RH_D_RH_H0_BB")
                              (Output::D_RESIDUAL_RH_D_RH_H0_NN, "D_RESIDUAL_RH_D_RH_H0_NN")
                              (Output::D_RESIDUAL_RH_D_WV_H0_BN, "D_RESIDUAL_RH_D_WV_H0_BN")
                              (Output::D_RESIDUAL_RH_D_WV_H0_NN, "D_RESIDUAL_RH_D_WV_H0_NN")
                              (Output::D_RESIDUAL_WV_D_RH_H0_NN, "D_RESIDUAL_WV_D_RH_H0_NN")
                              (Output::D_RESIDUAL_WV_D_WV_H0_BB, "D_RESIDUAL_WV_D_WV_H0_BB")
                              (Output::D_RESIDUAL_WV_D_WV_H0_BN, "D_RESIDUAL_WV_D_WV_H0_BN")
                              (Output::D_RESIDUAL_WV_D_WV_H0_NN, "D_RESIDUAL_WV_D_WV_H0_NN")
                              (Output::D_RESIDUAL_RH_D_RH_H1_NN, "D_RESIDUAL_RH_D_RH_H1_NN")
                              (Output::D_RESIDUAL_RH_D_WV_H1_NN, "D_RESIDUAL_RH_D_WV_H1_NN")
                              (Output::D_RESIDUAL_WV_D_WV_H1_NN, "D_RESIDUAL_WV_D_WV_H1_NN")
                              (Output::RESIDUAL_NORM_FACTOR_DISPLACEMENTS, "RESIDUAL_NORM_FACTOR_DISPLACEMENTS")
                              (Output::RESIDUAL_NORM_FACTOR_RELATIVE_HUMIDITY, "RESIDUAL_NORM_FACTOR_RELATIVE_HUMIDITY")
                              (Output::RESIDUAL_NORM_FACTOR_WATER_VOLUME_FRACTION, "RESIDUAL_NORM_FACTOR_WATER_VOLUME_FRACTION")
                              (Output::NONLOCAL_PARAMETER_XI, "NONLOCAL_PARAMETER_XI");

  std::map< Output::eOutput, std::string >::const_iterator it = lut.find( e );
  if ( lut.end() != it )
    return it->second;

  return std::string("undefined");
}
}//Constitutive
}//Nuto
#endif /* CONSTITUTIVEENUM_H_ */
