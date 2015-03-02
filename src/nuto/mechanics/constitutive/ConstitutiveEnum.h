// $Id$
#ifndef CONSTITUTIVEENUM_H_
#define CONSTITUTIVEENUM_H_

#include <map>
#include <boost/assign/list_of.hpp>

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
    MOISTURE_TRANSPORT                             //!< moisture transport model
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

namespace Input
{
enum eInput
{
    DEFORMATION_GRADIENT_1D,            //!<
    DEFORMATION_GRADIENT_2D,            //!<
    DEFORMATION_GRADIENT_3D,            //!<
    ENGINEERING_STRAIN_1D,              //!<
    ENGINEERING_STRAIN_2D,              //!<
	ENGINEERING_STRAIN_3D,
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
    WATER_PHASE_FRACTION                //!<
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
	ENGINEERING_STRESS_ELASTIC_3D,	 //!<
	D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_1D,
	D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_2D,
	D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_3D,
	D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_ELASTIC_3D,
    D_ENGINEERING_STRESS_D_NONLOCAL_EQ_PLASTIC_STRAIN_1D,
    D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_1D,
    D_ENGINEERING_STRESS_D_NONLOCAL_TOTAL_STRAIN_1D,
	D_ENGINEERING_STRESS_D_TEMPERATURE_1D,
	D_ENGINEERING_STRESS_D_TEMPERATURE_2D,
	D_ENGINEERING_STRESS_D_TEMPERATURE_3D,
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
    ENGINEERING_STRESS_REAL_1D,
	ENGINEERING_STRAIN_VIRT_1D,
	D_ENGINEERING_STRESS_REAL_D_ENGINEERING_STRAIN_REAL_1D,
	D_ENGINEERING_STRESS_REAL_D_NONLOCAL_TOTAL_STRAIN_REAL_1D,
	D_ENGINEERING_STRAIN_VIRT_D_STRESS_REAL_1D,
    D_ENGINEERING_STRAIN_VIRT_D_NONLOCAL_TOTAL_STRAIN_VIRT_1D,
    PHASE_MASS_EXCHANGE_RATE,
    PHASE_MASS_EXCHANGE_RATE_TIMES_EQUILIBRIUM_SORPTION_CURVE,
    VAPOR_PHASE_DIFFUSION_COEFFICIENT,
    VAPOR_PHASE_SATURATION_DENSITY_TIMES_VAPOR_PHASE_VOLUME_FRACTION,
    VAPOR_PHASE_SATURATION_DENSITY_TIMES_RELATIVE_HUMIDITY,
    WATER_PHASE_DIFFUSION_COEFFICIENT,
    WATER_PHASE_DENSITY,
	FATIGUE_SAVE_STATIC_DATA,
	FATIGUE_RESTORE_STATIC_DATA,
	FATIGUE_EXTRAPOLATE_STATIC_DATA,
	VARIABLE_NONLOCAL_RADIUS
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
                              (Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_1D,"D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_1D")
                              (Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_2D,"D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_2D")
                              (Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_3D,"D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_3D")
                              (Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_PLASTIC_STRAIN_1D,"D_ENGINEERING_STRESS_D_NONLOCAL_EQ_PLASTIC_STRAIN_1D")
                              (Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_1D,"D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_1D")
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
                              (Output::ENGINEERING_STRESS_REAL_1D,"ENGINEERING_STRESS_REAL_1D")
                              (Output::ENGINEERING_STRAIN_VIRT_1D,"ENGINEERING_STRAIN_VIRT_1D")
                              (Output::D_ENGINEERING_STRESS_REAL_D_ENGINEERING_STRAIN_REAL_1D,"D_ENGINEERING_STRESS_REAL_D_ENGINEERING_STRAIN_REAL_1D")
                              (Output::D_ENGINEERING_STRESS_REAL_D_NONLOCAL_TOTAL_STRAIN_REAL_1D,"D_ENGINEERING_STRESS_REAL_D_NONLOCAL_TOTAL_STRAIN_REAL_1D")
                              (Output::D_ENGINEERING_STRAIN_VIRT_D_STRESS_REAL_1D,"D_ENGINEERING_STRAIN_VIRT_D_STRESS_REAL_1D")
                              (Output::D_ENGINEERING_STRAIN_VIRT_D_NONLOCAL_TOTAL_STRAIN_VIRT_1D,"D_ENGINEERING_STRAIN_VIRT_D_NONLOCAL_TOTAL_STRAIN_VIRT_1D")
                              (Output::VAPOR_PHASE_DIFFUSION_COEFFICIENT, "VAPOR_PHASE_DIFFUSION_COEFFICIENT")
                              (Output::WATER_PHASE_DIFFUSION_COEFFICIENT, "WATER_PHASE_DIFFUSION_COEFFICIENT")
                              (Output::VARIABLE_NONLOCAL_RADIUS, "VARIABLE_NONLOCAL_RADIUS");

  std::map< Output::eOutput, std::string >::const_iterator it = lut.find( e );
  if ( lut.end() != it )
    return it->second;

  return std::string("undefined");
}
}//Constitutive
}//Nuto
#endif /* CONSTITUTIVEENUM_H_ */
