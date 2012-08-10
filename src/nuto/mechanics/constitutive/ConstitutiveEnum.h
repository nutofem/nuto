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
    LINEAR_ELASTIC,            //!< linear-elastic behavior
    LINEAR_ELASTIC_ENGINEERING_STRESS,            //!< linear-elastic behavior
    MISES_PLASTICITY_ENGINEERING_STRESS,          //!< mises plasticity with isotropic and kinematic hardening
    NONLOCAL_DAMAGE_PLASTICITY_ENGINEERING_STRESS,//!< nonlocal damage model with plasticity in the effective stress space
    MULTISCALE,                //!< multiscale model, where the average stress is calculated from a full fine scale model
    LATTICE_CONCRETE,          //!< material law for lattice model
    LINEAR_HEAT_FLUX           //!< material law for lattice model
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

enum eInput
{
	DEFORMATION_GRADIENT_1D,           //!<
	DEFORMATION_GRADIENT_2D,           //!<
	DEFORMATION_GRADIENT_3D,           //!<
	TEMPERATURE,                       //!<
	TEMPERATURE_GRADIENT_1D,           //!<
	TEMPERATURE_GRADIENT_2D,           //!<
	TEMPERATURE_GRADIENT_3D,           //!<
};

static inline std::string InputToString ( const eInput& e )
{
	const std::map< eInput, std::string > lut =
    boost::assign::map_list_of(DEFORMATION_GRADIENT_1D, "DEFORMATION_GRADIENT_1D")
                              (DEFORMATION_GRADIENT_2D, "DEFORMATION_GRADIENT_2D")
                              (DEFORMATION_GRADIENT_3D,"DEFORMATION_GRADIENT_3D")
                              (TEMPERATURE,"TEMPERATURE")
                              (TEMPERATURE_GRADIENT_1D,"TEMPERATURE_GRADIENT_1D")
                              (TEMPERATURE_GRADIENT_2D,"TEMPERATURE_GRADIENT_2D")
                              (TEMPERATURE_GRADIENT_3D,"TEMPERATURE_GRADIENT_3D");
  std::map< eInput, std::string >::const_iterator it = lut.find( e );
  if ( lut.end() != it )
    return it->second;

  return std::string("undefined");
}


enum eOutput
{
	ENGINEERING_STRAIN_1D,           //!<
	ENGINEERING_STRAIN_2D,           //!<
	ENGINEERING_STRAIN_3D,           //!<
	ENGINEERING_PLASTIC_STRAIN_3D,   //!<
	ENGINEERING_STRESS_1D,           //!<
	ENGINEERING_STRESS_2D,           //!<
	ENGINEERING_STRESS_3D,           //!<
	D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_1D,
	D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_2D,
	D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_3D,
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
	UPDATE_TMP_STATIC_DATA
};

static inline std::string OutputToString( const eOutput& e )
{
	const std::map< eOutput, std::string > lut =
    boost::assign::map_list_of(ENGINEERING_STRAIN_1D, "ENGINEERING_STRAIN_1D" )
                              (ENGINEERING_STRAIN_2D, "ENGINEERING_STRAIN_2D" )
                              (ENGINEERING_STRAIN_3D, "ENGINEERING_STRAIN_3D" )
                              (ENGINEERING_PLASTIC_STRAIN_3D, "ENGINEERING_PLASTIC_STRAIN_3D" )
                              (ENGINEERING_STRESS_1D,"ENGINEERING_STRESS_1D")
                              (ENGINEERING_STRESS_2D,"ENGINEERING_STRESS_2D")
                              (ENGINEERING_STRESS_3D,"ENGINEERING_STRESS_3D")
                              (D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_1D,"D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_1D")
                              (D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_2D,"D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_2D")
                              (D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_3D,"D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_3D")
                              (D_ENGINEERING_STRESS_D_TEMPERATURE_1D,"D_ENGINEERING_STRESS_D_TEMPERATURE_1D")
                              (D_ENGINEERING_STRESS_D_TEMPERATURE_2D,"D_ENGINEERING_STRESS_D_TEMPERATURE_2D")
                              (D_ENGINEERING_STRESS_D_TEMPERATURE_3D,"D_ENGINEERING_STRESS_D_TEMPERATURE_3D")
                              (D_HEAT_FLUX_D_TEMPERATURE_RATE_1D,"D_HEAT_FLUX_D_TEMPERATURE_RATE_1D")
                              (D_HEAT_FLUX_D_TEMPERATURE_RATE_2D,"D_HEAT_FLUX_D_TEMPERATURE_RATE_2D")
                              (D_HEAT_FLUX_D_TEMPERATURE_RATE_3D,"D_HEAT_FLUX_D_TEMPERATURE_RATE_3D")
                              (HEAT_FLUX_1D,"HEAT_FLUX_1D")
                              (HEAT_FLUX_2D,"HEAT_FLUX_2D")
                              (HEAT_FLUX_3D,"HEAT_FLUX_3D")
                              (D_HEAT_FLUX_D_TEMPERATURE_GRADIENT_1D,"D_HEAT_FLUX_D_TEMPERATURE_GRADIENT_1D")
                              (D_HEAT_FLUX_D_TEMPERATURE_GRADIENT_2D,"D_HEAT_FLUX_D_TEMPERATURE_GRADIENT_2D")
                              (D_HEAT_FLUX_D_TEMPERATURE_GRADIENT_3D,"D_HEAT_FLUX_D_TEMPERATURE_GRADIENT_3D")
                              (DAMAGE,"DAMAGE")
                              (UPDATE_STATIC_DATA,"UPDATE_STATIC_DATA")
                              (UPDATE_TMP_STATIC_DATA,"UPDATE_TMP_STATIC_DATA");
  std::map< eOutput, std::string >::const_iterator it = lut.find( e );
  if ( lut.end() != it )
    return it->second;

  return std::string("undefined");
}
}//Constitutive
}//Nuto
#endif /* CONSTITUTIVEENUM_H_ */
