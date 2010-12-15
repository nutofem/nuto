// $Id$ 
#ifndef CONSTITUTIVEENUM_H_
#define CONSTITUTIVEENUM_H_

namespace NuTo
{
namespace Constitutive
{
enum eConstitutiveType
{
    LINEAR_ELASTIC,            //!< linear-elastic behavior
    MISES_PLASTICITY,          //!< mises plasticity with isotropic and kinematic hardening
    NONLOCAL_DAMAGE_PLASTICITY,//!< nonlocal damage model with plasticity in the effective stress space
    MULTISCALE                 //!< multiscale model, where the average stress is calculated from a full fine scale model
};

enum eNonlocalDamageYieldSurface
{
	RANKINE_ROUNDED,            //!< Rankine yield surface (rounded in tension
	ONLY_DP,                    //!< only Drucker-Prager yield surface
	COMBINED_ROUNDED,           //!< Drucker-Prager and rounded Rankine
	COMBINED_SHARP,             //!< Drucker-Prager and Rankine (no rounding)
	TWO_DP                      //!< two Drucker Prager cones
};

}
}
#endif /* CONSTITUTIVEENUM_H_ */
