// $ld: $ 
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
    NONLOCAL_DAMAGE_PLASTICITY //!< nonlocal damage model with plasticity in the effective stress space
};
}
}
#endif /* CONSTITUTIVEENUM_H_ */
