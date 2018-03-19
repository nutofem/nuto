#pragma once

namespace NuTo
{
namespace Constitutive
{

//! @brief interface for common damage laws
class DamageLaw
{
public:
    //! @brief calculates the damage for a given history variable kappa
    //! @param kappa history variable
    //! @return damage
    virtual double Damage(double kappa) const = 0;

    //! @brief calculates the derivative of the damage with respect to the history variable kappa
    //! @param kappa history variable
    //! @return damage derivative
    virtual double Derivative(double kappa) const = 0;
};
} /* Constitutive */
} /* NuTo */
