#pragma once

namespace NuTo
{
namespace Constitutive
{

//! @brief base class for common damage laws
class DamageLaw
{
public:

    //! @brief ctor with common member variable for all damage laws
    //! @param kappa0 initial kappa, values below kappa0 do not cause damage
    DamageLaw(double kappa0) : mKappa0(kappa0) {}

    //! @brief dtor
    virtual ~DamageLaw() = default;

    //! @brief public nonvirtual interface for the damage calculation
    //! @param kappa history variable
    //! @return 0 for kappa <= kappa0, damage otherwise
    double CalculateDamage(const double kappa) const
    {
        if (kappa <= mKappa0)
            return 0;
        return Damage(kappa);
    }

    //! @brief public nonvirtual interface for the damage derivative calculation
    //! @param kappa history variable
    //! @return 0 for kappa <= kappa0, damage derivative otherwise
    double CalculateDerivative(const double kappa) const
    {
        if (kappa <= mKappa0)
            return 0;
        return Derivative(kappa);
    }

    //! @brief trivial getter for kappa0
    //! @return kappa0 
    double GetKappa0() const
    {
        return mKappa0;
    }

protected:
    //! @brief protected virtual method for the damage calculation
    //!        the case kappa <= kappa0 is already covered in the public interface
    //! @param kappa history variable
    //! @return damage
    virtual double Damage(const double kappa) const = 0;
    
    //! @brief protected virtual method for the damage derivative calculation
    //!        the case kappa <= kappa0 is already covered in the public interface
    //! @param kappa history variable
    //! @return damage derivative
    virtual double Derivative(const double kappa) const = 0;

    //! initial kappa, values below kappa0 do not cause damage
    const double mKappa0;
};
} /* Constitutive */
} /* NuTo */
