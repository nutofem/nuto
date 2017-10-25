#pragma once

#include <cmath>
#include <eigen3/Eigen/Core>

//! Porous medium class.
class PorousMedium
{
public:
    //! Porosity of the medium.
    virtual double Porosity() const = 0;

    //! \f[ S_w = \left[1 + \left(\frac{p^c}{a} \right)^{b/(b-1)} \right]^{-1/b} \f]
    //! @param capillaryPressure Capillary pressure.
    //! @return Saturation of the pores with liquid water.
    virtual double Saturation(const double capillaryPressure) const = 0;

    //! Derivative of the saturation w.r.t. to the capillary pressure.
    //! @param capillaryPressure Capillary pressure.
    //! @return Derivative of saturation of the pores with liquid water w.r.t. the capillary pressure.
    virtual double DerivativeSaturation(const double capillaryPressure) const = 0;

    //! Intrinsic permeability in m².
    //! @param gasPressure Gas pressure in MPa.
    virtual double IntrinsicPermeability(const double gasPressure) const = 0;

    virtual double EffectiveDiffusivity(const double capillaryPressure, const double gasPressure,
                                        const double temperature) const = 0;

    virtual double GasRelativePermeability(const double capillaryPressure) const = 0;

    virtual double WaterRelativePermeability(const double capillaryPressure) const = 0;

    virtual ~PorousMedium()
    {
    }
};
