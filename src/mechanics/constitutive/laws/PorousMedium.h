#pragma once

#include <cmath>
#include <eigen3/Eigen/Core>

namespace NuTo
{

class PorousMedium
{
public:
    //! Porous medium class.
    //! @param porosity Porosity of the medium.
    //! @param a Parameter for Saturation(), common values range from 18.6MPa for ordinary concrete to 46.9MPa for HPC.
    //! @param b Parameter for Saturation(), common values range from 2.27 for ordinary concrete to 2.06 for HPC.
    // clang-format off
    //! @param fs Structure factor, see [Gawin et al., 1999](https://dx.doi.org/10.1002%2F(SICI)1099-1484(199901)4%3A1%3C37%3A%3AAID-CFM58%3E3.0.CO%3B2-S)
    // clang-format on
    PorousMedium(double porosity, double a, double b, double fs)
        : mPorosity(porosity)
        , mA(a)
        , mB(b)
        , mFs(fs)
    {
    }

    double Porosity() const
    {
        return mPorosity;
    }

    //! \f[ S_w = \left[1 + \left(\frac{p^c}{a} \right)^{b/(b-1)} \right]^{-1/b} \f]
    //! @param capillaryPressure Capillary pressure.
    //! @return Saturation of the pores with liquid water.
    double Saturation(const double capillaryPressure) const
    {
        return std::pow(1 + std::pow(capillaryPressure / mA, mB / (mB - 1.0)), -1.0 / mB);
    }

    // the equation is too long; reformating would introduce artifacts in the output
    // clang-format off
    //! \f[ \frac{∂S_w}{∂p^c} = - \frac{\left(\frac{p_{c}}{a}\right)^{b/(b - 1)}}{p_{c} \left(b - 1\right)} \left[\left(\frac{p_{c}}{a}\right)^{b/(b-1)} + 1\right]^{-(b+1)/b} \f]
    // clang-format on
    //! @param capillaryPressure Capillary pressure.
    //! @return Derivative of saturation of the pores with liquid water w.r.t. the capillary pressure.
    double DerivativeSaturation(const double capillaryPressure) const
    {
        const double ratio = capillaryPressure / mA;
        const double exp = mB / (mB - 1.0);
        const double nom = -std::pow(ratio, exp) * std::pow(std::pow(ratio, exp) + 1.0, -(mB + 1.0) / mB);
        const double den = capillaryPressure * (mB - 1.0);
        return nom / den;
    }

    //! Intrinsic permeability in m².
    //! @param gasPressure Gas pressure in MPa.
    //! @remark Source: Gawin et al. "What physical phenomena can be neglected when modeling concrete at high
    //!                 temperature? A compartive study. Part 1", 2011, DOI:
    //!                 [10.1016/j.ijsolstr.2011.03.004](https://dx.doi.org/10.1016/j.ijsolstr.2011.03.004)
    double IntrinsicPermeability(const double gasPressure) const
    {
        // TODO: add effects of dehydration and damage
        const double k0 = 2e-19;
        const double Ap = 0.36848;
        const double pg0 = 0.1; // MPa
        return k0 * std::pow(gasPressure / pg0, Ap);
    }

    double EffectiveDiffusivity(const double capillaryPressure, const double gasPressure,
                                const double temperature) const
    {
        const double T0 = 273.15;
        const double D_v0 = 2.58e-5;
        const double B_v = 1.667;
        const double A_v = 1.0;
        const double p_0 = 0.1;
        const double S_w = Saturation(capillaryPressure);
        return mPorosity * std::pow(1.0 - S_w, A_v) * mFs * D_v0 * std::pow(temperature / T0, B_v) * p_0 / gasPressure;
    }

    double GasRelativePermeability(const double capillaryPressure) const
    {
        // TODO: this is an empirical relation, could be very different
        const double S_w = Saturation(capillaryPressure);
        return 1.0 - S_w;
    }


    double WaterRelativePermeability(const double capillaryPressure) const
    {
        // TODO: this is an empirical relation, could be very different
        const double S_w = Saturation(capillaryPressure);
        return std::pow(S_w, 6);
    }

private:
    double mPorosity = 0.0;
    double mA = 0.0;
    double mB = 0.0;
    double mFs = 0.0;
};

} // namespace NuTo
