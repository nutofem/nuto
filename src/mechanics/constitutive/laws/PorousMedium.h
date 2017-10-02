#pragma once

#include <cmath>

namespace NuTo
{

class PorousMedium
{
public:
    //! Porous medium class.
    //! @param porosity Porosity of the medium.
    //! @param a Parameter for Saturation(), common values range from 18.6MPa for ordinary concrete to 46.9MPa for HPC.
    //! @param b Parameter for Saturation(), common values range from 2.27 for ordinary concrete to 2.06 for HPC.
    //! @remark Values
    PorousMedium(double porosity, double a, double b)
        : mPorosity(porosity)
        , mA(a)
        , mB(b)
    {
    }

    double GetPorosity() const
    {
        return mPorosity;
    }

    //! \f[ S_w = \left[1 + \left(\frac{p^c}{a} \right)^{b/b-1} \right]^{-1/b} \f]
    //! @param capillaryPressure Capillary pressure.
    //! @return Saturation of the pores with liquid water.
    double Saturation(const double capillaryPressure) const
    {
        return std::pow(1 + std::pow(capillaryPressure / mA, mB / (mB - 1.0)), -1.0 / mB);
    }

private:
    double mPorosity = 0.0;
    double mA = 0.0;
    double mB = 0.0;
};

} // namespace NuTo
