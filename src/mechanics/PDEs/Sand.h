#pragma once

#include <cmath>
#include <eigen3/Eigen/Core>
#include "PorousMedium.h"

//! Sand as desribed by Lewis & Schrefler, Section 5.7
//! @remark [ISBN 978-0-47-192809-6](https://en.wikipedia.org/wiki/Special:BookSources?isbn=9780471928096)
class Sand : public PorousMedium
{
public:
    double Porosity() const override
    {
        return 0.2975;
    }

    double Saturation(const double capillaryPressure) const override
    {
        if (capillaryPressure < 1e-12)
            return 1.0;
        const double p_c = capillaryPressure * 1e+6;
        const double exponent = 2.4279;
        const double factor = 1.9722e-11;
        return std::max(0.0, 1.0 - factor * std::pow(p_c, exponent));
    }

    double DerivativeSaturation(const double capillaryPressure) const override
    {
        if (capillaryPressure < 1e-12)
            return 0.0;
        const double p_c = capillaryPressure * 1e+6;
        const double exponent = 2.4279;
        const double factor = 1.9722e-11;
        return -exponent * factor * std::pow(p_c, exponent) / p_c;
    }

    double IntrinsicPermeability(const double gasPressure) const override
    {
        return 4.5e-13;
    }


    double EffectiveDiffusivity(const double capillaryPressure, const double gasPressure,
                                const double temperature) const override
    {
        return 0;
    }


    double GasRelativePermeability(const double capillaryPressure) const override
    {
        const double lower_limit = 0.0001;
        const double S = Saturation(capillaryPressure);
        const double S_e = std::max(0.0, (S - 0.2) / 0.8);
        const double K_rg = std::pow(1.0 - S_e, 2.0) * (1.0 - std::pow(S_e, 5.0/3.0));
        return std::max(K_rg, lower_limit);
    }


    double WaterRelativePermeability(const double capillaryPressure) const override
    {
        const double saturation = Saturation(capillaryPressure);
        const double exponent = 1.0121;
        const double factor = 2207;
        return 1.0 - factor * std::pow(1.0 - saturation, exponent);
    }
};
