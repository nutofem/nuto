#pragma once

#include "base/Exception.h"
#include "mechanics/constitutive/ConstitutivePlaneStateEnum.h"
#include "mechanics/constitutive/EngineeringTangent.h"
#include "mechanics/constitutive/EngineeringStrain.h"
#include "mechanics/constitutive/EngineeringStress.h"

namespace NuTo
{
namespace Laws
{

template <int TDim>
class LinearElasticDamage
{
    using EigenVDim = Eigen::Matrix<double, Voigt::Dim(TDim), 1>;

public:
    LinearElasticDamage(double E, double nu, bool isUnilateral = false,
                        ePlaneState planeState = ePlaneState::PLANE_STRAIN)
        : m3K(3 * E / (3 - 6 * nu))
        , m2G(2 * E / (2 + 2 * nu))
        , mIsUnilateral(isUnilateral)
        , mPlaneState(planeState)
    {
        if (mPlaneState == ePlaneState::PLANE_STRESS)
            throw Exception(__PRETTY_FUNCTION__, "Plane stress is not implemented. ");

        if (TDim == 1)
        {
            m3K = E;
            m2G = E;
        }
    }

    EngineeringStress<TDim> Stress(EngineeringStrain<TDim> strain, double omega) const
    {
        const double eV = 1. / 3. * D().transpose() * strain;
        const EngineeringStrain<TDim> e = strain - D() * eV;

        const double sV = (1. - omega * H(eV)) * m3K * eV;
        return D() * sV + (1. - omega) * m2G * PinvDiag().asDiagonal() * e;
    }

    EngineeringTangent<TDim> DstressDstrain(EngineeringStrain<TDim> strain, double omega) const
    {
        const double eV = 1. / 3. * D().transpose() * strain;

        const EngineeringTangent<TDim> dSigma_deV = m3K * (1. - omega * H(eV)) * mIv;
        const EngineeringTangent<TDim> dSigma_de = (1. - omega) * m2G * mPinvId;

        return dSigma_deV + dSigma_de;
    }

    EngineeringStress<TDim> DstressDomega(EngineeringStrain<TDim> strain, double) const
    {
        const double eV = 1. / 3. * D().transpose() * strain;
        const EngineeringStrain<TDim> e = strain - D() * eV;
        const double dsigmaV_dOmega = -m3K * eV * H(eV);

        return D() * dsigmaV_dOmega - m2G * PinvDiag().asDiagonal() * e;
    }

    static EigenVDim PinvDiag()
    {
        EigenVDim p = EigenVDim::Constant(0.5);
        p.template head<TDim>() = Eigen::Matrix<double, TDim, 1>::Ones();
        return p;
    }

    static EigenVDim D()
    {
        EigenVDim d = EigenVDim::Zero();
        d.template head<TDim>() = Eigen::Matrix<double, TDim, 1>::Ones();
        return d;
    }

    double H(double eV) const
    {
        if (mIsUnilateral)
            return eV < 0 ? 0 : 1;
        else
            return 1;
    }

public:
    double m3K;
    double m2G;
    bool mIsUnilateral;
    ePlaneState mPlaneState;

    const EngineeringTangent<TDim> mIv = 1. / 3. * D() * D().transpose();
    const EngineeringTangent<TDim> mPinvId = PinvDiag().asDiagonal() * (EngineeringTangent<TDim>::Identity() - mIv);
};

} /* Laws */
} /* NuTo */
