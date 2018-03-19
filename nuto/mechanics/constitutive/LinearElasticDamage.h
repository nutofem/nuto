#pragma once

#include "nuto/base/Exception.h"
#include "nuto/mechanics/constitutive/ConstitutivePlaneStateEnum.h"
#include "nuto/mechanics/constitutive/EngineeringTangent.h"
#include "nuto/mechanics/constitutive/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/EngineeringStress.h"

namespace NuTo
{
namespace Laws
{
enum eDamageApplication
{
    FULL,
    UNILATERAL
};

//! Applies an isotropic damage variable to the linear elastic hookes law. A hydrostatic/deviatoric split is performed
//! that allows applying damage selectivly:
//! `eDamageApplication::FULL` for the full application
//! `eDamageApplication::UNILATERAL` for excluding the compressive hydrostatic part
//!
//! So it corresponds to `NuTo::Laws::LinearElastic` when you pass `omega = 0`.
template <int TDim>
class LinearElasticDamage
{
    using EigenVDim = Eigen::Matrix<double, Voigt::Dim(TDim), 1>;

public:
    LinearElasticDamage(double E, double nu, eDamageApplication damageApplication = FULL,
                        ePlaneState planeState = ePlaneState::PLANE_STRAIN)
        : m3K(3 * E / (3 - 6 * nu))
        , m2G(2 * E / (2 + 2 * nu))
        , mDamageApplication(damageApplication)
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
        if (mDamageApplication == UNILATERAL)
            return eV < 0 ? 0 : 1;
        else
            return 1;
    }

public:
    double m3K;
    double m2G;
    eDamageApplication mDamageApplication;
    ePlaneState mPlaneState;

    const EngineeringTangent<TDim> mIv = 1. / 3. * D() * D().transpose();
    const EngineeringTangent<TDim> mPinvId = PinvDiag().asDiagonal() * (EngineeringTangent<TDim>::Identity() - mIv);
};

} /* Laws */
} /* NuTo */
