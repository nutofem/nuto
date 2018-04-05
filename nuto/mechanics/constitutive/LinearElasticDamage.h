#pragma once

#include "nuto/base/Exception.h"
#include "ConstitutivePlaneStateEnum.h"
#include "EngineeringTangent.h"
#include "EngineeringStrain.h"
#include "EngineeringStress.h"
#include "damageLaws/SofteningMaterial.h"

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
        if (TDim == 1)
        {
            m3K = E;
            m2G = E;
        }

        mD.setZero();
        mD.head(TDim) = Eigen::VectorXd::Ones(TDim);

        mD2 = mD;
        if (TDim == 2 && mPlaneState == ePlaneState::PLANE_STRESS)
            mD2 *= (1 + nu / (nu - 1));

        mPinvDiag = EigenVDim::Constant(0.5);
        mPinvDiag.head(TDim) = Eigen::VectorXd::Ones(TDim);

        mIv = 1. / 3. * mD * mD2.transpose();
        mPinvId = mPinvDiag.asDiagonal() * (EngineeringTangent<TDim>::Identity() - mIv);
    }

    LinearElasticDamage(Material::Softening m, eDamageApplication damageApplication = FULL)
        : LinearElasticDamage::LinearElasticDamage(m.E, m.nu, damageApplication, m.planeState)
    {
    }


    EngineeringStress<TDim> Stress(EngineeringStrain<TDim> strain, double omega) const
    {
        const double eV = 1. / 3. * mD2.transpose() * strain;
        const EngineeringStrain<TDim> e = strain - mD * eV;

        const double sV = (1. - omega * H(eV)) * m3K * eV;
        return mD * sV + (1. - omega) * m2G * mPinvDiag.asDiagonal() * e;
    }

    EngineeringTangent<TDim> DstressDstrain(EngineeringStrain<TDim> strain, double omega) const
    {
        const double eV = 1. / 3. * mD2.transpose() * strain;

        const EngineeringTangent<TDim> dSigma_deV = m3K * (1. - omega * H(eV)) * mIv;
        const EngineeringTangent<TDim> dSigma_de = (1. - omega) * m2G * mPinvId;

        return dSigma_deV + dSigma_de;
    }

    EngineeringStress<TDim> DstressDomega(EngineeringStrain<TDim> strain, double) const
    {
        const double eV = 1. / 3. * mD2.transpose() * strain;
        const EngineeringStrain<TDim> e = strain - mD * eV;
        const double dsigmaV_dOmega = -m3K * eV * H(eV);

        return mD * dsigmaV_dOmega - m2G * mPinvDiag.asDiagonal() * e;
    }

    double H(double eV) const
    {
        if (mDamageApplication == UNILATERAL)
            return eV < 0 ? 0 : 1;
        else
            return 1;
    }

public:
    //! three times the bulk modulus K
    double m3K;

    //! two times the shear modulus G
    double m2G;

    eDamageApplication mDamageApplication;
    ePlaneState mPlaneState;

    //! Kronecker delta in engineering (voigt) notation
    EigenVDim mD;

    //! modified Kronecker delta that gives the volumetric strain
    //! @remark this differs from `mD` for 2D PLANE_STRESS conditions
    EigenVDim mD2;

    //! diagonal scaling matrix to account for (voigt) gamma_xy == (tensor) 2 * epsilon_xy
    EigenVDim mPinvDiag;

    //! volumetric projection tensor I_v in engineering (voigt) notation
    EngineeringTangent<TDim> mIv;

    //! product of mPinv and the deviatoric projection tensor I_d in engineering (voigt) notation
    EngineeringTangent<TDim> mPinvId;
};

} /* Laws */
} /* NuTo */
