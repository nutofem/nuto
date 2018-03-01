#pragma once

#include "mechanics/integrands/MomentumBalance.h"
#include "mechanics/constitutive/LinearElastic.h"
#include "mechanics/constitutive/ModifiedMisesStrainNorm.h"

namespace NuTo
{
namespace Integrands
{

//! Implicit gradient enhanced damage model
//! Peerlings RHJ et al.
//! https://dx.doi.org/10.1002/(SICI)1097-0207(19961015)39:19<3391::AID-NME7>3.0.CO;2-D
//! @param TDim global dimension
//! @tparam TDamageLaw damage law that provides .Damage(double) and .Derivative(double)
template <int TDim, typename TDamageLaw>
class GradientDamage
{
public:
    //! ctor
    //! @param disp dof type associated with displacements
    //! @param eeq scalar dof type associated with nonlocal equivalent strains
    //! @param c nonlocal parameter unit length squared
    //! @param linearElasticLaw linear elastic law
    //! @param damageLaw damage law that provides .Damage(double) and .Derivative(double)
    //! @param strainNorm modified mises strain norm
    GradientDamage(DofType disp, ScalarDofType eeq, double c, Laws::LinearElastic<TDim> linearElasticLaw,
                   TDamageLaw damageLaw, Constitutive::ModifiedMisesStrainNorm<TDim> strainNorm)
        : mDisp(disp)
        , mEeq(eeq)
        , mC(c)
        , mElasticLaw(linearElasticLaw)
        , mDamageLaw(damageLaw)
        , mNorm(strainNorm)
    {
    }

    virtual ~GradientDamage() = default;

    DofVector<double> Gradient(const CellIpData& data)
    {
        DofVector<double> gradient;

        double eeq = data.Value(mEeq);
        double omega = mDamageLaw.Damage(Kappa(data));

        auto eeqGradient = data.Apply(mEeq, Nabla::Gradient());
        NuTo::EngineeringStrain<TDim> strain = data.Apply(mDisp, Nabla::Strain());

        NMatrix Neeq = data.N(mEeq);
        NMatrix Beeq = data.B(mEeq, Nabla::Gradient());
        BMatrixStrain Bdisp = data.B(mDisp, Nabla::Strain());

        gradient[mDisp] = Bdisp.transpose() * ((1 - omega) * mElasticLaw.Stress(strain));
        gradient[mEeq] = Neeq.transpose() * (eeq - mNorm.Value(strain)) + Beeq.transpose() * (mC * eeqGradient);

        return gradient;
    }

    DofMatrix<double> Hessian0(const CellIpData& data)
    {
        DofMatrix<double> hessian0;

        double kappa = Kappa(data);
        double omega = mDamageLaw.Damage(kappa);
        double dKappa_dEeq = DkappaDeeq(data);

        NuTo::EngineeringStrain<TDim> strain = data.Apply(mDisp, Nabla::Strain());

        NMatrix Neeq = data.N(mEeq);
        BMatrixGradient Beeq = data.B(mEeq, Nabla::Gradient());
        BMatrixStrain Bdisp = data.B(mDisp, Nabla::Strain());

        hessian0(mDisp, mDisp) = Bdisp.transpose() * ((1 - omega) * mElasticLaw.Tangent(strain)) * Bdisp;
        hessian0(mEeq, mDisp) = -Neeq.transpose() * mNorm.Derivative(strain).transpose() * Bdisp;
        hessian0(mEeq, mEeq) = Neeq.transpose() * Neeq + mC * Beeq.transpose() * Beeq;
        hessian0(mDisp, mEeq) =
                Bdisp.transpose() * ((-mDamageLaw.Derivative(kappa) * dKappa_dEeq) * mElasticLaw.Stress(strain)) * Neeq;

        return hessian0;
    }

    virtual void Update(const CellIpData& data)
    {
        mKappas(data.Ids().cellId, data.Ids().ipId) = Kappa(data);
    }

    virtual double Kappa(const CellIpData& data) const
    {
        return std::max(mKappas(data.Ids().cellId, data.Ids().ipId), data.Value(mEeq));
    }

    virtual double DkappaDeeq(const CellIpData& data) const
    {
        return data.Value(mEeq) >= mKappas(data.Ids().cellId, data.Ids().ipId) ? 1 : 0;
    }

    Eigen::MatrixXd mKappas;

    DofType mDisp;
    ScalarDofType mEeq;
    double mC;
    Laws::LinearElastic<TDim> mElasticLaw;
    TDamageLaw mDamageLaw;
    Constitutive::ModifiedMisesStrainNorm<TDim> mNorm;
};
} /* Integrand */
} /* NuTo */
