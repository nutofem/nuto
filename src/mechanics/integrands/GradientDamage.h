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
    //! @param eeq dof type associated with nonlocal equivalent strains
    //! @param c nonlocal parameter unit length squared
    //! @param linearElasticLaw linear elastic law
    //! @param damageLaw damage law that provides .Damage(double) and .Derivative(double)
    //! @param strainNorm modified mises strain norm
    GradientDamage(DofType disp, DofType eeq, double c, Laws::LinearElastic<TDim> linearElasticLaw,
                   TDamageLaw damageLaw, Constitutive::ModifiedMisesStrainNorm<TDim> strainNorm)
        : mDisp(disp)
        , mEeq(eeq)
        , mC(c)
        , mElasticLaw(linearElasticLaw)
        , mDamageLaw(damageLaw)
        , mNorm(strainNorm)
    {
    }

    DofVector<double> Gradient(const CellData& cellData, const CellIpData& cellIpData)
    {
        DofVector<double> gradient;

        // shape functions and their derivatives
        NMatrix Neeq = cellIpData.GetNMatrix(mEeq);
        BMatrixGradient Beeq = cellIpData.GetBMatrixGradient(mEeq);
        BMatrixStrain Bdisp = cellIpData.GetBMatrixStrain(mDisp);

        // node values
        double eeq = (Neeq * cellData.GetNodeValues(mEeq))[0];
        double kappa = std::max(mKappas(cellData.GetCellId(), cellIpData.GetIpId()), eeq);
        double omega = mDamageLaw.Damage(kappa);

        Eigen::Matrix<double, TDim, 1> eeqGradient = Beeq * cellData.GetNodeValues(mEeq);
        NuTo::EngineeringStrain<TDim> strain = Bdisp * cellData.GetNodeValues(mDisp);

        // build terms
        gradient[mDisp] = Bdisp.transpose() * (1 - omega) * mElasticLaw.Stress(strain, 0, 0);
        gradient[mEeq] = Neeq.transpose() * (eeq - mNorm.Value(strain)) + Beeq.transpose() * mC * eeqGradient;

        return gradient;
    }

    DofMatrix<double> Hessian0(const CellData& cellData, const CellIpData& cellIpData)
    {
        DofMatrix<double> hessian0;

        // shape functions and their derivatives
        NMatrix Neeq = cellIpData.GetNMatrix(mEeq);
        BMatrixGradient Beeq = cellIpData.GetBMatrixGradient(mEeq);
        BMatrixStrain Bdisp = cellIpData.GetBMatrixStrain(mDisp);

        // node values
        double eeq = (Neeq * cellData.GetNodeValues(mEeq))[0];
        NuTo::EngineeringStrain<TDim> strain = Bdisp * cellData.GetNodeValues(mDisp);

        // evaluate new kappa
        double kappa = std::max(mKappas(cellData.GetCellId(), cellIpData.GetIpId()), eeq);
        double omega = mDamageLaw.Damage(kappa);

        double dKappa_dEeq = kappa == eeq ? 1 : 0;

        hessian0(mDisp, mDisp) = Bdisp.transpose() * (1 - omega) * mElasticLaw.Tangent(strain, 0, 0) * Bdisp;

        hessian0(mDisp, mEeq) = Bdisp.transpose() * (-mDamageLaw.Derivative(kappa) * dKappa_dEeq) *
                                mElasticLaw.Stress(strain, 0, 0) * Neeq;

        hessian0(mEeq, mDisp) = -Neeq.transpose() * mNorm.Derivative(strain).transpose() * Bdisp;

        hessian0(mEeq, mEeq) = Neeq.transpose() * Neeq + mC * Beeq.transpose() * Beeq;

        return hessian0;
    }

    void Update(const CellData& cellData, const CellIpData& cellIpData)
    {
        NMatrix Neeq = cellIpData.GetNMatrix(mEeq);
        double eeq = (Neeq * cellData.GetNodeValues(mEeq))[0];

        double& oldKappa = mKappas(cellData.GetCellId(), cellIpData.GetIpId());
        oldKappa = std::max(oldKappa, eeq);
    }

    Eigen::MatrixXd mKappas;

private:
    DofType mDisp;
    DofType mEeq;
    double mC;
    Laws::LinearElastic<TDim> mElasticLaw;
    TDamageLaw mDamageLaw;
    Constitutive::ModifiedMisesStrainNorm<TDim> mNorm;
};
} /* Integrand */
} /* NuTo */
