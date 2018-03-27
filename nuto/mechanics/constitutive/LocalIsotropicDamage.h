#pragma once

#include "MechanicsInterface.h"
#include "LinearElasticDamage.h"
#include "ModifiedMisesStrainNorm.h"
#include "damageLaws/DamageLawExponential.h"

namespace NuTo
{
namespace Laws
{
template <int TDim>
class EvolutionImplicit;

/**
  * Local damage law with an isotropic damage variable
  *
  * \f[
  *  \boldsymbol \sigma = \left(1 - \omega(\kappa (\boldsymbol \varepsilon)\right) \boldsymbol
  * \sigma_\text{elastic}(\boldsymbol \varepsilon)
  * \f]
  *
  * following a policy based design (hopefully applied correctly...), where ...
  *
  * @tparam TDamageLaw the damage law provides the \f$\omega(\kappa)\f$. This requires the methods `.Damage(double)` and
  * `.Derivative(double)` to be implemented, e.g.  NuTo::Constitutive::DamageLaw.
  *
  * @tparam TEvolution the evolution equation \f$ \kappa(\boldsymbol \varepsilon) \f$. This requires the methods
  * `.Kappa(strain)` and `.DkappaDstrain(strain)` to be implemented, e.g. NuTo::Laws::EvolutionImplicit.
  *
  * @tparam TDim dimension
  */
template <int TDim, typename TDamageLaw = Constitutive::DamageLawExponential,
          typename TEvolution = EvolutionImplicit<TDim>>
class LocalIsotropicDamage : public MechanicsInterface<TDim>
{
public:
    LocalIsotropicDamage(LinearElasticDamage<TDim> elasticDamage, TDamageLaw damageLaw, TEvolution evolution)
        : mElasticDamage(elasticDamage)
        , mDamageLaw(damageLaw)
        , mEvolution(evolution)
    {
    }

    LocalIsotropicDamage(Material::Softening m, eDamageApplication damageApplication = eDamageApplication::FULL)
        : mElasticDamage(m, damageApplication)
        , mDamageLaw(m)
        , mEvolution(m)
    {
    }

    EngineeringStress<TDim> Stress(EngineeringStrain<TDim> strain, double deltaT, CellIds ids) const override
    {
        double kappa = mEvolution.Kappa(strain, deltaT, ids);
        double omega = mDamageLaw.Damage(kappa);
        return mElasticDamage.Stress(strain, omega);
    }

    EngineeringTangent<TDim> Tangent(EngineeringStrain<TDim> strain, double deltaT, CellIds ids) const override
    {
        double kappa = mEvolution.Kappa(strain, deltaT, ids);
        double omega = mDamageLaw.Damage(kappa);
        double dOmegadKappa = mDamageLaw.Derivative(kappa);

        return mElasticDamage.DstressDstrain(strain, omega) +
               mElasticDamage.DstressDomega(strain, omega) * dOmegadKappa *
                       mEvolution.DkappaDstrain(strain, deltaT, ids);
    }

    void Update(EngineeringStrain<TDim> strain, double deltaT, CellIds ids)
    {
        mEvolution.Update(strain, deltaT, ids);
    }

public:
    // Intentionally made public. I assume that we know what we do. And this avoids lots of code, either getters
    // (possibly nonconst) or forwarding functions like `double Kappa(...), double Damage(...)`
    LinearElasticDamage<TDim> mElasticDamage;
    TDamageLaw mDamageLaw;
    TEvolution mEvolution;
};

/** Explicit evolution equation for the NuTo::LocalIsotropicDamageLaw that implements
 * \f[
 *  \kappa(\boldsymbol \varepsilon) = \max \left(\kappa, \varepsilon_\text{eq}(\boldsymbol \varepsilon) \right)
 * \f]
 *
 * @tparam TDim dimension
 *
 */
template <int TDim>
class EvolutionImplicit
{

public:
    //! Constructor. As this evolution equation requires history data, they are also allocated.
    //! @param strainNorm strain norm, see class documentation
    //! @param numCells number of cells for the history data allocation
    //! @param numIpsPerCell nummer of integraiton points per cell for the history data allocation
    EvolutionImplicit(Constitutive::ModifiedMisesStrainNorm<TDim> strainNorm, size_t numCells = 1,
                      size_t numIpsPerCell = 1)
        : mStrainNorm(strainNorm)
    {
        ResizeHistoryData(numCells, numIpsPerCell);
    }

    //! Constructor. As this evolution equation requires history data, they are also allocated.
    //! @param m softening material parameters
    //! @param numCells number of cells for the history data allocation
    //! @param numIpsPerCell nummer of integraiton points per cell for the history data allocation
    EvolutionImplicit(Material::Softening m, size_t numCells = 1, size_t numIpsPerCell = 1)
        : mStrainNorm(m)
    {
        ResizeHistoryData(numCells, numIpsPerCell);
    }


    double Kappa(EngineeringStrain<TDim> strain, double, CellIds ids) const
    {
        return std::max(mStrainNorm.Value(strain), mKappas(ids.cellId, ids.ipId));
    }

    Eigen::Matrix<double, 1, Voigt::Dim(TDim)> DkappaDstrain(EngineeringStrain<TDim> strain, double, CellIds ids) const
    {
        if (mStrainNorm.Value(strain) >= mKappas(ids.cellId, ids.ipId))
            return mStrainNorm.Derivative(strain).transpose();
        return Eigen::Matrix<double, 1, Voigt::Dim(TDim)>::Zero();
    }

    void Update(EngineeringStrain<TDim> strain, double deltaT, CellIds ids)
    {
        mKappas(ids.cellId, ids.ipId) = Kappa(strain, deltaT, ids);
    }

    void ResizeHistoryData(size_t numCells, size_t numIpsPerCell)
    {
        mKappas.setZero(numCells, numIpsPerCell);
    }

public:
    Constitutive::ModifiedMisesStrainNorm<TDim> mStrainNorm;
    Eigen::MatrixXd mKappas;
};

} /* Laws */
} /* NuTo */
