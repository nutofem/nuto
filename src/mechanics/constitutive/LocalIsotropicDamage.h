#pragma once

#include <vector>
#include <mechanics/constitutive/MechanicsInterface.h>

namespace NuTo
{
namespace Laws
{

/**
  * Local damage law with an isotropic damage variable
  *
  * \f[
  *  \boldsymbol \sigma = \left(1 - \omega(\kappa (\boldsymbol \varepsilon)\right) \boldsymbol
  * \sigma_\text{elastic}(\boldsymbol \varepsilon)
  * \f]
  *
  * ... following a policy based design (hopefully applied correctly...), where ...
  *
  * @tparam TDamageLaw the damage law provides the \f$\omega(\kappa)\f$. This requires the methods `.Damage(double)` and
  * `.Derivative(double)` to be implemented, e.g.  NuTo::Constitutive::DamageLaw.
  *
  * @tparam TEvolution the evolution equation \f$ \kappa(\boldsymbol \varepsilon) \f$. This requires the methods
  * `.Kappa(strain)` and `.DkappaDstrain(strain)` to be implemented, e.g. NuTo::Laws::EvolutionImplicit.
  *
  * @tparam TElasticLaw the elastic constitutive law \f$ \sigma_\text{elastic}(\boldsymbol \varepsilon) \f$. This
  * requires the methods `.Stress(strain, ...)` and `.Tangent(strain, ...)` to be implemented, e.g.
  * NuTo::Laws::LinearElastic.
  *
  * @tparam TDim dimension
  */
template <int TDim, typename TDamageLaw, typename TEvolution, typename TElasticLaw>
class LocalIsotropicDamage : public MechanicsInterface<TDim>
{
public:
    LocalIsotropicDamage(TDamageLaw damageLaw, TEvolution evolution, TElasticLaw elasticLaw)
        : mDamageLaw(damageLaw)
        , mEvolution(evolution)
        , mElasticLaw(elasticLaw)
    {
    }

    virtual ~LocalIsotropicDamage() = default;

    EngineeringStress<TDim> Stress(EngineeringStrain<TDim> strain, double deltaT, CellIds ids) const override
    {
        auto kappa = mEvolution.Kappa(strain, deltaT, ids);
        return (1 - mDamageLaw.Damage(kappa)) * mElasticLaw.Stress(strain);
    }

    typename MechanicsInterface<TDim>::MechanicsTangent Tangent(EngineeringStrain<TDim> strain, double deltaT,
                                                                CellIds ids) const override
    {
        auto C = mElasticLaw.Tangent(strain, deltaT, ids);
        auto sigma = mElasticLaw.Stress(strain, deltaT, ids);

        auto kappa = mEvolution.Kappa(strain, deltaT, ids);
        auto omega = mDamageLaw.Damage(kappa);
        auto dOmegadKappa = mDamageLaw.Derivative(kappa);

        return (1 - omega) * C - sigma * dOmegadKappa * mEvolution.DkappaDstrain(strain, deltaT, ids);
    }

    void Update(EngineeringStrain<TDim> strain, double deltaT, CellIds ids)
    {
        mEvolution.Update(strain, deltaT, ids);
    }

public:
    // Intentionally made public. I assume that we know what we do. And this avoids lots of code, either getters
    // (possibly nonconst) or forwarding functions like `double Kappa(...), double Damage(...)`
    TDamageLaw mDamageLaw;
    TEvolution mEvolution;
    TElasticLaw mElasticLaw;
};

/** Explicit evolution equation for the NuTo::LocalIsotropicDamageLaw that implements
 * \f[
 *  \kappa(\boldsymbol \varepsilon) = \max \left(\kappa, \varepsilon_\text{eq}(\boldsymbol \varepsilon) \right)
 * \f]
 *
 * @tparam TDim dimension
 *
 * @tparam TStrainNorm strain norm \f$ \varepsilon_\text{eq}(\boldsymbol \varepsilon) \f$. This requires the methods
 * `Kappa(strain)` and `DkappaDstrain(strain)` to be implemented, e.g. NuTo::Constitutive::ModifiedMisesStrainNorm.
 */
template <int TDim, typename TStrainNorm>
class EvolutionImplicit
{

public:
    //! Constructor. As this evolution equation requires history data, they are also allocated.
    //! @param strainNorm strain norm, see class documentation
    //! @param numCells number of cells for the history data allocation
    //! @param numIpsPerCell nummer of integraiton points per cell for the history data allocation
    EvolutionImplicit(TStrainNorm strainNorm, size_t numCells = 0, size_t numIpsPerCell = 0)
        : mStrainNorm(strainNorm)
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
    TStrainNorm mStrainNorm;
    Eigen::MatrixXd mKappas;
};

} /* Laws */
} /* NuTo */
