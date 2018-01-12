#pragma once

#include <vector>
#include <mechanics/constitutive/MechanicsInterface.h>

namespace NuTo
{
namespace Laws
{

//! Local damage law with an isotropic damage variable
//! \f[
//!  \bm \sigma = (1 - \omega(\kappa\(\bm \varepsilon)) \sigma_\text{elastic}(\bm \varepsilon)
//! \f]
//! ... following a policy based design (hopefully applied correctly...), where ...
//! @tparam TDamageLaw the damage law provides the \f$\omega(\kappa)\f$. This requires the methods `.Damage(double)` and
//! `.Derivative(double)` to be implemented.
//! @tparam TEvolution the evolution equation \f$ \kappa(\bm \varepsilon) \f$. This requires the methods
//! `.Kappa(strain)` and `.dKappadStrain(strain)` to be implemented.
//! @tparam TElasticLaw the elastic constitutive law \f$ \sigma_\text{elastic}(\bm \varepsilon) \f$. This requires the
//! methods `.Stress(strain, ...)` and `.Tangent(strain, ...)` to be implemented.
//! @tparam TDim dimension
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

    EngineeringStress<TDim> Stress(EngineeringStrain<TDim> strain, double deltaT, int cellId, int ipId) const override
    {
        auto kappa = mEvolution.Kappa(strain, deltaT, cellId, ipId);
        return (1 - mDamageLaw.Damage(kappa)) * mElasticLaw.Stress(strain);
    }

    typename MechanicsInterface<TDim>::MechanicsTangent Tangent(EngineeringStrain<TDim> strain, double deltaT,
                                                                int cellId, int ipId) const override
    {
        auto C = mElasticLaw.Tangent(strain, deltaT, cellId, ipId);
        auto sigma = mElasticLaw.Stress(strain, deltaT, cellId, ipId);

        auto kappa = mEvolution.Kappa(strain, deltaT, cellId, ipId);
        auto omega = mDamageLaw.Damage(kappa);
        auto dOmegadKappa = mDamageLaw.Derivative(kappa);

        return (1 - omega) * C - sigma * dOmegadKappa * mEvolution.dKappadStrain(strain, deltaT, cellId, ipId);
    }

    void Update(EngineeringStrain<TDim> strain, double deltaT, int cellId, int ipId)
    {
        mEvolution.Update(strain, deltaT, cellId, ipId);
    }

private:
    TDamageLaw mDamageLaw;
    TEvolution mEvolution;
    TElasticLaw mElasticLaw;
};

template <int TDim, typename TDamageLaw, typename TStrainNorm, typename TElasticLaw>
static auto CreateLocalIsotropicDamage(TDamageLaw damageLaw, TStrainNorm strainNorm, TElasticLaw elasticLaw)
{
    return LocalIsotropicDamage<TDim, TDamageLaw, TStrainNorm, TElasticLaw>(damageLaw, strainNorm, elasticLaw);
}

template <int TDim, typename TStrainNorm>
class EvolutionImplicit
{
public:
    EvolutionImplicit(TStrainNorm strainNorm, size_t numCells, size_t numIpsPerCell)
        : mStrainNorm(strainNorm)
        , mNumIpsPerCell(numIpsPerCell)
        , mKappas(std::vector<double>(numCells * numIpsPerCell, 0.))
    {
    }

    double Kappa(EngineeringStrain<TDim> strain, double, int cellId, int ipId) const
    {
        return std::max(mStrainNorm.Value(strain), mKappas[Ip(cellId, ipId)]);
    }

    double dKappadStrain(EngineeringStrain<TDim> strain, double, int cellId, int ipId) const
    {
        if (mStrainNorm.Value(strain) > mKappas[Ip(cellId, ipId)])
            return 1.;
        return 0;
    }

    void Update(EngineeringStrain<TDim> strain, double deltaT, int cellId, int ipId)
    {
        mKappas[Ip(cellId, ipId)] = Kappa(strain, deltaT, cellId, ipId);
    }

protected:
    size_t Ip(int cellId, int ipId) const
    {
        return cellId * mNumIpsPerCell + ipId;
    }

private:
    TStrainNorm mStrainNorm;
    size_t mNumIpsPerCell;
    std::vector<double> mKappas;
};

} /* Laws */
} /* NuTo */
