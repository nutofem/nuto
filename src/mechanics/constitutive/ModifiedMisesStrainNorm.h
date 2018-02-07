#pragma once
#include "base/Exception.h"
#include "mechanics/constitutive/ConstitutivePlaneStateEnum.h"
#include "mechanics/constitutive/EngineeringStrain.h"
#include "mechanics/constitutive/EngineeringStrainInvariants.h"

namespace NuTo
{
namespace Constitutive
{
/**
 * Equivalent strain (strain norm) based on the modified mises norm.
 * De Vree et al. "Comparision of nonlocal approaches in continuum damage mechanics", Computers \& Structures,
 * 1995.
 * https://doi.org/10.1016/0045-7949(94)00501-S
 *
 * It includes the value `k` in a way that expresses the ratio of the materials compressive to tensile strength:
 * A
 * uniaxial tensile strain and a k-times higher uniaxial compressive strain both lead to the same equivalent
 * strain.
 * With the strain tensor invariant I1 and the deviatoric strain tensor invariant J2, it reads
 * \f[/
 * \varepsilon_\text{eq}(\boldsymbol \varepsilon) = \frac{k-1}{2k(1-2\nu)}I_1 +
 * \frac{1}{2k}\sqrt{\left(\frac{k-1}{1-2\nu}I_1\right)^2 + \frac{2k}{(1+\nu)^2}J_2}
 * \f]
 *
 * @tparam TDim dimension 1,2,3
 */
template <int TDim>
class ModifiedMisesStrainNorm
{
public:
    //! Constructor.
    //! @param nu Poissons ratio
    //! @param k ratio of compressive strength to tensile strength, ~10 for concrete
    //! @param planeState PLANE_STRAIN or PLANE_STRESS
    ModifiedMisesStrainNorm(double nu, double k, ePlaneState planeState = ePlaneState::PLANE_STRAIN);

    //! @param strain strain to evaluate
    //! @return the value of the modified Mises strain norm
    double Value(const NuTo::EngineeringStrain<TDim>& strain) const;

    //! @param strain strain to evaluate
    //! @return the derivative of the modified Mises strain norm with respect to the strain
    Eigen::Matrix<double, Voigt::Dim(TDim), 1> Derivative(const NuTo::EngineeringStrain<TDim>& strain) const;

private:
    //! Transforms `strain` of dimension TDim to dimension 3, depending on `nu` and `ePlaneState`
    //! @param nu poissons ratio
    //! @param planeState PLANE_STRAIN or PLANE_STRESS
    //! @remark This method may be useful for other classes. However, it cannot be a member of
    //!         EngineeringStrain since it is only valid for isotropic strains
    EngineeringStrain<3> Strain3D(const EngineeringStrain<TDim>& strain, double nu, ePlaneState) const;

    double mNu;
    double mK;
    ePlaneState mPlaneState;

    double mK1;
    double mK2;
};


template <int TDim>
inline ModifiedMisesStrainNorm<TDim>::ModifiedMisesStrainNorm(double nu, double k, ePlaneState planeState)
    : mNu(nu)
    , mK(k)
    , mPlaneState(planeState)
    , mK1((mK - 1.) / (2. * mK * (1 - 2 * mNu)))
    , mK2(3. / (mK * (1. + mNu) * (1. + mNu)))
{
}


template <>
inline EngineeringStrain<3> ModifiedMisesStrainNorm<1>::Strain3D(const EngineeringStrain<1>& strain, double nu,
                                                                 ePlaneState) const
{
    EngineeringStrain<3> v;
    v[0] = strain[0];
    v[1] = -nu * strain[0];
    v[2] = -nu * strain[0];
    v[3] = 0.;
    v[4] = 0.;
    v[5] = 0;
    return v;
}

template <>
inline EngineeringStrain<3> ModifiedMisesStrainNorm<2>::Strain3D(const EngineeringStrain<2>& strain, double nu,
                                                                 ePlaneState planeState) const
{
    EngineeringStrain<3> v;
    if (planeState == ePlaneState::PLANE_STRAIN)
    {
        v[0] = strain[0];
        v[1] = strain[1];
        v[2] = 0;
        v[3] = 0.;
        v[4] = 0.;
        v[5] = strain[2];
        return v;
    }
    if (planeState == ePlaneState::PLANE_STRESS)
    {
        v[0] = strain[0];
        v[1] = strain[1];
        v[2] = nu / (nu - 1.) * (strain[0] + strain[1]);
        v[3] = 0.;
        v[4] = 0.;
        v[5] = strain[2];
        return v;
    }
    throw Exception(__PRETTY_FUNCTION__, "Specify section behavior, either PLANE_STRAIN or PLANE_STRESS");
}

template <>
inline EngineeringStrain<3> ModifiedMisesStrainNorm<3>::Strain3D(const EngineeringStrain<3>& strain, double,
                                                                 ePlaneState) const
{
    return strain;
}


template <int TDim>
inline double ModifiedMisesStrainNorm<TDim>::Value(const NuTo::EngineeringStrain<TDim>& strain) const
{
    EngineeringStrain<3> strain3D = Strain3D(strain, mNu, mPlaneState);
    double I1 = NuTo::EngineeringStrainInvariants::I1(strain3D);
    double J2 = NuTo::EngineeringStrainInvariants::J2(strain3D);

    double A = std::sqrt(mK1 * mK1 * I1 * I1 + mK2 * J2);

    return mK1 * I1 + A;
}


template <>
inline Eigen::Matrix<double, 1, 1>
ModifiedMisesStrainNorm<1>::Derivative(const NuTo::EngineeringStrain<1>& strain) const
{
    Eigen::Matrix<double, 1, 1> derivative;
    double dJ2dexx = 2. / 3. * strain[0] * (1 + mNu) * (1 + mNu);
    double dI1dexx = (1 - 2 * mNu);

    EngineeringStrain<3> strain3D = Strain3D(strain, mNu, mPlaneState);
    double I1 = NuTo::EngineeringStrainInvariants::I1(strain3D);
    double J2 = NuTo::EngineeringStrainInvariants::J2(strain3D);

    double A = std::sqrt(mK1 * mK1 * I1 * I1 + mK2 * J2);

    if (A == 0)
        derivative[0] = mK1 * dI1dexx;
    else
        derivative[0] = mK1 * dI1dexx + 1. / (2 * A) * (2 * mK1 * mK1 * I1 * dI1dexx + mK2 * dJ2dexx);

    return derivative;
}

template <>
inline Eigen::Matrix<double, 3, 1>
ModifiedMisesStrainNorm<2>::Derivative(const NuTo::EngineeringStrain<2>& strain) const
{
    Eigen::Matrix<double, 3, 1> derivative;

    EngineeringStrain<3> strain3D = Strain3D(strain, mNu, mPlaneState);
    double I1 = NuTo::EngineeringStrainInvariants::I1(strain3D);
    double J2 = NuTo::EngineeringStrainInvariants::J2(strain3D);

    double A = std::sqrt(mK1 * mK1 * I1 * I1 + mK2 * J2);

    switch (mPlaneState)
    {
    case ePlaneState::PLANE_STRAIN:
    {
        if (A == 0)
        {
            derivative[0] = mK1;
            derivative[1] = mK1;
            derivative[2] = 0;
        }
        else
        {
            double dJ2dexx = 1. / 3. * (2 * strain[0] - strain[1]);
            double dJ2deyy = 1. / 3. * (2 * strain[1] - strain[0]);
            double dJ2dgxy = 0.5 * strain[2];

            derivative[0] = mK1 + 1. / (2 * A) * (2 * mK1 * mK1 * I1 + mK2 * dJ2dexx);
            derivative[1] = mK1 + 1. / (2 * A) * (2 * mK1 * mK1 * I1 + mK2 * dJ2deyy);
            derivative[2] = 1. / (2 * A) * (mK2 * dJ2dgxy);
        }
        break;
    }
    case ePlaneState::PLANE_STRESS:
    {
        double dI1dexxeyy = (1 + mNu / (mNu - 1));
        if (A == 0)
        {
            derivative[0] = dI1dexxeyy * mK1;
            derivative[1] = dI1dexxeyy * mK1;
            derivative[2] = 0;
        }
        else
        {
            double strainxy = mNu / (mNu - 1.) * (strain[0] + strain[1]);
            double dJ2dexx = 1. / 3. * (2. * strain[0] - strain[1] - 2. * strainxy + 2 * mNu / (mNu - 1) * strainxy);
            double dJ2deyy = 1. / 3. * (2. * strain[1] - strain[0] - 2. * strainxy + 2 * mNu / (mNu - 1) * strainxy);
            double dJ2dgxy = .5 * strain[2];

            derivative[0] = dI1dexxeyy * mK1 + 1. / (2. * A) * (2. * dI1dexxeyy * mK1 * mK1 * I1 + mK2 * dJ2dexx);
            derivative[1] = dI1dexxeyy * mK1 + 1. / (2. * A) * (2. * dI1dexxeyy * mK1 * mK1 * I1 + mK2 * dJ2deyy);
            derivative[2] = 1. / (2. * A) * (mK2 * dJ2dgxy);
        }
        break;
    }
    }
    return derivative;
}

template <>
inline Eigen::Matrix<double, 6, 1>
ModifiedMisesStrainNorm<3>::Derivative(const NuTo::EngineeringStrain<3>& strain) const
{
    Eigen::Matrix<double, 6, 1> derivative;

    double I1 = NuTo::EngineeringStrainInvariants::I1(strain);
    double J2 = NuTo::EngineeringStrainInvariants::J2(strain);

    double A = std::sqrt(mK1 * mK1 * I1 * I1 + mK2 * J2);

    if (A == 0)
    {
        derivative[0] = mK1;
        derivative[1] = mK1;
        derivative[2] = mK1;
        derivative[3] = 0;
        derivative[4] = 0;
        derivative[5] = 0;
    }
    else
    {
        double dJ2dexx = 1. / 3. * (2 * strain[0] - strain[1] - strain[2]);
        double dJ2deyy = 1. / 3. * (2 * strain[1] - strain[0] - strain[2]);
        double dJ2dezz = 1. / 3. * (2 * strain[2] - strain[0] - strain[1]);
        double dJ2dgyz = 0.5 * strain[3];
        double dJ2dgzx = 0.5 * strain[4];
        double dJ2dgxy = 0.5 * strain[5];

        derivative[0] = mK1 + 1. / (2 * A) * (2 * mK1 * mK1 * I1 + mK2 * dJ2dexx);
        derivative[1] = mK1 + 1. / (2 * A) * (2 * mK1 * mK1 * I1 + mK2 * dJ2deyy);
        derivative[2] = mK1 + 1. / (2 * A) * (2 * mK1 * mK1 * I1 + mK2 * dJ2dezz);
        derivative[3] = 1. / (2 * A) * mK2 * dJ2dgyz;
        derivative[4] = 1. / (2 * A) * mK2 * dJ2dgzx;
        derivative[5] = 1. / (2 * A) * mK2 * dJ2dgxy;
    }
    return derivative;
}

} /* Constitutive */
} /* NuTo */
