/*
 * EquivalentStrain.cpp
 *
 *  Created on: 10 Mar 2016
 *      Author: ttitsche
 */

#include "mechanics/constitutive/inputoutput/EquivalentStrain.h"


namespace NuTo
{
template <int TDim>
EquivalentStrainModifiedMises<TDim>::EquivalentStrainModifiedMises(const EngineeringStrain<TDim>& rStrain, double rK,
                                                                   double rNu, ePlaneState planeState)
    : mK1((rK - 1.) / (2. * rK * (1 - 2 * rNu)))
    , mK2(3 / (rK * (1 + rNu) * (1 + rNu)))
    , mStrain3D(rStrain.As3D(rNu, planeState))
    , mI1(mStrain3D.InvariantI1())
    , mJ2(mStrain3D.InvariantJ2())
    , mA(std::sqrt(mK1 * mK1 * mI1 * mI1 + mK2 * mJ2))
    , mNu(rNu)
    , mPlaneState(planeState)
{
}

template <int TDim>
double EquivalentStrainModifiedMises<TDim>::Get() const
{
    return mK1 * mI1 + mA;
}


template <>
ConstitutiveVector<1> EquivalentStrainModifiedMises<1>::GetDerivative() const
{
    ConstitutiveVector<1> tangent;
    double dJ2dexx = 2. / 3. * mStrain3D[0] * (1 + mNu) * (1 + mNu);
    double dI1dexx = (1 - 2 * mNu);

    if (mA == 0)
        tangent[0] = mK1 * dI1dexx;
    else
        tangent[0] = mK1 * dI1dexx + 1. / (2 * mA) * (2 * mK1 * mK1 * mI1 * dI1dexx + mK2 * dJ2dexx);

    return tangent;
}

template <>
ConstitutiveVector<3> EquivalentStrainModifiedMises<2>::GetDerivative() const
{

    ConstitutiveVector<3> tangent;
    if (mPlaneState == ePlaneState::PLANE_STRAIN)
    {
        if (mA == 0)
        {
            tangent[0] = mK1;
            tangent[1] = mK1;
            tangent[2] = 0;
        }
        else
        {
            double dJ2dexx = 1. / 3. * (2 * mStrain3D[0] - mStrain3D[1]);
            double dJ2deyy = 1. / 3. * (2 * mStrain3D[1] - mStrain3D[0]);
            double dJ2dgxy = 0.5 * mStrain3D[5];

            tangent[0] = mK1 + 1. / (2 * mA) * (2 * mK1 * mK1 * mI1 + mK2 * dJ2dexx);
            tangent[1] = mK1 + 1. / (2 * mA) * (2 * mK1 * mK1 * mI1 + mK2 * dJ2deyy);
            tangent[2] = 1. / (2 * mA) * (mK2 * dJ2dgxy);
        }
        return tangent;
    }

    if (mPlaneState == ePlaneState::PLANE_STRESS)
    {
        double dI1dexxeyy = (1 + mNu / (mNu - 1));
        if (mA == 0)
        {
            tangent[0] = dI1dexxeyy * mK1;
            tangent[1] = dI1dexxeyy * mK1;
            tangent[2] = 0;
        }
        else
        {
            double dJ2dexx = 1. / 3. * (2. * mStrain3D[0] - mStrain3D[1] - 2. * mStrain3D[2] +
                                        2 * mNu / (mNu - 1) * mStrain3D[2]);
            double dJ2deyy = 1. / 3. * (2. * mStrain3D[1] - mStrain3D[0] - 2. * mStrain3D[2] +
                                        2 * mNu / (mNu - 1) * mStrain3D[2]);
            double dJ2dgxy = .5 * mStrain3D[5];

            tangent[0] = dI1dexxeyy * mK1 + 1. / (2. * mA) * (2. * dI1dexxeyy * mK1 * mK1 * mI1 + mK2 * dJ2dexx);
            tangent[1] = dI1dexxeyy * mK1 + 1. / (2. * mA) * (2. * dI1dexxeyy * mK1 * mK1 * mI1 + mK2 * dJ2deyy);
            tangent[2] = 1. / (2. * mA) * (mK2 * dJ2dgxy);
        }
        return tangent;
    }

    throw Exception(__PRETTY_FUNCTION__,
                             "Section type undefined. Choose either PLANE_STRAIN or PLANE_STRESS.");
}

template <>
ConstitutiveVector<6> EquivalentStrainModifiedMises<3>::GetDerivative() const
{
    ConstitutiveVector<6> tangent;

    if (mA == 0)
    {
        tangent.SetZero();
        tangent[0] = mK1;
        tangent[1] = mK1;
        tangent[2] = mK1;
    }
    else
    {
        double dJ2dexx = 1. / 3. * (2 * mStrain3D[0] - mStrain3D[1] - mStrain3D[2]);
        double dJ2deyy = 1. / 3. * (2 * mStrain3D[1] - mStrain3D[0] - mStrain3D[2]);
        double dJ2dezz = 1. / 3. * (2 * mStrain3D[2] - mStrain3D[0] - mStrain3D[1]);
        double dJ2dgyz = 0.5 * mStrain3D[3];
        double dJ2dgzx = 0.5 * mStrain3D[4];
        double dJ2dgxy = 0.5 * mStrain3D[5];

        tangent[0] = mK1 + 1. / (2 * mA) * (2 * mK1 * mK1 * mI1 + mK2 * dJ2dexx);
        tangent[1] = mK1 + 1. / (2 * mA) * (2 * mK1 * mK1 * mI1 + mK2 * dJ2deyy);
        tangent[2] = mK1 + 1. / (2 * mA) * (2 * mK1 * mK1 * mI1 + mK2 * dJ2dezz);
        tangent[3] = 1. / (2 * mA) * mK2 * dJ2dgyz;
        tangent[4] = 1. / (2 * mA) * mK2 * dJ2dgzx;
        tangent[5] = 1. / (2 * mA) * mK2 * dJ2dgxy;
    }
    return tangent;
}

} // namespace NuTo


// NuTo::LocalEqStrain NuTo::GradientDamageEngineeringStress::CalculateLocalEqStrain2D(const NuTo::EngineeringStrain2D&
// rStrain2D) const
//{
//    NuTo::LocalEqStrain localEqStrain;
//    // calculate principal strains e1 and e2
//
//    double A = (rStrain2D(0) + rStrain2D(1)) / 2.;
//    double B = std::sqrt(std::pow(rStrain2D(0) - rStrain2D(1), 2.) / 4. + std::pow(rStrain2D(2), 2.) / 4.);
//
//    // macaulay brackets of principal strains
//    double e1 = std::max(A + B, 0.);
//    double e2 = std::max(A - B, 0.);
//
//    localEqStrain(0) = std::sqrt(e1 * e1 + e2 * e2);
//    return localEqStrain;
//}
//
// NuTo::LocalEqStrain NuTo::GradientDamageEngineeringStress::CalculateLocalEqStrain3D(const NuTo::EngineeringStrain3D&
// rStrain3D) const
//{
//    throw NuTo::Exception("[NuTo::GradientDamageEngineeringStress::CalculateLocalEqStrain3D] not
//    implemented");
//}
//
// NuTo::ConstitutiveTangentLocal<3, 1> NuTo::GradientDamageEngineeringStress::CalculateLocalEqStrainTangent2D(const
// NuTo::EngineeringStrain2D& rStrain2D) const
//{
//    NuTo::ConstitutiveTangentLocal<3, 1> tangent;
//    double A = (rStrain2D(0) + rStrain2D(1)) / 2.;
//    double B = 0.5 * std::sqrt(std::pow(rStrain2D(0) - rStrain2D(1), 2.) + std::pow(rStrain2D(2), 2.));
//
//    // macaulay brackets of principal strains
//    double e1 = std::max(A + B, 0.);
//    double e2 = std::max(A - B, 0.);
//
//    if (e1 + e2 == 0)
//    {
//        tangent.setZero();
//        return tangent;
//    }
//    if (e1 != 0 and e2 != 0)
//    {
//        double eMazar = std::sqrt(e1 * e1 + e2 * e2);
//        tangent(0) = rStrain2D(0) / eMazar;
//        tangent(1) = rStrain2D(1) / eMazar;
//        tangent(2) = rStrain2D(2) / eMazar * 0.5;
//        return tangent;
//    }
//    if (e1 != 0)
//    {
//        tangent(0) = (e1 - rStrain2D[1]) / (2 * B);
//        tangent(1) = (e1 - rStrain2D[0]) / (2 * B);
//        tangent(2) = 0.5 * rStrain2D[2] / (2 * B);
//        return tangent;
//    }
//    if (e2 != 0)
//    {
//        tangent(0) = (rStrain2D[0] - e2) / (2 * B);
//        tangent(1) = (rStrain2D[1] - e2) / (2 * B);
//        tangent(2) = -0.5 * rStrain2D[2] / (2 * B);
//        return tangent;
//    }
//
////    std::cout << e1 << std::endl;
////    std::cout << e2 << std::endl;
//
//    throw NuTo::Exception("[NuTo::GradientDamageEngineeringStress::CalculateLocalEqStrainTangent2D] error");
//}
//
// NuTo::ConstitutiveTangentLocal<6, 1> NuTo::GradientDamageEngineeringStress::CalculateLocalEqStrainTangent3D(const
// NuTo::EngineeringStrain3D& rStrain3D) const
//{
//    throw NuTo::Exception("[NuTo::GradientDamageEngineeringStress::CalculateLocalEqStrainTangent3D] not
//    implemented");
//}

template class NuTo::EquivalentStrainModifiedMises<1>;
template class NuTo::EquivalentStrainModifiedMises<2>;
template class NuTo::EquivalentStrainModifiedMises<3>;
