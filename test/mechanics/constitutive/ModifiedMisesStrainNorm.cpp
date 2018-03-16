#include "BoostUnitTest.h"
#include "nuto/mechanics/constitutive/ModifiedMisesStrainNorm.h"

using namespace NuTo;
using namespace NuTo::Constitutive;

template <int TDim>
std::vector<EngineeringStrain<TDim>> StrainTestCases()
{
    constexpr int VDim = Voigt::Dim(TDim);
    std::vector<EngineeringStrain<TDim>> strainCases(VDim * 2, EngineeringStrain<TDim>::Zero());

    // define test cases with exactly one value =! 0
    for (int i = 0; i < VDim; ++i)
    {
        strainCases[i][i] = M_PI;
        strainCases[i + VDim][i] = -M_PI;
    }

    // Plus some random cases
    for (int i = 0; i < 10; ++i)
        strainCases.push_back(Eigen::Matrix<double, VDim, 1>::Random());

    return strainCases;
}

//! param planeState Defaults to plane stress; for 1D/3D irrelevant
template <int TDim>
void CheckLocalEqStrainDerivativesMises(ePlaneState planeState = ePlaneState::PLANE_STRESS)
{
    double nu = 0.25;
    double k = 9.81;
    ModifiedMisesStrainNorm<TDim> modMises(nu, k, planeState);

    // check derivatives
    for (auto strain : StrainTestCases<TDim>())
    {
        double localEqStrain0 = modMises.Value(strain);

        constexpr int VDim = Voigt::Dim(TDim);
        Eigen::Matrix<double, VDim, 1> tangent = modMises.Derivative(strain);
        Eigen::Matrix<double, VDim, 1> tangent_CDF;

        // calculate derivative numerically
        for (int i = 0; i < VDim; ++i)
        {
            double delta = 1.e-8;
            strain[i] += delta;
            tangent_CDF[i] = (modMises.Value(strain) - localEqStrain0) / delta;
            strain[i] -= delta;
        }
        BoostUnitTest::CheckEigenMatrix(tangent, tangent_CDF, 1.e-6);
    }
}


BOOST_AUTO_TEST_CASE(CheckEqStrainDerivatives1d)
{
    CheckLocalEqStrainDerivativesMises<1>();
}

BOOST_AUTO_TEST_CASE(CheckEqStrainDerivatives2dPlaneStrain)
{
    CheckLocalEqStrainDerivativesMises<2>(ePlaneState::PLANE_STRAIN);
}

BOOST_AUTO_TEST_CASE(CheckEqStrainDerivatives2dPlaneStress)
{
    CheckLocalEqStrainDerivativesMises<2>(ePlaneState::PLANE_STRESS);
}

BOOST_AUTO_TEST_CASE(CheckEqStrainDerivatives3d)
{
    CheckLocalEqStrainDerivativesMises<3>();
}

BOOST_AUTO_TEST_CASE(EqStrainZero)
{
    const double k = 10;
    const double nu = 0.;
    BOOST_CHECK_SMALL(ModifiedMisesStrainNorm<1>(nu, k).Value(EngineeringStrain<1>::Zero()), 1.e-6);
    BOOST_CHECK_SMALL(ModifiedMisesStrainNorm<2>(nu, k).Value(EngineeringStrain<2>::Zero()), 1.e-6);
    BOOST_CHECK_SMALL(ModifiedMisesStrainNorm<3>(nu, k).Value(EngineeringStrain<3>::Zero()), 1.e-6);
}

BOOST_AUTO_TEST_CASE(EqStrain)
{
    const double k = 10;
    const double nu = 0.;
    EngineeringStrain<3> strain = EngineeringStrain<3>::Zero();
    strain[0] = 2.;
    BOOST_CHECK_CLOSE(ModifiedMisesStrainNorm<3>(nu, k).Value(strain), 2., 1.e-6);
    strain[0] = 0.;
    strain[2] = -2 * k;
    BOOST_CHECK_CLOSE(ModifiedMisesStrainNorm<3>(nu, k).Value(strain), 2., 1.e-6);
}
