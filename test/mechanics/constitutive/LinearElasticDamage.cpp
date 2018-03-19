#include "BoostUnitTest.h"
#include "nuto/math/EigenCompanion.h"
#include "nuto/mechanics/constitutive/LinearElasticDamage.h"
#include "nuto/mechanics/constitutive/LinearElastic.h"

constexpr double E = 6174;
constexpr double nu = 0.36;

// https://en.wikipedia.org/wiki/Hooke%27s_law

template <int TDim>
void CheckIsotropic()
{
    NuTo::Laws::LinearElasticDamage<3> linearElasticDamage(E, nu);
    NuTo::Laws::LinearElastic<3> linearElastic(E, nu);
    NuTo::EngineeringStrain<3> strain = NuTo::EngineeringStrain<3>::Random();
    double omega = 0.4;
    BoostUnitTest::CheckEigenMatrix((1. - omega) * linearElastic.Stress(strain),
                                    linearElasticDamage.Stress(strain, omega));

    BoostUnitTest::CheckEigenMatrix((1. - omega) * linearElastic.Tangent(strain),
                                    linearElasticDamage.DstressDstrain(strain, omega));

    BoostUnitTest::CheckEigenMatrix(-linearElastic.Stress(strain), linearElasticDamage.DstressDomega(strain, omega));
}

BOOST_AUTO_TEST_CASE(VsLinearElastic)
{
    CheckIsotropic<1>();
    CheckIsotropic<2>();
    CheckIsotropic<3>();
}

void CheckTangents(std::initializer_list<double> values, double omega)
{
    //
    // check DstressDstrain
    //
    NuTo::EngineeringStrain<2> strain = NuTo::EigenCompanion::ToEigen(values);
    BOOST_TEST_MESSAGE("Checking tangent for e = " << strain.transpose() << " and w = " << omega);
    NuTo::Laws::LinearElasticDamage<2> law(20000, 0.2);

    NuTo::EngineeringTangent<2> tangent = law.DstressDstrain(strain, omega);
    NuTo::EngineeringTangent<2> tangent_diff = tangent * 0.;

    const double delta = 1.e-8;

    for (int i = 0; i < 3; ++i)
    {
        NuTo::EngineeringStress<2> s0 = law.Stress(strain, omega);
        strain[i] += delta;
        NuTo::EngineeringStress<2> s1 = law.Stress(strain, omega);
        strain[i] -= delta;
        tangent_diff.col(i) = (s1 - s0) / delta;
    }

    BoostUnitTest::CheckEigenMatrix(tangent, tangent_diff, tangent.maxCoeff() * 1.e-3);


    //
    // check DstressDomega
    //
    NuTo::EngineeringStress<2> dStressdOmega = law.DstressDomega(strain, omega);
    NuTo::EngineeringStress<2> dStressdOmega_diff =
            (law.Stress(strain, omega + delta) - law.Stress(strain, omega)) / delta;
    double tolerance = std::max(dStressdOmega.cwiseAbs().maxCoeff() * 1.e-6, 1.e-6);
    BoostUnitTest::CheckEigenMatrix(dStressdOmega, dStressdOmega_diff, tolerance);
}

BOOST_AUTO_TEST_CASE(Tangent2D)
{
    for (double omega : {0., 0.5, 0.8})
    {
        CheckTangents({0, 0, 0}, omega);

        CheckTangents({2, 0, 0}, omega);
        CheckTangents({-2, 0, 0}, omega);
        CheckTangents({0, 2, 0}, omega);
        CheckTangents({0, -2, 0}, omega);
        CheckTangents({0, 0, 2}, omega);
        CheckTangents({0, 0, -2}, omega);

        CheckTangents({1, 1, 0}, omega);
        CheckTangents({0, 1, 1}, omega);
        CheckTangents({1, 0, 1}, omega);

        CheckTangents({-1, 1, 0}, omega);
        CheckTangents({0, -1, 1}, omega);
        CheckTangents({1, 0, -1}, omega);

        CheckTangents({1, -1, 0}, omega);
        CheckTangents({0, 1, -1}, omega);
        CheckTangents({-1, 0, 1}, omega);
    }
}
