#include "BoostUnitTest.h"
#include "mechanics/constitutive/laws/LinearElastic.h"

constexpr double E = 6174;
constexpr double nu = 0.36;

// https://en.wikipedia.org/wiki/Hooke%27s_law

BOOST_AUTO_TEST_CASE(LinearElastic1D)
{
    NuTo::Laws::LinearElastic<1> law(E, nu);
    auto C = law.Tangent(NuTo::EngineeringStrainPDE<1>::Zero());
    BoostUnitTest::CheckEigenMatrix(C, Eigen::Matrix<double, 1, 1>::Constant(E));

    // I trust Eigen::operator* and check LinearElastic::Stress(...) only once here
    
    BOOST_CHECK_CLOSE(law.Stress(NuTo::EngineeringStrainPDE<1>::Constant(12))[0], E * 12, 1.e-10);
}

BOOST_AUTO_TEST_CASE(LinearElastic2DPlaneStress)
{
    NuTo::Laws::LinearElastic<2> law(E, nu, NuTo::ePlaneState::PLANE_STRESS);
    auto C = law.Tangent(NuTo::EngineeringStrainPDE<2>::Zero());

    Eigen::Matrix3d expected = Eigen::Matrix3d::Zero();
    expected(0, 0) = 1;
    expected(1, 1) = 1;
    expected(0, 1) = nu;
    expected(1, 0) = nu;
    expected(2, 2) = (1 - nu) / 2;
    expected *= E / (1 - nu * nu);
    BoostUnitTest::CheckEigenMatrix(C, expected);
}

BOOST_AUTO_TEST_CASE(LinearElastic2DPlaneStrain)
{
    NuTo::Laws::LinearElastic<2> law(E, nu, NuTo::ePlaneState::PLANE_STRAIN);
    auto C = law.Tangent(NuTo::EngineeringStrainPDE<2>::Zero());

    Eigen::Matrix3d expected = Eigen::Matrix3d::Zero();
    expected(0, 0) = 1 - nu;
    expected(1, 1) = 1 - nu;
    expected(0, 1) = nu;
    expected(1, 0) = nu;
    expected(2, 2) = (1 - 2 * nu) / 2;
    expected *= E / ((1 + nu) * (1 - 2 * nu));
    BoostUnitTest::CheckEigenMatrix(C, expected);
}

BOOST_AUTO_TEST_CASE(LinearElastic3D)
{
    NuTo::Laws::LinearElastic<3> law(E, nu);
    auto C = law.Tangent(NuTo::EngineeringStrainPDE<3>::Zero());

    Eigen::Matrix<double, 6, 6> expected = Eigen::Matrix<double, 6, 6>::Zero();
    expected(0, 0) = 1 - nu;
    expected(1, 1) = 1 - nu;
    expected(2, 2) = 1 - nu;

    expected(0, 1) = nu;
    expected(1, 0) = nu;
    expected(0, 2) = nu;
    expected(2, 0) = nu;
    expected(2, 1) = nu;
    expected(1, 2) = nu;

    expected(3, 3) = (1 - 2 * nu) / 2;
    expected(4, 4) = (1 - 2 * nu) / 2;
    expected(5, 5) = (1 - 2 * nu) / 2;
    expected *= E / ((1 + nu) * (1 - 2 * nu));
    BoostUnitTest::CheckEigenMatrix(C, expected);
}
