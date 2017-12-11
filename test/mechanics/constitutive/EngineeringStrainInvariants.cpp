#include "BoostUnitTest.h"
#include <eigen3/Eigen/Dense> // for Eigen::determinant()
#include "mechanics/constitutive/EngineeringStrainInvariants.h"

using namespace NuTo::EngineeringStrainInvariants;
using Strain = NuTo::EngineeringStrain<3>;

BOOST_AUTO_TEST_CASE(TensorToTensor)
{
    Strain v = Strain::Random();
    Eigen::Matrix3d m;
    m << v[0], .5 * v[5], .5 * v[4], .5 * v[5], v[1], .5 * v[3], .5 * v[4], .5 * v[3], v[2];
    BoostUnitTest::CheckEigenMatrix(ToTensor(v), m);
}

BOOST_AUTO_TEST_CASE(TensorInvariantI1)
{
    Strain v = Strain::Random();
    BOOST_CHECK_CLOSE(I1(v), ToTensor(v).trace(), 1.e-10);
}

BOOST_AUTO_TEST_CASE(TensorInvariantI2)
{
    Strain v = Strain::Random();
    auto m = ToTensor(v);
    BOOST_CHECK_CLOSE(I2(v), .5 * (m.trace() * m.trace() - (m * m).trace()), 1.e-10);
}

BOOST_AUTO_TEST_CASE(TensorInvariantI3)
{
    Strain v = Strain::Random();
    BOOST_CHECK_CLOSE(I3(v), ToTensor(v).determinant(), 1.e-10);
}

BOOST_AUTO_TEST_CASE(TensorInvariantJ2)
{
    Strain v = Strain::Random();
    BOOST_CHECK_CLOSE(J2(v), 1. / 3. * I1(v) * I1(v) - I2(v), 1.e-10);
}

BOOST_AUTO_TEST_CASE(TensorDeviatoric)
{
    Strain v = Strain::Random();
    BOOST_CHECK_CLOSE(J2(v), -I2(Deviatoric(v)), 1.e-10);
}
