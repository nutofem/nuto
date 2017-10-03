#include <iostream>

#include "BoostUnitTest.h"
#include "mechanics/PDEs/DryAirMassBalance.h"
#include "mechanics/constitutive/laws/PorousMedium.h"

BOOST_AUTO_TEST_CASE(DensityChange)
{
    Hygro::DryAirMassBalance massBalance;
    NuTo::PorousMedium material(0.5, 20.0, 2.0);

    Eigen::Matrix<double, 1, 2> N;
    N << 0.5, 0.5;

    const Eigen::MatrixXd result = massBalance.DensityChangeDueToGasPressure(material, N, 20.0);

    // manually computed result
    const Eigen::MatrixXd expected = N.transpose() * 1.7114883834233898e-06 * N;

    BoostUnitTest::CheckEigenMatrix(result, expected);
}

BOOST_AUTO_TEST_CASE(SaturationVariation)
{
    Hygro::DryAirMassBalance massBalance;
    NuTo::PorousMedium material(0.5, 20.0, 2.0);

    Eigen::Matrix<double, 1, 2> N;
    N << 0.5, 0.5;

    const Eigen::MatrixXd result = massBalance.VariationOfSaturation(material, N, 20.0, 5.0);

    // manually computed result
    const Eigen::MatrixXd expected = N.transpose() * 5.16204155559069e-7 * N;

    BoostUnitTest::CheckEigenMatrix(result, expected);
}
