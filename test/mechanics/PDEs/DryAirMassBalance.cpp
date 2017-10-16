#include <iostream>

#include "BoostUnitTest.h"
#include "mechanics/PDEs/DryAirMassBalance.h"
#include "mechanics/constitutive/laws/PorousMedium.h"

BOOST_AUTO_TEST_CASE(DensityChange)
{
    NuTo::PorousMedium material(0.5, 20.0, 2.0);

    Eigen::Matrix<double, 1, 2> N;
    N << 0.5, 0.5;

    const Eigen::MatrixXd result = Hygro::DryAirMassBalance::DensityChangeDueToGasPressure(material, N, 20.0);

    // manually computed result
    const Eigen::MatrixXd expected = N.transpose() * 1.86813202093239e-6 * N;

    BoostUnitTest::CheckEigenMatrix(result, expected);
}


BOOST_AUTO_TEST_CASE(SaturationVariation)
{
    NuTo::PorousMedium material(0.5, 20.0, 2.0);

    Eigen::Matrix<double, 1, 2> N;
    N << 0.5, 0.5;

    const Hygro::PoreState state(20.0, 5.0, 273.15);
    const Eigen::MatrixXd result = Hygro::DryAirMassBalance::VariationOfSaturation(material, N, 20.0, state);

    // manually computed result
    const Eigen::MatrixXd expected = N.transpose() * 5.63699763427398e-7 * N;

    BoostUnitTest::CheckEigenMatrix(result, expected);
}
