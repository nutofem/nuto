#include <iostream>

#include "BoostUnitTest.h"
#include "mechanics/PDEs/DryAirMassBalance.h"
#include "mechanics/PDEs/ConcreteMedium.h"

BOOST_AUTO_TEST_CASE(DensityChange)
{
    ConcreteMedium concrete(0.5, 20.0, 2.0, 0.1);

    Eigen::Matrix<double, 1, 2> N;
    N << 0.5, 0.5;

    const Hygro::PoreState state(20.0, 5.0, 273.15);
    const Eigen::MatrixXd result = Hygro::DryAirMassBalance::DensityChangeDueToGasPressure(state, concrete, N);

    // manually computed result
    const Eigen::MatrixXd expected = N.transpose() * 1.86813202093239e-6 * N;

    BoostUnitTest::CheckEigenMatrix(result, expected);
}


BOOST_AUTO_TEST_CASE(SaturationVariation)
{
    ConcreteMedium concrete(0.5, 20.0, 2.0, 0.1);

    Eigen::Matrix<double, 1, 2> N;
    N << 0.5, 0.5;

    const Hygro::PoreState state(20.0, 5.0, 273.15);
    const Eigen::MatrixXd result = Hygro::DryAirMassBalance::VariationOfSaturation(state, concrete, N);

    // manually computed result
    const Eigen::MatrixXd expected = N.transpose() * 5.63699763427398e-7 * N;

    BoostUnitTest::CheckEigenMatrix(result, expected);
}
