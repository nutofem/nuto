#include <iostream>
#include "BoostUnitTest.h"
#include "mechanics/PDEs/HygroHelpers.h"

const double eps = 1e-9;

BOOST_AUTO_TEST_CASE(water_density)
{
    // Tested againts it's own result, seemed reasonable. Guards against regressions.
    BOOST_CHECK_CLOSE(Hygro::DensityOfWater(273.15), 1015.11, eps);
    BOOST_CHECK_CLOSE(Hygro::DensityOfWater(373.15), 962.707, eps);
    BOOST_CHECK_CLOSE(Hygro::DensityOfWater(473.15), 867.442, eps);
}


BOOST_AUTO_TEST_CASE(saturation_pressure)
{
    // Tested againts it's own result, seemed reasonable. Guards against regressions.
    BOOST_CHECK_CLOSE(Hygro::SaturationPressure(273.15), 0.0006112128459398556, eps);
    BOOST_CHECK_CLOSE(Hygro::SaturationPressure(373.15), 0.10141799381792783, eps);
    BOOST_CHECK_CLOSE(Hygro::SaturationPressure(473.15), 1.5549392220497635, eps);
    BOOST_CHECK_CLOSE(Hygro::SaturationPressure(573.15), 8.587867486373652, eps);
}


BOOST_AUTO_TEST_CASE(kelvin_eq)
{
    // Tested againts it's own result, seemed reasonable. Guards against regressions.
    BOOST_CHECK_CLOSE(Hygro::KelvinEquation::VapourPressure(0.0, 273.15), 0.0006112128459398556, eps);
    BOOST_CHECK_CLOSE(Hygro::KelvinEquation::VapourPressure(0.0, 373.15), 0.10141799381792783, eps);
    BOOST_CHECK_CLOSE(Hygro::KelvinEquation::VapourPressure(0.0, 473.15), 1.5549392220497635, eps);
    BOOST_CHECK_CLOSE(Hygro::KelvinEquation::VapourPressure(0.0, 573.15), 8.587867486373652, eps);

    BOOST_CHECK_CLOSE(Hygro::KelvinEquation::VapourPressure(1.0, 573.15), 8.5425925026695655, eps);
    BOOST_CHECK_CLOSE(Hygro::KelvinEquation::VapourPressure(10.0, 573.15), 8.1457090076851948, eps);
    BOOST_CHECK_CLOSE(Hygro::KelvinEquation::VapourPressure(100.0, 573.15), 5.0619855955789381, eps);

    // TODO: check derivatives
}


BOOST_AUTO_TEST_CASE(air_viscosity)
{
    // vapour
    BOOST_CHECK_CLOSE(Hygro::DynamicViscosityOfAir(0.0, 0.1, 273.15), 8.85e-6, eps);
    BOOST_CHECK_CLOSE(Hygro::DynamicViscosityOfAir(0.0, 0.1, 373.15), 1.238e-5, eps);
    BOOST_CHECK_CLOSE(Hygro::DynamicViscosityOfAir(0.0, 0.1, 473.15), 1.591e-5, eps);

    // dry air
    BOOST_CHECK_CLOSE(Hygro::DynamicViscosityOfAir(0.1, 0.1, 273.15), 1.717e-5, eps);
    BOOST_CHECK_CLOSE(Hygro::DynamicViscosityOfAir(0.1, 0.1, 373.15), 2.2122e-5, eps);
    BOOST_CHECK_CLOSE(Hygro::DynamicViscosityOfAir(0.1, 0.1, 473.15), 2.7518e-5, eps);

    // moist air
    BOOST_CHECK_CLOSE(Hygro::DynamicViscosityOfAir(0.1, 0.2, 273.15), 1.4308798819102556e-05, eps);
    BOOST_CHECK_CLOSE(Hygro::DynamicViscosityOfAir(0.1, 0.2, 373.15), 1.8771781021117443e-05, eps);
    BOOST_CHECK_CLOSE(Hygro::DynamicViscosityOfAir(0.1, 0.2, 473.15), 2.352607412165174e-05, eps);
}
