#include <iostream>
#include "BoostUnitTest.h"
#include "mechanics/PDEs/HygroHelpers.h"

const double eps = 1e-9;

BOOST_AUTO_TEST_CASE(water_density)
{
    // Tested againts it's own result, seemed reasonable. Guards against regressions.
    BOOST_CHECK_CLOSE(Hygro::WaterDensity(273.15), 1015.11, eps);
    BOOST_CHECK_CLOSE(Hygro::WaterDensity(373.15), 962.707, eps);
    BOOST_CHECK_CLOSE(Hygro::WaterDensity(473.15), 867.442, eps);
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
    BOOST_CHECK_CLOSE(Hygro::KelvinEquation(0.0, 273.15), 0.0006112128459398556, eps);
    BOOST_CHECK_CLOSE(Hygro::KelvinEquation(0.0, 373.15), 0.10141799381792783, eps);
    BOOST_CHECK_CLOSE(Hygro::KelvinEquation(0.0, 473.15), 1.5549392220497635, eps);
    BOOST_CHECK_CLOSE(Hygro::KelvinEquation(0.0, 573.15), 8.587867486373652, eps);

    BOOST_CHECK_CLOSE(Hygro::KelvinEquation(1.0, 573.15), 8.542592568927079, eps);
    BOOST_CHECK_CLOSE(Hygro::KelvinEquation(10.0, 573.15), 8.145709639477523, eps);
    BOOST_CHECK_CLOSE(Hygro::KelvinEquation(100.0, 573.15), 5.061989521725555, eps);
}
