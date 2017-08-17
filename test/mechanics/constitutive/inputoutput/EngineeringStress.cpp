#include "BoostUnitTest.h"
#include "mechanics/constitutive/inputoutput/EngineeringStress.h"
#include <cfenv> // catching NaN

using namespace NuTo;

BOOST_AUTO_TEST_CASE(ListConstructor)
{
    EngineeringStress<1> a({1.});
    EngineeringStress<2> b({1., 2., 3.});
    EngineeringStress<3> c({1., 2., 3., 4., 5., 6.});

    BOOST_CHECK_EQUAL(a[0], 1.0);

    BOOST_CHECK_EQUAL(b[0], 1.0);
    BOOST_CHECK_EQUAL(b[1], 2.0);
    BOOST_CHECK_EQUAL(b[2], 3.0);

    BOOST_CHECK_EQUAL(c[0], 1.0);
    BOOST_CHECK_EQUAL(c[1], 2.0);
    BOOST_CHECK_EQUAL(c[2], 3.0);
    BOOST_CHECK_EQUAL(c[3], 4.0);
    BOOST_CHECK_EQUAL(c[4], 5.0);
    BOOST_CHECK_EQUAL(c[5], 6.0);
}

BOOST_AUTO_TEST_CASE(As3D)
{
    EngineeringStress<1> a({1.});
    EngineeringStress<2> b({1., 2., 3.});
    EngineeringStress<3> c({1., 2., 3., 4., 5., 6.});

    BoostUnitTest::CheckVector(a.As3D(), std::vector<double>{1., 0., 0., 0., 0., 0.}, 6);
    BoostUnitTest::CheckVector(b.As3D(), std::vector<double>{1., 2., 0., 0., 0., 3.}, 6);
    BoostUnitTest::CheckVector(c.As3D(), std::vector<double>{1., 2., 3., 4., 5., 6.}, 6);
}

//! @remark This is rather poor. Better than nothing though.
BOOST_AUTO_TEST_CASE(Mises)
{
    feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT); // Enable all floating point exceptions but FE_INEXACT

    EngineeringStress<1> a({5.});
    EngineeringStress<2> b({5., 0., 0.});
    EngineeringStress<3> c({5., 0., 0., 0., 0., 0.});

    BOOST_CHECK_CLOSE(a.VonMisesStress(), 5., 1.e-10);
    BOOST_CHECK_CLOSE(b.VonMisesStress(), 5., 1.e-10);
    BOOST_CHECK_CLOSE(c.VonMisesStress(), 5., 1.e-10);

    BOOST_CHECK_CLOSE(EngineeringStress<2>({-2., -2., 0}).VonMisesStress(), 2., 1.e-10);


    BOOST_CHECK_SMALL(EngineeringStress<1>({0}).VonMisesStress(), 1.e-10);
    BOOST_CHECK_SMALL(EngineeringStress<2>({0, 0, 0}).VonMisesStress(), 1.e-10);
    BOOST_CHECK_SMALL(EngineeringStress<3>({0, 0, 0, 0, 0, 0}).VonMisesStress(), 1.e-10);
}

//! @remark This is rather poor. Better than nothing though.
BOOST_AUTO_TEST_CASE(Rankine)
{
    feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT); // Enable all floating point exceptions but FE_INEXACT

    EngineeringStress<1> a({5.});
    EngineeringStress<2> b({5., 0., 0.});
    EngineeringStress<3> c({5., 0., 0., 0., 0., 0.});

    BOOST_CHECK_CLOSE(a.SmoothRankine(), 5., 1.e-10);
    BOOST_CHECK_CLOSE(b.SmoothRankine(), 5., 1.e-10);
    BOOST_CHECK_CLOSE(c.SmoothRankine(), 5., 1.e-10);

    BOOST_CHECK_CLOSE(EngineeringStress<2>({-2., -2., 0}).SmoothRankine(), 0., 1.e-10);

    BOOST_CHECK_SMALL(EngineeringStress<1>({0}).SmoothRankine(), 1.e-10);
    BOOST_CHECK_SMALL(EngineeringStress<2>({0, 0, 0}).SmoothRankine(), 1.e-10);
    BOOST_CHECK_SMALL(EngineeringStress<3>({0, 0, 0, 0, 0, 0}).SmoothRankine(), 1.e-10);
}
