#include "BoostUnitTest.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include <eigen3/Eigen/Dense> // for ::determinant()

using namespace NuTo;


void CheckInvariants(EngineeringStrain<3> r)
{

    Eigen::Matrix3d m;
    m << r[0], .5 * r[5], .5 * r[4], .5 * r[5], r[1], .5 * r[3], .5 * r[4], .5 * r[3], r[2];

    BOOST_CHECK_CLOSE(r.InvariantI1(), m.trace(), 1.e-10);

    BOOST_CHECK_CLOSE(r.InvariantI2(), .5 * (m.trace() * m.trace() - (m * m).trace()), 1.e-10);

    BOOST_CHECK_CLOSE(r.InvariantI3(), m.determinant(), 1.e-10);

    // note: the invariants are not equal but with opposite sign
    BOOST_CHECK_CLOSE(r.InvariantJ2(), -r.Deviatoric().InvariantI2(), 1.e-10);

    BOOST_CHECK_CLOSE(r.InvariantJ2(), 1. / 3. * r.InvariantI1() * r.InvariantI1() - r.InvariantI2(), 1.e-10);
}


BOOST_AUTO_TEST_CASE(strainInvariants)
{
    EngineeringStrain<1> e1D;
    EngineeringStrain<2> e2D;
    EngineeringStrain<3> e3D;

    e1D.AsVector() = Eigen::Matrix<double, 1, 1>::Random();
    e2D.AsVector() = Eigen::Matrix<double, 3, 1>::Random();
    e3D.AsVector() = Eigen::Matrix<double, 6, 1>::Random();

    CheckInvariants(e1D.As3D(0.3));

    CheckInvariants(e2D.As3D(0.3, ePlaneState::PLANE_STRAIN));
    CheckInvariants(e2D.As3D(0.3, ePlaneState::PLANE_STRESS));

    CheckInvariants(e3D);
}


BOOST_AUTO_TEST_CASE(ListConstructor)
{
    EngineeringStrain<1> a({1.});
    EngineeringStrain<2> b({1., 2., 3.});
    EngineeringStrain<3> c({1., 2., 3., 4., 5., 6.});

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
