#include "BoostUnitTest.h"
#include <type_traits>

#include "mechanics/elements/Element.h"
#include "mechanics/interpolation/InterpolationTriangle.h"

BOOST_AUTO_TEST_CASE(ElementCopyMove)
{
    BOOST_CHECK(std::is_copy_constructible<NuTo::Element>::value);
    BOOST_CHECK(std::is_move_constructible<NuTo::Element>::value);
}


struct TestElement : public NuTo::Element
{
    TestElement()
        : NuTo::Element({&n0, &n1, &n2}, interpolation)
    {
    }

    NuTo::NodeSimple n0 = NuTo::NodeSimple(Eigen::Vector2d({1, 1}));
    NuTo::NodeSimple n1 = NuTo::NodeSimple(Eigen::Vector2d({5, 1}));
    NuTo::NodeSimple n2 = NuTo::NodeSimple(Eigen::Vector2d({1, 7}));
    NuTo::InterpolationTriangle interpolation = NuTo::InterpolationTriangle(NuTo::eInterpolation::GAUSS, 1);
};

BOOST_AUTO_TEST_CASE(ElementExtractNodeValues)
{
    auto nodeValues = TestElement().ExtractNodeValues();
    BOOST_CHECK_CLOSE(nodeValues[0], 1, 1.e-10);
    BOOST_CHECK_CLOSE(nodeValues[1], 1, 1.e-10);
    BOOST_CHECK_CLOSE(nodeValues[2], 5, 1.e-10);
    BOOST_CHECK_CLOSE(nodeValues[3], 1, 1.e-10);
    BOOST_CHECK_CLOSE(nodeValues[4], 1, 1.e-10);
    BOOST_CHECK_CLOSE(nodeValues[5], 7, 1.e-10);
}


BOOST_AUTO_TEST_CASE(ElementInterpolate)
{
    auto interpolatedValues = TestElement().Interpolate(Eigen::Vector2d({0.5, 0.5}));
    BOOST_CHECK_CLOSE(interpolatedValues[0], 3., 1.e-10);
    BOOST_CHECK_CLOSE(interpolatedValues[1], 4., 1.e-10);
}
