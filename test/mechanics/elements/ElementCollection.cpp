#include "BoostUnitTest.h"

#include "nuto/mechanics/elements/ElementCollection.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleLinear.h"

BOOST_AUTO_TEST_CASE(ElementCollectionAccess)
{
    NuTo::NodeSimple n0(Eigen::Vector2d(0, 0));
    NuTo::NodeSimple n1(Eigen::Vector2d(0, 0));
    NuTo::NodeSimple n2(Eigen::Vector2d(0, 0));
    NuTo::InterpolationTriangleLinear interpolation;


    NuTo::ElementCollectionFem elements(NuTo::DofElementFem({n0, n1, n2}, interpolation));
    BOOST_CHECK_NO_THROW(elements.CoordinateElement().Interpolation());

    NuTo::ElementCollection& cellElements = elements;
    BOOST_CHECK_EQUAL(cellElements.CoordinateElement().GetDofDimension(), 2);
    BOOST_CHECK_EQUAL(cellElements.CoordinateElement().GetNumNodes(), 3);
}
