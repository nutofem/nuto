#include "BoostUnitTest.h"

#include "nuto/mechanics/elements/ElementCollection.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleLinear.h"

BOOST_AUTO_TEST_CASE(ElementCollectionAccess)
{
    NuTo::CoordinateNode n0(Eigen::Vector2d(0, 0));
    NuTo::CoordinateNode n1(Eigen::Vector2d(0, 0));
    NuTo::CoordinateNode n2(Eigen::Vector2d(0, 0));
    NuTo::InterpolationTriangleLinear interpolation;


    NuTo::CoordinateElementFem cElm({n0, n1, n2}, interpolation);
    NuTo::ElementCollectionFem elements(cElm);
    BOOST_CHECK_NO_THROW(elements.CoordinateElement().Interpolation());

    NuTo::ElementCollection& cellElements = elements;
    BOOST_CHECK_EQUAL(cellElements.CoordinateElement().GetDofDimension(), 2);
    BOOST_CHECK_EQUAL(cellElements.CoordinateElement().GetNumNodes(), 3);
}
