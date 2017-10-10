#include "BoostUnitTest.h"

#include "mechanics/elements/ElementCollection.h"
#include "mechanics/interpolation/InterpolationTriangleLinear.h"

BOOST_AUTO_TEST_CASE(ElementCollectionAccess)
{
    NuTo::NodeSimple n0(Eigen::Vector2d(0,0));
    NuTo::NodeSimple n1(Eigen::Vector2d(0,0));
    NuTo::NodeSimple n2(Eigen::Vector2d(0,0));
    NuTo::InterpolationTriangleLinear interpolation(2);


    NuTo::ElementCollectionFem elements(NuTo::ElementFem({&n0, &n1, &n2}, interpolation));
    BOOST_CHECK_NO_THROW(elements.CoordinateElement().Interpolation());

    NuTo::ElementCollection& cellElements = elements;
    BOOST_CHECK_EQUAL(cellElements.CoordinateElement().GetDofDimension(), 2);
    BOOST_CHECK_EQUAL(cellElements.CoordinateElement().GetNumNodes(), 3);
}

