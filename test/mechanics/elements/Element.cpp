#include "BoostUnitTest.h"
#include "TypeTraits.h"

#include "mechanics/elements/Element.h"


BOOST_AUTO_TEST_CASE(ElementExtractNodeValues)
{
    NuTo::NodeSimple n0(Eigen::Vector2d{1, 1});
    NuTo::NodeSimple n1(Eigen::Vector2d{5, 1});
    NuTo::NodeSimple n2(Eigen::Vector2d{1, 7});

    NuTo::Element e({&n0, &n1, &n2});

    Eigen::VectorXd nodeValues = e.ExtractNodeValues();
    BOOST_CHECK_CLOSE(nodeValues[0], 1, 1.e-10);
    BOOST_CHECK_CLOSE(nodeValues[1], 1, 1.e-10);
    BOOST_CHECK_CLOSE(nodeValues[2], 5, 1.e-10);
    BOOST_CHECK_CLOSE(nodeValues[3], 1, 1.e-10);
    BOOST_CHECK_CLOSE(nodeValues[4], 1, 1.e-10);
    BOOST_CHECK_CLOSE(nodeValues[5], 7, 1.e-10);
}
