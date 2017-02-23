#include "BoostUnitTest.h"
#include "mechanics/nodes/NodeSimple.h"

BOOST_AUTO_TEST_CASE(NodeSimpleValues)
{
    NuTo::NodeSimple node(Eigen::Vector2d{1., 2.});
    auto nodalValue = node.GetValues();
    BOOST_CHECK_CLOSE(nodalValue[0], 1., 1.e-10);
    BOOST_CHECK_CLOSE(nodalValue[1], 2., 1.e-10);

    BOOST_CHECK_EQUAL(node.GetNumValues(), 2);
}

BOOST_AUTO_TEST_CASE(NodeSimpleDofNumbers)
{
    NuTo::NodeSimple node(Eigen::Vector3d{1., 2., 3.});
    BOOST_CHECK_EQUAL(node.GetDofNumber(0), 0);
}
