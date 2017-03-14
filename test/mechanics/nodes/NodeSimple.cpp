#include "BoostUnitTest.h"
#include "TypeTraits.h"
#include "mechanics/nodes/NodeSimple.h"


BOOST_AUTO_TEST_CASE(NodeSimpleCopyMove)
{
    NuTo::Test::Copy<NuTo::NodeSimple>();
    NuTo::Test::Move<NuTo::NodeSimple>();
}

BOOST_AUTO_TEST_CASE(NodeSimpleValues)
{
    NuTo::NodeSimple node(Eigen::Vector2d{1., 2.});
    BOOST_CHECK_EQUAL(node.GetNumValues(), 2);

    BoostUnitTest::CheckVector(node.GetValues(), std::vector<double>({1., 2.}), 2);

    node.SetValue(0, 11.);
    node.SetValue(1, 12.);
    BoostUnitTest::CheckVector(node.GetValues(), std::vector<double>({11., 12.}), 2);
}

BOOST_AUTO_TEST_CASE(NodeSimpleDofNumbers)
{
    NuTo::NodeSimple node(Eigen::Vector3d{1., 2., 3.});
    BOOST_CHECK_EQUAL(node.GetDofNumber(0), 0);

    node.SetDofNumber(1, 42);
    BOOST_CHECK_EQUAL(node.GetDofNumber(1), 42);
}
