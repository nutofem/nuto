#include "BoostUnitTest.h"
#include "mechanics/nodes/DofContainer.h"

BOOST_AUTO_TEST_CASE(DofContainerAccess)
{
    NuTo::DofContainer<int> container;

    NuTo::DofType dof0("0", 1, 0);
    NuTo::DofType dof1("1", 1, 1);

    BOOST_CHECK_NO_THROW(container[dof0] = 42);
    BOOST_CHECK_NO_THROW(container[dof1] = 6174);

    BOOST_CHECK_EQUAL(container[dof0], 42);
    BOOST_CHECK_EQUAL(container[dof1], 6174);
}
