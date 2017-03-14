#include "BoostUnitTest.h"
#include "mechanics/nodes/DofContainer.h"

//! @brief check access with dof index > 10 (magic number)
BOOST_AUTO_TEST_CASE(DofContainerAccess)
{
    NuTo::DofContainer<int> container;
    for (int i = 0; i < 100; ++i)
    {
        NuTo::DofType dof(" ", 1, i);
        BOOST_CHECK_NO_THROW(container[dof]);
    }
}
