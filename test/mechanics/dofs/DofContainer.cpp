#include "BoostUnitTest.h"
#include "TypeTraits.h"
#include "nuto/mechanics/dofs/DofContainer.h"

BOOST_AUTO_TEST_CASE(DofContainerCopyMove)
{
    NuTo::Test::Copy<NuTo::DofContainer<int>>();
    NuTo::Test::Move<NuTo::DofContainer<int>>();
}

BOOST_AUTO_TEST_CASE(DofContainerAccess)
{
    NuTo::DofContainer<int> container;

    NuTo::DofType dof0("0", 1);
    NuTo::DofType dof1("1", 1);

    BOOST_CHECK_NO_THROW(container[dof0] = 42);
    BOOST_CHECK_NO_THROW(container[dof1] = 6174);

    BOOST_CHECK_EQUAL(container[dof0], 42);
    BOOST_CHECK_EQUAL(container[dof1], 6174);
}

BOOST_AUTO_TEST_CASE(DofContainerInsert)
{
    NuTo::DofContainer<int> container;
    NuTo::DofType dof("0", 1);
    container.Insert(dof, 42);
    BOOST_CHECK_EQUAL(container[dof], 42);
    BOOST_CHECK_THROW(container.Insert(dof, 6174), NuTo::Exception);
    container[dof] = 6174;
    BOOST_CHECK_EQUAL(container[dof], 6174);
}

BOOST_AUTO_TEST_CASE(DofContainerAccessMany)
{
    NuTo::DofContainer<int> container;
    for (int i = 0; i < 100; ++i)
    {
        NuTo::DofType dof(" ", 1);
        BOOST_CHECK_NO_THROW(container[dof]);
    }
}
