#include "BoostUnitTest.h"
#include "TypeTraits.h"
#include "nuto/mechanics/dofs/DofType.h"

BOOST_AUTO_TEST_CASE(DofTypeCopy)
{
    NuTo::Test::Copy<NuTo::DofType>();
}

BOOST_AUTO_TEST_CASE(DofTypeMove)
{
    NuTo::Test::Move<NuTo::DofType>();
}

BOOST_AUTO_TEST_CASE(DofTypeMembers)
{
    NuTo::DofType dof("Special :)", 42);

    BOOST_CHECK(dof.GetName() == std::string("Special :)"));
    BOOST_CHECK(dof.GetNum() == 42);
}

BOOST_AUTO_TEST_CASE(DofTypeEmpty)
{
    BOOST_CHECK_THROW(NuTo::DofType("", 1), NuTo::Exception);
    BOOST_CHECK_THROW(NuTo::DofType("Something", 0), NuTo::Exception);
}


BOOST_AUTO_TEST_CASE(ScalarDofType)
{
    NuTo::ScalarDofType scalarDof("Special :)");

    BOOST_CHECK(scalarDof.GetName() == std::string("Special :)"));
    BOOST_CHECK(scalarDof.GetNum() == 1);

    NuTo::DofType dof(scalarDof);
    BOOST_CHECK(dof.GetName() == std::string("Special :)"));
    BOOST_CHECK(dof.GetNum() == 1);
}
