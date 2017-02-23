#include "BoostUnitTest.h"
#include "TypeTraits.h"
#include "mechanics/nodes/DofType.h"

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
    NuTo::DofType dof("SpacialDof :)", 42, 0);

    BOOST_CHECK(dof.GetName() == std::string("SpacialDof :)"));
    BOOST_CHECK(dof.GetNum() == 42);
}


BOOST_AUTO_TEST_CASE(DofTypeEmpty)
{
    BOOST_CHECK_THROW(NuTo::DofType("", 1, 0), NuTo::MechanicsException);
    BOOST_CHECK_THROW(NuTo::DofType("Something", 0, 0), NuTo::MechanicsException);
}
