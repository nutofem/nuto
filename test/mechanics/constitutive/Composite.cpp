#include "nuto/mechanics/constitutive/staticData/Composite.h"
#include "nuto/mechanics/constitutive/staticData/Leaf.h"

#define BOOST_TEST_MODULE CompositeTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

// necessary to build with clang when boost has been compiled by gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size() ));
}

#include <iostream>

using namespace NuTo::Constitutive::StaticData;

void foo(Component* bla) {}

BOOST_AUTO_TEST_CASE(build_a_simple_tree)
{
    auto& composite = *Composite::Create();
    composite.AddComponent(Leaf<double>::Create(0.0));
    composite.AddComponent(Composite::Create());
    auto& subComposite = dynamic_cast<Composite&>(composite.GetComponent(1));
    Composite* subCompositeCopy = dynamic_cast<Composite&>(composite.GetComponent(1)).Clone();
    subComposite.AddComponent(Leaf<double>::Create(0.0));
    subComposite.AddComponent(Leaf<double>::Create(0.0));
    foo(&subComposite.GetComponent(0));
}
