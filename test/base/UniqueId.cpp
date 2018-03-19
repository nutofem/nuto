#include "BoostUnitTest.h"
#include "nuto/base/UniqueId.h"

BOOST_AUTO_TEST_CASE(IdsCreate)
{
    struct Foo : NuTo::UniqueId<Foo> {};
    for (int i = 0; i < 10; ++i)
        BOOST_CHECK_EQUAL(Foo().Id(), i);
    struct Bar : NuTo::UniqueId<Bar> {};
    for (int i = 0; i < 10; ++i)
        BOOST_CHECK_EQUAL(Bar().Id(), i);
}

BOOST_AUTO_TEST_CASE(IdsCopyMove)
{
    struct FooBase : NuTo::UniqueId<FooBase> {};
    struct Foo : FooBase {};
    Foo foo;
    int id = foo.Id();
    Foo fooCopyCtor(foo); // copy ctor
    BOOST_CHECK_EQUAL(fooCopyCtor.Id(), id);

    Foo fooCopyAssign; // should create a new unique id
    BOOST_CHECK(fooCopyAssign.Id() != id);

    fooCopyAssign = foo; // copy assign copies the id
    BOOST_CHECK_EQUAL(fooCopyAssign.Id(), id);

    Foo fooMoveCtor(std::move(foo));
    BOOST_CHECK_EQUAL(fooMoveCtor.Id(), id);
    
    Foo fooMoveAssign;
    BOOST_CHECK(fooMoveAssign.Id() != id);
    fooMoveAssign = std::move(fooMoveCtor); // foo is invalid after move 
    BOOST_CHECK_EQUAL(fooMoveAssign.Id(), id);
}
