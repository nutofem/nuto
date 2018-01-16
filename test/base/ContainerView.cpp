#include "BoostUnitTest.h"
#include "base/ContainerView.h"
#include <vector>
#include <set>
#include <algorithm>

BOOST_AUTO_TEST_CASE(FromVector)
{
    std::vector<int> v{1, 2, 3, 4};
    NuTo::ContainerView<int> cv(v);

    BOOST_CHECK(std::equal(cv.begin(), cv.end(), v.begin()));
}
BOOST_AUTO_TEST_CASE(FromSet)
{
    std::set<int> v{1, 2, 3, 4};
    NuTo::ContainerView<int> cv(v);

    BOOST_CHECK(std::equal(cv.begin(), cv.end(), v.begin()));
}

struct Interface
{
    virtual ~Interface() = default;
    virtual int Value() const = 0;
    bool operator==(const Interface& other)
    {
        return Value() == other.Value();
    }
};
struct Impl : Interface
{
    Impl(int value)
        : mValue(value)
    {
    }
    int Value() const override
    {
        return mValue;
    }
    int mValue;
};

BOOST_AUTO_TEST_CASE(FromBase)
{
    std::vector<Impl> v{1, 2, 3};
    NuTo::ContainerView<Interface> cv(v);
    BOOST_CHECK(std::equal(cv.begin(), cv.end(), v.begin()));
}

void AcceptingInitList(NuTo::ContainerView<Interface> view)
{
    BOOST_CHECK_EQUAL(view.begin()->Value(), 1);
}

struct WithInitList
{
    WithInitList(NuTo::ContainerView<Interface> view)
        : mView(view)
    {
    }

    void Check()
    {
        BOOST_CHECK_EQUAL(mView.begin()->Value(), 1);
    }

    NuTo::ContainerView<Interface> mView;
};

BOOST_AUTO_TEST_CASE(FromInitializerList)
{
    Impl a(1), b(2), c(3);
    AcceptingInitList({a, b, c});

    // WithInitList w({a, b, c});
    // w.Check();

    // NuTo::ContainerView<Interface> view({a, b, c});
    // BOOST_CHECK_EQUAL(view.begin()->Value(), 1);
}
