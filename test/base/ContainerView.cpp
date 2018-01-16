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
