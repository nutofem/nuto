#include "BoostUnitTest.h"

#include <iostream>
#include "math/NaturalCoordinateMemoizer.h"

int counter = 0;
Eigen::Vector3d CountedFunction(const Eigen::Vector3d& v)
{
    counter++;
    return v;
}

struct Lobatto
{
    std::vector<Eigen::Vector3d> ips;
    Lobatto()
    {
        std::vector<double> ip1D = {-1., -0.654653670707977087, 0., +0.654653670707977087, +1.};
        ips.reserve(125);
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 5; j++)
                for (int k = 0; k < 5; k++)
                    ips.push_back({ip1D[i], ip1D[j], ip1D[k]});
    }
};

template <typename TMemoizer>
void CheckMemoizer()
{
    counter = 0;
    TMemoizer memo(CountedFunction);

    Lobatto lobatto;
    for (int i = 0; i < 12; ++i)
        for (const auto& ip : lobatto.ips)
            BoostUnitTest::CheckVector(memo.Get(ip), ip, 3);

    // even though the method is called 12 times, it should only
    // be evaluated once per integration point
    int expectedCounter = lobatto.ips.size();
    BOOST_CHECK_EQUAL(counter, expectedCounter);

    memo.ClearCache();
    for (int i = 0; i < 12; ++i)
        for (const auto& ip : lobatto.ips)
            BoostUnitTest::CheckVector(memo.Get(ip), ip, 3);

    expectedCounter += lobatto.ips.size();
    BOOST_CHECK_EQUAL(counter, expectedCounter);
}

BOOST_AUTO_TEST_CASE(memoizationMap)
{
    using Memoizer = NuTo::NaturalCoordinateMemoizerMap<Eigen::Vector3d, Eigen::Vector3d>;
    CheckMemoizer<Memoizer>();
}
