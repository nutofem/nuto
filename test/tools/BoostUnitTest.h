#pragma once
/*
 * collect boost::unittest boilerplate here
 */

#define BOOST_TEST_MODULE DummyTestName
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

// necessary to build with clang when boost has been compiled by gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin() + 1, name.size() - 1) : std::string(name.begin(), name.size()));
}

namespace BoostUnitTest
{
template <typename T1, typename T2>
void CheckVector(const T1& r1, const T2& r2, int rSize, double rTolerance = 1.e-10)
{
    for (int i = 0; i < rSize; ++i)
        BOOST_CHECK_CLOSE(r1[i], r2[i], rTolerance);
}

template <typename T1, typename T2>
void CheckEigenMatrix(const T1& r1, const T2& r2, double rTolerance = 1.e-10)
{
    BOOST_CHECK_EQUAL(r1.rows(), r2.rows());
    BOOST_CHECK_EQUAL(r1.cols(), r2.cols());

    bool isApprox = r1.isApprox(r2);
    if (not isApprox)
    {
        BOOST_TEST_MESSAGE("1st matrix: \n" << r1);
        BOOST_TEST_MESSAGE("2nd matrix: \n" << r2);
    }
    BOOST_CHECK(isApprox);
}
} /* BoostUnitTest */
