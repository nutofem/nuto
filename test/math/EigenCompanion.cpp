#define BOOST_TEST_MODULE EigenCompanion
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "math/EigenCompanion.h"
#include "math/MathException.h"

// necessary to build with clang when boost has been compiled by gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size() ));
}

BOOST_AUTO_TEST_CASE(AppendRows)
{
    Eigen::MatrixXd top = Eigen::Matrix<double, 2, 2>::Ones();
    Eigen::MatrixXd bottom = 2*Eigen::Matrix<double, 1, 2>::Ones();
    NuTo::EigenCompanion::AppendRows(top, bottom);

    Eigen::Matrix<double, 3, 2> expected;
    expected << 1, 1, 1, 1, 2, 2;
    BOOST_CHECK_EQUAL(top, expected);

    Eigen::MatrixXd wrongBottom = Eigen::Matrix<double, 1, 3>::Ones();
    BOOST_CHECK_THROW(NuTo::EigenCompanion::AppendRows(top, wrongBottom), NuTo::MathException);
}
