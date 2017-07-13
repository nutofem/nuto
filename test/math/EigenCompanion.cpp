#include "BoostUnitTest.h"

#include "math/EigenCompanion.h"
#include "base/Exception.h"

BOOST_AUTO_TEST_CASE(AppendRows)
{
    Eigen::MatrixXd top = Eigen::Matrix<double, 2, 2>::Ones();
    Eigen::MatrixXd bottom = 2 * Eigen::Matrix<double, 1, 2>::Ones();
    NuTo::EigenCompanion::AppendRows(top, bottom);

    Eigen::Matrix<double, 3, 2> expected;
    expected << 1, 1, 1, 1, 2, 2;
    BOOST_CHECK_EQUAL(top, expected);

    Eigen::MatrixXd wrongBottom = Eigen::Matrix<double, 1, 3>::Ones();
    BOOST_CHECK_THROW(NuTo::EigenCompanion::AppendRows(top, wrongBottom), NuTo::Exception);
}

BOOST_AUTO_TEST_CASE(ReadWriteFile)
{
    Eigen::MatrixXd toFile(2, 2);
    toFile << 1, 2, 3, 4;
    NuTo::EigenCompanion::WriteToFile(toFile, "EigenCompanionFile.dat");

    Eigen::MatrixXd fromFile = NuTo::EigenCompanion::ReadFromFile("EigenCompanionFile.dat");

    BOOST_CHECK_CLOSE((toFile - fromFile).norm(), 0, 1.e-10);
}
