#include "BoostUnitTest.h"
#include "nuto/math/EigenIO.h"
#include "nuto/base/Exception.h"

BOOST_AUTO_TEST_CASE(ReadWriteFile)
{
    Eigen::MatrixXd toFile(2, 2);
    toFile << 1, 2, 3, 4;
    NuTo::EigenIO::WriteToFile(toFile, "EigenCompanionFile.dat");

    Eigen::MatrixXd fromFile = NuTo::EigenIO::ReadFromFile("EigenCompanionFile.dat");

    BOOST_CHECK_CLOSE((toFile - fromFile).norm(), 0, 1.e-10);

    BOOST_CHECK_THROW(NuTo::EigenIO::WriteToFile(toFile, "/foo/bar.dat"), NuTo::Exception);
    BOOST_CHECK_THROW(NuTo::EigenIO::ReadFromFile("/foo/bar.dat"), NuTo::Exception);
}
