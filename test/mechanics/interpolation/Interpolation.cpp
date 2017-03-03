#include "BoostUnitTest.h"
#include "mechanics/interpolation/Interpolation.h"
#include "mechanics/MechanicsException.h"


class MockInterpolation : public NuTo::Interpolation
{
public:
    MockInterpolation(int rDofDimension)
        : NuTo::Interpolation(NuTo::eInterpolation::GAUSS, 1, rDofDimension)
    {
    }
    Eigen::VectorXd GetShapeFunctions(const Eigen::VectorXd& rLocalIPCoords) const override
    {
        return Eigen::Vector3d({42, 6174, 12});
    }
    Eigen::VectorXd GetLocalCoords(int rNodeId) const override
    {
        return Eigen::Vector3d::Random();
    }
    int GetNumNodes() const override
    {
        return 3;
    }
};

class MockInterpolation1D : public MockInterpolation
{
public:
    MockInterpolation1D(int rDofDimension)
        : MockInterpolation(rDofDimension)
    {
    }
    Eigen::MatrixXd GetDerivativeShapeFunctions(const Eigen::VectorXd& rLocalIPCoords) const override
    {
        Eigen::MatrixXd m(3, 1);
        m(0, 0) = 0;
        m(1, 0) = 10;
        m(2, 0) = 20;
        return m;
    }
};

class MockInterpolation2D : public MockInterpolation
{
public:
    MockInterpolation2D(int rDofDimension)
        : MockInterpolation(rDofDimension)
    {
    }
    Eigen::MatrixXd GetDerivativeShapeFunctions(const Eigen::VectorXd& rLocalIPCoords) const override
    {
        Eigen::MatrixXd m(3, 2);
        m(0, 0) = 0;
        m(0, 1) = 1;
        m(1, 0) = 10;
        m(1, 1) = 11;
        m(2, 0) = 20;
        m(2, 1) = 21;
        return m;
    }
};

class MockInterpolation3D : public MockInterpolation
{
public:
    MockInterpolation3D(int rDofDimension)
        : MockInterpolation(rDofDimension)
    {
    }
    Eigen::MatrixXd GetDerivativeShapeFunctions(const Eigen::VectorXd& rLocalIPCoords) const override
    {
        Eigen::MatrixXd m(3, 3);
        m(0, 0) = 0;
        m(0, 1) = 1;
        m(0, 2) = 2;
        m(1, 0) = 10;
        m(1, 1) = 11;
        m(1, 2) = 12;
        m(2, 0) = 20;
        m(2, 1) = 21;
        m(2, 2) = 22;
        return m;
    }
};


void TestN(int rDim)
{
    MockInterpolation1D p(rDim);
    Eigen::MatrixXd expected = Eigen::MatrixXd::Zero(rDim, 3 * rDim);
    for (int iDim = 0; iDim < rDim; ++iDim)
    {
        expected(iDim, iDim)               = 42;
        expected(iDim, iDim + rDim)        = 6174;
        expected(iDim, iDim + rDim + rDim) = 12;
    }
    BoostUnitTest::CheckEigenMatrix(p.GetN(Eigen::VectorXd()), expected);
}

BOOST_AUTO_TEST_CASE(InterpolationN)
{
    for (int dimension = 0; dimension < 3; ++dimension)
        TestN(dimension);
}

BOOST_AUTO_TEST_CASE(InterpolationBGradient)
{
    // dof dimension == 2, only defined for scalars
    BOOST_CHECK_THROW(MockInterpolation2D(2).GetBGradient(Eigen::VectorXd()), NuTo::MechanicsException);


    Eigen::MatrixXd expected(2, 3);
    expected(0, 0) = 0;
    expected(0, 1) = 10;
    expected(0, 2) = 20;
    expected(1, 0) = 1;
    expected(1, 1) = 11;
    expected(1, 2) = 21;
    MockInterpolation2D p(1);
    BoostUnitTest::CheckEigenMatrix(p.GetBGradient(Eigen::VectorXd()), expected);
}

BOOST_AUTO_TEST_CASE(InterpolationBStrain1D)
{
    MockInterpolation1D p(1);
    BoostUnitTest::CheckEigenMatrix(p.GetBStrain(Eigen::VectorXd()), p.GetBGradient(Eigen::VectorXd()));
}

BOOST_AUTO_TEST_CASE(InterpolationBStrain2D)
{
    Eigen::MatrixXd expected = Eigen::MatrixXd::Zero(3, 6);
    expected(0, 0) = 0;
    expected(1, 1) = 1;
    expected(2, 0) = 1;
    expected(2, 1) = 0;
    
    expected(0, 2) = 10;
    expected(1, 3) = 11;
    expected(2, 2) = 11;
    expected(2, 3) = 10;
   
    expected(0, 4) = 20;
    expected(1, 5) = 21;
    expected(2, 4) = 21;
    expected(2, 5) = 20;
    MockInterpolation2D p(2);
    BoostUnitTest::CheckEigenMatrix(p.GetBStrain(Eigen::VectorXd()), expected);
}

BOOST_AUTO_TEST_CASE(InterpolationBStrain3D)
{
    Eigen::MatrixXd expected = Eigen::MatrixXd::Zero(6, 9);
    expected(0, 0) = 0;
    expected(1, 1) = 1;
    expected(2, 2) = 2;
    expected(3, 1) = 2;
    expected(3, 2) = 1;
    expected(4, 0) = 2;
    expected(4, 2) = 0;
    expected(5, 0) = 1;
    expected(5, 1) = 0;
    
    expected(0, 3) = 10;
    expected(1, 4) = 11;
    expected(2, 5) = 12;
    expected(3, 4) = 12;
    expected(3, 5) = 11;
    expected(4, 3) = 12;
    expected(4, 5) = 10;
    expected(5, 3) = 11;
    expected(5, 4) = 10;
    
    expected(0, 6) = 20;
    expected(1, 7) = 21;
    expected(2, 8) = 22;
    expected(3, 7) = 22;
    expected(3, 8) = 21;
    expected(4, 6) = 22;
    expected(4, 8) = 20;
    expected(5, 6) = 21;
    expected(5, 7) = 20;
    
    MockInterpolation3D p(3);
    BoostUnitTest::CheckEigenMatrix(p.GetBStrain(Eigen::VectorXd()), expected);
}

