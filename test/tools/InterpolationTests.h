#pragma once
#include "BoostUnitTest.h"
#include <vector>
#include "mechanics/interpolation/InterpolationSimple.h"
#include "TypeTraits.h"

namespace NuTo
{
namespace Test
{

template <typename T>
void CheckCopyMove()
{
    NuTo::Test::Copy<T>();
    NuTo::Test::Move<T>();
}

void CheckPartitionOfUnity(const NuTo::InterpolationSimple& r, const std::vector<Eigen::VectorXd>& rPoints)
{
    for (const auto& point : rPoints)
        BOOST_CHECK_CLOSE(r.GetShapeFunctions(point).sum(), 1, 1.e-10);
}

void CheckDerivativeShapeFunctionsCDF(const NuTo::InterpolationSimple& r, const std::vector<Eigen::VectorXd>& rPoints,
                                      double rDelta = 1.e-6)
{
    for (const auto& point : rPoints)
    {
        auto shapeFunctions = r.GetShapeFunctions(point);
        auto derivativeShapeFunctions = r.GetDerivativeShapeFunctions(point);

        auto cdf = derivativeShapeFunctions;
        cdf.setZero();

        for (int iDim = 0; iDim < point.rows(); ++iDim)
        {
            auto pointDelta = point;
            pointDelta[iDim] += rDelta;
            cdf.col(iDim) = (r.GetShapeFunctions(pointDelta) - shapeFunctions) / rDelta;
        }
        BoostUnitTest::CheckEigenMatrix(cdf, derivativeShapeFunctions, 1.e-6);
    }
}

void CheckShapeFunctionsAndNodePositions(const NuTo::InterpolationSimple& r)
{

    for (int iNode = 0; iNode < r.GetNumNodes(); ++iNode)
    {
        auto localNodeCoordinates = r.GetLocalCoords(iNode);
        auto N = r.GetShapeFunctions(localNodeCoordinates);
        for (int i = 0; i < r.GetNumNodes(); ++i)
        {
            if (i == iNode)
                BOOST_CHECK_CLOSE(N[i], 1, 1.e-10);
            else
                BOOST_CHECK_SMALL(N[i], 1.e-10);
        }
    }
}

template <typename T>
void RunTests(const std::vector<Eigen::VectorXd>& rPoints)
{
    T interpolation(1);

    BOOST_TEST_MESSAGE("Checking copy move...");
    CheckCopyMove<T>();
    BOOST_TEST_MESSAGE("Checking match of shape functions and node positions...");
    CheckShapeFunctionsAndNodePositions(interpolation);
    BOOST_TEST_MESSAGE("Checking patrition of unity...");
    CheckPartitionOfUnity(interpolation, rPoints);
    BOOST_TEST_MESSAGE("Checking shape function derivatives via CDF...");
    CheckDerivativeShapeFunctionsCDF(interpolation, rPoints);
}
} /* Test */
} /* NuTo */
