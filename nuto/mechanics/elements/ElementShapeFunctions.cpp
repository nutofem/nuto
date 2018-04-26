/*
 * ElementShapeFunctions.cpp
 *
 *  Created on: 30 Mar 2015
 *      Author: ttitsche
 */
#include <cassert>
#include "nuto/mechanics/elements/ElementShapeFunctions.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleQuadratic.h"
#include "nuto/base/Exception.h"

namespace NuTo
{

namespace ShapeFunctionsInterface2D // interval -1 to 1
{
////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::MatrixXd NodeCoordinatesInterface2dOrder1(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0:
        return Eigen::Vector2d(-1, -1);
    case 1:
        return Eigen::Vector2d(+1, -1);
    case 2:
        return Eigen::Vector2d(+1, +1);
    case 3:
        return Eigen::Vector2d(-1, +1);
    default:
        throw NuTo::Exception(std::string(__PRETTY_FUNCTION__) + ":\t node index out of range (0..3)");
    }
}

Eigen::MatrixXd ShapeFunctionsInterface2dOrder1(const Eigen::VectorXd& rCoordinates)
{
    const double N00 = 0.5 * (1. - rCoordinates[0]);
    const double N01 = 0.5 * (1. + rCoordinates[0]);

    return Eigen::Vector4d(-N00, -N01, N01, N00);
}

Eigen::MatrixXd DerivativeShapeFunctionsInterface2dOrder1()
{
    // this interface element does not need any shape function derivatives
    return Eigen::Matrix2d::Zero();
}

Eigen::MatrixXd NodeCoordinatesInterface2dOrder2(int rNodeIndex)
{
    switch (rNodeIndex)
    {
    case 0:
        return Eigen::Vector2d(-1., -1.);
    case 1:
        return Eigen::Vector2d(+0., -1);
    case 2:
        return Eigen::Vector2d(+1, -1);
    case 3:
        return Eigen::Vector2d(+1, +1);
    case 4:
        return Eigen::Vector2d(0, +1);
    case 5:
        return Eigen::Vector2d(-1, +1);
    default:
        throw NuTo::Exception(std::string(__PRETTY_FUNCTION__) + ":\t node index out of range (0..5)");
    }
}

Eigen::MatrixXd ShapeFunctionsInterface2dOrder2(const Eigen::VectorXd& rCoordinates)
{
    const double xi = rCoordinates(0, 0);

    const double N00 = 0.5 * xi * (xi - 1.);
    const double N01 = (1. - xi * xi);
    const double N02 = 0.5 * xi * (xi + 1.);

    return (Eigen::MatrixXd(6, 1) << N00, N01, N02, -N02, -N01, -N00).finished();
}

Eigen::MatrixXd DerivativeShapeFunctionsInterface2dOrder2(const Eigen::VectorXd&)
{
    // this interface element does not need any shape function derivatives
    return Eigen::Matrix2d::Zero();
}
}

// In order to maintain equal shape functions on each element the Bézier extraction together with Bernstein polynomials
// is used (see. Borden et. al. 2011)
namespace ShapeFunctionsIGA
{

/////////////////////////////// BSPLINE ////////////////////////////////////////////////////////////////////////

int FindSpan(double rParameter, int rDegree, const Eigen::VectorXd& rKnots)
{
    int size = rKnots.rows();
    assert(rParameter >= rKnots(0) && rParameter <= rKnots(size - 1));
    int numBasisFuns = size - rDegree - 1;
    if (rParameter == rKnots[numBasisFuns])
        return numBasisFuns - 1;

    int low = rDegree;
    int high = numBasisFuns;
    int mid = (low + high) / 2;

    while (rParameter < rKnots[mid] || rParameter >= rKnots[mid + 1])
    {
        if (rParameter < rKnots[mid])
            high = mid;
        else
            low = mid;

        mid = (low + high) / 2;
    }
    return mid;
}

Eigen::VectorXd BasisFunctions(double rParameter, int spanIdx, int rDegree, const Eigen::VectorXd& rKnots)
{
    Eigen::VectorXd rBasisFunctions(rDegree + 1);

    rBasisFunctions[0] = 1.;

    Eigen::VectorXd left(rDegree + 1);
    Eigen::VectorXd right(rDegree + 1);

    for (int j = 1; j <= rDegree; j++)
    {
        left[j] = rParameter - rKnots[spanIdx + 1 - j];
        right[j] = rKnots[spanIdx + j] - rParameter;
        double saved = 0.;
        for (int r = 0; r < j; r++)
        {
            double temp = rBasisFunctions[r] / (right[r + 1] + left[j - r]);
            rBasisFunctions[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        rBasisFunctions[j] = saved;
    }

    return rBasisFunctions;
}

Eigen::MatrixXd BasisFunctionsAndDerivatives(int der, double rParameter, int spanIdx, int rDegree,
                                             const Eigen::VectorXd& rKnots)
{
    Eigen::MatrixXd ndu(rDegree + 1, rDegree + 1);

    ndu(0, 0) = 1.0;

    Eigen::VectorXd left(rDegree + 1);
    Eigen::VectorXd right(rDegree + 1);

    for (int j = 1; j <= rDegree; j++)
    {
        left[j] = rParameter - rKnots[spanIdx + 1 - j];
        right[j] = rKnots[spanIdx + j] - rParameter;
        double saved = 0.;
        for (int r = 0; r < j; r++)
        {
            ndu(j, r) = right[r + 1] + left[j - r]; // lower triange (knot differences)
            double temp = ndu(r, j - 1) / ndu(j, r);
            ndu(r, j) = saved + right[r + 1] * temp; // upper triangle
            saved = left[j - r] * temp;
        }
        ndu(j, j) = saved;
    }

    Eigen::MatrixXd ders(der + 1, rDegree + 1);

    for (int j = 0; j <= rDegree; j++)
        ders(0, j) = ndu(j, rDegree);

    Eigen::MatrixXd a(2, rDegree + 1);
    for (int r = 0; r <= rDegree; r++)
    {
        int s1 = 0;
        int s2 = 1;

        a(0, 0) = 1.0;

        int j = 0;
        for (int k = 1; k <= der; k++)
        {
            double d = 0.;
            int rk = r - k;
            int pk = rDegree - k;
            if (r >= k)
            {
                a(s2, 0) = a(s1, 0) / ndu(pk + 1, rk);
                d = a(s2, 0) * ndu(rk, pk);
            }
            int j1 = 0;
            if (rk >= -1)
                j1 = 1;
            else
                j1 = -rk;

            int j2 = 0;
            if (r - 1 <= pk)
                j2 = k - 1;
            else
                j2 = rDegree - r;

            for (j = j1; j <= j2; j++)
            {
                a(s2, j) = (a(s1, j) - a(s1, j - 1)) / ndu(pk + 1, rk + j);
                d += a(s2, j) * ndu(rk + j, pk);
            }

            if (r <= pk)
            {
                a(s2, k) = -a(s1, k - 1) / ndu(pk + 1, r);
                d += a(s2, k) * ndu(r, pk);
            }

            ders(k, r) = d;
            j = s1;
            s1 = s2;
            s2 = j;
        }
    }

    int r = rDegree;
    for (int k = 1; k <= der; k++)
    {
        for (int j = 0; j <= rDegree; j++)
            ders(k, j) *= r;
        r *= (rDegree - k);
    }


    return ders;
}

Eigen::VectorXd BasisFunctionsRat(double rParameter, int spanIdx, int rDegree, const Eigen::VectorXd& rKnots,
                                  const Eigen::VectorXd& rWeights)
{
    Eigen::VectorXd rBasisFunctions = BasisFunctions(rParameter, spanIdx, rDegree, rKnots);

    // NURBS specific ...
    double sum = 0.;
    for (int i = 0; i <= rDegree; i++)
        sum += rBasisFunctions(i) * rWeights(spanIdx - rDegree + i);
    for (int i = 0; i <= rDegree; i++)
        rBasisFunctions(i) *= (rWeights(spanIdx - rDegree + i) / sum);

    return rBasisFunctions;
}

Eigen::VectorXd BasisFunctionsAndDerivativesRat(int der, double rParameter, int spanIdx, int rDegree,
                                                const Eigen::VectorXd& rKnots, const Eigen::VectorXd& rWeights)
{
    if (der < 0 || der > 2)
        throw NuTo::Exception(std::string(__PRETTY_FUNCTION__) +
                              ":\t der greater than 2 not implemented, possible values 0,1,2!");

    Eigen::MatrixXd ders = BasisFunctionsAndDerivatives(der, rParameter, spanIdx, rDegree, rKnots);

    // NURBS specific ...
    Eigen::VectorXd sum(der + 1);
    sum.setZero(der + 1);
    for (int numDer = 0; numDer < der; numDer++)
    {
        for (int i = 0; i <= rDegree; i++)
        {
            sum(numDer) += ders(numDer, i) * rWeights(spanIdx - rDegree + i);
        }
    }


    Eigen::VectorXd dersRat(rDegree + 1);
    dersRat.setZero(rDegree + 1);

    for (int i = 0; i <= rDegree; i++)
    {
        double weight = rWeights(spanIdx - rDegree + i);
        if (der == 0)
        {
            dersRat(i) = ders(der, i) * weight / sum(0);
        }
        else if (der == 1)
        {
            dersRat(i) = (ders(der, i) * sum(0) - ders(0, i) * sum(1)) * weight / (sum(0) * sum(0));
        }
        else
        {
            double sum2 = sum(0) * sum(0);
            dersRat(i) =
                    weight * (ders(der, i) / sum(0) - 2 * ders(1, i) * sum(1) / (sum2)-ders(0, i) * sum(2) / (sum2) +
                              2 * ders(0, i) * sum(1) * sum(1) / (sum2 * sum(0)));
        }
    }

    return dersRat;
}

Eigen::VectorXd BasisFunctions2DRat(const Eigen::VectorXd& rCoordinates, const Eigen::Vector2i& rSpanIdx,
                                    const Eigen::Vector2i& rDegree, const Eigen::VectorXd& rKnotsX,
                                    const Eigen::VectorXd& rKnotsY, const Eigen::MatrixXd& rWeights)
{
    Eigen::VectorXd xBasis = BasisFunctions(rCoordinates(0), rSpanIdx(0), rDegree(0), rKnotsX);
    Eigen::VectorXd yBasis = BasisFunctions(rCoordinates(1), rSpanIdx(1), rDegree(1), rKnotsY);

    Eigen::VectorXd basis((rDegree(0) + 1) * (rDegree(1) + 1));
    basis.setZero((rDegree(0) + 1) * (rDegree(1) + 1));

    double sum = 0.;
    for (int i = 0; i <= rDegree(1); i++)
    {
        for (int j = 0; j <= rDegree(0); j++)
        {
            double weight = rWeights(rSpanIdx(1) - rDegree(1) + i, rSpanIdx(0) - rDegree(0) + j);
            sum += xBasis(j) * yBasis(i) * weight;
        }
    }

    int count = 0;
    for (int i = 0; i <= rDegree(1); i++)
    {
        for (int j = 0; j <= rDegree(0); j++)
        {
            double weight = rWeights(rSpanIdx(1) - rDegree(1) + i, rSpanIdx(0) - rDegree(0) + j);
            basis(count) = xBasis(j) * yBasis(i) * weight / sum;
            count++;
        }
    }

    return basis;
}

Eigen::MatrixXd BasisFunctionsAndDerivatives2DRat(int der, const Eigen::VectorXd& rCoordinates,
                                                  const Eigen::Vector2i& rSpanIdx, const Eigen::Vector2i& rDegree,
                                                  const Eigen::VectorXd& rKnotsX, const Eigen::VectorXd& rKnotsY,
                                                  const Eigen::MatrixXd& rWeights)
{
    if (der < 0 || der > 2)
        throw NuTo::Exception(std::string(__PRETTY_FUNCTION__) +
                              ":\t 'der' greater than 2 not implemented, possible values 0,1,2!");

    Eigen::MatrixXd xBasisDer = BasisFunctionsAndDerivatives(der, rCoordinates(0), rSpanIdx(0), rDegree(0), rKnotsX);
    Eigen::MatrixXd yBasisDer = BasisFunctionsAndDerivatives(der, rCoordinates(1), rSpanIdx(1), rDegree(1), rKnotsY);

    Eigen::MatrixXd ders((rDegree(0) + 1) * (rDegree(1) + 1), der + 1);
    ders.setZero((rDegree(0) + 1) * (rDegree(1) + 1), der + 1);


    int num = ((der + 1) * (der + 1) + (der + 1)) / 2;
    Eigen::VectorXd sum(num);
    sum.setZero(num);

    for (int i = 0; i <= rDegree(1); i++)
    {
        for (int j = 0; j <= rDegree(0); j++)
        {
            double weight = rWeights(rSpanIdx(1) - rDegree(1) + i, rSpanIdx(0) - rDegree(0) + j);

            if (der == 0)
            {
                sum(0) += xBasisDer(0, j) * yBasisDer(0, i) * weight;
            }
            else if (der == 1)
            {
                sum(0) += xBasisDer(0, j) * yBasisDer(0, i) * weight;
                sum(1) += xBasisDer(1, j) * yBasisDer(0, i) * weight;
                sum(2) += xBasisDer(0, j) * yBasisDer(1, i) * weight;
            }
            else
            {
                sum(0) += xBasisDer(0, j) * yBasisDer(0, i) * weight;
                sum(1) += xBasisDer(1, j) * yBasisDer(0, i) * weight;
                sum(2) += xBasisDer(2, j) * yBasisDer(0, i) * weight;

                sum(3) += xBasisDer(0, j) * yBasisDer(1, i) * weight;
                sum(4) += xBasisDer(0, j) * yBasisDer(2, i) * weight;

                sum(5) += xBasisDer(1, j) * yBasisDer(1, i) * weight;
            }
        }
    }

    int count = 0;
    for (int i = 0; i <= rDegree(1); i++)
    {
        for (int j = 0; j <= rDegree(0); j++)
        {
            double weight = rWeights(rSpanIdx(1) - rDegree(1) + i, rSpanIdx(0) - rDegree(0) + j);

            if (der == 0)
            {
                ders(count, 0) = xBasisDer(0, j) * yBasisDer(0, i) * weight / sum(0);
                count++;
            }
            else if (der == 1)
            {
                ders(count, 0) =
                        (xBasisDer(1, j) * yBasisDer(0, i) * sum(0) - xBasisDer(0, j) * yBasisDer(0, i) * sum(1)) *
                        weight / (sum(0) * sum(0));
                ders(count, 1) =
                        (xBasisDer(0, j) * yBasisDer(1, i) * sum(0) - xBasisDer(0, j) * yBasisDer(0, i) * sum(2)) *
                        weight / (sum(0) * sum(0));
                count++;
            }
            else
            {
                ders(count, 0) =
                        (xBasisDer(2, j) * yBasisDer(0, i) / sum(0) -
                         2 * xBasisDer(1, j) * yBasisDer(0, i) * sum(1) / (sum(0) * sum(0)) -
                         xBasisDer(0, j) * yBasisDer(0, i) * sum(2) / (sum(0) * sum(0)) +
                         2 * xBasisDer(0, j) * yBasisDer(0, i) * sum(1) * sum(1) / (sum(0) * sum(0) * sum(0))) *
                        weight;

                ders(count, 1) =
                        (xBasisDer(0, j) * yBasisDer(2, i) / sum(0) -
                         2 * xBasisDer(0, j) * yBasisDer(1, i) * sum(3) / (sum(0) * sum(0)) -
                         xBasisDer(0, j) * yBasisDer(0, i) * sum(4) / (sum(0) * sum(0)) +
                         2 * xBasisDer(0, j) * xBasisDer(0, i) * sum(3) * sum(3) / (sum(0) * sum(0) * sum(0))) *
                        weight;

                ders(count, 2) =
                        (xBasisDer(1, j) * yBasisDer(1, i) / sum(0) -
                         xBasisDer(1, j) * yBasisDer(0, i) * sum(3) / (sum(0) * sum(0)) -
                         xBasisDer(0, j) * yBasisDer(1, i) * sum(1) / (sum(0) * sum(0)) -
                         xBasisDer(1, j) * yBasisDer(0, i) * sum(5) * sum(5) / (sum(0) * sum(0)) +
                         2 * xBasisDer(0, j) * yBasisDer(0, i) * sum(1) * sum(3) / (sum(0) * sum(0) * sum(0))) *
                        weight;
                count++;
            }
        }
    }

    return ders;
}

///////////////////////////////////// BERNSTEIN ////////////////////////////////////////////////////////////////

Eigen::VectorXd Bernstein1DOrder1(double rParameter)
{
    return Bernstein1D(rParameter, 1);
}

Eigen::VectorXd Bernstein1DOrder2(double rParameter)
{
    return Bernstein1D(rParameter, 2);
}

Eigen::VectorXd Bernstein1DOrder3(double rParameter)
{
    return Bernstein1D(rParameter, 3);
}

Eigen::VectorXd Bernstein1DOrder4(double rParameter)
{
    return Bernstein1D(rParameter, 4);
}

Eigen::VectorXd Bernstein1D(double rParameter, int rOrder)
{
    Eigen::VectorXd values(rOrder);

    double paramUnivariate = (rParameter + 1.) / 2;
    double u1 = 1. - paramUnivariate;

    for (int j = 1; j < rOrder; j++)
    {
        double saved = 0.;
        for (int k = 0; k < j; k++)
        {
            double temp = values(k);
            values(k) = saved + u1 * temp;
            saved = paramUnivariate * temp;
        }
        values(j) = saved;
    }
    return values;
}


Eigen::VectorXd DerivativeBernstein1DOrder1(double rParameter)
{
    return DerivativeBernstein1D(rParameter, 1);
}

Eigen::VectorXd DerivativeBernstein1DOrder2(double rParameter)
{
    return DerivativeBernstein1D(rParameter, 2);
}

Eigen::VectorXd DerivativeBernstein1DOrder3(double rParameter)
{
    return DerivativeBernstein1D(rParameter, 3);
}

Eigen::VectorXd DerivativeBernstein1DOrder4(double rParameter)
{
    return DerivativeBernstein1D(rParameter, 4);
}

Eigen::VectorXd DerivativeBernstein1D(double, int rOrder)
{
    Eigen::VectorXd values(rOrder);
    double u1 = rOrder - 1;

    for (int j = 1; j < rOrder; j++)
    {
        double saved = 0.;
        for (int k = 0; k < j; k++)
        {
            double temp = values(k);
            values(k) = saved + u1 * temp;
            saved = u1 * temp;
        }
        values(j) = saved;
    }
    return values;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
} // namespace ShapeFunctionsIGA1D
}
