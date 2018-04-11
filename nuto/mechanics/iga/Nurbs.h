#pragma once

#include "nuto/mechanics/nodes/NodeSimple.h"
#include "nuto/base/Exception.h"
#include <array>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <vector>
#include "nuto/base/ValueVector.h"

namespace NuTo
{
//! @author Peter Otto, BAM
//! @date September, 2017
//! @brief ... class for NURBS curves, with IGA specific functions
//! @brief ... NURBS specific algorithms taken from Piegl, Tiller 'The NURBS book' 1996
//! @brief ... TDimParameter is the dimension of the parametric space: curve is 1D and surface is 2D
template <int TDimParameter>
class Nurbs
{
public:
    enum mParametrization
    {
        chord,
        centripetal
    };

    /** Constructors **/

    //! @brief ... constructor
    //! @param rDegree ... degree of the polynomial
    //! @param rKnots ... knot vector
    //! @param rControlPoints ... control points
    Nurbs(const std::array<std::vector<double>, TDimParameter>& knots, const std::vector<NodeSimple*>& controlPoints,
          const std::vector<double>& weights, const std::array<int, TDimParameter>& degree)
        : mKnots(knots)
        , mControlPoints(controlPoints)
        , mWeights(weights)
        , mDegree(degree)
    {
        assert(degree.size() == knots.size());
        assert(controlPoints.size() == weights.size());

        size_t numControlPoints = 1;
        for (size_t i = 0; i < mKnots.size(); i++)
        {
            int numControlPointsDir = knots[i].size() - degree[i] - 1;
            mNumControlPointsInDirection[i] = numControlPointsDir;
            numControlPoints *= numControlPointsDir;
        }

        assert(numControlPoints == controlPoints.size());
    }

    //! @brief ... constructor
    //! @param knots ... knot vector
    //! @param controlPoints ... the coordinates of the control points (to be refined)
    //! @param weights ... weights to the control points
    //! @param degree ... degrees of the polynomial in each parameter direction
    //! @param refinementLevel ... the number of refinements for each direction (0 = no refinement)
    //! @param mesh ... after the refinement the control points are stored in the mesh and the pointers in this
    Nurbs(const std::array<std::vector<double>, TDimParameter>& knots,
          const std::vector<Eigen::VectorXd>& controlPoints, const std::vector<double>& weights,
          const std::array<int, TDimParameter>& degree, const std::array<int, TDimParameter>& refinementLevel,
          std::vector<NodeSimple>& nodes)
    {
        assert(degree.size() == knots.size());
        assert(controlPoints.size() == weights.size());

        size_t numControlPoints = 1;
        for (size_t i = 0; i < mKnots.size(); i++)
        {
            int numControlPointsDir = knots[i].size() - degree[i] - 1;
            mNumControlPointsInDirection[i] = numControlPointsDir;
            numControlPoints *= numControlPointsDir;
        }

        assert(numControlPoints == controlPoints.size());

        // Refine - call routine ot this class

        // save the control points
    }

    /** Getter **/

    //! @brief ... get the dimension of the curve = dimension of each control point
    //! @return ... dimension of the curve
    int GetDimension() const
    {
        return mControlPoints[0]->GetNumValues();
    }

    //! @brief ... get the number of control points for each IGA element (one parametric span) in a specific direction
    //! @return ... degree + 1
    int GetNumControlPointsElement(int dir) const
    {
        assert(dir <= TDimParameter && dir >= 0);
        return mDegree[dir] + 1;
    }

    //! @brief ... get the number of control points for each IGA element (one parametric span)
    //! @return ... product of the number of control points in each direction
    int GetNumControlPointsElement() const
    {
        int numCPs = 1;
        for (int degree : mDegree)
            numCPs *= degree + 1;

        return numCPs;
    }

    //! @brief ... get the knots to given knot ids
    //! @return ... knots
    std::array<Eigen::Vector2d, TDimParameter> GetKnotVectorElement(std::array<int, TDimParameter> knotIDs) const
    {
        std::array<Eigen::Vector2d, TDimParameter> knots;
        Eigen::Vector2d knotCoordinates;
        for (int i = 0; i < TDimParameter; i++)
        {
            knotCoordinates(0) = mKnots[i][knotIDs[i]];
            knotCoordinates(1) = mKnots[i][knotIDs[i] + 1];
            knots[i] = knotCoordinates;
        }
        return knots;
    }

    Eigen::VectorXi GetControlPointCoordinatesElementHelper(const std::array<int, TDimParameter>& knotID) const
    {
        for (int i = 0; i < TDimParameter; i++)
            assert(knotID[i] >= mDegree[i] && knotID[i] < mKnots[i].size());

        Eigen::VectorXi indexBegin(TDimParameter);

        // tensor index
        for (int i = 0; i < TDimParameter; i++)
            indexBegin[i] = knotID[i] - mDegree[i];

        return indexBegin;
    }

    Eigen::VectorXd GetControlPointCoordinatesElement(const std::array<int, TDimParameter>& knotID,
                                                      int instance = 0) const;

    const NodeSimple* GetControlPointElement(const std::array<int, TDimParameter>& knotID, int i) const;

    NodeSimple* GetControlPoint(const Eigen::VectorXi& ids)
    {
        int index = ids[0];

        for (int i = 1; i < TDimParameter; i++)
            for (int j = 0; j < i; j++)
                index += mNumControlPointsInDirection[j] * ids[i];

        return mControlPoints[index];
    }

    /** Evaluation **/

    Eigen::VectorXd Evaluate(const Eigen::Matrix<double, TDimParameter, 1>& parameter, int derivativeOrder = 0) const
    {
        const std::array<int, TDimParameter> spanIdx = FindSpan(parameter);
        Eigen::MatrixXd shapeFunctions = BasisFunctionsAndDerivativesRational(derivativeOrder, parameter);

        Eigen::VectorXd coordinates(GetDimension());
        coordinates.setZero(GetDimension());

        Eigen::VectorXd cpCoords = GetControlPointCoordinatesElement(spanIdx);

        int numCPs = GetNumControlPointsElement();

        assert(numCPs == shapeFunctions.rows());

        int dim = GetDimension();

        for (int i = 0; i < dim; i++)
            for (int j = 0; j < numCPs; j++)
                coordinates(i) += shapeFunctions(j, 0) * cpCoords(dim * j + i);

        return coordinates;
    }

    /** Basis Functions **/

    static int FindSpan(double parameter, int degree, const std::vector<double>& knots)
    {
        int size = knots.size();
        int iterations = 0;

        if (parameter < knots[0] || parameter > knots[size - 1])
            throw NuTo::Exception(__PRETTY_FUNCTION__, "The parameter is out of the range of the knot vector.");

        int numBasisFuns = size - degree - 1;
        if (parameter == knots[numBasisFuns])
            return numBasisFuns - 1;

        int low = degree;
        int high = numBasisFuns;
        int mid = (low + high) / 2;

        while (parameter < knots[mid] || parameter >= knots[mid + 1])
        {
            if (parameter < knots[mid])
                high = mid;
            else
                low = mid;

            mid = (low + high) / 2;
            iterations++;
            if (iterations > size)
                throw NuTo::Exception(__PRETTY_FUNCTION__,
                                      "The maximum number of iterations for finding the span exceeded.");
        }
        return mid;
    }

    //! @brief ... calculates the basisfunctions and derivatives, see Piegl/Tiller 'NURBS Book' 2nd ed., Page 72
    //! @param der ... up to derivative der
    //! @param parameter ... parameter at which the functions are calculated
    //! @param spanIdx ... spanIdx to the parameter
    //! @param degree ... degree of the given NURBS curve
    //! @param knots ... knot vector
    //! @return Eigen::MatrixXd ... matrix containing the shape functions up to derivative der, the row index represents
    //! the derivative order
    static Eigen::MatrixXd BasisFunctionsAndDerivatives(int derivativeOrder, double parameter, int spanIdx, int degree,
                                                        const std::vector<double>& knots)
    {
        // store shape functions and knot differences
        Eigen::MatrixXd ndu(degree + 1, degree + 1);

        ndu(0, 0) = 1.0;

        Eigen::VectorXd left(degree + 1);
        Eigen::VectorXd right(degree + 1);

        for (int j = 1; j <= degree; j++)
        {
            left[j] = parameter - knots[spanIdx + 1 - j];
            right[j] = knots[spanIdx + j] - parameter;
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

        Eigen::MatrixXd basisFctsDerivatives(derivativeOrder + 1, degree + 1);

        for (int j = 0; j <= degree; j++)
            basisFctsDerivatives(0, j) = ndu(j, degree);

        // the shape functions
        Eigen::MatrixXd a(2, degree + 1);
        for (int r = 0; r <= degree; r++)
        {
            int s1 = 0;
            int s2 = 1;

            a(0, 0) = 1.0;

            int j = 0;
            for (int k = 1; k <= derivativeOrder; k++)
            {
                double d = 0.;
                int rk = r - k;
                int pk = degree - k;
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
                    j2 = degree - r;

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

                basisFctsDerivatives(k, r) = d;
                j = s1;
                s1 = s2;
                s2 = j;
            }
        }

        int r = degree;
        for (int k = 1; k <= derivativeOrder; k++)
        {
            for (int j = 0; j <= degree; j++)
                basisFctsDerivatives(k, j) *= r;
            r *= (degree - k);
        }

        return basisFctsDerivatives;
    }

    //! @brief ... calculates the knot span to given parameters
    //! @param parameter ... parameter(s) at which the functions are calculated
    //! @return std::array ... knot span(s) to given parameter(s)
    const std::array<int, TDimParameter> FindSpan(const Eigen::Matrix<double, TDimParameter, 1>& parameter) const
    {
        std::array<int, TDimParameter> parameterIDs;
        for (int i = 0; i < parameter.rows(); i++)
        {
            parameterIDs[i] = FindSpan(parameter[i], mDegree[i], mKnots[i]);
        }
        return parameterIDs;
    }

    Eigen::MatrixXd
    BasisFunctionsAndDerivativesRational(int der, const Eigen::Matrix<double, TDimParameter, 1>& parameter) const;

    /** Change the discretization (e.g. refinement) **/

    void InsertKnot(const Eigen::Matrix<double, TDimParameter, 1>& knotToInsert, int rMultiplicity)
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Iga - Not implemented yet!");
    }

    void RefineKnots(const std::vector<Eigen::Matrix<double, TDimParameter, 1>>& knotsToInsert)
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Iga - Not implemented yet!");
    }

    void DuplicateKnots()
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Iga - Not implemented yet!");
    }

    /** Projection (minimumDistance) **/

    Eigen::Matrix<double, TDimParameter, 1> minimumDistance(const Eigen::VectorXd& coordinatesSlave) const
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Iga - Not implemented yet!");
    }

private:
    //! @brief Knot vector (in isogeometric framework each segment between two
    //! knots is an element)
    std::array<std::vector<double>, TDimParameter> mKnots;

    //! @brief Control points of the BSpline curve (# rows = num control points)
    std::vector<NodeSimple*> mControlPoints;

    //! @brief Weights to NURBS
    std::vector<double> mWeights;

    //! @brief Degree of the polynomials (order = mDegree+1)
    std::array<int, TDimParameter> mDegree;

    std::array<int, TDimParameter> mNumControlPointsInDirection;
};

template <>
inline Eigen::VectorXd Nurbs<1>::GetControlPointCoordinatesElement(const std::array<int, 1>& knotID, int instance) const
{
    Eigen::VectorXi indexBegin = GetControlPointCoordinatesElementHelper(knotID);

    int dim = GetDimension();
    int numCPs = GetNumControlPointsElement();

    Eigen::VectorXd nodeValues(numCPs * dim);

    int count = 0;
    for (int i = indexBegin[0]; i <= indexBegin[0] + mDegree[0]; i++)
    {
        nodeValues.segment(count, dim) = mControlPoints[i]->GetValues(instance);
        count += dim;
    }

    return nodeValues;
}

template <>
inline Eigen::VectorXd Nurbs<2>::GetControlPointCoordinatesElement(const std::array<int, 2>& knotID, int instance) const
{
    int dim = GetDimension();
    int numCPs = GetNumControlPointsElement();

    Eigen::VectorXd nodeValues(numCPs * dim);

    Eigen::VectorXi indexBegin = GetControlPointCoordinatesElementHelper(knotID);

    int count = 0;
    for (int j = indexBegin[1]; j <= indexBegin[1] + mDegree[1]; j++)
        for (int i = indexBegin[0]; i <= indexBegin[0] + mDegree[0]; i++)
        {
            int index = i + j * mNumControlPointsInDirection[0];
            nodeValues.segment(count, dim) = mControlPoints[index]->GetValues(instance);
            count += dim;
        }

    return nodeValues;
}

template <>
inline const NodeSimple* Nurbs<1>::GetControlPointElement(const std::array<int, 1>& knotID, int i) const
{
    assert(i >= GetNumControlPointsElement());
    assert(knotID[0] >= mDegree[0]);

    return mControlPoints[knotID[0] - mDegree[0] + i];
}

template <>
inline const NodeSimple* Nurbs<2>::GetControlPointElement(const std::array<int, 2>& knotID, int i) const
{
    assert(i < GetNumControlPointsElement());
    assert(knotID[0] >= mDegree[0] && knotID[1] >= mDegree[1]);

    return mControlPoints[i];
}


template <>
inline Eigen::MatrixXd
Nurbs<1>::BasisFunctionsAndDerivativesRational(int der, const Eigen::Matrix<double, 1, 1>& parameter) const
{
    assert(der >= 0 && der <= 2);

    int spanIdx = FindSpan(parameter)[0];
    Eigen::MatrixXd ders = BasisFunctionsAndDerivatives(der, parameter[0], spanIdx, mDegree[0], mKnots[0]);

    // NURBS specific ...
    Eigen::VectorXd sum = Eigen::VectorXd::Zero(der + 1);

    for (int i = 0; i <= mDegree[0]; i++)
    {
        double weight = mWeights[spanIdx - mDegree[0] + i];
        if (der == 0)
            sum(0) += ders(0, i) * weight;
        else if (der == 1)
        {
            sum(0) += ders(0, i) * weight;
            sum(1) += ders(1, i) * weight;
        }
        else
        {
            sum(0) += ders(0, i) * weight;
            sum(1) += ders(1, i) * weight;
            sum(2) += ders(2, i) * weight;
        }
    }

    Eigen::VectorXd dersRat(mDegree[0] + 1);
    dersRat.setZero(mDegree[0] + 1);

    for (int i = 0; i <= mDegree[0]; i++)
    {
        double weight = mWeights[spanIdx - mDegree[0] + i];
        if (der == 0)
            dersRat(i) = ders(0, i) * weight / sum(0);
        else if (der == 1)
            dersRat(i) = (ders(1, i) * sum(0) - ders(0, i) * sum(1)) * weight / (sum(0) * sum(0));
        else
        {
            double sum2 = sum(0) * sum(0);
            dersRat(i) = weight * (ders(2, i) / sum(0) - 2 * ders(1, i) * sum(1) / (sum2)-ders(0, i) * sum(2) / (sum2) +
                                   2 * ders(0, i) * sum(1) * sum(1) / (sum2 * sum(0)));
        }
    }
    return dersRat;
}

template <>
inline Eigen::MatrixXd
Nurbs<2>::BasisFunctionsAndDerivativesRational(int der, const Eigen::Matrix<double, 2, 1>& parameter) const
{
    assert(der >= 0 && der <= 2);

    const std::array<int, 2> spanIdx = FindSpan(parameter);
    std::array<Eigen::MatrixXd, 2> shapeFunctions;

    int numDers = 1;
    for (int i = 0; i < 2; i++)
    {
        shapeFunctions[i] = BasisFunctionsAndDerivatives(der, parameter[i], spanIdx[i], mDegree[i], mKnots[i]);
        numDers *= mDegree[i] + 1;
    }

    Eigen::MatrixXd ders(numDers, der + 1);
    ders.setZero(numDers, der + 1);

    int num = ((der + 1) * (der + 1) + (der + 1)) / 2;
    Eigen::VectorXd sum(num);
    sum.setZero(num);

    for (int i = 0; i <= mDegree[1]; i++)
    {
        for (int j = 0; j <= mDegree[0]; j++)
        {
            int index = (spanIdx[0] - mDegree[0] + j) + (spanIdx[1] - mDegree[1] + i) * mNumControlPointsInDirection[0];
            double weight = mWeights[index];

            if (der == 0)
            {
                sum(0) += shapeFunctions[0](0, j) * shapeFunctions[1](0, i) * weight;
            }
            else if (der == 1)
            {
                sum(0) += shapeFunctions[0](0, j) * shapeFunctions[1](0, i) * weight;
                sum(1) += shapeFunctions[0](1, j) * shapeFunctions[1](0, i) * weight;
                sum(2) += shapeFunctions[0](0, j) * shapeFunctions[1](1, i) * weight;
            }
            else
            {
                sum(0) += shapeFunctions[0](0, j) * shapeFunctions[1](0, i) * weight;
                sum(1) += shapeFunctions[0](1, j) * shapeFunctions[1](0, i) * weight;
                sum(2) += shapeFunctions[0](2, j) * shapeFunctions[1](0, i) * weight;

                sum(3) += shapeFunctions[0](0, j) * shapeFunctions[1](1, i) * weight;
                sum(4) += shapeFunctions[0](0, j) * shapeFunctions[1](2, i) * weight;

                sum(5) += shapeFunctions[0](1, j) * shapeFunctions[1](1, i) * weight;
            }
        }
    }

    int count = 0;
    for (int i = 0; i <= mDegree[1]; i++)
    {
        for (int j = 0; j <= mDegree[0]; j++)
        {
            int index = (spanIdx[0] - mDegree[0] + j) + (spanIdx[1] - mDegree[1] + i) * mNumControlPointsInDirection[0];
            double weight = mWeights[index];

            if (der == 0)
            {
                ders(count, 0) = shapeFunctions[0](0, j) * shapeFunctions[1](0, i) * weight / sum(0);
                count++;
            }
            else if (der == 1)
            {
                ders(count, 0) = (shapeFunctions[0](1, j) * shapeFunctions[1](0, i) * sum(0) -
                                  shapeFunctions[0](0, j) * shapeFunctions[1](0, i) * sum(1)) *
                                 weight / (sum(0) * sum(0));
                ders(count, 1) = (shapeFunctions[0](0, j) * shapeFunctions[1](1, i) * sum(0) -
                                  shapeFunctions[0](0, j) * shapeFunctions[1](0, i) * sum(2)) *
                                 weight / (sum(0) * sum(0));
                count++;
            }
            else
            {
                ders(count, 0) = (shapeFunctions[0](2, j) * shapeFunctions[1](0, i) / sum(0) -
                                  2 * shapeFunctions[0](1, j) * shapeFunctions[1](0, i) * sum(1) / (sum(0) * sum(0)) -
                                  shapeFunctions[0](0, j) * shapeFunctions[1](0, i) * sum(2) / (sum(0) * sum(0)) +
                                  2 * shapeFunctions[0](0, j) * shapeFunctions[1](0, i) * sum(1) * sum(1) /
                                          (sum(0) * sum(0) * sum(0))) *
                                 weight;

                ders(count, 1) = (shapeFunctions[0](0, j) * shapeFunctions[1](2, i) / sum(0) -
                                  2 * shapeFunctions[0](0, j) * shapeFunctions[1](1, i) * sum(3) / (sum(0) * sum(0)) -
                                  shapeFunctions[0](0, j) * shapeFunctions[1](0, i) * sum(4) / (sum(0) * sum(0)) +
                                  2 * shapeFunctions[0](0, j) * shapeFunctions[0](0, i) * sum(3) * sum(3) /
                                          (sum(0) * sum(0) * sum(0))) *
                                 weight;

                ders(count, 2) =
                        (shapeFunctions[0](1, j) * shapeFunctions[1](1, i) / sum(0) -
                         shapeFunctions[0](1, j) * shapeFunctions[1](0, i) * sum(3) / (sum(0) * sum(0)) -
                         shapeFunctions[0](0, j) * shapeFunctions[1](1, i) * sum(1) / (sum(0) * sum(0)) -
                         shapeFunctions[0](1, j) * shapeFunctions[1](0, i) * sum(5) * sum(5) / (sum(0) * sum(0)) +
                         2 * shapeFunctions[0](0, j) * shapeFunctions[1](0, i) * sum(1) * sum(3) /
                                 (sum(0) * sum(0) * sum(0))) *
                        weight;
                count++;
            }
        }
    }
    return ders;
}
}
