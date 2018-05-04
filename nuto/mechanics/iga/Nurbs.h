#pragma once

#include "nuto/mechanics/nodes/NodeSimple.h"
#include "nuto/base/Exception.h"
#include <array>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <vector>
#include "nuto/base/ValueVector.h"
#include "nuto/mechanics/DirectionEnum.h"

namespace NuTo
{
//! @author Peter Otto, BAM
//! @date September, 2017
//! @brief ... Class for NURBS curves, with IGA specific functions. NURBS specific algorithms taken from Piegl, Tiller
//! 'The NURBS book' 1996. TDimParameter is the dimension of the parametric space: curve is 1D and surface is 2D ...
template <int TDimParameter>
class Nurbs
{
public:
    enum eParametrization
    {
        chord,
        centripetal
    };

    /** Constructors **/

    //! @brief ... constructor
    //! @param degree ... degree of the polynomial
    //! @param knots ... knot vector
    //! @param controlPoints ... control points
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

    std::array<int, TDimParameter> VectorToTensor(int id)
    {
        std::array<int, TDimParameter> ind;

        // TODO.. remove
        std::array<int, TDimParameter + 1> numCPInDir;
        numCPInDir[0] = 1;
        for (int i = 1; i < TDimParameter + 1; i++)
            numCPInDir[i] = mNumControlPointsInDirection[i - 1];
        //

        ind[0] = (id / numCPInDir[0]) % numCPInDir[1];
        for (int i = 1; i < TDimParameter; i++)
            ind[i] = (id / (numCPInDir[i - 1] * numCPInDir[i])) % numCPInDir[i];

        return ind;
    }

    static void GetCPIndicesForRefinement(int accessIdLeft, int accessIdRight,
                                          const std::array<int, TDimParameter + 1>& numCPinDirRight,
                                          const std::array<int, TDimParameter + 1>& numCPinDirLeft, int dir,
                                          int numCPsRight,
                                          std::vector<int>& indicesLeft, // better std::vector<std::pair<int,int>>
                                          std::vector<int>& indicesRight)
    {
        indicesLeft.clear();
        indicesRight.clear();

        for (int id = 0; id < numCPsRight; id++)
        {
            std::array<int, TDimParameter> indexRight;
            indexRight[0] = (id / numCPinDirRight[0]) % numCPinDirRight[1];
            for (int i = 1; i < TDimParameter; i++)
                indexRight[i] = (id / (numCPinDirRight[i - 1] * numCPinDirRight[i])) % numCPinDirRight[i];

            if (indexRight[dir] == accessIdRight)
            {
                int indexRightVectorized = indexRight[0];
                for (int i = 1; i < TDimParameter; i++)
                    indexRightVectorized += indexRight[i] * numCPinDirRight[i];

                indicesRight.push_back(indexRightVectorized);

                indexRight[dir] = accessIdLeft;
                int indexLeftVectorized = indexRight[0];
                for (int i = 1; i < TDimParameter; i++)
                    indexLeftVectorized += indexRight[i] * numCPinDirLeft[i];

                indicesLeft.push_back(indexLeftVectorized);
            }
        }
    }

    //! @brief ... constructor
    //! @param knots ... knot vector
    //! @param controlPoints ... the coordinates of the control points (to be refined)
    //! @param weights ... weights to the control points
    //! @param degree ... degrees of the polynomial in each parameter direction
    //! @param refinementLevel ... the number of refinements for each direction (0 = no refinement)
    static void Refinement(const std::array<std::vector<double>, TDimParameter>& knots,
                           const std::vector<Eigen::VectorXd>& controlPoints, const std::vector<double>& weights,
                           const std::array<int, TDimParameter>& degree,
                           const std::array<std::vector<double>, TDimParameter>& knotsInserted,
                           std::vector<Eigen::VectorXd>& rControlPoints, std::vector<double>& rWeights,
                           std::array<std::vector<double>, TDimParameter>& rKnots)
    {
        assert(weights.size() == controlPoints.size());

        size_t numControlPoints = 1;
        for (size_t i = 0; i < knots.size(); i++)
            numControlPoints *= (knots[i].size() - degree[i] - 1);
        assert(numControlPoints == controlPoints.size());

        int dim = controlPoints[0].rows() + 1;

        // define actual knots, cps, weigts
        std::array<std::vector<double>, TDimParameter> actualKnots(knots);
        std::vector<Eigen::VectorXd> actualControlPoints(controlPoints);
        std::vector<double> actualWeights(weights);

        for (int dir = 0; dir < TDimParameter; dir++)
        {
            int numControlPointsDir = actualKnots[dir].size() - degree[dir] - 1;

            int numInsert = knotsInserted[dir].size();
            int begin = Nurbs<TDimParameter>::FindSpan(knotsInserted[dir][0], degree[dir], actualKnots[dir]);
            int end = Nurbs<TDimParameter>::FindSpan(knotsInserted[dir][numInsert - 1], degree[dir], actualKnots[dir]);
            end++;

            // new knot vector
            std::vector<double> newKnots(actualKnots[dir].size() + numInsert);
            for (int i = 0; i <= begin; i++)
                newKnots[i] = actualKnots[dir][i];
            for (int i = end + degree[0]; i < actualKnots[dir].size(); i++)
                newKnots[numInsert + i] = actualKnots[dir][i];

            // since we're dealing with NURBS and not B-Splines, a projection is needed ...
            std::vector<Eigen::VectorXd> controlPointsProjected;
            controlPointsProjected.resize(actualControlPoints.size());

            int count = 0;
            for (auto& it : actualControlPoints)
            {
                Eigen::VectorXd point(dim);
                point.block(0, 0, dim - 1, 1) = actualWeights[count] * it;
                point(dim - 1) = actualWeights[count];
                controlPointsProjected[count] = point;
                count++;
            }

            int numControlPointsNew = 1;
            for (int i = 0; i < TDimParameter; i++)
            {
                if (i == dir)
                    numControlPointsNew *= actualKnots[i].size() - degree[i] - 1 + numInsert;
                else
                    numControlPointsNew *= actualKnots[i].size() - degree[i] - 1;
            }

            std::vector<Eigen::VectorXd> newControlPoints;
            newControlPoints.resize(numControlPointsNew);

            std::array<int, TDimParameter + 1> numCPinDirOld;
            numCPinDirOld[0] = 1;
            for (int i = 1; i < TDimParameter + 1; i++)
                numCPinDirOld[i] = actualKnots[i - 1].size() - degree[i - 1] - 1;

            std::array<int, TDimParameter + 1> numCPinDirNew;
            numCPinDirNew[0] = 1;
            for (int i = 1; i < TDimParameter + 1; i++)
            {
                if (i - 1 == dir)
                    numCPinDirNew[i] = newKnots.size() - degree[i - 1] - 1;
                else
                    numCPinDirNew[i] = actualKnots[i - 1].size() - degree[i - 1] - 1;
            }

            std::vector<int> indicesLeft, indicesRight;

            for (int j = 0; j <= begin - degree[dir]; j++)
            {
                GetCPIndicesForRefinement(j, j, numCPinDirOld, numCPinDirNew, dir, actualControlPoints.size(),
                                          indicesLeft, indicesRight);
                for (int i = 0; i < indicesLeft.size(); i++)
                    newControlPoints[indicesLeft[i]] = controlPointsProjected[indicesRight[i]];
            }

            for (int j = end - 1; j < numControlPointsDir; j++)
            {
                GetCPIndicesForRefinement(j + numInsert, j, numCPinDirOld, numCPinDirNew, dir,
                                          actualControlPoints.size(), indicesLeft, indicesRight);
                for (int i = 0; i < indicesLeft.size(); i++)
                    newControlPoints[indicesLeft[i]] = controlPointsProjected[indicesRight[i]];
            }

            int i = end + degree[dir] - 1;
            int k = end + degree[dir] + numInsert - 1;

            for (int j = numInsert - 1; j >= 0; j--)
            {
                while (knotsInserted[dir][j] <= actualKnots[dir][i] && i > begin)
                {
                    newKnots[k] = actualKnots[dir][i];

                    GetCPIndicesForRefinement(k - degree[dir] - 1, i - degree[dir] - 1, numCPinDirOld, numCPinDirNew,
                                              dir, actualControlPoints.size(), indicesLeft, indicesRight);
                    for (int i = 0; i < indicesLeft.size(); i++)
                        newControlPoints[indicesLeft[i]] = controlPointsProjected[indicesRight[i]];

                    k--;
                    i--;
                }

                GetCPIndicesForRefinement(k - degree[dir] - 1, k - degree[dir], numCPinDirNew, numCPinDirNew, dir,
                                          newControlPoints.size(), indicesLeft, indicesRight);
                for (int i = 0; i < indicesLeft.size(); i++)
                    newControlPoints[indicesLeft[i]] = newControlPoints[indicesRight[i]];

                for (int l = 1; l <= degree[dir]; l++)
                {
                    int ind = k - degree[dir] + l;
                    double alpha = newKnots[k + l] - knotsInserted[dir][j];
                    if (std::fabs(alpha) == 0.0)
                    {
                        GetCPIndicesForRefinement(ind - 1, ind, numCPinDirNew, numCPinDirNew, dir,
                                                  newControlPoints.size(), indicesLeft, indicesRight);
                        for (int i = 0; i < indicesLeft.size(); i++)
                            newControlPoints[indicesLeft[i]] = newControlPoints[indicesRight[i]];
                    }
                    else
                    {
                        alpha /= (newKnots[k + l] - actualKnots[dir][i - degree[dir] + l]);

                        GetCPIndicesForRefinement(ind - 1, ind, numCPinDirNew, numCPinDirNew, dir,
                                                  newControlPoints.size(), indicesLeft, indicesRight);
                        for (int i = 0; i < indicesLeft.size(); i++)
                            newControlPoints[indicesLeft[i]] = alpha * newControlPoints[indicesLeft[i]] +
                                                               (1. - alpha) * newControlPoints[indicesRight[i]];
                    }
                }
                newKnots[k] = knotsInserted[dir][j];
                k--;
            }

            actualKnots[dir] = newKnots;
            actualWeights.clear();
            actualControlPoints.clear();

            for (Eigen::VectorXd& it : newControlPoints)
            {
                Eigen::VectorXd point(dim - 1);
                actualWeights.push_back(it(dim - 1));
                point = (1. / it(dim - 1)) * it.block(0, 0, dim - 1, 1);
                actualControlPoints.push_back(point);
            }
        }

        rControlPoints = actualControlPoints;
        rWeights = actualWeights;
        rKnots = actualKnots;
    }

    /** Getter **/

    //! @brief ... get the dimension of the curve = dimension of each control point
    //! @return ... dimension of the curve
    int GetDimension() const
    {
        return mControlPoints[0]->GetNumValues();
    }

    //! @brief ... get the number of control points in total
    int GetNumControlPoints() const
    {
        return mControlPoints.size();
    }

    //! @brief ... get the number of control points for each IGA element (one parametric span) in a specific
    //! direction
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

    //! @brief ... get the number of elements = # of non-vanishing knot spans
    //! @return ... number of isogeometric elements in this nurbs
    int GetNumElementsTotal()
    {
        int numElements = 1;
        for (auto& knots : mKnots)
        {
            int numElementsDir = 0;
            for (int i = 1; i < knots.size(); i++)
            {
                if (knots[i - 1] < knots[i])
                    numElementsDir++;
            }
            numElements *= numElementsDir;
        }
        return numElements;
    }

    //! @brief ... get the number of elements = # of non-vanishing knot spans
    //! @return ... number of isogeometric elements in this nurbs
    int GetNumElementsInDirection(int dir)
    {
        assert(dir < TDimParameter);

        std::vector<double> knots = mKnots[dir];

        int numElementsDir = 0;
        for (int i = 1; i < knots.size(); i++)
        {
            if (knots[i - 1] < knots[i])
                numElementsDir++;
        }

        return numElementsDir;
    }

    //! @brief  ... get the knot ids of the element with the local id
    //! @return ... knot ids of the element with local id
    std::array<int, TDimParameter> GetElement(const std::array<int, TDimParameter>& elementId)
    {
        std::array<int, TDimParameter> knotIDs;

        int count = 0;
        for (auto& knots : mKnots)
        {
            int countElementsDir = 0;
            for (int i = 1; i < knots.size(); i++)
            {
                if (knots[i - 1] < knots[i])
                    countElementsDir++;

                if (countElementsDir == elementId[count])
                {
                    knotIDs[count] = countElementsDir;
                    count++;
                }
            }
        }

        assert(count == TDimParameter);

        return knotIDs;
    }

    //! @brief  ... get the knot ids of the element with the local id
    //! @return ... knot ids of the element with local id
    std::array<int, TDimParameter> GetElement(int id)
    {
        std::array<int, TDimParameter> knotIDs;

        std::array<int, TDimParameter> index;
        int power = 1;
        for (int dim = 0; dim < TDimParameter; dim++)
        {
            int numElementsInOneDir = GetNumElementsInDirection(dim);
            index[dim] = (id / power) % numElementsInOneDir;
            power *= numElementsInOneDir;
        }

        int count = 0, countCheck = 0;
        for (auto& knots : mKnots)
        {
            int countElementsDir = -1;
            for (int i = 1; i < knots.size(); i++)
            {
                if (knots[i - 1] < knots[i])
                    countElementsDir++;

                if (countElementsDir == index[count])
                {
                    knotIDs[count] = i - 1;
                    countCheck++;
                    break;
                }
            }
            count++;
        }

        assert(TDimParameter == countCheck);

        return knotIDs;
    }

    const NodeSimple* GetControlPointElement(const std::array<int, TDimParameter>& knotID, int id) const
    {
        assert(id < GetNumControlPointsElement());
        int index = GetElementControlPointIndex(knotID, id);
        return mControlPoints[index];
    }

    NodeSimple* GetControlPointElement(const std::array<int, TDimParameter>& knotID, int id)
    {
        assert(id < GetNumControlPointsElement());
        int index = GetElementControlPointIndex(knotID, id);
        return mControlPoints[index];
    }

    //! @brief ... get the knot vector in the direction
    //! @return ... knots
    std::vector<double> GetKnotVectorDirection(eDirection dir) const
    {
        const int directionComponent = ToComponentIndex(dir);
        assert(directionComponent < TDimParameter);
        return mKnots[directionComponent];
    }

    int GetDegreeDirection(eDirection dir) const
    {
        const int directionComponent = ToComponentIndex(dir);
        assert(directionComponent < TDimParameter);
        return mDegree[directionComponent];
    }

    void GetControlPointsAnsWeightsBoundary(eDirection dirNotIncluded, int locationBoundary,
                                            std::vector<NodeSimple*>& rControlPoints, std::vector<double>& rWeights)
    {
        rControlPoints.clear();
        rWeights.clear();

        const int dirNotIncludedComponent = ToComponentIndex(dirNotIncluded);
        assert(dirNotIncludedComponent < TDimParameter);

        int idTensorCheck = -1;
        if (locationBoundary == 0)
            idTensorCheck = 0;
        else
            idTensorCheck = mNumControlPointsInDirection[dirNotIncludedComponent] - 1;

        for (int id = 0; id < mControlPoints.size(); id++)
        {
            std::array<int, TDimParameter> idTensor = VectorToTensor(id);
            if (idTensor[dirNotIncludedComponent] == idTensorCheck)
            {
                rControlPoints.push_back(mControlPoints[id]);
                rWeights.push_back(mWeights[id]);
            }
        }
    }


    //! @brief ... get the knots to given knot ids
    //! @return ... knots
    std::array<Eigen::Vector2d, TDimParameter> GetKnotVectorElement(const std::array<int, TDimParameter>& knotIDs) const
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

    std::vector<std::array<double, 2>> GetKnots(const std::array<int, TDimParameter>& knotID) const
    {
        std::vector<std::array<double, 2>> ret;
        for (int i = 0; i < TDimParameter; i++)
        {
            assert(knotID[i] < mKnots[i].size() - mDegree[i] - 1);
            std::array<double, 2> knots;
            knots[0] = mKnots[i][knotID[i]];
            knots[1] = mKnots[i][knotID[i] + 1];
            ret.push_back(knots);
        }

        return ret;
    }

    Eigen::Matrix<double, TDimParameter, 1> Transformation(Eigen::VectorXd ipCoords,
                                                           const std::array<int, TDimParameter>& knotIDs) const
    {
        std::array<Eigen::Vector2d, TDimParameter> knots = GetKnotVectorElement(knotIDs);
        Eigen::Matrix<double, TDimParameter, 1> coordinateTransformed;
        for (int i = 0; i < TDimParameter; i++)
        {
            coordinateTransformed[i] = (knots[i](0) + 0.5 * (ipCoords(i) + 1) * (knots[i](1) - knots[i](0)));
        }
        return coordinateTransformed;
    }

    std::array<int, TDimParameter> GetElementControlPointStartIndex(const std::array<int, TDimParameter>& knotID) const
    {
        for (int i = 0; i < TDimParameter; i++)
            assert(knotID[i] >= mDegree[i] && knotID[i] < mKnots[i].size());

        std::array<int, TDimParameter> indexBegin;

        // tensor index
        for (int i = 0; i < TDimParameter; i++)
            indexBegin[i] = knotID[i] - mDegree[i];

        return indexBegin;
    }

    int GetElementControlPointIndex(const std::array<int, TDimParameter>& knotID, int id) const
    {
        assert(id < GetNumControlPointsElement());

        std::array<int, TDimParameter> indexBegin = GetElementControlPointStartIndex(knotID);

        std::array<int, TDimParameter> indexAdd;
        int power = 1;
        for (int dim = 0; dim < TDimParameter; dim++)
        {
            int numCPInElementDir = mDegree[dim] + 1;
            indexAdd[dim] = (id / power) % numCPInElementDir;
            power *= numCPInElementDir;
        }

        int index = indexBegin[0] + indexAdd[0];

        for (int i = 1; i < TDimParameter; i++)
        {
            int indexTemp = indexBegin[i] + indexAdd[i];
            for (int j = 0; j < i; j++)
                indexTemp *= mNumControlPointsInDirection[j];

            index += indexTemp;
        }

        return index;
    }


    Eigen::VectorXd GetControlPointCoordinatesElement(const std::array<int, TDimParameter>& knotID,
                                                      int instance = 0) const;

    NodeSimple* GetControlPoint(const Eigen::VectorXi& ids)
    {
        int index = ids[0];

        for (int i = 1; i < TDimParameter; i++)
        {
            int indexTemp = ids[i];
            for (int j = 0; j < i; j++)
                indexTemp *= mNumControlPointsInDirection[j];

            index += indexTemp;
        }

        return mControlPoints[index];
    }

    NodeSimple* GetControlPoint(int id)
    {
        assert(id < mControlPoints.size());
        return mControlPoints[id];
    }

    Eigen::VectorXd GetControlPointCoordinates(int id, int instance = 0)
    {
        assert(id < mControlPoints.size());
        return mControlPoints[id]->GetValues(instance);
    }

    /** Evaluation **/

    Eigen::VectorXd Evaluate(const Eigen::Matrix<double, TDimParameter, 1>& parameter, int derivativeOrder = 0) const
    {
        const std::array<int, TDimParameter> spanIdx = FindSpan(parameter);
        Eigen::MatrixXd shapeFunctions = BasisFunctionsAndDerivativesRational(derivativeOrder, parameter, spanIdx);

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
    //! @return Eigen::MatrixXd ... matrix containing the shape functions up to derivative der, the row index
    //! represents
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

    Eigen::MatrixXd BasisFunctionsAndDerivativesRational(int der,
                                                         const Eigen::Matrix<double, TDimParameter, 1>& parameter,
                                                         const std::array<int, TDimParameter>& knotIDs) const;

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
    std::array<int, 1> indexBegin = GetElementControlPointStartIndex(knotID);

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

    std::array<int, 2> indexBegin = GetElementControlPointStartIndex(knotID);

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
inline Eigen::MatrixXd Nurbs<1>::BasisFunctionsAndDerivativesRational(int der,
                                                                      const Eigen::Matrix<double, 1, 1>& parameter,
                                                                      const std::array<int, 1>& knotIDs) const
{
    assert(der >= 0 && der <= 2);

    // int spanIdx = FindSpan(parameter)[0];
    int spanIdx = knotIDs[0];
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
inline Eigen::MatrixXd Nurbs<2>::BasisFunctionsAndDerivativesRational(int der,
                                                                      const Eigen::Matrix<double, 2, 1>& parameter,
                                                                      const std::array<int, 2>& spanIdx) const
{
    assert(der >= 0 && der <= 2);

    //    const std::array<int, 2> spanIdx = FindSpan(parameter);
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

} // namespace NuTo
