#pragma once

#include<eigen3/Eigen/Dense>
#include <array>
#include <vector>

#include "mechanics/nodes/NodeSimple.h"

namespace NuTo
{
template <int TDimParameter>
class NURBS
{
public:

    enum mParametrization {chord, centripetal};

    /** Constructors **/

    //! @brief ... default constructor
    NURBS(){}

    //! @brief ... constructor
    //! @param rDegree ... degree of the polynomial
    //! @param rKnots ... knot vector
    //! @param rControlPoints ... control points
    NURBS(const std::array<std::vector<double>, TDimParameter> &rKnots,
          const std::vector<std::vector<NodeSimple>>           &rControlPoints,
          const std::vector<std::vector<double>>               &rWeights,
          const std::array<int,TDimParameter>                  &rDegree)
        : mKnots(rKnots), mControlPoints(rControlPoints), mWeights(rWeights), mDegree(rDegree)
    {}

    /** Getter **/

    //! @brief ... get the dimension of the curve
    //! @return ... dimension of the curve
    int GetDimension() const
    {
        return mControlPoints[0][0].GetNumValues();
    }

    int GetNumControlPointsElement(int dir) const
    {
        assert(dir <= TDimParameter && dir >= 0);
        return mDegree[dir]+1;
    }

    int GetNumControlPointsElement() const
    {
        int numCPs = 1;
        for(int i = 0; i < TDimParameter; i++)
            numCPs *= mDegree[i]+1;

        return numCPs;
    }

    std::array<Eigen::Vector2d, TDimParameter> GetKnotVectorElement(std::array<int, TDimParameter> knotIDs) const
    {
        std::array<Eigen::Vector2d, 2> knots;
        Eigen::Vector2d knotCoordinates;
        for(int i = 0; i < TDimParameter; i++)
        {
            knotCoordinates(0) = mKnots[i][knotIDs[i]];
            knotCoordinates(1) = mKnots[i][knotIDs[i]+1];
            knots[i] = knotCoordinates;
        }
        return knots;
    }


    Eigen::VectorXd GetControlPointsElement(std::array<int, TDimParameter> knotID) const
    {
        throw new std::string("Not implemented for an arbitrary dimension!");
    }

    /** Evaluation **/

    Eigen::VectorXd Evaluate(const Eigen::Matrix<double, TDimParameter, 1> &parameter) const
    {
        const std::array<int, TDimParameter> spanIdx = FindSpan(parameter);
        Eigen::MatrixXd shapeFunctions = BasisFunctionsAndDerivativesRational(0, parameter);

        Eigen::VectorXd coordinates(GetDimension());
        coordinates.setZero(GetDimension());

        Eigen::VectorXd cpCoords = GetControlPointsElement(spanIdx);

        int numCPs = GetNumControlPointsElement();

        assert(numCPs == shapeFunctions.rows());

        int dim = GetDimension();

        for(int i = 0; i < dim; i++)
            for(int j = 0; j < numCPs; j++)
                coordinates(i) += shapeFunctions(j,0) * cpCoords(dim*j+i);

        return coordinates;
    }

    /** Basis Functions **/

    static int FindSpan(double rParameter, int rDegree, const std::vector<double> &rKnots)
    {
        int size = rKnots.size();
        int iterations = 0;

        if(rParameter < rKnots[0] || rParameter > rKnots[size - 1])
            throw new std::string("The parameter is out of the range of the knot vector.");

        int numBasisFuns = size - rDegree - 1;
        if(rParameter == rKnots[numBasisFuns]) return numBasisFuns-1;

        int low  = rDegree;
        int high = numBasisFuns;
        int mid  = (low + high)/2;

        while(rParameter < rKnots[mid] || rParameter >= rKnots[mid+1])
        {
            if(rParameter < rKnots[mid]) high = mid;
            else                         low = mid;

            mid = (low + high)/2;
            iterations++;
            if(iterations > size)
                throw new std::string("The maximum number of iterations for finding the span exceeded.");
        }
        return mid;
    }

    static Eigen::MatrixXd BasisFunctionsAndDerivatives(int der, double parameter, int spanIdx, int rDegree, const std::vector<double> &rKnots)
    {
        Eigen::MatrixXd ndu(rDegree+1, rDegree+1);

        ndu(0,0) = 1.0;

        Eigen::VectorXd left(rDegree + 1);
        Eigen::VectorXd right(rDegree + 1);

        for (int j = 1; j <= rDegree; j++)
        {
            left[j]  = parameter - rKnots[spanIdx + 1 - j];
            right[j] = rKnots[spanIdx + j] - parameter;
            double saved = 0.;
            for(int r = 0; r < j; r++)
            {
                ndu(j,r) = right[r+1] + left[j-r]; // lower triange (knot differences)
                double temp = ndu(r, j-1)/ndu(j,r);
                ndu(r,j) = saved + right[r+1]*temp; // upper triangle
                saved = left[j-r]*temp;
            }
            ndu(j,j) = saved;
        }

        Eigen::MatrixXd ders(der+1, rDegree+1);

        for(int j = 0; j <= rDegree; j++) ders(0,j) = ndu(j,rDegree);

        Eigen::MatrixXd a(2, rDegree+1);
        for(int r = 0; r<=rDegree; r++)
        {
            int s1 = 0;
            int s2 = 1;

            a(0,0) = 1.0;

            int j = 0;
            for(int k = 1; k<=der; k++)
            {
                double d = 0.;
                int rk = r-k;
                int pk = rDegree-k;
                if(r >= k)
                {
                    a(s2,0) = a(s1,0)/ndu(pk+1, rk);
                    d = a(s2,0)*ndu(rk,pk);
                }
                int j1 = 0;
                if(rk >= -1) j1 = 1;
                else         j1 = -rk;

                int j2 = 0;
                if(r-1 <= pk) j2 = k-1;
                else          j2 = rDegree-r;

                for(j = j1; j <= j2; j++)
                {
                    a(s2, j) = (a(s1,j) - a(s1, j-1))/ndu(pk+1, rk+j);
                    d += a(s2,j)*ndu(rk+j, pk);
                }

                if(r <= pk)
                {
                    a(s2,k) = -a(s1,k-1)/ndu(pk+1,r);
                    d += a(s2,k)*ndu(r,pk);
                }

                ders(k,r) = d;
                j = s1;
                s1 = s2;
                s2 = j;
            }
        }

        int r = rDegree;
        for(int k = 1; k <= der; k++)
        {
            for (int j = 0; j <= rDegree; j++) ders(k,j) *=r;
            r*= (rDegree-k);
        }


        return ders;
    }

    static Eigen::VectorXd BasisFunctionsAndDerivativesRational(int der, double rParameter, int spanIdx, int rDegree, const std::vector<double> &rKnots, const std::vector<double> &rWeights)
    {
        assert(der >= 0 || der <= 2 );

        Eigen::MatrixXd ders = BasisFunctionsAndDerivatives(der, rParameter, spanIdx, rDegree, rKnots);

        // NURBS specific ...
        Eigen::VectorXd sum(der+1);
        sum.setZero(der+1);

        for(int i = 0; i <= rDegree; i++)
        {
            double weight = rWeights[spanIdx - rDegree + i];
            if     (der == 0)
                sum(0) += ders(0,i)*weight;
            else if(der == 1)
            {
                sum(0) += ders(0,i)*weight;
                sum(1) += ders(1,i)*weight;
            }
            else
            {
                sum(0) += ders(0,i)*weight;
                sum(1) += ders(1,i)*weight;
                sum(2) += ders(2,i)*weight;
            }
        }

        Eigen::VectorXd dersRat(rDegree+1);
        dersRat.setZero(rDegree+1);

        for(int i = 0; i <= rDegree; i++)
        {
            double weight = rWeights[spanIdx - rDegree + i];
            if     (der == 0)
                dersRat(i) = ders(0, i)*weight/sum(0);
            else if(der == 1)
                dersRat(i) = (ders(1, i)*sum(0) - ders(0,i)*sum(1))* weight/(sum(0)*sum(0));
            else
            {
                double sum2 = sum(0)*sum(0);
                dersRat(i) = weight*(ders(2, i)/sum(0) - 2*ders(1,i)*sum(1)/(sum2) - ders(0,i)*sum(2)/(sum2) + 2*ders(0,i)*sum(1)*sum(1)/(sum2*sum(0)));
            }
        }

        return dersRat;
    }

    const std::array<int, TDimParameter> FindSpan(const Eigen::Matrix<double, TDimParameter, 1> &parameter) const
    {
       std::array<int, TDimParameter> parameterIDs;
       for(int i = 0; i < parameter.rows(); i++)
       {
           parameterIDs[i] = FindSpan(parameter[i], mDegree[i], mKnots[i]);
       }
       return parameterIDs;
    }

    Eigen::MatrixXd BasisFunctionsAndDerivativesRational(int der, const Eigen::Matrix<double, TDimParameter, 1> &parameter) const
    {
        throw new std::string("Not implemented for an arbitrary dimension!");
    }

private:
    //! @brief Knot vector (in isogeometric framework each segment between two knots is an element)
    std::array<std::vector<double>, TDimParameter> mKnots;

    //! @brief Control points of the BSpline curve (# rows = num control points)
    std::vector<std::vector<NuTo::NodeSimple>> mControlPoints;

    //! @brief Weights to NURBS
    std::vector<std::vector<double>> mWeights;

    //! @brief Degree of the polynomials (order = mDegree+1)
    std::array<int,TDimParameter> mDegree;
};

template<>
Eigen::VectorXd NURBS<1>::GetControlPointsElement(std::array<int, 1> knotID) const
{
    int dim = GetDimension();

    int numCPs = GetNumControlPointsElement();

    Eigen::VectorXd nodeValues(numCPs*dim);

    for (int i = 0; i < numCPs; i++)
        nodeValues.segment(dim * i, dim) = mControlPoints[0][knotID[0] - mDegree[0] + i].GetValues();

    return nodeValues;
}

template<>
Eigen::VectorXd NURBS<2>::GetControlPointsElement(std::array<int, 2> knotID) const
{
    int dim = GetDimension();

    int numCPs = GetNumControlPointsElement();

    Eigen::VectorXd nodeValues(numCPs*dim);

    int count = 0;
    for (int i = 0; i < mDegree[1]+1; i++)
    {
        for (int j = 0; j < mDegree[0]+1; j++)
        {
            nodeValues.segment(count, dim) = mControlPoints[knotID[1] - mDegree[1] + i][knotID[0] - mDegree[0] + j].GetValues();
            count+=dim;
        }
    }
    return nodeValues;
}

template<>
Eigen::MatrixXd NURBS<1>::BasisFunctionsAndDerivativesRational(int der, const Eigen::Matrix<double, 1, 1> &parameter) const
{
    assert(der >= 0 || der <= 2);

    int spanIdx = FindSpan(parameter)[0];
    Eigen::MatrixXd ders = BasisFunctionsAndDerivatives(der, parameter[0], spanIdx, mDegree[0], mKnots[0]);

    // NURBS specific ...
    Eigen::VectorXd sum(der+1);
    sum.setZero(der+1);

    for(int i = 0; i <= mDegree[0]; i++)
    {
        double weight = mWeights[0][spanIdx - mDegree[0] + i];
        if     (der == 0)
            sum(0) += ders(0,i)*weight;
        else if(der == 1)
        {
            sum(0) += ders(0,i)*weight;
            sum(1) += ders(1,i)*weight;
        }
        else
        {
            sum(0) += ders(0,i)*weight;
            sum(1) += ders(1,i)*weight;
            sum(2) += ders(2,i)*weight;
        }
    }

    Eigen::VectorXd dersRat(mDegree[0]+1);
    dersRat.setZero(mDegree[0]+1);

    for(int i = 0; i <= mDegree[0]; i++)
    {
        double weight = mWeights[0][spanIdx - mDegree[0] + i];
        if     (der == 0)
            dersRat(i) = ders(0, i)*weight/sum(0);
        else if(der == 1)
            dersRat(i) = (ders(1, i)*sum(0) - ders(0,i)*sum(1))* weight/(sum(0)*sum(0));
        else
        {
            double sum2 = sum(0)*sum(0);
            dersRat(i) = weight*(ders(2, i)/sum(0) - 2*ders(1,i)*sum(1)/(sum2) - ders(0,i)*sum(2)/(sum2) + 2*ders(0,i)*sum(1)*sum(1)/(sum2*sum(0)));
        }
    }
    return dersRat;
}

template<>
Eigen::MatrixXd NURBS<2>::BasisFunctionsAndDerivativesRational(int der, const Eigen::Matrix<double, 2, 1> &parameter) const
{
    assert(der >= 0 || der <= 2 );

    const std::array<int, 2> spanIdx = FindSpan(parameter);
    std::array<Eigen::MatrixXd, 2> shapeFunctions;

    int numDers = 1;
    for(int i = 0; i < 2; i++)
    {
         shapeFunctions[i] = BasisFunctionsAndDerivatives(der, parameter[i], spanIdx[i], mDegree[i], mKnots[i]);
         numDers *= mDegree[i]+1;
    }

    Eigen::MatrixXd ders(numDers, der+1);
    ders.setZero(numDers, der+1);

    int num = ( (der+1)*(der+1) + (der+1) )/2;
    Eigen::VectorXd sum(num);
    sum.setZero(num);

    for(int i = 0; i <= mDegree[1]; i++)
    {
        for(int j = 0; j <= mDegree[0]; j++)
        {
            double weight = mWeights[spanIdx[1] - mDegree[1] + i][spanIdx[0] - mDegree[0] + j];

            if     (der == 0)
            {
                sum(0) += shapeFunctions[0](0, j)*shapeFunctions[1](0, i)*weight;
            }
            else if(der == 1)
            {
                sum(0) += shapeFunctions[0](0, j)*shapeFunctions[1](0, i)*weight;
                sum(1) += shapeFunctions[0](1, j)*shapeFunctions[1](0, i)*weight;
                sum(2) += shapeFunctions[0](0, j)*shapeFunctions[1](1, i)*weight;
            }
            else
            {
                sum(0) += shapeFunctions[0](0, j)*shapeFunctions[1](0, i)*weight;
                sum(1) += shapeFunctions[0](1, j)*shapeFunctions[1](0, i)*weight;
                sum(2) += shapeFunctions[0](2, j)*shapeFunctions[1](0, i)*weight;

                sum(3) += shapeFunctions[0](0, j)*shapeFunctions[1](1, i)*weight;
                sum(4) += shapeFunctions[0](0, j)*shapeFunctions[1](2, i)*weight;

                sum(5) += shapeFunctions[0](1, j)*shapeFunctions[1](1, i)*weight;
            }
        }
    }

    int count = 0;
    for(int i = 0; i <= mDegree[1]; i++)
    {
        for(int j = 0; j <= mDegree[0]; j++)
        {
            double weight = mWeights[spanIdx[1] - mDegree[1] + i][spanIdx[0] - mDegree[0] + j];

            if     (der == 0)
            {
                ders(count,0) = shapeFunctions[0](0,j)*shapeFunctions[1](0,i)*weight/sum(0);
                count++;
            }
            else if(der == 1)
            {
                ders(count,0) = (shapeFunctions[0](1,j)*shapeFunctions[1](0,i)*sum(0)
                               - shapeFunctions[0](0,j)*shapeFunctions[1](0,i)*sum(1)) * weight/(sum(0)*sum(0));
                ders(count,1) = (shapeFunctions[0](0,j)*shapeFunctions[1](1,i)*sum(0)
                               - shapeFunctions[0](0,j)*shapeFunctions[1](0,i)*sum(2)) * weight/(sum(0)*sum(0));
                count++;
            }
            else
            {
                ders(count,0) =(shapeFunctions[0](2,j)*shapeFunctions[1](0,i)/sum(0)
                            - 2*shapeFunctions[0](1,j)*shapeFunctions[1](0,i)*sum(1)/(sum(0)*sum(0))
                              - shapeFunctions[0](0,j)*shapeFunctions[1](0,i)*sum(2)/(sum(0)*sum(0))
                            + 2*shapeFunctions[0](0,j)*shapeFunctions[1](0,i)*sum(1)*sum(1)/(sum(0)*sum(0)*sum(0)) ) * weight;

                ders(count,1) = (shapeFunctions[0](0,j)*shapeFunctions[1](2,i)/sum(0)
                             - 2*shapeFunctions[0](0,j)*shapeFunctions[1](1,i)*sum(3)/(sum(0)*sum(0))
                               - shapeFunctions[0](0,j)*shapeFunctions[1](0,i)*sum(4)/(sum(0)*sum(0))
                             + 2*shapeFunctions[0](0,j)*shapeFunctions[0](0,i)*sum(3)*sum(3)/(sum(0)*sum(0)*sum(0)) ) * weight;

                ders(count,2) = (shapeFunctions[0](1,j)*shapeFunctions[1](1,i)/sum(0)
                               - shapeFunctions[0](1,j)*shapeFunctions[1](0,i)*sum(3)/(sum(0)*sum(0))
                               - shapeFunctions[0](0,j)*shapeFunctions[1](1,i)*sum(1)/(sum(0)*sum(0))
                               - shapeFunctions[0](1,j)*shapeFunctions[1](0,i)*sum(5)*sum(5)/(sum(0)*sum(0))
                             + 2*shapeFunctions[0](0,j)*shapeFunctions[1](0,i)*sum(1)*sum(3)/(sum(0)*sum(0)*sum(0)) ) * weight;
                count++;
            }
        }
    }

    return ders;
}

}
