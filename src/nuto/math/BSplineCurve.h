#pragma once

#include "nuto/math/FullMatrix.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/export.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

namespace NuTo
{
class BSplineCurve
{
public:
    /** Constructors **/

    //! @brief ... default constructor
    BSplineCurve(){}

    //! @brief ... constructor
    //! @param rDegree ... degree of the polynomial
    //! @param rKnots ... knot vector
    //! @param rControlPoints ... control points
    BSplineCurve(int rDegree,
                 const FullVector<double, Eigen::Dynamic>& rKnots,
                 const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rControlPoints);

    //! @brief ... constructor
    //! @param rDegree ... degree of the polynomial
    //! @param rPoints ... points to interpolate
    BSplineCurve(int rDegree,
                 const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rPoints);

    /** Getter **/

    //! @brief ... get the degree of the underlying polynomials
    //! @return ... degree of the polynomial
    int GetDegree(){return mDegree;}

    //! @brief ... get the dimension of the curve
    //! @return ... dimension of the curve
    int GetDimension(){return mControlPoints.GetNumColumns();}

    //! @brief ... get the number of control points
    //! @return ... number of control
    int GetNumControlPoints(){return mControlPoints.GetNumRows();}

    //! @brief ... find span of a parameter (rParameter) in the knot vector
    //! @param rParameter ... parameter to find the span of
    //! @return ... the span index
    int FindSpan(double rParameter);

    /** Basis Functions **/

    //! @brief ... find span of a parameter (rParameter) in the knot vector
    //! @param rParameter ... parameter to find the span of
    //! @return ... the span index
    void BasisFuns(double rParameter, int rSpan, NuTo::FullVector<double, Eigen::Dynamic>& rBasisFunctions);

    //! @brief ... find span of a parameter (rParameter) in the knot vector
    //! @param rParameter ... parameter to find the span of
    //! @return ... the span index
    void DerivativeBasisFuns(double rParameter, int rSpan, NuTo::FullVector<double, Eigen::Dynamic>& rDersBasisFunctions);

    /** Parametrization **/

    //! @brief ... parametrization for given points to interpolate (chord length method)
    //! @param rPoints ... points to interpolate
    void ParametrizationChordLengthMethod(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rPoints,
                                          FullVector<double, Eigen::Dynamic>& rParameters);

    //! @brief ... parametrization for given points to interpolate (centripetal method)
    //! @param rPoints ... points to interpolate
    void ParametrizationCentripetalMethod(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rPoints);

private:
    //! @brief Knot vector (in isogeometric framework each segment between two knots is a element)
    FullVector<double, Eigen::Dynamic> mKnots;
    //! @brief Control points of the BSpline curve
    FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> mControlPoints;
    //! @brief Degree of the polynomials (order = mDegree+1)
    int mDegree;
};
}
