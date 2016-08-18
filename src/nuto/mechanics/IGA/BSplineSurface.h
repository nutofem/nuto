#pragma once

#include "nuto/mechanics/IGA/BSplineCurve.h"

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
class BSplineSurface
{
public:
    /** Constructors **/

    //! @brief ... default constructor
    BSplineSurface(){}

    //! @brief ... constructor
    //! @param rDegree ... degree of the polynomial
    //! @param rKnots ... knot vector
    //! @param rControlPoints ... control points
    BSplineSurface(const Eigen::Vector2i &rDegree,
                   const Eigen::VectorXd &rKnotsX,
                   const Eigen::VectorXd &rKnotsY,
                   const Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> &rControlPoints);

    //! @brief ... constructor (interpolation of a point cloud)
    //! @param rDegree ... degree of the polynomial
    //! @param rPoints ... points to interpolate
    BSplineSurface(const Eigen::Vector2i &rDegree,
                   const Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> &rPoints,
                   Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &AInv);

    /** Getter **/

    // Control points //
    int GetNumControlPoints() const;
    int GetNumControlPoints(int dir) const;
    Eigen::VectorXd GetControlPoint(int rControlPointIDX, int rControlPointIDY) const;


    // Knot vector //
    const Eigen::VectorXd& GetKnotVector(int dir) const;
    int GetElementFirstKnotID(int rElementIDinDir, int dir) const;
    Eigen::MatrixXd GetElementKnots(int rElementIDX, int rElementIDY) const;
    Eigen::VectorXi GetElementKnotIDs(int rElementIDX, int rElementIDY) const;

    // Elements //
    int GetNumIGAElements(int dir) const;
    Eigen::VectorXi GetElementControlPointIDs(int rElementIDX, int rElementIDY) const;

    /** Parametrization **/

    //! @brief ... parametrization for given points to interpolate (chord length method)
    //! @param rPoints ... points to interpolate
    //! @return rParameters ... parameters to the given points
    void ParametrizationChordLengthMethod(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rPoints,
                                          FullVector<double, Eigen::Dynamic>& rParametersX,
                                          FullVector<double, Eigen::Dynamic>& rParametersY);

private:

//! @brief Knot vector (in isogeometric framework each segment between two knots is an element)
Eigen::VectorXd mKnotsX;

//! @brief Knot vector (in isogeometric framework each segment between two knots is an element)
Eigen::VectorXd mKnotsY;

//! @brief Degree of the polynomials
Eigen::Vector2i mDegree;

//! @brief Control points of the BSpline curve
Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> mControlPoints;



};
}
