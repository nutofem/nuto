#pragma once

#include "mechanics/iga/BSplineCurve.h"
#include "mechanics/nodes/NodeEnum.h"
#include <set>

namespace NuTo
{
//! @author Peter Otto, BAM
//! @date July, 2016
//! @brief ... class for B spline surfaces, with IGA specific functions
//! @brief ... B spline/NURBS specific algorithms taken from Piegl, Tiller 'The NURBS book' 1996

class Structure;

class BSplineSurface
{
public:
    /** Constructors **/

    //! @brief ... default constructor
    BSplineSurface()
    {
    }

    //! @brief ... constructor
    //! @param rDegree ... degree of the polynomial
    //! @param rKnots ... knot vector
    //! @param rControlPoints ... control points
    BSplineSurface(const Eigen::Vector2i& rDegree, const Eigen::VectorXd& rKnotsX, const Eigen::VectorXd& rKnotsY,
                   const Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic>& rControlPoints,
                   const Eigen::MatrixXd& rWeights);

    //! @brief ... constructor (interpolation of a point cloud)
    //! @param rDegree ... degree of the polynomial
    //! @param rPoints ... points to interpolate
    BSplineSurface(const Eigen::Vector2i& rDegree,
                   const Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic>& rPoints,
                   Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& AInv);

    /** Getter **/

    // --- Control points --- //

    int GetDimension() const
    {
        return mControlPoints(0).rows();
    }

    int GetNumControlPoints() const;

    int GetNumControlPoints(int dir) const;

    Eigen::VectorXd GetControlPoint(int rControlPointIDY, int rControlPointIDX) const;

    const Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic>& GetControlPoints() const
    {
        return mControlPoints;
    }

    const Eigen::MatrixXd& GetWeights() const
    {
        return mWeights;
    }

    // --- Knot vector --- //

    int GetNumKnots(int dir) const;

    const Eigen::VectorXd& GetKnotVector(int dir) const;

    int GetElementFirstKnotID(int rElementIDinDir, int dir) const;

    Eigen::MatrixXd GetElementKnots(int rElementIDX, int rElementIDY) const;

    Eigen::VectorXi GetElementKnotIDs(int rElementIDX, int rElementIDY) const;

    int GetMultiplicityOfKnot(double rKnot, int dir) const;

    // --- Elements --- //

    int GetNumIGAElements(int dir) const;

    int GetNumIGAElements() const;

    Eigen::VectorXi GetElementControlPointIDs(int rElementIDX, int rElementIDY) const;

    Eigen::VectorXi GetElementControlPointIDsGlobal(int rElementIDX, int rElementIDY,
                                                    const Eigen::MatrixXi& rNodeIDs) const;

    Eigen::MatrixXd GetKnotIDControlPoints(const Eigen::Vector2i& rKnotIDs) const;

    // --- Points --- //

    Eigen::VectorXd SurfacePoint(const Eigen::Vector2d& rParameter) const;

    Eigen::MatrixXd SurfacePoints(const Eigen::MatrixXd& rParameter) const;

    // --- Degree --- //

    const Eigen::Vector2i& GetDegree() const
    {
        return mDegree;
    }

    /** Parametrization **/

    /** Knot refinement **/
    void DuplicateKnots(int dir);

    void InsertKnot(double rKnotToInsert, int rMultiplicity, int dir);

    void RefineKnots(const Eigen::VectorXd& rKnotsToInsert, int dir);

    /** build structure **/

    void buildIGAStructure(NuTo::Structure& rStructure, const std::set<NuTo::Node::eDof>& setOfDOFS, int rGroupElements,
                           int rGroupNodes);

private:
    //! @brief Knot vector (in isogeometric framework each segment between two knots is an element)
    Eigen::VectorXd mKnotsX;

    //! @brief Knot vector (in isogeometric framework each segment between two knots is an element)
    Eigen::VectorXd mKnotsY;

    //! @brief Degree of the polynomials
    Eigen::Vector2i mDegree;

    //! @brief weights for nurbs
    Eigen::MatrixXd mWeights;

    //! @brief Control points of the BSpline curve
    Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> mControlPoints;
};
}
