/*
 * ElementShapeFunctions.h
 *
 *  Created on: 30 Mar 2015
 *      Author: ttitsche
 */

#pragma once

#include <Eigen/Core>

namespace NuTo
{

namespace ShapeFunctionsInterface2D
{
////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::MatrixXd NodeCoordinatesInterface2dOrder1(int rNodeIndex);

Eigen::MatrixXd ShapeFunctionsInterface2dOrder1(const Eigen::VectorXd& rCoordinates);

Eigen::MatrixXd DerivativeShapeFunctionsInterface2dOrder1();

Eigen::MatrixXd NodeCoordinatesInterface2dOrder2(int rNodeIndex);

Eigen::MatrixXd ShapeFunctionsInterface2dOrder2(const Eigen::VectorXd& rCoordinates);

Eigen::MatrixXd DerivativeShapeFunctionsInterface2dOrder2(const Eigen::VectorXd& rCoordinates);

////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

// In order to maintain equal shape functions on each element the Bézier extraction together with Bernstein polynomials
// is used (see. Borden et. al. 2011)
namespace ShapeFunctionsIGA
{
/////////////////////////////// BSPLINE ////////////////////////////////////////////////////////////////////////
//! @brief find span of a parameter (rParameter) in the knot vector
//! @param rParameter parameter to find the span of
//!
//! @return the span index
int FindSpan(double rParameter, int rDegree, const Eigen::VectorXd& rKnots);

Eigen::VectorXd BasisFunctions(double rParameter, int spanIdx, int rDegree, const Eigen::VectorXd& rKnots);
Eigen::MatrixXd BasisFunctionsAndDerivatives(int der, double rParameter, int spanIdx, int rDegree,
                                             const Eigen::VectorXd& rKnots);

Eigen::VectorXd BasisFunctionsRat(double rParameter, int spanIdx, int rDegree, const Eigen::VectorXd& rKnots,
                                  const Eigen::VectorXd& rWeights);
Eigen::VectorXd BasisFunctionsAndDerivativesRat(int der, double rParameter, int spanIdx, int rDegree,
                                                const Eigen::VectorXd& rKnots, const Eigen::VectorXd& rWeights);

//! @brief calculate the derivative of NURBS basis functions
//! @param rCoordinates parameters (u,v)
//! @param rSpanIdx span index of u and v
//! @param rDegree degree in x and y directions
//! @param rKnotsX knot vector in x direction
//! @param rKnotsY knot vector in y direction
//! @param rWeights weights for NURBS (if every entry == 1 => B spline)
//! @return the vector containing the derivatives
Eigen::VectorXd BasisFunctions2DRat(const Eigen::VectorXd& rCoordinates, const Eigen::Vector2i& rSpanIdx,
                                    const Eigen::Vector2i& rDegree, const Eigen::VectorXd& rKnotsX,
                                    const Eigen::VectorXd& rKnotsY, const Eigen::MatrixXd& rWeights);

//! @brief calculate the derivative of NURBS basis functions
//! @param der the derivative to calculate (possibilities: 0,1,2 -- no more implemented)
//! @param rCoordinates parameters (u,v)
//! @param rSpanIdx span index of u and v
//! @param rDegree degree in x and y directions
//! @param rKnotsX knot vector in x direction
//! @param rKnotsY knot vector in y direction
//! @param rWeights weights for NURBS (if every entry == 1 => B spline)
//! @return the matrix containing the derivatives
Eigen::MatrixXd BasisFunctionsAndDerivatives2DRat(int der, const Eigen::VectorXd& rCoordinates,
                                                  const Eigen::Vector2i& rSpanIdx, const Eigen::Vector2i& rDegree,
                                                  const Eigen::VectorXd& rKnotsX, const Eigen::VectorXd& rKnotsY,
                                                  const Eigen::MatrixXd& rWeights);


/////////////////////////////// BERNSTEIN //////////////////////////////////////////////////////////////////////

Eigen::VectorXd Bernstein1DOrder1(double rParameter);
Eigen::VectorXd Bernstein1DOrder2(double rParameter);
Eigen::VectorXd Bernstein1DOrder3(double rParameter);
Eigen::VectorXd Bernstein1DOrder4(double rParameter);
Eigen::VectorXd Bernstein1D(double rParameter, int rOrder);


Eigen::VectorXd DerivativeBernstein1DOrder1(double rParameter);
Eigen::VectorXd DerivativeBernstein1DOrder2(double rParameter);
Eigen::VectorXd DerivativeBernstein1DOrder3(double rParameter);
Eigen::VectorXd DerivativeBernstein1DOrder4(double rParameter);
Eigen::VectorXd DerivativeBernstein1D(double rParameter, int rOrder);


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
} // namespace ShapeFunctionsIGA


} /* namespace NuTo */
