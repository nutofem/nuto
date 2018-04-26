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
namespace ShapeFunctions1D
{
////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 1, 1> NodeCoordinatesTrussOrder1(int rNodeIndex);

Eigen::Matrix<double, 2, 1> ShapeFunctionsTrussOrder1(const Eigen::VectorXd& rCoordinates);

Eigen::Matrix<double, 2, 1> DerivativeShapeFunctionsTrussOrder1();

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 1, 1> NodeCoordinatesTrussOrder2(int rNodeIndex);

Eigen::Matrix<double, 3, 1> ShapeFunctionsTrussOrder2(const Eigen::VectorXd& rCoordinates);

Eigen::Matrix<double, 3, 1> DerivativeShapeFunctionsTrussOrder2(const Eigen::VectorXd& rCoordinates);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 1, 1> NodeCoordinatesTrussOrder3(int rNodeIndex);

Eigen::Matrix<double, 4, 1> ShapeFunctionsTrussOrder3(const Eigen::VectorXd& rCoordinates);

Eigen::Matrix<double, 4, 1> DerivativeShapeFunctionsTrussOrder3(const Eigen::VectorXd& rCoordinates);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 1, 1> NodeCoordinatesTrussOrder4(int rNodeIndex);

Eigen::Matrix<double, 5, 1> ShapeFunctionsTrussOrder4(const Eigen::VectorXd& rCoordinates);

Eigen::Matrix<double, 5, 1> DerivativeShapeFunctionsTrussOrder4(const Eigen::VectorXd& rCoordinates);
}

namespace ShapeFunctions2D
{


////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 2, 1> NodeCoordinatesTriangleOrder2(int rNodeIndex);

Eigen::Matrix<double, 6, 1> ShapeFunctionsTriangleOrder2(const Eigen::VectorXd& rCoordinates);

Eigen::Matrix<double, 6, 2> DerivativeShapeFunctionsTriangleOrder2(const Eigen::VectorXd& rCoordinates);

////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 2, 1> NodeCoordinatesTriangleOrder3(int rNodeIndex);

Eigen::Matrix<double, 10, 1> ShapeFunctionsTriangleOrder3(const Eigen::VectorXd& rCoordinates);

Eigen::Matrix<double, 10, 2> DerivativeShapeFunctionsTriangleOrder3(const Eigen::VectorXd& rCoordinates);

////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 2, 1> NodeCoordinatesTriangleOrder4(int rNodeIndex);

Eigen::Matrix<double, 15, 1> ShapeFunctionsTriangleOrder4(const Eigen::VectorXd& rCoordinates);

Eigen::Matrix<double, 15, 2> DerivativeShapeFunctionsTriangleOrder4(const Eigen::VectorXd& rCoordinates);

////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 2, 1> NodeCoordinatesQuadOrder1(int rNodeIndex);

Eigen::Matrix<double, 4, 1> ShapeFunctionsQuadOrder1(const Eigen::VectorXd& rCoordinates);

Eigen::Matrix<double, 4, 2> DerivativeShapeFunctionsQuadOrder1(const Eigen::VectorXd& rCoordinates);

////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 2, 1> NodeCoordinatesQuadOrder2(int rNodeIndex);

Eigen::Matrix<double, 8, 1> ShapeFunctionsQuadOrder2(const Eigen::VectorXd& rCoordinates);

Eigen::Matrix<double, 8, 2> DerivativeShapeFunctionsQuadOrder2(const Eigen::VectorXd& rCoordinates);

////////////////////////////////////////////////////////////////////////////////////////////////////////////
}


namespace ShapeFunctions3D
{

////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 3, 1> NodeCoordinatesTetrahedronOrder1(int rNodeIndex);

Eigen::Matrix<double, 4, 1> ShapeFunctionsTetrahedronOrder1(const Eigen::VectorXd& rCoordinates);

Eigen::Matrix<double, 4, 3> DerivativeShapeFunctionsTetrahedronOrder1();

////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 3, 1> NodeCoordinatesTetrahedronOrder2(int rNodeIndex);

Eigen::Matrix<double, 10, 1> ShapeFunctionsTetrahedronOrder2(const Eigen::VectorXd& rCoordinates);

Eigen::Matrix<double, 10, 3> DerivativeShapeFunctionsTetrahedronOrder2(const Eigen::VectorXd& rCoordinates);

////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 3, 1> NodeCoordinatesBrickOrder1(int rNodeIndex);

Eigen::Matrix<double, 8, 1> ShapeFunctionsBrickOrder1(const Eigen::VectorXd& rCoordinates);

Eigen::Matrix<double, 8, 3> DerivativeShapeFunctionsBrickOrder1(const Eigen::VectorXd& rCoordinates);

////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 3, 1> NodeCoordinatesBrickOrder2(int rNodeIndex);

Eigen::Matrix<double, 20, 1> ShapeFunctionsBrickOrder2(const Eigen::VectorXd& rCoordinates);

Eigen::Matrix<double, 20, 3> DerivativeShapeFunctionsBrickOrder2(const Eigen::VectorXd& rCoordinates);

////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 3, 1> NodeCoordinatesPrismOrder1(int rNodeIndex);

Eigen::Matrix<double, 6, 1> ShapeFunctionsPrismOrder1(const Eigen::VectorXd& rCoordinates);

Eigen::Matrix<double, 6, 3> DerivativeShapeFunctionsPrismOrder1(const Eigen::VectorXd& rCoordinates);

////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 3, 1> NodeCoordinatesPrismOrder2(int rNodeIndex);

Eigen::Matrix<double, 18, 1> ShapeFunctionsPrismOrder2(const Eigen::VectorXd& rCoordinates);

Eigen::Matrix<double, 18, 3> DerivativeShapeFunctionsPrismOrder2(const Eigen::VectorXd& rCoordinates);

////////////////////////////////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix<double, 3, 1> NodeCoordinatesPyramidOrder1(int rNodeIndex);

Eigen::Matrix<double, 5, 1> ShapeFunctionsPyramidOrder1(const Eigen::VectorXd& rCoordinates);

Eigen::Matrix<double, 5, 3> DerivativeShapeFunctionsPyramidOrder1(const Eigen::VectorXd& rCoordinates);

////////////////////////////////////////////////////////////////////////////////////////////////////////////
}


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
