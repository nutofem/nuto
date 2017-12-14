#pragma once
#include <Eigen/Core>


namespace NuTo
{

typedef Eigen::VectorXd NodeValues;
typedef Eigen::VectorXd NaturalCoords;
typedef Eigen::VectorXd ShapeFunctions;
typedef Eigen::MatrixXd DerivativeShapeFunctionsNatural;
typedef Eigen::MatrixXd DerivativeShapeFunctionsGlobal;
typedef Eigen::MatrixXd BMatrixStrain;
typedef Eigen::MatrixXd BMatrixGradient;
typedef Eigen::MatrixXd NMatrix;

} /* NuTo */
