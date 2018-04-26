#pragma once
#include <Eigen/Core>


namespace NuTo
{

typedef Eigen::VectorXd NodeValues;
// typedef Eigen::VectorXd NaturalCoords;
typedef Eigen::MatrixXd DerivativeShapeFunctionsNatural;
typedef Eigen::MatrixXd DerivativeShapeFunctionsGlobal;
typedef Eigen::MatrixXd BMatrixStrain;
typedef Eigen::MatrixXd BMatrixGradient;
typedef Eigen::MatrixXd NMatrix;
//
// constexpr int maxNumNodes = 125;
constexpr int maxDim = 3;
// typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor, NuTo::maxNumNodes * NuTo::maxDim, 1> NodeValues;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor, NuTo::maxDim, 1> NaturalCoords;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor, NuTo::maxDim, 1> GlobalCoords;
// typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor, NuTo::maxNumNodes, 1> ShapeFunctions;
// typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, NuTo::maxNumNodes, NuTo::maxDim>
//        DerivativeShapeFunctionsNatural;
// typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, NuTo::maxNumNodes, NuTo::maxDim>
//        DerivativeShapeFunctionsGlobal;
// typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, 6, NuTo::maxDim * NuTo::maxNumNodes>
//        BMatrixStrain;
// typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, NuTo::maxDim, NuTo::maxNumNodes>
//        BMatrixGradient;
// typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, NuTo::maxDim, NuTo::maxNumNodes>
// NMatrix;

} /* NuTo */
