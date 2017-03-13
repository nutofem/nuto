#pragma once
#include <eigen3/Eigen/Core>


#define maxNumNodes 125
#define maxDim 3

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


/*
typedef Eigen::Matrix<double, Eigen::Dynamic, 1,Eigen::ColMajor,maxNumNodes*maxDim,1> NodeValues;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1,Eigen::ColMajor,maxDim,1> NaturalCoords;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1,Eigen::ColMajor,maxNumNodes,1>ShapeFunctions;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,Eigen::ColMajor,maxNumNodes,maxDim> DerivativeShapeFunctionsNatural;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,Eigen::ColMajor,maxNumNodes,maxDim> DerivativeShapeFunctionsGlobal;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,Eigen::ColMajor,6,maxDim*maxNumNodes> BMatrixStrain;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,Eigen::ColMajor,maxDim,maxNumNodes> BMatrixGradient;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,Eigen::ColMajor,maxDim,maxNumNodes> NMatrix;
*/

} /* NuTo */ 
