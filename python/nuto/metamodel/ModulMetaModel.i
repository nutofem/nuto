%module(package="nuto") ModulMetaModel
%feature("autodoc","1");
%{
//Put headers and other declarations here to be added in the wrapper files
#include "metamodel/NeuralNetwork.h"
%}

%include "typemaps.i"
%include "eigen.i"
%eigen_typemaps(Eigen::VectorXd)
%eigen_typemaps(Eigen::MatrixXd)
%eigen_typemaps(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>)
%eigen_typemaps(Eigen::Matrix<double, Eigen::Dynamic, 1>)
%eigen_typemaps(Eigen::Matrix<double, 1, 1>)
%eigen_typemaps(Eigen::Matrix<double, 2, 2>)
%eigen_typemaps(Eigen::Matrix<double, 3, 3>)
%eigen_typemaps(Eigen::Matrix<double, 2, 1>)
%eigen_typemaps(Eigen::Matrix<double, 3, 1>)
%eigen_typemaps(Eigen::Matrix<double, 6, 1>)
// convert python string to std::string
%include "std_string.i"
// convert python tuple to std::vector
%include "std_vector.i"
namespace std {
   %template(IntVector) vector<int>;
   %template(DoubleVector) vector<double>;
}

%ignore Exception;
%include "base/ModulNuToBase.i"

%include "metamodel/Metamodel.h"
%include "optimize/CallbackHandler.h"
%include "metamodel/NeuralNetwork.h"
