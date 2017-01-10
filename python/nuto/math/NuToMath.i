// $Id$
%module NuToMath
%feature("autodoc","1");
%{
//Put headers and other declarations here to be added in the wrapper files
#include "math/NuToMath.h"
%}

//declare inputs of the functions to be used as output on python level
%include "typemaps.i"
// typemaps.i is a built-in swig interface that lets us map c++ types to other
// types python. We'll use it to map Eigen matrices to Numpy arrays.
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

%eigen_typemaps(Eigen::VectorXi)
%eigen_typemaps(Eigen::MatrixXi)
%eigen_typemaps(Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>)
%eigen_typemaps(Eigen::Matrix<int, Eigen::Dynamic, 1>)

// convert python string to std::string
%include "std_string.i"
// convert python tuple to std::vector
%include "std_vector.i"
namespace std {
   %template(IntVector) vector<int>;
   %template(DoubleVector) vector<double>;
}

// use exceptions, but build no interface for NUTO::Exception
%ignore Exception;
%include "base/ModulNuToBase.i"

// provide the public interface
%include "math/NuToMath.h"

