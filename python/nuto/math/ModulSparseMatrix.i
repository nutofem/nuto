%module ModulSparseMatrix
%feature("autodoc","1");
//remove the warning for not exposing the base class features of eigen to python via swig
#pragma SWIG nowarn=401
%{
//Put headers and other declarations here to be added in the wrapper files
#include "math/Operator.h"
#include "math/SparseMatrix.h"
#include "math/SparseMatrixCSR.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "math/SparseMatrixCSRSymmetric.h"
#include "math/SparseMatrixCSRVector2.h"
#include "math/SparseMatrixCSRVector2General.h"
#include "math/SparseMatrixCSRVector2Symmetric.h"
%}

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

// convert python string to std::string
%include "std_string.i"
// convert python tuple to std::vector
%include "std_vector.i"
// use exceptions, but build no interface for NUTO::Exception
%ignore Exception;
%include "base/ModulNuToBase.i"

%import "math/NuToMath.i"
%import "math/ModulMatrix.i"

%include "math/SparseMatrix.h"
%include "math/SparseMatrixCSR.h"
%include "math/SparseMatrixCSRGeneral_Def.h"
%include "math/SparseMatrixCSRSymmetric_Def.h"
%include "math/SparseMatrixCSRVector2.h"
%include "math/SparseMatrixCSRVector2General_Def.h"
%include "math/SparseMatrixCSRVector2Symmetric_Def.h"
%template(DoubleSparseMatrix) NuTo::SparseMatrix<double>;
%template(IntSparseMatrix) NuTo::SparseMatrix<int>;
%template(DoubleSparseMatrixCSR) NuTo::SparseMatrixCSR<double>;
%template(IntSparseMatrixCSR) NuTo::SparseMatrixCSR<int>;
%template(DoubleSparseMatrixCSRGeneral) NuTo::SparseMatrixCSRGeneral<double>;
%template(IntSparseMatrixCSRGeneral) NuTo::SparseMatrixCSRGeneral<int>;
%template(DoubleSparseMatrixCSRSymmetric) NuTo::SparseMatrixCSRSymmetric<double>;
%template(IntSparseMatrixCSRSymmetric) NuTo::SparseMatrixCSRSymmetric<int>;
%template(DoubleSparseMatrixCSRVector2) NuTo::SparseMatrixCSRVector2<double>;
%template(IntSparseMatrixCSRVector2) NuTo::SparseMatrixCSRVector2<int>;
%template(DoubleSparseMatrixCSRVector2General) NuTo::SparseMatrixCSRVector2General<double>;
%template(IntSparseMatrixCSRVector2General) NuTo::SparseMatrixCSRVector2General<int>;
%template(DoubleSparseMatrixCSRVector2Symmetric) NuTo::SparseMatrixCSRVector2Symmetric<double>;
%template(IntSparseMatrixCSRVector2Symmetric) NuTo::SparseMatrixCSRVector2Symmetric<int>;
