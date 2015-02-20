%module ModulSparseMatrix
%feature("autodoc","1");
//remove the warning for not exposing the base class features of eigen to python via swig
#pragma SWIG nowarn=401
%{
//Put headers and other declarations here to be added in the wrapper files
#include "nuto/math/FullMatrix.h"
#include "nuto/math/Operator.h"
#include "nuto/math/SparseMatrix.h"
#include "nuto/math/SparseMatrixCSR.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseMatrixCSRSymmetric.h"
#include "nuto/math/SparseMatrixCSRVector2.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"
#include "nuto/math/SparseMatrixCSRVector2Symmetric.h"
%}

// convert python string to std::string
%include "std_string.i"
// convert python tuple to std::vector
%include "std_vector.i"
// use exceptions, but build no interface for NUTO::Exception
%ignore Exception;
%include "nuto/base/ModulNuToBase.i"

%import "nuto/math/NuToMath.i"
%import "nuto/math/ModulMatrix.i"

%include "nuto/math/SparseMatrix.h"
%include "nuto/math/SparseMatrixCSR.h"
%include "nuto/math/SparseMatrixCSRGeneral_Def.h"
%include "nuto/math/SparseMatrixCSRSymmetric_Def.h"
%include "nuto/math/SparseMatrixCSRVector2.h"
%include "nuto/math/SparseMatrixCSRVector2General_Def.h"
%include "nuto/math/SparseMatrixCSRVector2Symmetric_Def.h"
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
