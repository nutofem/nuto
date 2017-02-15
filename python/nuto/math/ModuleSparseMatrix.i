%module ModuleSparseMatrix
%feature("autodoc","1");
%{
//Put headers and other declarations here to be added in the wrapper files
#include "math/SparseMatrix.h"
#include "math/SparseMatrixCSR.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "math/SparseMatrixCSRSymmetric.h"
#include "math/SparseMatrixCSRVector2.h"
#include "math/SparseMatrixCSRVector2General.h"
#include "math/SparseMatrixCSRVector2Symmetric.h"
%}


// use exceptions, but build no interface for NUTO::Exception
%ignore Exception;
%include "base/ModuleNuToBase.i"

%include "math/NuToMath.i" // defines typenames for std::vector and Eigen::Matrix

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
