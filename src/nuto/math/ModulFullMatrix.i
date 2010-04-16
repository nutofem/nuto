%module(package="nuto") ModulFullMatrix
%feature("autodoc","1");
%{
//Put headers and other declarations here to be added in the wrapper files
#include "nuto/math/FullMatrix.h"
#include "nuto/math/Operator.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseMatrixCSRSymmetric.h"
#include <eigen2/Eigen/Array>
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
%import "nuto/math/ModulSparseMatrix.i"

%include "nuto/math/FullMatrix.h"
%template(DoubleFullMatrix) NuTo::FullMatrix<double>;
%template(IntFullMatrix) NuTo::FullMatrix<int>;
