%module ModulSparseSolver
%feature("autodoc","1");
//remove the warning for not exposing the base class features of eigen to python via swig
#pragma SWIG nowarn=401
%{
//Put headers and other declarations here to be added in the wrapper files
#include "math/FullVector.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "math/SparseMatrixCSRSymmetric.h"
#include "math/SparseMatrixCSRVector2General.h"
#include "math/SparseDirectSolver.h"
#include "math/SparseDirectSolverMKLDSS.h"
#include "math/SparseDirectSolverMKLPardiso.h"
#include "math/SparseDirectSolverMUMPS.h"
%}

// convert python string to std::string
%include "std_string.i"
// convert python tuple to std::vector
%include "std_vector.i"
// use exceptions, but build no interface for NUTO::Exception
%ignore Exception;
%include "base/ModulNuToBase.i"

%import "math/NuToMath.i"
%import "math/ModulMatrix.i"
%import "math/ModulFullMatrix.i"
%import "math/ModulFullVector.i"
%import "math/ModulSparseMatrix.i"

/* solver */
%include "math/SparseDirectSolver.h"
%include "math/SparseDirectSolverMKLDSS.h"
%include "math/SparseDirectSolverMKLPardiso.h"
%include "math/SparseDirectSolverMUMPS.h"
