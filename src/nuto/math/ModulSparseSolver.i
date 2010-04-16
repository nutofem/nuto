%module(package="nuto") ModulSparseSolver
%feature("autodoc","1");
%{
//Put headers and other declarations here to be added in the wrapper files
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseMatrixCSRSymmetric.h"
#include "nuto/math/SparseDirectSolver.h"
#include "nuto/math/SparseDirectSolverMKLDSS.h"
#include "nuto/math/SparseDirectSolverMKLPardiso.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
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
%import "nuto/math/ModulFullMatrix.i"
%import "nuto/math/ModulSparseMatrix.i"

/* solver */
%include "nuto/math/SparseDirectSolver.h"
%include "nuto/math/SparseDirectSolverMKLDSS.h"
%include "nuto/math/SparseDirectSolverMKLPardiso.h"
%include "nuto/math/SparseDirectSolverMUMPS.h"
