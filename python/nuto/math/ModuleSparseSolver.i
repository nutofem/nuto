%module ModuleSparseSolver
%feature("autodoc","1");
//remove the warning for not exposing the base class features of eigen to python via swig
#pragma SWIG nowarn=401
%{
//Put headers and other declarations here to be added in the wrapper files
#include "math/SparseMatrixCSRGeneral.h"
#include "math/SparseMatrixCSRSymmetric.h"
#include "math/SparseMatrixCSRVector2General.h"
#include "math/SparseDirectSolver.h"
#include "math/SparseDirectSolverMKLPardiso.h"
#include "math/SparseDirectSolverMUMPS.h"
%}


%include "math/ModuleSparseMatrix.i"

/* solver */
%include "math/SparseDirectSolver.h"
%include "math/SparseDirectSolverMKLPardiso.h"
%include "math/SparseDirectSolverMUMPS.h"
