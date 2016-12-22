%module(package="nuto") ModulOptimizer
%feature("autodoc","1");
%{
//Put headers and other declarations here to be added in the wrapper files
#include "math/FullMatrix.h"
#include "optimize/Optimizer.h"
#include "optimize/ConjugateGradientGrid.h"
#include "optimize/ConjugateGradientNonLinear.h"
#include "optimize/CallbackHandler.h"
#include "optimize/CallbackHandlerGrid.h"
#include "optimize/CallbackHandlerPython.h"
#include "optimize/Jacobi.h"
#include "optimize/MisesWielandt.h"
#include "optimize/OptimizeException.h"
%}

// convert python string to std::string
%include "std_string.i"
// convert python tuple to std::vector
%include "std_vector.i"
// use exceptions, but build no interface for NUTO::Exception
%ignore Exception;
%include "base/ModulNuToBase.i"

%include "optimize/CallbackHandler.h"
%include "optimize/CallbackHandlerGrid.h"
%include "optimize/CallbackHandlerPython.h"
%include "optimize/Optimizer.h"
%include "optimize/ConjugateGradientGrid.h"
%include "optimize/ConjugateGradientNonLinear.h"
%include "optimize/Jacobi.h"
%include "optimize/MisesWielandt.h"
