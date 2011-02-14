%module(package="nuto") ModulOptimizer
%feature("autodoc","1");
%{
//Put headers and other declarations here to be added in the wrapper files
#include "nuto/math/FullMatrix.h"
#include "nuto/optimize/Optimizer.h"
#include "nuto/optimize/ConjugateGradientGrid.h"
#include "nuto/optimize/ConjugateGradientNonLinear.h"
#include "nuto/optimize/CallbackHandler.h"
#include "nuto/optimize/CallbackHandlerGrid.h"
#include "nuto/optimize/CallbackHandlerPython.h"
#include "nuto/optimize/OptimizeException.h"
%}

// convert python string to std::string
%include "std_string.i"
// convert python tuple to std::vector
%include "std_vector.i"
// use exceptions, but build no interface for NUTO::Exception
%ignore Exception;
%include "nuto/base/ModulNuToBase.i"

%include "nuto/optimize/CallbackHandler.h"
%include "nuto/optimize/CallbackHandlerGrid.h"
%include "nuto/optimize/CallbackHandlerPython.h"
%include "nuto/optimize/Optimizer.h"
%include "nuto/optimize/ConjugateGradientGrid.h"
%include "nuto/optimize/ConjugateGradientNonLinear.h"
