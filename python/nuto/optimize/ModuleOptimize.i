%module(package="nuto") ModuleOptimize
%feature("autodoc","1");
%{
//Put headers and other declarations here to be added in the wrapper files
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

%include "math/NuToMath.i" // defines typenames for std::vector and Eigen::Matrix

%include "optimize/CallbackHandler.h"
%include "optimize/CallbackHandlerGrid.h"
%include "optimize/CallbackHandlerPython.h"
%include "optimize/Optimizer.h"
%include "optimize/ConjugateGradientGrid.h"
%include "optimize/ConjugateGradientNonLinear.h"
%include "optimize/Jacobi.h"
%include "optimize/MisesWielandt.h"
