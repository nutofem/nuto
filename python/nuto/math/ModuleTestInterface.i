%module ModuleTestInterface
%feature("autodoc","1");
%{
//Put headers and other declarations here to be added in the wrapper files
#include "math/TestNumpy.h"
%}

%include "math/NuToMath.i"

%import "math/NuToMath.i"
%include "math/TestNumpy.h"
