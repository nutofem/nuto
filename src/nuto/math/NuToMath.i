// $Id$
%module NuToMath
%feature("autodoc","1");
//remove the warning for not exposing the base class features of eigen to python via swig
#pragma SWIG nowarn=401
%{
//Put headers and other declarations here to be added in the wrapper files
#include "nuto/math/NuToMath.h"
#include <cstddef> // normally not necessary  
%}

// convert python string to std::string
%include "std_string.i"
// convert python tuple to std::vector
%include "std_vector.i"
// use exceptions, but build no interface for NUTO::Exception
%ignore Exception;
%include "nuto/base/ModulNuToBase.i"

// provide the public interface
%include "nuto/math/NuToMath.h"

// define some functions
double sin(double);
double cos(double);
double tan(double);

#define I_CONST 5             // An integer constant
#define Pi      3.14159       // A Floating point constant
#define S_CONST "hello world" // A string constant
#define NEWLINE '\n'          // Character constant
