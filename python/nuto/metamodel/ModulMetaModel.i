%module(package="nuto") ModulMetaModel
%feature("autodoc","1");
%{
//Put headers and other declarations here to be added in the wrapper files
#include "math/FullMatrix.h"
#include "metamodel/NeuralNetwork.h"
#include "metamodel/MultipleLinearRegression.h"
%}

// convert python string to std::string
%include "std_string.i"
// convert python tuple to std::vector
%include "std_vector.i"
%ignore Exception;
%include "base/ModulNuToBase.i"

%include "metamodel/Metamodel.h"
%include "metamodel/NeuralNetwork.h"
%include "metamodel/MultipleLinearRegression.h"
