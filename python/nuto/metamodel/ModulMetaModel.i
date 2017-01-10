%module(package="nuto") ModulMetaModel
%feature("autodoc","1");
%{
//Put headers and other declarations here to be added in the wrapper files
#include "metamodel/NeuralNetwork.h"
%}

%include "math/NuToMath.i" // defines typenames for std::vector and Eigen::Matrix

%include "metamodel/Metamodel.h"
%include "optimize/CallbackHandler.h"
%include "metamodel/NeuralNetwork.h"

