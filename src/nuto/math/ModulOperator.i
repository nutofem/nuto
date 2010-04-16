%module(package="nuto") ModulOperator
%feature("autodoc","1");
%{
//Put headers and other declarations here to be added in the wrapper files
#include "nuto/math/Operator.h"
%}

// convert python string to std::string
%include "std_string.i"
// convert python tuple to std::vector
%include "std_vector.i"
// use exceptions, but build no interface for NUTO::Exception
%ignore Exception;
%include "nuto/base/ModulNuToBase.i"

%include "nuto/math/Operator.h"

%template(DoubleMOperator)    NuTo::MonadicOperator<double>;
%template(IntMOperator)       NuTo::MonadicOperator<int>;
%template(DoubleDOperator)    NuTo::DyadicOperator<double>;
%template(IntDOperator)       NuTo::DyadicOperator<int>;

%template(DoubleMOperatorAbs) NuTo::MOperatorAbs<double>;
%template(IntMOperatorAbs)    NuTo::MOperatorAbs<int>;
%template(DoubleMOperatorMin) NuTo::MOperatorMin<double>;
%template(IntMOperatorMin)    NuTo::MOperatorMin<int>;
%template(DoubleMOperatorMax) NuTo::MOperatorMax<double>;
%template(IntMOperatorMax)    NuTo::MOperatorMax<int>;

%template(DoubleDOperatorMin) NuTo::DOperatorMin<double>;
%template(IntDOperatorMin)    NuTo::DOperatorMin<int>;
%template(DoubleDOperatorMax) NuTo::DOperatorMax<double>;
%template(IntDOperatorMax)    NuTo::DOperatorMax<int>;
%template(DoubleDOperatorAdd) NuTo::DOperatorAdd<double>;
%template(IntDOperatorAdd)    NuTo::DOperatorAdd<int>;
%template(DoubleDOperatorSub) NuTo::DOperatorSub<double>;
%template(IntDOperatorSub) NuTo::DOperatorSub<int>;
%template(DoubleDOperatorMul) NuTo::DOperatorMul<double>;
%template(IntDOperatorMul) NuTo::DOperatorMul<int>;
