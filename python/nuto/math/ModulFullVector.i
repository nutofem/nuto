%module ModulFullVector
%feature("autodoc","1");
//remove the warning for not exposing the base class features of eigen to python via swig
#pragma SWIG nowarn=401
%{
//Put headers and other declarations here to be added in the wrapper files
#include "math/FullMatrix.h"
#include "math/FullVector.h"
#include "math/Operator.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "math/SparseMatrixCSRSymmetric.h"
#include "math/SparseMatrixCSRVector2General.h"
%}

// convert python string to std::string
%include "std_string.i"
// convert python tuple to std::vector
%include "std_vector.i"
// use exceptions, but build no interface for NUTO::Exception
%ignore Exception;
%include "base/ModulNuToBase.i"

%import "math/NuToMath.i"
%import "math/ModulMatrix.i"
%import "math/ModulFullMatrix.i"

%include "math/FullVector_Def.h"

//this is for vectors
// extend the python-interface with FullVectoroperators, since the c++ operators are only defined in the base class (which is not exposed to python) 
%extend NuTo::FullVector<double,Eigen::Dynamic> 
{
  FullVector<double,Eigen::Dynamic> operator*(NuTo::FullVector<double,Eigen::Dynamic> rOther)
  {
      return NuTo::FullVector<double,Eigen::Dynamic>($self->operator*(rOther));
  }
}

%extend NuTo::FullVector<double,Eigen::Dynamic>
{
  FullVector<double,Eigen::Dynamic> operator+(NuTo::FullVector<double,Eigen::Dynamic> rOther)
  {
      return NuTo::FullVector<double,Eigen::Dynamic>($self->operator+(rOther));
  }
}

%extend NuTo::FullVector<double,Eigen::Dynamic>
{
  FullVector<double,Eigen::Dynamic> operator+ (NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> rOther)
  {
      return NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>($self->operator+(rOther));
  }
}

%extend NuTo::FullVector<double,Eigen::Dynamic> 
{
  FullVector<double,Eigen::Dynamic> operator-(NuTo::FullVector<double,Eigen::Dynamic> rOther)
  {
      return NuTo::FullVector<double,Eigen::Dynamic>($self->operator-(rOther));
  }
}

%extend NuTo::FullVector<double,Eigen::Dynamic> 
{
  FullVector<double,Eigen::Dynamic>& operator+=(const NuTo::FullVector<double,Eigen::Dynamic> rOther)
  {
      $self->operator+=(rOther);
      return *$self;
  }
}

%extend NuTo::FullVector<double,Eigen::Dynamic> 
{
  FullVector<double,Eigen::Dynamic> operator*(double rOther)
  {
      return NuTo::FullVector<double,Eigen::Dynamic>($self->operator*(rOther));
  }
}

%extend NuTo::FullVector<int,Eigen::Dynamic> 
{
  FullVector<int,Eigen::Dynamic> operator*(NuTo::FullVector<int,Eigen::Dynamic> rOther)
  {
      return NuTo::FullVector<int,Eigen::Dynamic>($self->operator*(rOther));
  }
}

%extend NuTo::FullVector<int,Eigen::Dynamic> 
{
  FullVector<int,Eigen::Dynamic> operator+(NuTo::FullVector<int,Eigen::Dynamic> rOther)
  {
      return NuTo::FullVector<int,Eigen::Dynamic>($self->operator+(rOther));
  }
}

%extend NuTo::FullVector<int,Eigen::Dynamic> 
{
  FullVector<int,Eigen::Dynamic> operator-(NuTo::FullVector<int,Eigen::Dynamic> rOther)
  {
      return NuTo::FullVector<int,Eigen::Dynamic>($self->operator-(rOther));
  }
}

%extend NuTo::FullVector<int,Eigen::Dynamic> 
{
  FullVector<int,Eigen::Dynamic> operator*(int rOther)
  {
      return NuTo::FullVector<int,Eigen::Dynamic>($self->operator*(rOther));
  }
}

%extend NuTo::FullVector<int,Eigen::Dynamic> 
{
  FullVector<int,Eigen::Dynamic>& operator+=(const NuTo::FullVector<int,Eigen::Dynamic> rOther)
  {
      $self->operator+=(rOther);
      return *$self;
  }
}


%template(DoubleFullVector) NuTo::FullVector<double,Eigen::Dynamic>;
%template(IntFullVector) NuTo::FullVector<int,Eigen::Dynamic>;
