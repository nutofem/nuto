%module(package="nuto") ModulFullMatrix
%feature("autodoc","1");
//remove the warning for not exposing the base class features of eigen to python via swig
#pragma SWIG nowarn=401
%{
//Put headers and other declarations here to be added in the wrapper files
#include "nuto/math/FullMatrix.h"
#include "nuto/math/Operator.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseMatrixCSRSymmetric.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"
%}

// convert python string to std::string
%include "std_string.i"
// convert python tuple to std::vector
%include "std_vector.i"
// use exceptions, but build no interface for NUTO::Exception
%ignore Exception;
%include "nuto/base/ModulNuToBase.i"

%import "nuto/math/NuToMath.i"
%import "nuto/math/ModulMatrix.i"
%import "nuto/math/ModulSparseMatrix.i"

%include "nuto/math/FullMatrix_Def.h"

// extend the python-interface with FullMatrixoperators, since the c++ operators are only defined in the base class (which is not exposed to pyhton) 
%extend NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> 
{
  FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> operator*(NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> rOther)
  {
      return NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>($self->operator*(rOther));
  }
}

%extend NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> 
{
  FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> operator+(NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> rOther)
  {
      return NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>($self->operator+(rOther));
  }
}

%extend NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> 
{
  FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> operator-(NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> rOther)
  {
      return NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>($self->operator-(rOther));
  }
}

%extend NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> 
{
  FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& operator+=(const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> rOther)
  {
      $self->operator+=(rOther);
      return *$self;
  }
}

%extend NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> 
{
  FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> operator*(double rOther)
  {
      return NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>($self->operator*(rOther));
  }
}

%extend NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> 
{
  FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> operator*(NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> rOther)
  {
      return NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic>($self->operator*(rOther));
  }
}

%extend NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> 
{
  FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> operator+(NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> rOther)
  {
      return NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic>($self->operator+(rOther));
  }
}

%extend NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> 
{
  FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> operator-(NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> rOther)
  {
      return NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic>($self->operator-(rOther));
  }
}

%extend NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> 
{
  FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> operator*(int rOther)
  {
      return NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic>($self->operator*(rOther));
  }
}

%extend NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> 
{
  FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic>& operator+=(const NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> rOther)
  {
      $self->operator+=(rOther);
      return *$self;
  }
}
%template(DoubleFullMatrix) NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>;
%template(IntFullMatrix) NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic>;
