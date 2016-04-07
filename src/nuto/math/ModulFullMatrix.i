%module ModulFullMatrix
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

//// parts added
#ifdef ENABLE_NUMPY
%include "numpy.i"


%init %{
import_array();
%}


%apply (double * IN_ARRAY2, int DIM1, int DIM2){(double * inData,int rRow, int rCol)};
%apply (double* INPLACE_ARRAY2, int DIM1,int DIM2){(double * indata,int rRow,int rCol)};

%apply (int* IN_ARRAY2, int DIM1, int DIM2){(int * inData,int rRow, int rCol)};
%apply (int* INPLACE_ARRAY2, int DIM1,int DIM2){(int * indata,int rRow,int rCol)};
#endif


//////////////////////////


// convert python string to std::string
%include "std_string.i"
// convert python tuple to std::vector
// %include "std_vector.i"

// use exceptions, but build no interface for NUTO::Exception
%ignore Exception;
%include "nuto/base/ModulNuToBase.i"

%import "nuto/math/NuToMath.i"
%import "nuto/math/ModulMatrix.i"
%import "nuto/math/ModulSparseMatrix.i"

%include "nuto/math/FullMatrix_Def.h"

//this is for matrices
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
  FullVector<double,Eigen::Dynamic> operator*(NuTo::FullVector<double,Eigen::Dynamic> rOther)
  {
      return NuTo::FullVector<double,Eigen::Dynamic>($self->operator*(rOther));
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

// this is necessary to be able for the derived vector class call the base class methods 
%template(DoubleFullVectorBase) NuTo::FullMatrix<double,Eigen::Dynamic,1>;
%template(IntFullVectorBase) NuTo::FullMatrix<int,Eigen::Dynamic,1>;
