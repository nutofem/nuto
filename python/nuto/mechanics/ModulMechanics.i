%module(package="nuto") ModulMechanics
%feature("autodoc","1");
%{
//Put headers and other declarations here to be added in the wrapper files
#include "math/SparseMatrix.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "mechanics/dofSubMatrixStorage/BlockFullMatrix.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/MechanicsException.h"
#include "mechanics/nodes/NodeEnum.h"
#include "base/Logger.h"
#include "mechanics/elements/ElementEnum.h"
#include "base/CallbackInterface.h"
#include "mechanics/timeIntegration/TimeIntegrationBase.h"
#include "mechanics/timeIntegration/NewmarkBase.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
%}

// typemaps.i is a built-in swig interface that lets us map c++ types to other
// types python. We'll use it to map Eigen matrices to Numpy arrays.
%include "typemaps.i"
%include "eigen.i"
%eigen_typemaps(Eigen::VectorXd)
%eigen_typemaps(Eigen::MatrixXd)
%eigen_typemaps(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>)
%eigen_typemaps(Eigen::Matrix<double, Eigen::Dynamic, 1>)
%eigen_typemaps(Eigen::Matrix<double, 1, 1>)
%eigen_typemaps(Eigen::Matrix<double, 2, 2>)
%eigen_typemaps(Eigen::Matrix<double, 3, 3>)
%eigen_typemaps(Eigen::Matrix<double, 2, 1>)
%eigen_typemaps(Eigen::Matrix<double, 3, 1>)
%eigen_typemaps(Eigen::Matrix<double, 6, 1>)

%eigen_typemaps(Eigen::VectorXi)
%eigen_typemaps(Eigen::MatrixXi)
%eigen_typemaps(Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>)
%eigen_typemaps(Eigen::Matrix<int, Eigen::Dynamic, 1>)
// convert python string to std::string
%include "std_string.i"
// convert python tuple to std::vector
%include "std_vector.i"
namespace std {
   %template(IntVector) vector<int>;
   %template(DoubleVector) vector<double>;
}

%ignore Exception;
%include "base/ModulNuToBase.i"
%include "base/CallbackInterface.h"

%include "mechanics/structures/StructureBase.h"
%include "mechanics/structures/unstructured/Structure.h"
%include "mechanics/elements/ElementEnum.h"

%include "mechanics/dofSubMatrixStorage/BlockStorageBase.h"
%include "mechanics/dofSubMatrixStorage/DofStatus.h"
%include "mechanics/dofSubMatrixStorage/BlockScalar.h"
%include "mechanics/dofSubMatrixStorage/BlockFullMatrix.h"
%include "mechanics/dofSubMatrixStorage/BlockFullVector.h"
%include "mechanics/dofSubMatrixStorage/BlockSparseMatrix.h"

%include "mechanics/nodes/NodeEnum.h"

%include "mechanics/structures/StructureOutputBase.h"
%include "mechanics/structures/StructureOutputBlockMatrix.h"
%include "mechanics/structures/StructureOutputBlockVector.h"

%include "mechanics/timeIntegration/TimeIntegrationBase.h"
%include "mechanics/timeIntegration/NewmarkBase.h"
%include "mechanics/timeIntegration/NewmarkDirect.h"

%template(DoubleBlockFullMatrix) NuTo::BlockFullMatrix<double>;
%template(DoubleBlockFullVector) NuTo::BlockFullVector<double>;
%template(IntBlockFullVector) NuTo::BlockFullVector<int>;



