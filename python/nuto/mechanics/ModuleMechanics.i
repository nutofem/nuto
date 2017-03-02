%module(package="nuto") ModuleMechanics
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
#include "mechanics/timeIntegration/RungeKuttaBase.h"
#include "mechanics/timeIntegration/RungeKutta4.h"
#include "mechanics/mesh/MeshGenerator.h"
%}


%include "math/NuToMath.i" // defines typenames for std::vector and Eigen::Matrix
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
%include "mechanics/timeIntegration/RungeKuttaBase.h"
%include "mechanics/timeIntegration/RungeKutta4.h"

%include "mechanics/mesh/MeshGenerator.h"

%template(DoubleBlockFullMatrix) NuTo::BlockFullMatrix<double>;
%template(DoubleBlockFullVector) NuTo::BlockFullVector<double>;
%template(IntBlockFullVector) NuTo::BlockFullVector<int>;




