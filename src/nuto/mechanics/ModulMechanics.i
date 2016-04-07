%module(package="nuto") ModulMechanics
%feature("autodoc","1");
%{
//Put headers and other declarations here to be added in the wrapper files
#include "nuto/math/FullMatrix.h"
#include "nuto/math/FullVector.h"
#include "nuto/math/SparseMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/base/Logger.h"
#include "nuto/mechanics/elements/ElementEnum.h"
%}

// convert python string to std::string
%include "std_string.i"
// convert python tuple to std::vector
%include "std_vector.i"
%ignore Exception;
%include "nuto/base/ModulNuToBase.i"


%include "nuto/mechanics/structures/StructureBase.h"
%include "nuto/mechanics/structures/unstructured/Structure.h"
%include "nuto/mechanics/elements/ElementEnum.h"

%include "nuto/mechanics/dofSubMatrixStorage/BlockStorageBase.h"
%include "nuto/mechanics/dofSubMatrixStorage/DofStatus.h"
%include "nuto/mechanics/dofSubMatrixStorage/BlockScalar.h"
%include "nuto/mechanics/dofSubMatrixStorage/BlockFullMatrix.h"
%include "nuto/mechanics/dofSubMatrixStorage/BlockFullVector.h"
%include "nuto/mechanics/dofSubMatrixStorage/BlockSparseMatrix.h"

%include "nuto/mechanics/structures/StructureOutputBase.h"
%include "nuto/mechanics/structures/StructureOutputBlockMatrix.h"
%include "nuto/mechanics/structures/StructureOutputBlockVector.h"

%template(DoubleBlockFullMatrix) NuTo::BlockFullMatrix<double>;
%template(DoubleBlockFullVector) NuTo::BlockFullVector<double>;
%template(IntBlockFullVector) NuTo::BlockFullVector<int>;




