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
#include "mechanics/DirectionEnum.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeDof.h"
#include "base/Logger.h"
#include "mechanics/elements/ElementEnum.h"
#include "base/CallbackInterface.h"
#include "mechanics/timeIntegration/TimeIntegrationBase.h"
#include "mechanics/timeIntegration/NewmarkBase.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/timeIntegration/RungeKuttaBase.h"
#include "mechanics/timeIntegration/RungeKutta4.h"
#include "mechanics/mesh/MeshGenerator.h"
#include "mechanics/constraints/Term.h"
#include "mechanics/constraints/Equation.h"
#include "mechanics/constraints/Constraints.h"
#include "mechanics/constraints/ConstraintCompanion.h" 
#include "mechanics/structures/Assembler.h"

#include <stdexcept>

class PyCallback
{
    PyObject *func;
    PyCallback& operator=(const PyCallback&) = delete;
public:
    PyCallback(const PyCallback& o) : func(o.func)
    {
        Py_XINCREF(func);
    }

    PyCallback(PyObject *func) : func(func)
    {
        Py_XINCREF(this->func);
        assert(PyCallable_Check(this->func) and
                "Not a Python function.");
    }

    ~PyCallback()
    {
        Py_XDECREF(func);
    }

    double operator()(double x)
    {
        if (!func || Py_None == func || !PyCallable_Check(func))
            throw std::runtime_error("Not a Python function.");

        PyObject *args = Py_BuildValue("(d)", x);
        PyObject *pyResult = PyObject_CallObject(func, args);
        double result = PyFloat_AsDouble(pyResult);

        Py_DECREF(args);
        Py_XDECREF(pyResult);
        return result;
    }
};

using namespace NuTo;
using namespace NuTo::Constraint;
%}


%include "math/NuToMath.i" // defines typenames for std::vector and Eigen::Matrix

%include "mechanics/Sections.i"

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

%include "mechanics/DirectionEnum.h"

%include "mechanics/nodes/NodeEnum.h"
%include "mechanics/nodes/NodeBase.h"
%include "mechanics/nodes/NodeDof.h"

%include "mechanics/structures/StructureOutputBase.h"
%include "mechanics/structures/StructureOutputBlockMatrix.h"
%include "mechanics/structures/StructureOutputBlockVector.h"

%typemap(in) std::vector<NuTo::eDirection>
%{
    {
    // Hello World
    std::vector<NuTo::eDirection> dirVec;
    PyObject* seq = PySequence_Fast($input, "Argument not a list.");
    int len = PySequence_Size($input);
    for (int i = 0; i < len; i++)
    {
        PyObject* item  = PySequence_Fast_GET_ITEM(seq, i);
        dirVec.push_back(static_cast<eDirection>(PyInt_AsLong(item)));
    }

    $1 = dirVec;
    }
%}

%typemap(typecheck, precedence=141) std::vector<NuTo::eDirection>
%{
    $1 = PyList_Check($input) ? 1 : 0;
%}

%include "mechanics/constraints/Term.h"
%include "mechanics/constraints/Equation.h"
%include "mechanics/constraints/Constraints.h"
%include "mechanics/constraints/ConstraintCompanion.h" 
%include "mechanics/structures/Assembler.h"

%include "mechanics/timeIntegration/TimeIntegrationBase.h"
%include "mechanics/timeIntegration/NewmarkBase.h"
%include "mechanics/timeIntegration/NewmarkDirect.h"
%include "mechanics/timeIntegration/RungeKuttaBase.h"
%include "mechanics/timeIntegration/RungeKutta4.h"

%include "mechanics/mesh/MeshGenerator.h"

%template(DoubleBlockFullMatrix) NuTo::BlockFullMatrix<double>;
%template(DoubleBlockFullVector) NuTo::BlockFullVector<double>;
%template(IntBlockFullVector) NuTo::BlockFullVector<int>;

%extend NuTo::Constraint::Equation
{
    Equation(PyObject* callback)
    {
        return new Equation(std::function<double(double)>(PyCallback(callback)));
    }
};

namespace std {
    // these are necessary for vectors of types with no default constructor
    %ignore vector<NuTo::Constraint::Equation>::vector(size_type);
    %ignore vector<NuTo::Constraint::Equation>::resize;
    %template(EquationVector) vector<NuTo::Constraint::Equation>;
}
