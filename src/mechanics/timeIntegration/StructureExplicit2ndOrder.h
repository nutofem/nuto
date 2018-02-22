#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "boost/operators.hpp"

namespace NuTo
{
namespace TimeIntegration
{
class StructureStateExplicit2ndOrder : boost::addable<StructureStateExplicit2ndOrder>,
                                       boost::multipliable<StructureStateExplicit2ndOrder, double>
{
public:
    NuTo::BlockFullVector<double> dof0;
    NuTo::BlockFullVector<double> dof1;

    StructureStateExplicit2ndOrder(const StructureStateExplicit2ndOrder& m)
        : dof0(m.dof0)
        , dof1(m.dof1)
    {
    }

    StructureStateExplicit2ndOrder(NuTo::BlockFullVector<double> a, NuTo::BlockFullVector<double> b)
        : dof0(a)
        , dof1(b)
    {
    }

    StructureStateExplicit2ndOrder& operator+=(const StructureStateExplicit2ndOrder& p)
    {
        dof0 += p.dof0;
        dof1 += p.dof1;
        return *this;
    }

    StructureStateExplicit2ndOrder& operator*=(const double& a)
    {
        dof0 *= a;
        dof1 *= a;
        return *this;
    }
};

class StructureRhsExplicit2ndOrder
{

    NuTo::Structure& mS;
    NuTo::StructureOutputBlockMatrix mHessian2;

public:
    StructureRhsExplicit2ndOrder(NuTo::Structure& s);

    void operator()(const StructureStateExplicit2ndOrder& x, StructureStateExplicit2ndOrder& dxdt, const double t);
};
}
}
