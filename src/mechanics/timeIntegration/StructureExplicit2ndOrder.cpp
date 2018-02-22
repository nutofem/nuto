#include "mechanics/timeIntegration/StructureExplicit2ndOrder.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/Assembler.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"

NuTo::TimeIntegration::StructureRhsExplicit2ndOrder::StructureRhsExplicit2ndOrder(NuTo::Structure& s)
    : mS(s)
    , mHessian2(s.BuildGlobalHessian2Lumped())
{
    mHessian2.ApplyCMatrix(mS.GetAssembler().GetConstraintMatrix());
    mHessian2.CwiseInvert();
}

void NuTo::TimeIntegration::StructureRhsExplicit2ndOrder::
operator()(const StructureStateExplicit2ndOrder& x, StructureStateExplicit2ndOrder& dxdt, const double t)
{
    // ------------------------------------------------------------------
    // merge x with structure, update dependent dofs including velocities
    // ------------------------------------------------------------------
    auto cmat = mS.GetAssembler().GetConstraintMatrix();
    mS.GetAssembler().ConstraintUpdateRhs(t);
    auto dof0K = mS.NodeCalculateDependentDofValues(x.dof0);
    auto dof1K = -1. * (cmat * x.dof1);
    mS.NodeMergeDofValues(0, x.dof0, dof0K);
    mS.NodeMergeDofValues(1, x.dof1, dof1K);
    // -----------------------
    // compute right hand side
    // -----------------------
    NuTo::StructureOutputBlockVector Fext = mS.BuildGlobalExternalLoadVector();
    NuTo::StructureOutputBlockVector Fint = mS.BuildGlobalInternalGradient();
    NuTo::StructureOutputBlockVector Ftotal = Fext - Fint;
    NuTo::BlockFullVector<double> Fmod = NuTo::Assembler::ApplyCMatrix(Ftotal, cmat);

    dxdt.dof0 = x.dof1;
    dxdt.dof1 = mHessian2.JJ * Fmod;
}
