#include "nuto/mechanics/tools/TimeDependentProblem.h"
#include "nuto/mechanics/dofs/DofNumbering.h"

using namespace NuTo;

TimeDependentProblem::TimeDependentProblem(MeshFem* rMesh)
    : mMerger(rMesh)
{
}

GlobalDofVector TimeDependentProblem::RenumberDofs(Constraint::Constraints constraints, std::vector<DofType> dofTypes,
                                                   GlobalDofVector oldDofValues)
{
    DofInfo dofInfos;

    GlobalDofVector renumberedValues;

    for (auto dofType : dofTypes)
    {
        if (oldDofValues.J[dofType].rows() != 0) // initial case: No merge of an uninitialized vector
            mMerger.Merge(oldDofValues, {dofType});

        auto dofInfo = DofNumbering::Build(mMerger.Nodes(dofType), dofType, constraints);
        dofInfos.Merge(dofType, dofInfo);

        renumberedValues.J[dofType].resize(dofInfo.numIndependentDofs[dofType]);
        renumberedValues.K[dofType].resize(dofInfo.numDependentDofs[dofType]);

        mMerger.Extract(&renumberedValues, {dofType});
    }
    mAssembler.SetDofInfo(dofInfos);
    return renumberedValues;
}

void TimeDependentProblem::AddGradientFunction(Group<CellInterface> group, GradientFunction f)
{
    mGradientFunctions.push_back({group, f});
}

void TimeDependentProblem::AddHessian0Function(Group<CellInterface> group, HessianFunction f)
{
    mHessian0Functions.push_back({group, f});
}

void TimeDependentProblem::AddUpdateFunction(Group<CellInterface> group, UpdateFunction f)
{
    mUpdateFunctions.push_back({group, f});
}

template <typename TCellInterfaceFunction, typename TTimeDepFunction>
TCellInterfaceFunction Apply(TTimeDepFunction& f, double t, double dt)
{
    using namespace std::placeholders;
    return std::bind(f, _1, t, dt);
}

GlobalDofVector TimeDependentProblem::Gradient(const GlobalDofVector& dofValues, std::vector<DofType> dofs, double t,
                                               double dt)
{
    mMerger.Merge(dofValues, dofs);
    GlobalDofVector gradient;
    for (auto& gradientFunction : mGradientFunctions)
        gradient += mAssembler.BuildVector(gradientFunction.first, dofs,
                                           Apply<CellInterface::VectorFunction>(gradientFunction.second, t, dt));
    return gradient;
}

GlobalDofMatrixSparse TimeDependentProblem::Hessian0(const GlobalDofVector& dofValues, std::vector<DofType> dofs,
                                                     double t, double dt)
{
    mMerger.Merge(dofValues, dofs);
    GlobalDofMatrixSparse hessian0;
    for (auto& hessian0Function : mHessian0Functions)
        hessian0 += mAssembler.BuildMatrix(hessian0Function.first, dofs,
                                           Apply<CellInterface::MatrixFunction>(hessian0Function.second, t, dt));
    return hessian0;
}

void TimeDependentProblem::UpdateHistory(const GlobalDofVector& dofValues, std::vector<DofType> dofs, double t,
                                         double dt)
{
    mMerger.Merge(dofValues, dofs);
    for (auto& updateFunction : mUpdateFunctions)
        for (auto& cell : updateFunction.first)
            cell.Apply(Apply<CellInterface::VoidFunction>(updateFunction.second, t, dt));
}
