#include "mechanics/tools/EquationSystem.h"

using namespace NuTo;

EquationSystem::EquationSystem(MeshFem* rMesh)
    : mMerger(rMesh)
{
}

void EquationSystem::SetDofInfo(DofInfo dofInfo)
{
    mAssembler.SetDofInfo(dofInfo);
}

void EquationSystem::AddGradientFunction(Group<CellInterface> group, GradientFunction f)
{
    mGradientFunctions.push_back({group, f});
}

void EquationSystem::AddHessian0Function(Group<CellInterface> group, HessianFunction f)
{
    mHessian0Functions.push_back({group, f});
}

void EquationSystem::AddUpdateFunction(Group<CellInterface> group, UpdateFunction f)
{
    mUpdateFunctions.push_back({group, f});
}

template <typename TCellInterfaceFunction, typename TTimeDepFunction>
TCellInterfaceFunction Apply(TTimeDepFunction& f, double t, double dt)
{
    using namespace std::placeholders;
    return std::bind(f, _1, _2, t, dt);
};

GlobalDofVector EquationSystem::Gradient(GlobalDofVector dofValues, std::vector<DofType> dofs, double t, double dt)
{
    mMerger.Merge(dofValues, dofs);
    GlobalDofVector gradient;
    for (auto& gradientFunction : mGradientFunctions)
        gradient += mAssembler.BuildVector(gradientFunction.first, dofs,
                                           Apply<CellInterface::VectorFunction>(gradientFunction.second, t, dt));
    return gradient;
}

GlobalDofMatrixSparse EquationSystem::Hessian0(GlobalDofVector dofValues, std::vector<DofType> dofs, double t,
                                               double dt)
{
    mMerger.Merge(dofValues, dofs);
    GlobalDofMatrixSparse hessian0;
    for (auto& hessian0Function : mHessian0Functions)
        hessian0 += mAssembler.BuildMatrix(hessian0Function.first, dofs,
                                           Apply<CellInterface::MatrixFunction>(hessian0Function.second, t, dt));
    return hessian0;
}

void EquationSystem::UpdateHistory(GlobalDofVector dofValues, std::vector<DofType> dofs, double t, double dt)
{
    mMerger.Merge(dofValues, dofs);
    for (auto& updateFunction : mUpdateFunctions)
        for (auto& cell : updateFunction.first)
            cell.Apply(Apply<CellInterface::VoidFunction>(updateFunction.second, t, dt));
}
