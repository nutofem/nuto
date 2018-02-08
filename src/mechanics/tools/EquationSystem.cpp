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

void EquationSystem::AddGradientFunction(Group<CellInterface> group, CellInterface::VectorFunction f)
{
    mGradientFunctions.push_back({group, f});
}

void EquationSystem::AddHessian0Function(Group<CellInterface> group, CellInterface::MatrixFunction f)
{
    mHessian0Functions.push_back({group, f});
}

void EquationSystem::AddUpdateFunction(Group<CellInterface> group, CellInterface::VoidFunction f)
{
    mUpdateFunctions.push_back({group, f});
}

GlobalDofVector EquationSystem::Gradient(GlobalDofVector dofValues, std::vector<DofType> dofs)
{
    mMerger.Merge(dofValues, dofs);
    GlobalDofVector gradient;
    for (auto& gradientFunction : mGradientFunctions)
        gradient += mAssembler.BuildVector(gradientFunction.first, dofs, gradientFunction.second);
    return gradient;
}

GlobalDofMatrixSparse EquationSystem::Hessian0(GlobalDofVector dofValues, std::vector<DofType> dofs)
{
    mMerger.Merge(dofValues, dofs);
    GlobalDofMatrixSparse hessian0;
    for (auto& hessian0Function : mHessian0Functions)
        hessian0 += mAssembler.BuildMatrix(hessian0Function.first, dofs, hessian0Function.second);
    return hessian0;
}

void EquationSystem::UpdateHistory(GlobalDofVector dofValues, std::vector<DofType> dofs)
{
    mMerger.Merge(dofValues, dofs);
    for (auto& updateFunction : mUpdateFunctions)
        for (auto& cell : updateFunction.first)
            cell.Apply(updateFunction.second);
}
