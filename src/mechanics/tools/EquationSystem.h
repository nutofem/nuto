#pragma once
#include "mechanics/cell/CellInterface.h"
#include "mechanics/cell/SimpleAssembler.h"
#include "base/Group.h"
#include "mechanics/dofs/GlobalDofVector.h"
#include "mechanics/dofs/GlobalDofMatrixSparse.h"

#include "mechanics/tools/NodalValueMerger.h"

namespace NuTo
{

//! Equation system that contains Gradient(u) and its derivative with respect to u
class EquationSystem
{
public:
    EquationSystem(SimpleAssembler* rAssembler, MeshFem* rMesh);

    void AddGradientFunction(Group<CellInterface> group, CellInterface::VectorFunction f);
    void AddHessian0Function(Group<CellInterface> group, CellInterface::MatrixFunction f);
    void AddUpdateFunction(Group<CellInterface> group, CellInterface::VoidFunction f);

    GlobalDofVector Gradient(GlobalDofVector dofValues, std::vector<DofType> dofs);
    GlobalDofMatrixSparse Hessian0(GlobalDofVector dofValues, std::vector<DofType> dofs);

    void UpdateHistory(GlobalDofVector dofValues, std::vector<DofType> dofs);

private:
    SimpleAssembler& mAssembler;
    NodalValueMerger mMerger;

    using GradientPair = std::pair<Group<CellInterface>, CellInterface::VectorFunction>;
    using Hessian0Pair = std::pair<Group<CellInterface>, CellInterface::MatrixFunction>;
    using UpdatePair = std::pair<Group<CellInterface>, CellInterface::VoidFunction>;

    std::vector<GradientPair> mGradientFunctions;
    std::vector<Hessian0Pair> mHessian0Functions;
    std::vector<UpdatePair> mUpdateFunctions;
};
} /* NuTo */
