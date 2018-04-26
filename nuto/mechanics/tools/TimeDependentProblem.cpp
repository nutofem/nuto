#include "nuto/mechanics/tools/TimeDependentProblem.h"
#include "nuto/mechanics/dofs/DofNumbering.h"

namespace NuTo
{

template <unsigned int TNumTimeDer>
TimeDependentProblem<TNumTimeDer>::TimeDependentProblem(MeshFem* rMesh)
    : mMerger(rMesh)
{
}

template <unsigned int TNumTimeDer>
DofVector<double> TimeDependentProblem<TNumTimeDer>::RenumberDofs(Constraint::Constraints constraints,
                                                                  std::vector<DofType> dofTypes,
                                                                  DofVector<double> oldDofValues)
{
    DofInfo dofInfos;

    DofVector<double> renumberedValues;

    for (auto dofType : dofTypes)
    {
        if (oldDofValues[dofType].rows() != 0) // initial case: No merge of an uninitialized vector
            mMerger.Merge(oldDofValues, {dofType});

        auto dofInfo = DofNumbering::Build(mMerger.Nodes(dofType), dofType, constraints);
        dofInfos.Merge(dofType, dofInfo);

        renumberedValues[dofType].resize(dofInfo.numIndependentDofs[dofType] + dofInfo.numDependentDofs[dofType]);

        mMerger.Extract(&renumberedValues, {dofType});
    }
    mAssembler.SetDofInfo(dofInfos);
    return renumberedValues;
}

template <unsigned int TNumTimeDer>
void TimeDependentProblem<TNumTimeDer>::AddGradientFunction(Group<CellInterface> group, GradientFunction f)
{
    mGradientFunctions.push_back({group, f});
}


template <unsigned int TNumTimeDer>
void TimeDependentProblem<TNumTimeDer>::AddUpdateFunction(Group<CellInterface> group, UpdateFunction f)
{
    mUpdateFunctions.push_back({group, f});
}


template <unsigned int TNumTimeDer>
DofVector<double> TimeDependentProblem<TNumTimeDer>::Gradient(const DofVector<double>& dofValues,
                                                              std::vector<DofType> dofs, double t, double dt)
{
    mMerger.Merge(dofValues, dofs);
    DofVector<double> gradient;
    for (auto& gradientFunction : mGradientFunctions)
        gradient += mAssembler.BuildVector(gradientFunction.first, dofs,
                                           Apply<CellInterface::VectorFunction>(gradientFunction.second, t, dt));
    return gradient;
}


template <unsigned int TNumTimeDer>
void TimeDependentProblem<TNumTimeDer>::UpdateHistory(const DofVector<double>& dofValues, std::vector<DofType> dofs,
                                                      double t, double dt)
{
    mMerger.Merge(dofValues, dofs);
    for (auto& updateFunction : mUpdateFunctions)
        for (auto& cell : updateFunction.first)
            cell.Apply(Apply<CellInterface::VoidFunction>(updateFunction.second, t, dt));
}


template class TimeDependentProblem<0>;
template class TimeDependentProblem<1>;
template class TimeDependentProblem<2>;
} // namespace NuTo
