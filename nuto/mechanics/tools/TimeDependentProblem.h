#pragma once
#include "nuto/mechanics/cell/CellInterface.h"
#include "nuto/mechanics/cell/SimpleAssembler.h"
#include "nuto/base/Group.h"
#include "nuto/mechanics/dofs/DofVector.h"
#include "nuto/mechanics/dofs/DofMatrixSparse.h"

#include "nuto/mechanics/tools/NodalValueMerger.h"

#include <array>

namespace NuTo
{

template <typename TCellInterfaceFunction, typename TTimeDepFunction>
TCellInterfaceFunction Apply(TTimeDepFunction& f, double t, double dt)
{
    using namespace std::placeholders;
    return std::bind(f, _1, t, dt);
}

//! Equation system that contains R(u, u') + M u'' = 0 with
//! R = Gradient
//! dR/du = Hessian0
//! dR/du' = Hessian1
//! M = Hessian2

template <unsigned int TNumTimeDer = 0>
class TimeDependentProblem
{
public:
    // let us call these "time dependent functions"
    using GradientFunction = std::function<DofVector<double>(const CellIpData&, double t, double dt)>;
    using HessianFunction = std::function<DofMatrix<double>(const CellIpData&, double t, double dt)>;
    using UpdateFunction = std::function<void(const CellIpData&, double t, double dt)>;

    TimeDependentProblem(MeshFem* rMesh);

    DofVector<double> RenumberDofs(Constraint::Constraints constraints, std::vector<DofType> dofTypes,
                                   DofVector<double> oldDofValues);

    void AddGradientFunction(Group<CellInterface> group, GradientFunction f);
    template <int TTimeDer>
    void AddHessianFunction(Group<CellInterface> group, HessianFunction f);
    void AddUpdateFunction(Group<CellInterface> group, UpdateFunction f);

    DofVector<double> Gradient(std::vector<DofType> dofs, double t, double dt);
    DofVector<double> Gradient(const DofVector<double>& dofValues, std::vector<DofType> dofs, double t, double dt);

    // maybe use type traits to activate different versions depending on TNumTimeDer
    template <int TTimeDer>
    DofMatrixSparse<double> Hessian(std::vector<DofType> dofs, double t, double dt);
    template <int TTimeDer>
    DofMatrixSparse<double> Hessian(const DofVector<double>& dofValues, std::vector<DofType> dofs, double t, double dt);
    template <int TTimeDer>
    DofMatrixSparse<double> Hessian(const std::array<DofVector<double>, TNumTimeDer + 1>& dofValues,
                                    std::vector<DofType> dofs, double t, double dt);

    void UpdateHistory(const DofVector<double>& dofValues, std::vector<DofType> dofs, double t, double dt);

private:
    SimpleAssembler mAssembler;
    NodalValueMerger mMerger;

    using GradientPair = std::pair<Group<CellInterface>, GradientFunction>;
    using HessianPair = std::pair<Group<CellInterface>, HessianFunction>;
    using UpdatePair = std::pair<Group<CellInterface>, UpdateFunction>;

    std::vector<GradientPair> mGradientFunctions;
    std::array<std::vector<HessianPair>, TNumTimeDer + 1> mHessianFunctions;
    std::vector<UpdatePair> mUpdateFunctions;


    /*
     *
     *
     *
     *
     *
     * We better hide these two beauties in the basement of this class
     *
     *
     *
     *
     *
     *
     *
     *
     */

public:
    //! @brief binds a dt dependent integrand to a time dependent function
    template <typename TObject, typename TReturn>
    static auto Bind_dt(TObject& object, TReturn (TObject::*f)(const CellIpData&, double dt))
    {
        return [&object, f](const CellIpData& cellIpData, double, double dt) { return (object.*f)(cellIpData, dt); };
    }

    //! @brief binds a t dependent integrand to a time dependent function
    template <typename TObject, typename TReturn>
    static auto Bind_t(TObject& object, TReturn (TObject::*f)(const NuTo::CellIpData&, double t))
    {
        return [&object, f](const CellIpData& cellIpData, double t, double) { return (object.*f)(cellIpData, t); };
    }

    //! @brief binds nothing
    template <typename TObject, typename TReturn>
    static auto Bind(TObject& object, TReturn (TObject::*f)(const NuTo::CellIpData&))
    {
        return [&object, f](const CellIpData& cellIpData, double, double) { return (object.*f)(cellIpData); };
    }
};

template <unsigned int TNumTimeDer>
template <int TTimeDer>
void TimeDependentProblem<TNumTimeDer>::AddHessianFunction(Group<CellInterface> group,
                                                           TimeDependentProblem::HessianFunction f)
{
    static_assert(TTimeDer <= TNumTimeDer, "Hessian time derivative bigger than the problems time derivative.");
    mHessianFunctions[TTimeDer].push_back({group, f});
}


template <unsigned int TNumTimeDer>
template <int TTimeDer>
DofMatrixSparse<double> TimeDependentProblem<TNumTimeDer>::Hessian(std::vector<DofType> dofs, double t, double dt)
{
    static_assert(TTimeDer <= TNumTimeDer, "Hessian time derivative bigger than the problems time derivative.");
    DofMatrixSparse<double> hessian0;
    for (auto& hessianFunction : mHessianFunctions[TTimeDer])
        hessian0 += mAssembler.BuildMatrix(hessianFunction.first, dofs,
                                           Apply<CellInterface::MatrixFunction>(hessianFunction.second, t, dt));
    return hessian0;
}

template <unsigned int TNumTimeDer>
template <int TTimeDer>
DofMatrixSparse<double> TimeDependentProblem<TNumTimeDer>::Hessian(const DofVector<double>& dofValues,
                                                                   std::vector<DofType> dofs, double t, double dt)
{
    mMerger.Merge(dofValues, dofs);
    return Hessian<TTimeDer>(dofs, t, dt);
}

template <unsigned int TNumTimeDer>
template <int TTimeDer>
DofMatrixSparse<double>
TimeDependentProblem<TNumTimeDer>::Hessian(const std::array<DofVector<double>, TNumTimeDer + 1>& dofValues,
                                           std::vector<DofType> dofs, double t, double dt)
{
    for (size_t i = 0; i < dofValues.size(); ++i)
        mMerger.Merge(dofValues[i], dofs, i);
    return Hessian<TTimeDer>(dofs, t, dt);
}

} /* NuTo */
