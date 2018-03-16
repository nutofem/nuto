#pragma once
#include "nuto/mechanics/cell/CellInterface.h"
#include "nuto/mechanics/cell/SimpleAssembler.h"
#include "nuto/base/Group.h"
#include "nuto/mechanics/dofs/GlobalDofVector.h"
#include "nuto/mechanics/dofs/GlobalDofMatrixSparse.h"

#include "nuto/mechanics/tools/NodalValueMerger.h"

namespace NuTo
{

//! Equation system that contains R(u, u') + M u'' = 0 with
//! R = Gradient
//! dR/du = Hessian0
//! dR/du' = Hessian1
//! M = Hessian2
class TimeDependentProblem
{
public:
    // let us call these "time dependent functions"
    using GradientFunction = std::function<DofVector<double>(const CellIpData&, double t, double dt)>;
    using HessianFunction = std::function<DofMatrix<double>(const CellIpData&, double t, double dt)>;
    using UpdateFunction = std::function<void(const CellIpData&, double t, double dt)>;

    TimeDependentProblem(MeshFem* rMesh);

    GlobalDofVector RenumberDofs(Constraint::Constraints constraints, std::vector<DofType> dofTypes,
                                 GlobalDofVector oldDofValues);

    void AddGradientFunction(Group<CellInterface> group, GradientFunction f);
    void AddHessian0Function(Group<CellInterface> group, HessianFunction f);
    void AddUpdateFunction(Group<CellInterface> group, UpdateFunction f);

    GlobalDofVector Gradient(const GlobalDofVector& dofValues, std::vector<DofType> dofs, double t, double dt);
    GlobalDofMatrixSparse Hessian0(const GlobalDofVector& dofValues, std::vector<DofType> dofs, double t, double dt);

    void UpdateHistory(const GlobalDofVector& dofValues, std::vector<DofType> dofs, double t, double dt);

private:
    SimpleAssembler mAssembler;
    NodalValueMerger mMerger;

    using GradientPair = std::pair<Group<CellInterface>, GradientFunction>;
    using Hessian0Pair = std::pair<Group<CellInterface>, HessianFunction>;
    using UpdatePair = std::pair<Group<CellInterface>, UpdateFunction>;

    std::vector<GradientPair> mGradientFunctions;
    std::vector<Hessian0Pair> mHessian0Functions;
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
} /* NuTo */
