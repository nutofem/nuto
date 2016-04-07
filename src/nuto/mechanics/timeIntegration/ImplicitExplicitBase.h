/*
 * ImplicitExplicitBase.h
 *
 *  Created on: 17 Mar 2016
 *      Author: ttitsche
 */

#pragma once

#include "nuto/mechanics/timeIntegration/TimeIntegrationBase.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"

namespace NuTo
{

//! @brief Base class for implicit/explicit time integration schemes like ImplEx and CycleJump
class ImplicitExplicitBase : public TimeIntegrationBase
{
public:

    ImplicitExplicitBase(StructureBase* rStructure);

    virtual ~ImplicitExplicitBase() = default;

    //! @brief perform the time integration
    //! @param rTimeDelta ... length of the simulation
    NuTo::Error::eError Solve(double rTimeDelta) override;

    //! @brief returns true, if the method is only conditionally stable (for unconditional stable, this is false)
    bool HasCriticalTimeStep() const override {
        return false;
    }

    //! @brief calculate the critical time step for explicit routines
    //! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
    double CalculateCriticalTimeStep()const override;

    //! @brief ... Adds a calculation step to each timestep
    //! param rActiveDofs ... active Dofs of the calculation step
    void AddCalculationStep(const std::set<Node::eDof> &rActiveDofs);

    //! @brief ... Adds a dof type with a constant hessian0 matrix
    //! param rDofWithConstantHessian ... dof type
    void AddDofWithConstantHessian(Node::eDof rDofWithConstantHessian);

private:

    //! @brief Stores wich Dofs are active in which calculation step
    std::vector<std::set<Node::eDof>> mStepActiveDofs;

    std::set<Node::eDof> mDofsWithConstantHessian;

    void FactorizeConstantHessians(std::map<Node::eDof, SparseDirectSolverMUMPS>& rPreFactorizedHessians);
};

} /* namespace NuTo */

