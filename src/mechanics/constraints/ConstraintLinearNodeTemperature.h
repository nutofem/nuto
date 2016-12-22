#pragma once

#include "mechanics/constraints/ConstraintLinear.h"
#include "mechanics/constraints/ConstraintNode.h"

namespace NuTo
{
class ConstraintLinearNodeTemperature : public ConstraintLinear, public ConstraintNode
{
public:
    ConstraintLinearNodeTemperature(const NodeBase* rNode, double rValue);

    //! @brief Adds the constraint equations to the matrix.
    //! @param curConstraintEquation (is incremented during the function call)
    //! @param rConstraintMatrix (the first row where a constraint equation is
    //!     added is given by curConstraintEquation)
    virtual void AddToConstraintMatrix(int& curConstraintEquation,
            NuTo::SparseMatrix<double>& rConstraintMatrix) const override;

    //! @brief Returns the number of constraint equations.
    //! @return Number of constraints
    virtual int GetNumLinearConstraints() const override;

    //!@brief Returns the right hand side of the constraint equations.
    virtual double GetRHS() const override;

    //!@brief Sets/modifies the right hand side of the constraint equation.
    //!@param rRHS New right hand side
    virtual void SetRHS(double rRHS) override;

    //!@brief Writes for the current constraint equation(s) the rhs into the
    //!     vector. (in case of more than one equation per constraint,
    //!     curConstraintEquation is increased based on the number of constraint
    //!     equations per constraint)
    //! @param rCurConstraintEquation (is incremented during the function call)
    //! @param rConstraintMatrix (the first row where a constraint equation is
    //!     added is given by curConstraintEquation)
    virtual void GetRHS(int& rCurConstraintEquation,
            NuTo::FullVector<double,Eigen::Dynamic>& rRHS) const override;

    //! @brief Print information about the object.
    //! @param rVerboseLevel Verbosity of the information
    virtual void Info(unsigned short rVerboseLevel) const override;

    //! @brief Determines the dof type affected by the constraint.
    //! @return DOF type
    Node::eDof GetDofType() const override;

protected:
    double mRHS = 0.0;
};
}
