# include "ConstraintLinearNodeRelativeHumidity.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"




                            NuTo::ConstraintLinearNodeRelativeHumidity::ConstraintLinearNodeRelativeHumidity            (const NodeBase* rNode, double rValue)
    :   ConstraintNode(rNode),
        ConstraintLinear(),
        mRHS{rValue}
{}

//! @brief adds the constraint equations to the matrix
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void                        NuTo::ConstraintLinearNodeRelativeHumidity::AddToConstraintMatrix                           (int& curConstraintEquation, NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix) const
{
    // add constraint to constrain matrix
    if (mNode->GetNumRelativeHumidity()!=1)
    {
        throw MechanicsException("[NuTo::ConstraintLinearNodeRelativeHumidity::AddToConstraintMatrix] Node does not have a relative humidity component or has more than one relative humidity component.");
    }
    rConstraintMatrix.AddValue(curConstraintEquation,mNode->GetDofRelativeHumidity(), 1);

    // increase constraint equation number
    curConstraintEquation++;
}

//! @brief returns the number of constraint equations
//! @return number of constraints
int                         NuTo::ConstraintLinearNodeRelativeHumidity::GetNumLinearConstraints                         () const
{
    return 1;
}

//!@brief returns the right hand side of the constraint equations
//!@return rRHS
double                      NuTo::ConstraintLinearNodeRelativeHumidity::GetRHS                                          () const
{
    return mRHS;
}

//!@brief writes for the current constraint equation(s) the rhs into the vector
// (in case of more than one equation per constraint, curConstraintEquation is increased based on the number of constraint equations per constraint)
//! @param rCurConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void                        NuTo::ConstraintLinearNodeRelativeHumidity::GetRHS                                          (int& rCurConstraintEquation, NuTo::FullVector<double,Eigen::Dynamic>& rRHS) const
{
    // add constraint to constrain matrix
    if (mNode->GetNumRelativeHumidity()!=1)
    {
        throw MechanicsException("[NuTo::ConstraintLinearNodeRelativeHumidity::GetRHS] Node does not have a relative humidity component or has more than one relative humidity component.");
    }
    // set right hand side value
    rRHS(rCurConstraintEquation) = mRHS;

    // increase constraint equation number
    rCurConstraintEquation++;
}

//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void                        NuTo::ConstraintLinearNodeRelativeHumidity::Info                                            (unsigned short rVerboseLevel) const
{
    throw MechanicsException("[NuTo::ConstraintLinearNodeRelativeHumidity::Info] to be implemented.");
}
