#include "ConstraintLinearNodeWaterVolumeFraction.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"


NuTo::ConstraintLinearNodeWaterVolumeFraction::ConstraintLinearNodeWaterVolumeFraction(const NodeBase* rNode, double rValue)
    : ConstraintNode(rNode),
      ConstraintLinear(),
      mRHS{rValue}
{}

//! @brief adds the constraint equations to the matrix
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void NuTo::ConstraintLinearNodeWaterVolumeFraction::AddToConstraintMatrix(int& curConstraintEquation,
                                   NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix) const
{
    // add constraint to constrain matrix
    if (mNode->GetNumWaterVolumeFraction()!=1)
    {
        throw MechanicsException("[NuTo::ConstraintLinearNodeWaterVolumeFraction::AddToConstraintMatrix] Node does not have a water volume fraction component or has more than one water volume fraction component.");
    }
    rConstraintMatrix.AddValue(curConstraintEquation,mNode->GetDofWaterVolumeFraction(), 1);

    // increase constraint equation number
    curConstraintEquation++;
}


//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintLinearNodeWaterVolumeFraction::GetNumLinearConstraints() const
{
    return 1;
}

//!@brief returns the right hand side of the constraint equations
//!@return rRHS
double NuTo::ConstraintLinearNodeWaterVolumeFraction::GetRHS() const
{
    return mRHS;
}

//!@brief writes for the current constraint equation(s) the rhs into the vector
// (in case of more than one equation per constraint, curConstraintEquation is increased based on the number of constraint equations per constraint)
//! @param rCurConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void NuTo::ConstraintLinearNodeWaterVolumeFraction::GetRHS(int& rCurConstraintEquation, NuTo::FullVector<double,Eigen::Dynamic>& rRHS) const
{
    // add constraint to constrain matrix
    if (mNode->GetNumWaterVolumeFraction()!=1)
    {
        throw MechanicsException("[NuTo::ConstraintLinearNodeWaterVolumeFraction::GetRHS] Node does not have water volume fraction or has more than one water volume fraction component.");
    }
    // set right hand side value
    rRHS(rCurConstraintEquation) = mRHS;

    // increase constraint equation number
    rCurConstraintEquation++;
}


//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void NuTo::ConstraintLinearNodeWaterVolumeFraction::Info(unsigned short rVerboseLevel) const
{
    throw MechanicsException("[NuTo::ConstraintLinearNodeWaterVolumeFraction::Info] to be implemented.");
}
