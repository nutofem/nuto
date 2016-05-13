#include "ConstraintLinearNodeTemperature.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/math/SparseMatrix.h"

NuTo::ConstraintLinearNodeTemperature::ConstraintLinearNodeTemperature(const NodeBase* rNode, double rValue)
    :   ConstraintLinear(),
		ConstraintNode(rNode),
        mRHS{rValue}
{}

void  NuTo::ConstraintLinearNodeTemperature::AddToConstraintMatrix(
        int& curConstraintEquation, NuTo::SparseMatrix<double>& rConstraintMatrix) const
{
    // add constraint to constrain matrix
    if (mNode->GetNumTemperature()!=1)
    {
        throw MechanicsException(__PRETTY_FUNCTION__,
                "Node does not have a temperature component or has more than one temperature component.");
    }
    rConstraintMatrix.AddValue(curConstraintEquation,mNode->GetDofTemperature(), 1);

    // increase constraint equation number
    curConstraintEquation++;
}

int NuTo::ConstraintLinearNodeTemperature::GetNumLinearConstraints() const
{
    return 1;
}

double NuTo::ConstraintLinearNodeTemperature::GetRHS() const
{
    return mRHS;
}

void NuTo::ConstraintLinearNodeTemperature::SetRHS(double rRHS)
{
    mRHS = rRHS;
}

void NuTo::ConstraintLinearNodeTemperature::GetRHS(int& rCurConstraintEquation,
        NuTo::FullVector<double,Eigen::Dynamic>& rRHS) const
{
    // add constraint to constrain matrix
    if (mNode->GetNumTemperature() != 1)
    {
        throw MechanicsException(__PRETTY_FUNCTION__,
            "Node does not have a temperature component or has more than one temperature component.");
    }
    // set right hand side value
    rRHS(rCurConstraintEquation) = mRHS;

    // increase constraint equation number
    rCurConstraintEquation++;
}

void NuTo::ConstraintLinearNodeTemperature::Info(unsigned short rVerboseLevel) const
{
    throw MechanicsException(__PRETTY_FUNCTION__, " to be implemented.");
}
