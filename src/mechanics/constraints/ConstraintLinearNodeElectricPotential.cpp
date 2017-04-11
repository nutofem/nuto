#include "ConstraintLinearNodeElectricPotential.h"
#include "mechanics/MechanicsException.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "math/SparseMatrix.h"

NuTo::ConstraintLinearNodeElectricPotential::ConstraintLinearNodeElectricPotential(const NodeBase* rNode, double rValue)
    :   ConstraintLinear(),
		ConstraintNode(rNode),
        mRHS{rValue}
{}

void  NuTo::ConstraintLinearNodeElectricPotential::AddToConstraintMatrix(
        int& curConstraintEquation, NuTo::SparseMatrix<double>& rConstraintMatrix) const
{
    // add constraint to constrain matrix
    if (mNode->GetNum(Node::eDof::ELECTRICPOTENTIAL)!=1)
    {
        throw MechanicsException(__PRETTY_FUNCTION__,
                "Node does not have a ElectricPotential component or has more than one ElectricPotential component.");
    }
    rConstraintMatrix.AddValue(curConstraintEquation,mNode->GetDof(Node::eDof::ELECTRICPOTENTIAL), 1);

    // increase constraint equation number
    curConstraintEquation++;
}

int NuTo::ConstraintLinearNodeElectricPotential::GetNumLinearConstraints() const
{
    return 1;
}

double NuTo::ConstraintLinearNodeElectricPotential::GetRHS() const
{
    return mRHS;
}

void NuTo::ConstraintLinearNodeElectricPotential::SetRHS(double rRHS)
{
    mRHS = rRHS;
}

void NuTo::ConstraintLinearNodeElectricPotential::GetRHS(int& rCurConstraintEquation,
        Eigen::VectorXd& rRHS) const
{
    // add constraint to constrain matrix
    if (mNode->GetNum(Node::eDof::ELECTRICPOTENTIAL) != 1)
    {
        throw MechanicsException(__PRETTY_FUNCTION__,
            "Node does not have a ElectricPotential component or has more than one ElectricPotential component.");
    }
    // set right hand side value
    rRHS(rCurConstraintEquation) = mRHS;

    // increase constraint equation number
    rCurConstraintEquation++;
}

void NuTo::ConstraintLinearNodeElectricPotential::Info(unsigned short rVerboseLevel) const
{
    throw MechanicsException(__PRETTY_FUNCTION__, " to be implemented.");
}

NuTo::Node::eDof NuTo::ConstraintLinearNodeElectricPotential::GetDofType() const
{
    return Node::eDof::ELECTRICPOTENTIAL;
}
