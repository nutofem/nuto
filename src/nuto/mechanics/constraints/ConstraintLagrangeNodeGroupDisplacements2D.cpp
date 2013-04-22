// $Id: ConstraintLagrangeNodeDisplacements2D.cpp -1   $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include <iostream>

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/constraints/ConstraintLagrangeNodeGroupDisplacements2D.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRVector2Symmetric.h"

NuTo::ConstraintLagrangeNodeGroupDisplacements2D::ConstraintLagrangeNodeGroupDisplacements2D(const Group<NodeBase>* rGroup, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDirection, NuTo::Constraint::eEquationSign rEquationSign, double rRHS) :
        ConstraintNodeGroup(rGroup), ConstraintLagrange(rEquationSign)
{
    if (rDirection.GetNumColumns()!=1 || rDirection.GetNumRows()!=2)
        throw MechanicsException("[NuTo::ConstraintLagrangeNodeGroupDisplacements2D::ConstraintLagrangeNodeGroupDisplacements2D] Dimension of the direction matrix must be equal to the dimension of the structure.");

    memcpy(mDirection,rDirection.data(),2*sizeof(double));
    //normalize the direction
    double norm = sqrt(mDirection[0]*mDirection[0]+mDirection[1]*mDirection[1]);
    if (norm < 1e-14)
    {
        throw MechanicsException("[NuTo::ConstraintLagrangeNodeGroupDisplacements2D::ConstraintLagrangeNodeGroupDisplacements2D] direction vector has zero length.");
    }
    double invNorm = 1./norm;
    mDirection[0]*=invNorm;
    mDirection[1]*=invNorm;
    mRHS = rRHS;

    mLagrangeValue.resize(mGroup->GetNumMembers(),0);
    mLagrangeDOF.resize(mGroup->GetNumMembers());
}

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintLagrangeNodeGroupDisplacements2D::GetNumLagrangeMultipliers()const
{
    return mLagrangeValue.size();
}

//! @brief returns the Lagrange Multiplier
//! first col Lagrange, second column slack variables
void NuTo::ConstraintLagrangeNodeGroupDisplacements2D::GetLagrangeMultiplier(FullVector<double,Eigen::Dynamic>& rLagrangeMultiplier)const
{
    rLagrangeMultiplier.Resize(mGroup->GetNumMembers());
    for (unsigned int count=0; count<mLagrangeValue.size(); count++)
            rLagrangeMultiplier(count,0) = mLagrangeValue[count];
}

//! @brief returns the Lagrange Multiplier dofs
//! first col Lagrangedofs
void NuTo::ConstraintLagrangeNodeGroupDisplacements2D::GetDofsLagrangeMultiplier(FullVector<int,Eigen::Dynamic>& rLagrangeMultiplier)const
{
    rLagrangeMultiplier.Resize(mGroup->GetNumMembers());
    for (unsigned int count=0; count<mLagrangeDOF.size(); count++)
        rLagrangeMultiplier(count,0) = mLagrangeDOF[count];
}

//!@brief sets/modifies the right hand side of the constraint equations
//!@param rRHS new right hand side
void NuTo::ConstraintLagrangeNodeGroupDisplacements2D::SetRHS(double rRHS)
{
    mRHS = rRHS;
}

//! @brief sets the global dofs
//! @param rDOF current maximum DOF, this variable is increased within the routine
void NuTo::ConstraintLagrangeNodeGroupDisplacements2D::SetGlobalDofs(int& rDOF)
{
    assert((int)mLagrangeDOF.size()==mGroup->GetNumMembers());

    for (unsigned int count=0; count<mLagrangeDOF.size(); count++)
    {
        mLagrangeDOF[count] = rDOF++;
    }
}

//! @brief write dof values to constraints (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::ConstraintLagrangeNodeGroupDisplacements2D::SetGlobalDofValues(const FullVector<double,Eigen::Dynamic>& rActiveDofValues, const FullVector<double,Eigen::Dynamic>& rDependentDofValues)
{
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);

    for (unsigned int count=0; count<mLagrangeDOF.size(); count++)
    {
        int dof = this->mLagrangeDOF[count];
        double value;
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            value = rDependentDofValues(dof,0);
        }
        else
        {
            value = rActiveDofValues(dof,0);
        }
        this->mLagrangeValue[count] = value;
    }
}

//! @brief extract dof values from the constraints (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::ConstraintLagrangeNodeGroupDisplacements2D::GetGlobalDofValues(FullVector<double,Eigen::Dynamic>& rActiveDofValues, FullVector<double,Eigen::Dynamic>& rDependentDofValues) const
{
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);

    for (unsigned int count=0; count<mLagrangeDOF.size(); count++)
    {
        int dof = this->mLagrangeDOF[count];
        double value = this->mLagrangeValue[count];
        if (dof >= rActiveDofValues.GetNumRows())
        {
            dof -= rActiveDofValues.GetNumRows();
            assert(dof < rDependentDofValues.GetNumRows());
            rDependentDofValues(dof,0) = value;
        }
        else
        {
            rActiveDofValues(dof,0) = value;
        }
    }
}

//! @brief renumber the global dofs according to predefined ordering
//! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
void NuTo::ConstraintLagrangeNodeGroupDisplacements2D::RenumberGlobalDofs(const std::vector<int>& rMappingInitialToNewOrdering)
{
    for (unsigned int count=0; count<mLagrangeDOF.size(); count++)
    {
        mLagrangeDOF[count] = rMappingInitialToNewOrdering[mLagrangeDOF[count]];
    }
}

//! @brief cast to Lagrange constraint - the corresponding dofs are eliminated in the global system
NuTo::ConstraintLagrange* NuTo::ConstraintLagrangeNodeGroupDisplacements2D::AsConstraintLagrange()
{
    return this;
}

//! @brief cast to Lagrange constraint - the corresponding dofs are eliminated in the global system
const NuTo::ConstraintLagrange* NuTo::ConstraintLagrangeNodeGroupDisplacements2D::AsConstraintLagrange()const
{
    return this;
}

//! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the stiffness matrix (Hessian) of the total potential
//! @param rResult ... coefficient matrix
//! @param rGlobalDofs ... row and column numbers in global system
void NuTo::ConstraintLagrangeNodeGroupDisplacements2D::CalculateCoefficientMatrix_0(NuTo::SparseMatrixCSRVector2Symmetric<double>& rResult,
        std::vector<int>& rGlobalDofs)const
{
    //Lagrange mult + disp Dofs of all nodes
    int dof(3*mLagrangeDOF.size());
    rResult.Resize(dof,dof);
    rGlobalDofs.resize(dof,1);
    int curNodeEntry(0);
    int theNode(0);
    for (Group<NodeBase>::const_iterator itNode=mGroup->begin(); itNode!=mGroup->end();itNode++, curNodeEntry+=3, theNode++)
    {
        double disp[2];
        int dofx(itNode->second->GetDofDisplacement(0));
        int dofy(itNode->second->GetDofDisplacement(1));
        itNode->second->GetDisplacements2D(disp);

        rGlobalDofs[curNodeEntry] = mLagrangeDOF[theNode];
        rGlobalDofs[curNodeEntry+1] = dofx;
        rGlobalDofs[curNodeEntry+2] = dofy;

        switch(mEquationSign)
        {
        case NuTo::Constraint::EQUAL:
            //derivative with respect to ux and lambda
            rResult.AddValue(curNodeEntry,curNodeEntry+1,mDirection[0]);
            //derivative with respect to uy and lambda
            rResult.AddValue(curNodeEntry,curNodeEntry+2,mDirection[1]);
            //derivative with respect to ux^2
            rResult.AddValue(curNodeEntry+1,curNodeEntry+1,mPenalty*mDirection[0]*mDirection[0]);
            //derivative with respect to ux and uy
            rResult.AddValue(curNodeEntry+1,curNodeEntry+2,mPenalty*mDirection[0]*mDirection[1]);
            //derivative with respect to uy^2
            rResult.AddValue(curNodeEntry+2,curNodeEntry+2,mPenalty*mDirection[1]*mDirection[1]);
        break;
        case NuTo::Constraint::SMALLER:
            if (disp[0]*mDirection[0]+disp[1]*mDirection[1]-mRHS>-mLagrangeValue[theNode]/mPenalty)
            {
                //derivative with respect to ux and lambda
                rResult.AddValue(curNodeEntry,curNodeEntry+1,mDirection[0]);
                //derivative with respect to uy and lambda
                rResult.AddValue(curNodeEntry,curNodeEntry+2,mDirection[1]);
                //derivative with respect to ux^2
                rResult.AddValue(curNodeEntry+1,curNodeEntry+1,mPenalty*mDirection[0]*mDirection[0]);
                //derivative with respect to ux and uy
                rResult.AddValue(curNodeEntry+1,curNodeEntry+2,mPenalty*mDirection[0]*mDirection[1]);
                //derivative with respect to uy^2
                rResult.AddValue(curNodeEntry+2,curNodeEntry+2,mPenalty*mDirection[1]*mDirection[1]);
            }
            else
            {
                //derivative with respect to lambda^2
                rResult.AddValue(curNodeEntry,curNodeEntry,-1./mPenalty);
            }
        break;
        case NuTo::Constraint::GREATER:
            if (mRHS-disp[0]*mDirection[0]-disp[1]*mDirection[1]>-mLagrangeValue[theNode]/mPenalty)
            {
                //derivative with respect to ux and lambda
                rResult.AddValue(curNodeEntry,curNodeEntry+1,-mDirection[0]);
                //derivative with respect to uy and lambda
                rResult.AddValue(curNodeEntry,curNodeEntry+2,-mDirection[1]);
                //derivative with respect to ux^2
                rResult.AddValue(curNodeEntry+1,curNodeEntry+1,mPenalty*mDirection[0]*mDirection[0]);
                //derivative with respect to ux and uy
                rResult.AddValue(curNodeEntry+1,curNodeEntry+2,mPenalty*mDirection[0]*mDirection[1]);
                //derivative with respect to uy^2
                rResult.AddValue(curNodeEntry+2,curNodeEntry+2,mPenalty*mDirection[1]*mDirection[1]);
            }
            else
            {
                //derivative with respect to lambda^2
                rResult.AddValue(curNodeEntry,curNodeEntry,-1./mPenalty);
            }
        break;
        default:
            throw MechanicsException("[NuTo::ConstraintLagrangeNodeGroupDisplacements1D::CalculateCoefficientMatrix_0] equation sign should be either EQUAL, SMALLER or GREAER.");
        }
    }
}

//! @brief calculates the gradient of the internal potential
//! for a mechanical problem, this corresponds to the internal force vector
void NuTo::ConstraintLagrangeNodeGroupDisplacements2D::CalculateGradientInternalPotential(NuTo::FullVector<double,Eigen::Dynamic>& rResult,
        std::vector<int>& rGlobalDofs)const
{
    int dof(3*mLagrangeDOF.size());
    rResult.Resize(dof);
    rGlobalDofs.resize(dof);
    int curNodeEntry(0);
    int theNode(0);
    for (Group<NodeBase>::const_iterator itNode=mGroup->begin(); itNode!=mGroup->end();itNode++, curNodeEntry+=3, theNode++)
    {
        double disp[2];
        int dofx(itNode->second->GetDofDisplacement(0));
        int dofy(itNode->second->GetDofDisplacement(1));
        itNode->second->GetDisplacements2D(disp);

        rGlobalDofs[curNodeEntry] = mLagrangeDOF[theNode];
        rGlobalDofs[curNodeEntry+1] = dofx;
        rGlobalDofs[curNodeEntry+2] = dofy;

        switch(mEquationSign)
        {
        case NuTo::Constraint::EQUAL:
            //derivative with respect to lambda
            rResult(curNodeEntry,0)=disp[0]*mDirection[0]+disp[1]*mDirection[1]-mRHS;
            //derivative with respect to ux
            rResult(curNodeEntry+1,0)=mLagrangeValue[theNode]*mDirection[0]+mPenalty*mDirection[0]*(disp[0]*mDirection[0]+disp[1]*mDirection[1]-mRHS);
            //derivative with respect to uy
            rResult(curNodeEntry+2,0)=mLagrangeValue[theNode]*mDirection[1]+mPenalty*mDirection[1]*(disp[0]*mDirection[0]+disp[1]*mDirection[1]-mRHS);
        break;
        case NuTo::Constraint::SMALLER:
            if (disp[0]*mDirection[0]+disp[1]*mDirection[1]-mRHS>-mLagrangeValue[theNode]/mPenalty)
            {
                //derivative with respect to lambda
                rResult(curNodeEntry,0)=disp[0]*mDirection[0]+disp[1]*mDirection[1]-mRHS;
                //derivative with respect to ux
                rResult(curNodeEntry+1,0)=mLagrangeValue[theNode]*mDirection[0]+mPenalty*mDirection[0]*(disp[0]*mDirection[0]+disp[1]*mDirection[1]-mRHS);
                //derivative with respect to uy
                rResult(curNodeEntry+2,0)=mLagrangeValue[theNode]*mDirection[1]+mPenalty*mDirection[1]*(disp[0]*mDirection[0]+disp[1]*mDirection[1]-mRHS);
            }
            else
            {
                //derivative with respect to lambda
                rResult(curNodeEntry,0)=-mLagrangeValue[theNode]/mPenalty;
            }
        break;
        case NuTo::Constraint::GREATER:
            if (mRHS-disp[0]*mDirection[0]-disp[1]*mDirection[1]>-mLagrangeValue[theNode]/mPenalty)
            {
                //derivative with respect to lambda
                rResult(curNodeEntry,0)=mRHS - disp[0]*mDirection[0] - disp[1]*mDirection[1];
                //derivative with respect to ux
                rResult(curNodeEntry+1,0)=-mLagrangeValue[theNode]*mDirection[0]+mPenalty*mDirection[0]*(disp[0]*mDirection[0]+disp[1]*mDirection[1]-mRHS);
                //derivative with respect to uy
                rResult(curNodeEntry+2,0)=-mLagrangeValue[theNode]*mDirection[1]+mPenalty*mDirection[1]*(disp[0]*mDirection[0]+disp[1]*mDirection[1]-mRHS);
            }
            else
            {
                //derivative with respect to lambda
                rResult(curNodeEntry,0)=-mLagrangeValue[theNode]/mPenalty;
            }
        break;
        default:
            throw MechanicsException("[NuTo::ConstraintLagrangeNodeGroupDisplacements1D::CalculateCoefficientMatrix_0] equation sign should be either EQUAL, SMALLER or GREAER.");
        }
    }
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ConstraintLagrangeNodeGroupDisplacements2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrangeNodeGroupDisplacements2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrangeNodeGroupDisplacements2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrangeNodeGroupDisplacements2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrangeNodeGroupDisplacements2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrangeNodeGroupDisplacements2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintLagrangeNodeGroupDisplacements2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLagrangeNodeGroupDisplacements2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintNodeGroup)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintLagrange)
       & BOOST_SERIALIZATION_NVP(mRHS)
       & BOOST_SERIALIZATION_NVP(mDirection)
       & BOOST_SERIALIZATION_NVP(mLagrangeValue)
       & BOOST_SERIALIZATION_NVP(mLagrangeDOF);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLagrangeNodeGroupDisplacements2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintLagrangeNodeGroupDisplacements2D)
#endif // ENABLE_SERIALIZATION
