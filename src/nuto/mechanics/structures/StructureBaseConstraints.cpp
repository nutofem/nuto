// $Id$

#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseMatrixCSRVector2Symmetric.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/constraints/ConstraintEnum.h"
#include "nuto/mechanics/constraints/ConstraintLagrangeNodeGroupDisplacements1D.h"
#include "nuto/mechanics/constraints/ConstraintLagrangeNodeGroupDisplacements2D.h"
#include "nuto/mechanics/constraints/ConstraintLinearDisplacementsPeriodic2D.h"
#include "nuto/mechanics/constraints/ConstraintLinearEquation.h"
#include "nuto/mechanics/constraints/ConstraintLinearNodeDisplacements1D.h"
#include "nuto/mechanics/constraints/ConstraintLinearNodeDisplacements2D.h"
#include "nuto/mechanics/constraints/ConstraintLinearNodeDisplacements3D.h"
#include "nuto/mechanics/constraints/ConstraintLinearNodeFineScaleDisplacements2D.h"
#include "nuto/mechanics/constraints/ConstraintLinearNodeGroupDisplacements1D.h"
#include "nuto/mechanics/constraints/ConstraintLinearNodeGroupDisplacements2D.h"
#include "nuto/mechanics/constraints/ConstraintLinearNodeGroupDisplacements3D.h"
#include "nuto/mechanics/constraints/ConstraintLinearNodeGroupFineScaleDisplacements2D.h"
#include "nuto/mechanics/constraints/ConstraintLinearNodeGroupRotations2D.h"
#include "nuto/mechanics/constraints/ConstraintLinearNodeRotations2D.h"
#include "nuto/mechanics/constraints/ConstraintNonlinear.h"

//! @brief adds a displacement constraint equation for a node group solved using Lagrange multiplier
//! @param rGroupId group id
//! @param rDirection direction of the constraint (in 2D a point with 2 entries, in 3D 3 entries, in 1D not used)
//! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
//! @return integer id to delete or modify the constraint
int NuTo::StructureBase::ConstraintLagrangeSetDisplacementNodeGroup(int rGroupId, const NuTo::FullMatrix<double>& rDirection, const std::string& rEquationSign, double rValue)
{
    this->mNodeNumberingRequired = true;
    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::ConstraintLagrangeSetDisplacementNodeGroup] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Nodes)
        throw MechanicsException("[NuTo::Structure::ConstraintLagrangeSetDisplacementNodeGroup] Group is not a node group.");
    Group<NodeBase> *nodeGroup = dynamic_cast<Group<NodeBase>*>(itGroup->second);
    assert(nodeGroup!=0);
    //convert equation sign to upper case
    std::string equationSign;
    std::transform(rEquationSign.begin(), rEquationSign.end(), std::back_inserter(equationSign), (int(*)(int)) toupper);
    Constraint::eEquationSign eEquationSign;
    if(equationSign == "EQUAL")
    {
        eEquationSign = NuTo::Constraint::EQUAL;
    }
    else if(equationSign == "GREATER")
    {
        eEquationSign = NuTo::Constraint::GREATER;
    }
    else if(equationSign == "SMALLER")
    {
        eEquationSign = NuTo::Constraint::SMALLER;
    }
    else
    {
        throw MechanicsException("[NuTo::StructureBase::ConstraintLagrangeSetDisplacementNodeGroup] equation sign is either equal, smaller or greater.");
    }

    return ConstraintLagrangeSetDisplacementNodeGroup(nodeGroup,rDirection, eEquationSign, rValue);

}
//! @brief adds a displacement constraint equation for a node group solved using Lagrange multiplier
//! @param rGroup group pointer
//! @param rDirection direction of the constraint (in 2D a point with 2 entries, in 3D 3 entries, in 1D not used)
//! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
//! @return integer id to delete or modify the constraint
int NuTo::StructureBase::ConstraintLagrangeSetDisplacementNodeGroup(Group<NodeBase>* rGroup, const NuTo::FullMatrix<double>& rDirection, NuTo::Constraint::eEquationSign rEquationSign, double rValue)
{
    this->mNodeNumberingRequired = true;
    //find unused integer id
    int id(0);
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(id);
    while (it!=mConstraintMap.end())
    {
        id++;
        it = mConstraintMap.find(id);
    }

    switch (mDimension)
    {
    case 1:
    mConstraintMap.insert(id, new NuTo::ConstraintLagrangeNodeGroupDisplacements1D(rGroup,rDirection,rEquationSign,rValue));
    break;
    case 2:
        mConstraintMap.insert(id, new NuTo::ConstraintLagrangeNodeGroupDisplacements2D(rGroup,rDirection,rEquationSign,rValue));
        break;
    //case 3:
    //mConstraintMap.insert(id, new NuTo::ConstraintLagrangeNodeGroupDisplacements3D(rGroup,rDirection,rEquationSign,rValue));
    //break;
    default:
        throw MechanicsException("[NuTo::StructureBase::ConstraintSetDisplacementNodeGroup] Incorrect dimension of the structure.");
    }
    return id;

}

//! @brief adds a displacement constraint equation for a node
//! @param rNode pointer to node
//! @param rDirection direction of the constraint (in 2D a point with 2 entries, in 3D 3 entries, in 1D not used)
//! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
//! @return integer id to delete or modify the constraint
int NuTo::StructureBase::ConstraintLinearSetDisplacementNode(NodeBase* rNode, const NuTo::FullMatrix<double>& rDirection, double rValue)
{
	this->mNodeNumberingRequired = true;
	//find unused integer id
    int id(0);
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(id);
    while (it!=mConstraintMap.end())
    {
        id++;
        it = mConstraintMap.find(id);
    }

    switch (mDimension)
    {
    case 1:
        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeDisplacements1D(rNode,rDirection(0,0),rValue));
        break;
    case 2:
        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeDisplacements2D(rNode,rDirection,rValue));
        break;
    case 3:
        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeDisplacements3D(rNode,rDirection,rValue));
        break;
    default:
        throw MechanicsException("[NuTo::StructureBase::ConstraintSetDisplacementNode] Incorrect dimension of the structure.");
    }
    return id;
}

//! @brief adds a rotation constraint equation for a node
//! @param rNode pointer to node
//! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
//! @return integer id to delete or modify the constraint
int NuTo::StructureBase::ConstraintLinearSetRotationNode(NodeBase* rNode, double rValue)
{
	this->mNodeNumberingRequired = true;
	//find unused integer id
    int id(0);
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(id);
    while (it!=mConstraintMap.end())
    {
        id++;
        it = mConstraintMap.find(id);
    }

    switch (mDimension)
    {
    case 1:
        throw MechanicsException("[NuTo::StructureBase::ConstraintLinearSetDisplacementNode] not implemented for 1D.");
        break;
    case 2:
        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeRotations2D(rNode,rValue));
        break;
    case 3:
        throw MechanicsException("[NuTo::StructureBase::ConstraintLinearSetDisplacementNode] not implemented for 3D.");
        break;
    default:
        throw MechanicsException("[NuTo::StructureBase::ConstraintSetDisplacementNode] Incorrect dimension of the structure.");
    }
    return id;
}

//! @brief adds a fine scale displacement constraint equation for a node
//! @param rNode pointer to node
//! @param rDirection direction of the constraint (in 2D a point with 2 entries, in 3D 3 entries, in 1D not used)
//! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
//! @return integer id to delete or modify the constraint
int NuTo::StructureBase::ConstraintLinearSetFineScaleDisplacementNode(NodeBase* rNode, const NuTo::FullMatrix<double>& rDirection, double rValue)
{
    this->mNodeNumberingRequired = true;
    //find unused integer id
    int id(0);
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(id);
    while (it!=mConstraintMap.end())
    {
        id++;
        it = mConstraintMap.find(id);
    }

    switch (mDimension)
    {
    case 1:
        throw MechanicsException("[NuTo::StructureBase::ConstraintSetfineScaleDisplacementNode] not implemented for 1D.");
        break;
    case 2:
        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeFineScaleDisplacements2D(rNode,rDirection,rValue));
        break;
    case 3:
        throw MechanicsException("[NuTo::StructureBase::ConstraintSetfineScaleDisplacementNode] not implemented for 3D.");
        break;
    default:
        throw MechanicsException("[NuTo::StructureBase::ConstraintSetDisplacementNode] Incorrect dimension of the structure.");
    }
    return id;
}

//! @brief adds a displacement constraint equation for a node
//! @param rNode identifier for node
//! @param rComponent e.g. the first (count from zero) displacement component
//! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
int  NuTo::StructureBase::ConstraintLinearSetDisplacementNode(int rIdent, const NuTo::FullMatrix<double>& rDirection, double rValue)
{
	this->mNodeNumberingRequired = true;
    NodeBase* nodePtr;
    try
    {
        nodePtr = NodeGetNodePtr(rIdent);
    }
    catch (NuTo::MechanicsException &e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstraintSetDisplacementNode] Node with the given identifier could not be found.");
        throw e;
    }
    catch (...)
    {
        throw MechanicsException("[NuTo::StructureBase::ConstraintSetDisplacementNode] Node with the given identifier could not be found.");
    }

    return ConstraintLinearSetDisplacementNode(nodePtr,rDirection, rValue);
}

//! @brief adds a fine scale displacement constraint equation for a node
//! @param rNode identifier for node
//! @param rComponent e.g. the first (count from zero) displacement component
//! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
int  NuTo::StructureBase::ConstraintLinearSetFineScaleDisplacementNode(int rIdent, const NuTo::FullMatrix<double>& rDirection, double rValue)
{
    this->mNodeNumberingRequired = true;
    NodeBase* nodePtr;
    try
    {
        nodePtr = NodeGetNodePtr(rIdent);
    }
    catch (NuTo::MechanicsException &e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstraintSetFineScaleDisplacementNode] Node with the given identifier could not be found.");
        throw e;
    }
    catch (...)
    {
        throw MechanicsException("[NuTo::StructureBase::ConstraintSetFineScaleDisplacementNode] Node with the given identifier could not be found.");
    }

    return ConstraintLinearSetFineScaleDisplacementNode(nodePtr,rDirection, rValue);
}

//! @brief adds a rotation constraint equation for a node
//! @param rNode identifier for node
//! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
int  NuTo::StructureBase::ConstraintLinearSetRotationNode(int rIdent, double rValue)
{
	this->mNodeNumberingRequired = true;
    NodeBase* nodePtr;
    try
    {
        nodePtr = NodeGetNodePtr(rIdent);
    }
    catch (NuTo::MechanicsException &e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstraintLinearSetRotationNode] Node with the given identifier could not be found.");
        throw e;
    }
    catch (...)
    {
        throw MechanicsException("[NuTo::StructureBase::ConstraintLinearSetRotationNode] Node with the given identifier could not be found.");
    }

    return ConstraintLinearSetRotationNode(nodePtr, rValue);
}

//! @brief adds a displacement constraint equation for a group of node
//! @param rNode pointer to group of nodes
//! @param rDirection direction of the constraint (in 2D a point with 2 entries, in 3D 3 entries, in 1D not used)
//! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
//! @return integer id to delete or modify the constraint
int NuTo::StructureBase::ConstraintLinearSetDisplacementNodeGroup(Group<NodeBase>* rGroup, const NuTo::FullMatrix<double>& rDirection, double rValue)
{
	this->mNodeNumberingRequired = true;
    //find unused integer id
    int id(0);
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(id);
    while (it!=mConstraintMap.end())
    {
        id++;
        it = mConstraintMap.find(id);
    }

    switch (mDimension)
    {
    case 1:
        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeGroupDisplacements1D(rGroup,rDirection(0,0),rValue));
        break;
    case 2:
        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeGroupDisplacements2D(rGroup,rDirection,rValue));
        break;
    case 3:
        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeGroupDisplacements3D(rGroup,rDirection,rValue));
        break;
    default:
        throw MechanicsException("[NuTo::StructureBase::ConstraintSetDisplacementNodeGroup] Incorrect dimension of the structure.");
    }
    return id;
}

//! @brief adds a fine scale displacement constraint equation for a group of node
//! @param rNode pointer to group of nodes
//! @param rDirection direction of the constraint (in 2D a point with 2 entries, in 3D 3 entries, in 1D not used)
//! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
//! @return integer id to delete or modify the constraint
int NuTo::StructureBase::ConstraintLinearSetFineScaleDisplacementNodeGroup(Group<NodeBase>* rGroup, const NuTo::FullMatrix<double>& rDirection, double rValue)
{
    this->mNodeNumberingRequired = true;
    //find unused integer id
    int id(0);
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(id);
    while (it!=mConstraintMap.end())
    {
        id++;
        it = mConstraintMap.find(id);
    }

    switch (mDimension)
    {
    case 1:
        throw MechanicsException("[NuTo::StructureBase::ConstraintSetfineScaleDisplacementNode] not implemented for 1D.");
        break;
    case 2:
        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeGroupFineScaleDisplacements2D(rGroup,rDirection,rValue));
        break;
    case 3:
        throw MechanicsException("[NuTo::StructureBase::ConstraintSetfineScaleDisplacementNode] not implemented for 2D.");
        break;
    default:
        throw MechanicsException("[NuTo::StructureBase::ConstraintSetDisplacementNode] Incorrect dimension of the structure.");
    }
    return id;
}

//! @brief adds a rotation constraint equation for a group of node
//! @param rNode pointer to group of nodes
//! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
//! @return integer id to delete or modify the constraint
int NuTo::StructureBase::ConstraintLinearSetRotationNodeGroup(Group<NodeBase>* rGroup, double rValue)
{
	this->mNodeNumberingRequired = true;
    //find unused integer id
    int id(0);
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(id);
    while (it!=mConstraintMap.end())
    {
        id++;
        it = mConstraintMap.find(id);
    }

    switch (mDimension)
    {
    case 1:
        throw MechanicsException("[NuTo::StructureBase::ConstraintLinearSetRotationNodeGroup] not implemented for 1D.");
        break;
    case 2:
        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeGroupRotations2D(rGroup,rValue));
        break;
    case 3:
        throw MechanicsException("[NuTo::StructureBase::ConstraintLinearSetRotationNodeGroup] not implemented for 3D.");
        break;
    default:
        throw MechanicsException("[NuTo::StructureBase::ConstraintLinearSetRotationNodeGroup] Incorrect dimension of the structure.");
    }
    return id;
}


//! @brief adds a constraint equation for a group of nodes
//! @param rGroupIdent identifier for group of nodes
//! @param rAttribute displacements, rotations, temperatures
//! @param rComponent e.g. the first (count from zero) displacement component
//! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
int NuTo::StructureBase::ConstraintLinearSetDisplacementNodeGroup(int rGroupIdent, const NuTo::FullMatrix<double>& rDirection, double rValue)
{
	this->mNodeNumberingRequired = true;
    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::ConstraintSetDisplacementNodeGroup] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Nodes)
        throw MechanicsException("[NuTo::Structure::ConstraintSetDisplacementNodeGroup] Group is not a node group.");
    Group<NodeBase> *nodeGroup = dynamic_cast<Group<NodeBase>*>(itGroup->second);
    assert(nodeGroup!=0);

    return ConstraintLinearSetDisplacementNodeGroup(nodeGroup,rDirection, rValue);
}

//! @brief adds a constraint equation for a group of nodes
//! @param rGroupIdent identifier for group of nodes
//! @param rAttribute displacements, rotations, temperatures
//! @param rComponent e.g. the first (count from zero) displacement component
//! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
int NuTo::StructureBase::ConstraintLinearSetFineScaleDisplacementNodeGroup(int rGroupIdent, const NuTo::FullMatrix<double>& rDirection, double rValue)
{
    this->mNodeNumberingRequired = true;
    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::ConstraintSetFineScaleDisplacementNodeGroup] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Nodes)
        throw MechanicsException("[NuTo::Structure::ConstraintSetFineScaleDisplacementNodeGroup] Group is not a node group.");
    Group<NodeBase> *nodeGroup = dynamic_cast<Group<NodeBase>*>(itGroup->second);
    assert(nodeGroup!=0);

    return ConstraintLinearSetFineScaleDisplacementNodeGroup(nodeGroup,rDirection, rValue);
}

//! @brief adds a constraint equation for a group of nodes
//! @param rGroupIdent identifier for group of nodes
//! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
int NuTo::StructureBase::ConstraintLinearSetRotationNodeGroup(int rGroupIdent, double rValue)
{
	this->mNodeNumberingRequired = true;
    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::ConstraintLinearSetRotationNodeGroup] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Nodes)
        throw MechanicsException("[NuTo::Structure::ConstraintLinearSetRotationNodeGroup] Group is not a node group.");
    Group<NodeBase> *nodeGroup = dynamic_cast<Group<NodeBase>*>(itGroup->second);
    assert(nodeGroup!=0);

    return ConstraintLinearSetRotationNodeGroup(nodeGroup, rValue);
}


//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::StructureBase::ConstraintGetNumLinearConstraints()const
{
    int numLinearConstraints(0);
    for (boost::ptr_map<int,ConstraintBase>::const_iterator itConstraint = mConstraintMap.begin(); itConstraint != mConstraintMap.end(); itConstraint++ )
    {
        numLinearConstraints+=itConstraint->second-> GetNumLinearConstraints();
    }
    return numLinearConstraints;
}

//! @brief calculates the constraint matrix that builds relations between the nodal dagrees of freedom
//! rConstraintMatrix*DOFS = RHS
//! @param rConstraintMatrix constraint matrix
//! @param rRHS right hand side
void NuTo::StructureBase::ConstraintGetConstraintMatrix(NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix, NuTo::FullMatrix<double>& rRHS)
{
    if (mNodeNumberingRequired)
    {
        throw MechanicsException("[NuTo::StructureBase::ConstraintGetConstraintMatrix] build global numbering first");
    }
    int numLinearConstraints = ConstraintGetNumLinearConstraints();
    rConstraintMatrix.Resize(numLinearConstraints,mNumDofs);
    rRHS.Resize(numLinearConstraints,1);
    int curConstraintEquations(0);
    for (boost::ptr_map<int,ConstraintBase>::const_iterator itConstraint = mConstraintMap.begin(); itConstraint != mConstraintMap.end(); itConstraint++ )
    {
        if (itConstraint->second->GetNumLinearConstraints()>0)
        {
            try
            {
                itConstraint->second->AsConstraintLinear()->AddToConstraintMatrix(curConstraintEquations, rConstraintMatrix, rRHS);
            }
            catch (MechanicsException& e)
            {
                e.AddMessage("[NuTo::StructureBase::ConstraintGetConstraintMatrix] mechanics exception while building constraint matrix for constraint with nonzero number of linear components.");
                throw e;
            }
            catch (...)
            {
                throw MechanicsException("[NuTo::StructureBase::ConstraintGetConstraintMatrix] error building constraint matrix for constraint with nonzero number of linear components.");
            }
        }
    }

    if (curConstraintEquations!=numLinearConstraints)
    {
        std::cout << "curConstraintEquations " << curConstraintEquations << std::endl;
        std::cout << "numConstraintEquations " << numLinearConstraints << std::endl;
        throw MechanicsException("[NuTo::StructureBase::ConstraintGetConstraintMatrix] Internal error, there is something wrong with the constraint equations.");
    }
}

//!@brief sets/modifies the right hand side of the constraint equations
//!@param rRHS new right hand side
//!@param rRHS new right hand side
void NuTo::StructureBase::ConstraintSetRHS(int rConstraintEquation, double rRHS)
{
	this->mNodeNumberingRequired = true;
    //find unused integer id
    //int id(0);
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(rConstraintEquation);
    if (it==mConstraintMap.end())
    {
    	throw MechanicsException("[NuTo::StructureBase::ConstraintSetRHS] Constraint equation does not exist.");
    }
    it->second->SetRHS(rRHS);
}

//!@brief gets the right hand side of the constraint equations
//!@param rConstraintEquation constraint equation
//!@return rRHS
double NuTo::StructureBase::ConstraintGetRHS(int rConstraintEquation)const
{
    boost::ptr_map<int,ConstraintBase>::const_iterator it = mConstraintMap.find(rConstraintEquation);
    if (it==mConstraintMap.end())
    {
    	throw MechanicsException("[NuTo::StructureBase::ConstraintSetRHS] Constraint equation does not exist.");
    }
    return it->second->GetRHS();
}

//!@brief sets/modifies the crack opening of a constraint equation (works only for periodic bc)
//!@param rConstraintEquation id of the constraint equation
//!@param rCrackOpening new crack opening (x,y)
void NuTo::StructureBase::ConstraintPeriodicSetCrackOpening(int rConstraintEquation, NuTo::FullMatrix<double> rCrackOpening)
{
    this->mNodeNumberingRequired = true;
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(rConstraintEquation);
    if (it==mConstraintMap.end())
    {
        throw MechanicsException("[NuTo::StructureBase::ConstraintPeriodicSetCrackOpening] Constraint equation does not exist.");
    }
    it->second->SetCrackOpening(rCrackOpening);
}

//!@brief sets/modifies the strain of a constraint equation (works only for periodic bc)
//!@param rConstraintEquation id of the constraint equation
//!@param rStrain new strain
void NuTo::StructureBase::ConstraintPeriodicSetStrain(int rConstraintEquation, NuTo::FullMatrix<double> rStrain)
{
    if (rStrain.GetNumColumns()==1)
        throw MechanicsException("[NuTo::StructureBase::ConstraintPeriodicSetStrain] Matrix has to have exactly one column");
    if (rStrain.GetNumRows()==3)
    {
        EngineeringStrain2D engineeringStrain;
        engineeringStrain.SetData(rStrain.mEigenMatrix.data());
        ConstraintPeriodicSetStrain2D(rConstraintEquation,engineeringStrain);
    }
    else
    {
        throw MechanicsException("[NuTo::StructureBase::ConstraintPeriodicSetStrain] Matrix has to have exactly three rows");
    }
}
    //!@brief sets/modifies the strain of a constraint equation (works only for periodic bc)
    //!@param rConstraintEquation id of the constraint equation
    //!@param rStrain new strain
void NuTo::StructureBase::ConstraintPeriodicSetStrain2D(int rConstraintEquation, const NuTo::EngineeringStrain2D& rStrain)
{
    this->mNodeNumberingRequired = true;
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(rConstraintEquation);
    if (it==mConstraintMap.end())
    {
        throw MechanicsException("[NuTo::StructureBase::ConstraintPeriodicSetStrain] Constraint equation does not exist.");
    }
    it->second->SetStrain(rStrain);
}

// create a constraint equation
int NuTo::StructureBase::ConstraintLinearEquationCreate(int rNode, const std::string& rDof, double rCoefficient, double rRHS)
{
	this->mNodeNumberingRequired = true;
    //find unused integer id
    int id(0);
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(id);
    while (it != mConstraintMap.end())
    {
        id++;
        it = mConstraintMap.find(id);
    }

    // create constraint
    this->ConstraintLinearEquationCreate(id, rNode, rDof, rCoefficient, rRHS);

    // return integer id
    return id;
}

// create a constraint equation
void NuTo::StructureBase::ConstraintLinearEquationCreate(int rConstraint, int rNode, const std::string& rDof, double rCoefficient, double rRHS)
{
	this->mNodeNumberingRequired = true;
    try
    {
        // convert dof string
        NuTo::Node::eAttributes dofType;
        int dofComponent;
        this->ConstraintEquationGetDofInformationFromString(rDof, dofType, dofComponent);

        // create constraint equation
        this->ConstraintLinearEquationCreate(rConstraint, rNode, dofType, dofComponent, rCoefficient, rRHS);
    }
    catch(NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstraintEquationCreate] error creating constraint equation");
        throw e;
    }
}

// create a constraint equation
void NuTo::StructureBase::ConstraintLinearEquationCreate(int rConstraint, int rNode, NuTo::Node::eAttributes rDofType, int rDofComponent, double rCoefficient, double rRHS)
{
	this->mNodeNumberingRequired = true;
    // check if constraint equation already exists
    if(this->mConstraintMap.find(rConstraint) != this->mConstraintMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::ConstraintEquationCreate] constraint equation already exist.");
    }

    try
    {
        // get node pointer
        NodeBase* nodePtr = this->NodeGetNodePtr(rNode);

        // create new constraint equation term
        ConstraintBase* constraintPtr = new ConstraintLinearEquation(nodePtr, rDofType, rDofComponent, rCoefficient, rRHS);

        // insert constraint equation into map
        this->mConstraintMap.insert(rConstraint, constraintPtr);
    }
    catch(NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstraintEquationCreate] error creating constraint equation");
        throw e;
    }
}

// add a term to a constraint equation
void NuTo::StructureBase::ConstraintLinearEquationAddTerm(int rConstraint, int rNode, const std::string& rDof, double rCoefficient)
{
	this->mNodeNumberingRequired = true;
    try
    {
        // convert dof string
        NuTo::Node::eAttributes dofType;
        int dofComponent;
        this->ConstraintEquationGetDofInformationFromString(rDof, dofType, dofComponent);

        // add term to constraint equation
        this->ConstraintLinearEquationAddTerm(rConstraint, rNode, dofType, dofComponent, rCoefficient);
    }
    catch(NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstraintEquationAddTerm] error adding term to constraint equation");
        throw e;
    }
}

// add a term to a constraint equation
void NuTo::StructureBase::ConstraintLinearEquationAddTerm(int rConstraint, int rNode, NuTo::Node::eAttributes rDofType, int rDofComponent, double rCoefficient)
{
	this->mNodeNumberingRequired = true;
    // get iterator
    boost::ptr_map<int,ConstraintBase>::iterator it = this->mConstraintMap.find(rConstraint);
    if(it == this->mConstraintMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::ConstraintEquationCreate] constraint equation does not exist.");
    }
    try
    {
        // get node pointer
        NodeBase* nodePtr = this->NodeGetNodePtr(rNode);

        // add term
        ConstraintLinearEquation* equationPtr = dynamic_cast<ConstraintLinearEquation*>(it->second);
        equationPtr->AddTerm(nodePtr, rDofType, rDofComponent, rCoefficient);
    }
    catch(NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstraintEquationAddTerm] error adding term to constraint equation");
        throw e;
    }
}

void NuTo::StructureBase::ConstraintEquationGetDofInformationFromString(const std::string& rDof, NuTo::Node::eAttributes& rDofType, int& rDofComponent)
{
    // convert string to upper-case
    std::string dofString;
    std::transform(rDof.begin(), rDof.end(), std::back_inserter(dofString), (int(*)(int)) toupper);
    if(dofString == "X_DISPLACEMENT")
    {
        rDofType = NuTo::Node::DISPLACEMENTS;
        rDofComponent = 0;
    }
    else if(dofString == "Y_DISPLACEMENT")
    {
        if(this->mDimension < 2)
        {
            throw NuTo::MechanicsException("[NuTo::StructureBase::ConstraintEquationGetDofInformationFromString] y-displacement dofs are only available in 2D or 3D structures.");
        }
        rDofType = NuTo::Node::DISPLACEMENTS;
        rDofComponent = 1;
    }
    else if(dofString == "Z_DISPLACEMENT")
    {
        if(this->mDimension < 3)
        {
            throw NuTo::MechanicsException("[NuTo::StructureBase::ConstraintEquationGetDofInformationFromString] z-displacement dofs are only available in 3D structures.");
        }
        rDofType = NuTo::Node::DISPLACEMENTS;
        rDofComponent = 2;
    }
    else if(dofString == "X_ROTATION")
    {
        if(this->mDimension < 3)
        {
            throw NuTo::MechanicsException("[NuTo::StructureBase::ConstraintEquationGetDofInformationFromString] x-rotation dofs are only available in 3D structures.");
        }
        rDofType = NuTo::Node::ROTATIONS;
        rDofComponent = 0;
    }
    else if(dofString == "Y_ROTATION")
    {
        if(this->mDimension < 3)
        {
            throw NuTo::MechanicsException("[NuTo::StructureBase::ConstraintEquationGetDofInformationFromString] y-rotation dofs are only available in 3D structures.");
        }
        rDofType = NuTo::Node::ROTATIONS;
        rDofComponent = 1;
    }
    else if(dofString == "Z_ROTATION")
    {
        rDofType = NuTo::Node::ROTATIONS;
        if(this->mDimension == 2)
        {
            rDofComponent = 0; // in 2D only one rotation
        }
        else if(this->mDimension == 3)
        {
            rDofComponent = 2;
        }
        else
        {
            throw NuTo::MechanicsException("[NuTo::StructureBase::ConstraintEquationGetDofInformationFromString] z-rotation dofs are only available in 2D and 3D structures.");
        }
    }
    else if(dofString == "TEMPERATURE")
    {
        rDofType = NuTo::Node::TEMPERATURES;
        rDofComponent = 0;
    }
    else
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::ConstraintEquationGetDofInformationFromString] invalid dof string.");
    }
}
//! @brief ... set periodic boundary conditions according to a prescibed angle of a localization zone
//! @param  rAngle... angle in deg
//! @param  rStrain... average strain to be applied (epsilon_xx, epsilon_yy, gamma_xy)
//! @param  rNodeGroupUpper... all nodes on the upper boundary
//! @param  rNodeGrouplower... all nodes on the lower boundary
//! @param  rNodeGroupLeft... all nodes on the left boundary
//! @param  rNodeGroupRight...  all nodes on the right boundary
int NuTo::StructureBase::ConstraintLinearDisplacementsSetPeriodic2D(double rAngle, NuTo::FullMatrix<double> rStrain,
        NuTo::FullMatrix<double> rCrackOpening, double rRadiusToCrackWithoutConstraints,
        int rNodeGroupUpperId, int rNodeGroupLowerId, int rNodeGroupLeftId, int rNodeGroupRightId)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    //check dimension of the structure
    if (mDimension!=2)
        throw MechanicsException("[NuTo::StructureBase::ConstraintSetPeriodicBoundaryConditions2D] only implemented for 2D");
    this->mNodeNumberingRequired = true;

    //find unused integer id
    int id(0);
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(id);
    while (it != mConstraintMap.end())
    {
        id++;
        it = mConstraintMap.find(id);
    }

    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rNodeGroupUpperId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::ConstraintSetPeriodicBoundaryConditions2D] Group of upper nodes with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Nodes)
        throw MechanicsException("[NuTo::Structure::ConstraintSetPeriodicBoundaryConditions2D] Group of upper nodes is not a node group.");
    Group<NodeBase> *nodeGroupUpperPtr = dynamic_cast<Group<NodeBase>*>(itGroup->second);
    assert(nodeGroupUpperPtr!=0);

     itGroup = mGroupMap.find(rNodeGroupLowerId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::ConstraintSetPeriodicBoundaryConditions2D] Group of lower nodes with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Nodes)
        throw MechanicsException("[NuTo::Structure::ConstraintSetPeriodicBoundaryConditions2D] Group of lower nodes is not a node group.");
    Group<NodeBase> *nodeGroupLowerPtr = dynamic_cast<Group<NodeBase>*>(itGroup->second);
    assert(nodeGroupLowerPtr!=0);

    itGroup = mGroupMap.find(rNodeGroupLeftId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::ConstraintSetPeriodicBoundaryConditions2D] Group of left nodes with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Nodes)
        throw MechanicsException("[NuTo::Structure::ConstraintSetPeriodicBoundaryConditions2D] Group of left nodes is not a node group.");
    Group<NodeBase> *nodeGroupLeftPtr = dynamic_cast<Group<NodeBase>*>(itGroup->second);
    assert(nodeGroupLeftPtr!=0);

    itGroup = mGroupMap.find(rNodeGroupRightId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::ConstraintSetPeriodicBoundaryConditions2D] Group of right nodes with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Nodes)
        throw MechanicsException("[NuTo::Structure::ConstraintSetPeriodicBoundaryConditions2D] Group of right nodes is not a node group.");
    Group<NodeBase> *nodeGroupRightPtr = dynamic_cast<Group<NodeBase>*>(itGroup->second);
    assert(nodeGroupRightPtr!=0);

    try
    {
        // create new constraint equation term
        if (rStrain.GetNumRows()!=3 || rStrain.GetNumColumns()!=1)
            throw MechanicsException("[NuTo::ConstraintNodeDisplacementsPeriodic2D::ConstraintNodeDisplacementsPeriodic2D] the strain is matrix (3,1) with (e_xx, e_yy, gamma_xy)");

        EngineeringStrain2D engineeringStrain;
        engineeringStrain.SetData(rStrain.mEigenMatrix.data());
        ConstraintBase* constraintPtr = new NuTo::ConstraintLinearDisplacementsPeriodic2D(this, rAngle, engineeringStrain, rCrackOpening, rRadiusToCrackWithoutConstraints,
                   nodeGroupUpperPtr, nodeGroupLowerPtr, nodeGroupLeftPtr, nodeGroupRightPtr);

        // insert constraint equation into map
        this->mConstraintMap.insert(id, constraintPtr);
    }
    catch(NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstraintSetPeriodicBoundaryConditions2D] error creating periodic boundary conditions in 2D");
        throw e;
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ConstraintDisplacementsSetPeriodic2D] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif

    // return integer id
    return id;

}

//!@brief number the free DOFS in the constraints (Lagrange multipliers)
//!@param rDOF current maximum DOF number, increased in the number
void NuTo::StructureBase::ConstraintNumberGlobalDofs(int& rDOF)
{
    for(boost::ptr_map<int,ConstraintBase>::iterator it =mConstraintMap.begin(); it!=mConstraintMap.end();it++)
    {
        if (it->second->GetNumLagrangeMultipliers()>0)
        {
            ConstraintLagrange* constraintLagrangePtr(it->second->AsConstraintLagrange());
            constraintLagrangePtr->SetGlobalDofs(rDOF);
        }
    }
}

//! @brief renumber the dofs of the Lagrange multipliers according to predefined ordering
//! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
void NuTo::StructureBase::ConstraintRenumberGlobalDofs(const std::vector<int>& mappingInitialToNewOrdering)
{
    for(boost::ptr_map<int,ConstraintBase>::iterator it =mConstraintMap.begin(); it!=mConstraintMap.end();it++)
    {
        if (it->second->GetNumLagrangeMultipliers()>0)
        {
            ConstraintLagrange* constraintLagrangePtr(it->second->AsConstraintLagrange());
            constraintLagrangePtr->RenumberGlobalDofs(mappingInitialToNewOrdering);
        }
    }
}

//! @brief extract dof values from the node (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::StructureBase::ConstraintExtractGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues)const
{
    for(boost::ptr_map<int,ConstraintBase>::const_iterator it =mConstraintMap.begin(); it!=mConstraintMap.end();it++)
    {
        if (it->second->GetNumLagrangeMultipliers()>0)
        {
            const ConstraintLagrange* constraintLagrangePtr(it->second->AsConstraintLagrange());
            constraintLagrangePtr->GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
        }
    }
}

//! @brief write dof values to the Lagrange multipliers (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::StructureBase::ConstraintMergeGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& dependentDofValues)
{
    for(boost::ptr_map<int,ConstraintBase>::iterator it =mConstraintMap.begin(); it!=mConstraintMap.end();it++)
    {
        if (it->second->GetNumLagrangeMultipliers()>0)
        {
            ConstraintLagrange* constraintLagrangePtr(it->second->AsConstraintLagrange());
            constraintLagrangePtr->SetGlobalDofValues(rActiveDofValues, dependentDofValues);
        }
    }
}

//! @brief ... add the contribution of Lagrange multipliers to the global system of equations
//! @param rMatrixJJ ... matrix jj
//! @param rMatrixJK ... matrix jk
void NuTo::StructureBase::ConstraintsBuildGlobalCoefficientSubMatrices0General(SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK)const
{
    // define variables storing the element contribution outside the loop
    NuTo::SparseMatrixCSRVector2Symmetric<double> constraintMatrix;
    std::vector<int> constraintMatrixGlobalDofs;

    // loop over all constraints
    for(boost::ptr_map<int,ConstraintBase>::const_iterator constraintIter = this->mConstraintMap.begin(); constraintIter != this->mConstraintMap.end(); constraintIter++)
    {
        // calculate element contribution
        if (constraintIter->second->IsLinear())
            continue;
        const ConstraintNonlinear* nonlinearPtr(constraintIter->second->AsConstraintNonlinear());
        nonlinearPtr->CalculateCoefficientMatrix_0(constraintMatrix, constraintMatrixGlobalDofs);

        assert(static_cast<unsigned int>(constraintMatrix.GetNumRows()) == constraintMatrixGlobalDofs.size());
        assert(static_cast<unsigned int>(constraintMatrix.GetNumColumns()) == constraintMatrixGlobalDofs.size());

        const std::vector<std::vector<double> >& values(constraintMatrix.GetValues());
        const std::vector<std::vector<int> >& columns(constraintMatrix.GetColumns());

        // write constraint contribution to global matrix
        for (unsigned int rowCount = 0; rowCount < values.size(); rowCount++)
        {
            int globalRowDof = constraintMatrixGlobalDofs[rowCount];
            if (globalRowDof < this->mNumActiveDofs)
            {
                for (unsigned int colCount = 0; colCount < values[rowCount].size(); colCount++)
                {
                    if (fabs(values[rowCount][colCount])>mToleranceStiffnessEntries)
                    {
                        int globalColumnDof = constraintMatrixGlobalDofs[columns[rowCount][colCount]-constraintMatrix.HasOneBasedIndexing()];
                        if (globalColumnDof < this->mNumActiveDofs)
                        {
                            rMatrixJJ.AddEntry(globalRowDof, globalColumnDof, values[rowCount][colCount]);
                        }
                        else
                        {
                            rMatrixJK.AddEntry(globalRowDof, globalColumnDof - this->mNumActiveDofs, values[rowCount][colCount]);
                        }
                    }
                }
            }

            //same thing, but for the transpose
            int globalColumnDof = constraintMatrixGlobalDofs[rowCount];
			for (unsigned int colCount = 0; colCount < values[rowCount].size(); colCount++)
			{
				if (fabs(values[rowCount][colCount])>mToleranceStiffnessEntries)
				{
					int globalRowDof = constraintMatrixGlobalDofs[columns[rowCount][colCount]-constraintMatrix.HasOneBasedIndexing()];
					if (globalRowDof==globalColumnDof)
						continue;
					if (globalRowDof < this->mNumActiveDofs)
					{
                        if (globalColumnDof < this->mNumActiveDofs)
                        {
                            // add upper triangle and diagonal
                            rMatrixJJ.AddEntry(globalRowDof, globalColumnDof, values[rowCount][colCount]);
                        }
                        else
                        {
                            rMatrixJK.AddEntry(globalRowDof , globalColumnDof - this->mNumActiveDofs , values[rowCount][colCount]);
                        }
                    }
                }
            }
        }
    }
}

//! @brief ... add the contribution of Lagrange multipliers to the global system of equations
//! @param rMatrixJJ ... matrix jj
//! @param rMatrixJK ... matrix jk
//! @param rMatrixKJ ... matrix kj
//! @param rMatrixKK ... matrix kk
void NuTo::StructureBase::ConstraintBuildGlobalCoefficientSubMatrices0General(SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK, SparseMatrix<double>& rMatrixKJ, SparseMatrix<double>& rMatrixKK)const
{
    // define variables storing the element contribution outside the loop
    NuTo::SparseMatrixCSRVector2Symmetric<double> constraintMatrix;
    std::vector<int> constraintMatrixGlobalDofs;

    // loop over all constraints
    for(boost::ptr_map<int,ConstraintBase>::const_iterator constraintIter = this->mConstraintMap.begin(); constraintIter != this->mConstraintMap.end(); constraintIter++)
    {
        // calculate element contribution
        if (constraintIter->second->IsLinear())
            continue;
        const ConstraintNonlinear* nonlinearPtr(constraintIter->second->AsConstraintNonlinear());
        nonlinearPtr->CalculateCoefficientMatrix_0(constraintMatrix, constraintMatrixGlobalDofs);
        assert(static_cast<unsigned int>(constraintMatrix.GetNumRows()) == constraintMatrixGlobalDofs.size());
        assert(static_cast<unsigned int>(constraintMatrix.GetNumColumns()) == constraintMatrixGlobalDofs.size());

        const std::vector<std::vector<double> >& values(constraintMatrix.GetValues());
        const std::vector<std::vector<int> >& columns(constraintMatrix.GetColumns());

        // write constraint contribution to global matrix
        for (unsigned int rowCount = 0; rowCount < values.size(); rowCount++)
        {
            int globalRowDof = constraintMatrixGlobalDofs[rowCount];
            if (globalRowDof < this->mNumActiveDofs)
            {
                for (unsigned int colCount = 0; colCount < values[rowCount].size(); colCount++)
                {
                    if (fabs(values[rowCount][colCount])>mToleranceStiffnessEntries)
                    {
                        int globalColumnDof = constraintMatrixGlobalDofs[columns[rowCount][colCount]-constraintMatrix.HasOneBasedIndexing()];
                        if (globalColumnDof < this->mNumActiveDofs)
                        {
                            rMatrixJJ.AddEntry(globalRowDof, globalColumnDof, values[rowCount][colCount]);
                        }
                        else
                        {
                            rMatrixJK.AddEntry(globalRowDof, globalColumnDof - this->mNumActiveDofs, values[rowCount][colCount]);
                        }
                    }

                }
            }
            else
            {
                for (unsigned int colCount = 0; colCount < values[rowCount].size(); colCount++)
                {
                    if (fabs(values[rowCount][colCount])>mToleranceStiffnessEntries)
                    {
                        int globalColumnDof = constraintMatrixGlobalDofs[columns[rowCount][colCount]-constraintMatrix.HasOneBasedIndexing()];
                        if (globalColumnDof < this->mNumActiveDofs)
                        {
                             rMatrixKJ.AddEntry(globalRowDof - this->mNumActiveDofs, globalColumnDof, values[rowCount][colCount]);
                        }
                        else
                        {
                             rMatrixKK.AddEntry(globalRowDof - this->mNumActiveDofs, globalColumnDof - this->mNumActiveDofs, values[rowCount][colCount]);
                        }
                    }
                }
            }
        }

        // write constraint contribution to global matrix - transpose
        for (unsigned int rowCount = 0; rowCount < values.size(); rowCount++)
        {
            int globalColumnDof = constraintMatrixGlobalDofs[rowCount];
            if (globalColumnDof < this->mNumActiveDofs)
            {
                for (unsigned int colCount = 0; colCount < values[rowCount].size(); colCount++)
                {
                    if (fabs(values[rowCount][colCount])>mToleranceStiffnessEntries)
                    {
                        int globalRowDof = constraintMatrixGlobalDofs[columns[rowCount][colCount]-constraintMatrix.HasOneBasedIndexing()];
                        if (globalRowDof==globalColumnDof)
                            continue;
                        if (globalRowDof < this->mNumActiveDofs)
                        {
                            rMatrixJJ.AddEntry(globalRowDof, globalColumnDof, values[rowCount][colCount]);
                        }
                        else
                        {
                            rMatrixKJ.AddEntry(globalRowDof - this->mNumActiveDofs, globalColumnDof, values[rowCount][colCount]);
                        }
                    }

                }
            }
            else
            {
                for (unsigned int colCount = 0; colCount < values[rowCount].size(); colCount++)
                {
                    if (fabs(values[rowCount][colCount])>mToleranceStiffnessEntries)
                    {
                        int globalRowDof = constraintMatrixGlobalDofs[colCount];
                        if (globalRowDof==globalColumnDof)
                            continue;
                        if (globalRowDof < this->mNumActiveDofs)
                        {
                            rMatrixJK.AddEntry(globalRowDof, globalColumnDof - this->mNumActiveDofs, values[rowCount][colCount]);
                        }
                        else
                        {
                            rMatrixKK.AddEntry(globalRowDof - this->mNumActiveDofs, globalColumnDof - this->mNumActiveDofs, values[rowCount][colCount]);
                        }
                    }
                }
            }
        }
    }
}

//! @brief ... add the contribution of Lagrange multipliers to the global system of equations
//! @param rMatrixJJ ... matrix jj
//! @param rMatrixJK ... matrix jk
void NuTo::StructureBase::ConstraintBuildGlobalCoefficientSubMatrices0Symmetric(SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK)const
{
    // define variables storing the element contribution outside the loop
    NuTo::SparseMatrixCSRVector2Symmetric<double> constraintMatrix;
    std::vector<int> constraintMatrixGlobalDofs;

    // loop over all elements
    for(boost::ptr_map<int,ConstraintBase>::const_iterator constraintIter = this->mConstraintMap.begin(); constraintIter != this->mConstraintMap.end(); constraintIter++)
    {
        // calculate element contribution
        if (constraintIter->second->IsLinear())
            continue;
        const ConstraintNonlinear* nonlinearPtr(constraintIter->second->AsConstraintNonlinear());
        nonlinearPtr->CalculateCoefficientMatrix_0(constraintMatrix, constraintMatrixGlobalDofs);
        assert(static_cast<unsigned int>(constraintMatrix.GetNumRows()) == constraintMatrixGlobalDofs.size());
        assert(static_cast<unsigned int>(constraintMatrix.GetNumColumns()) == constraintMatrixGlobalDofs.size());

        const std::vector<std::vector<double> >& values(constraintMatrix.GetValues());
        const std::vector<std::vector<int> >& columns(constraintMatrix.GetColumns());

        // write constraint contribution to global matrix
        for (unsigned int rowCount = 0; rowCount < values.size(); rowCount++)
        {
            int globalRowDof = constraintMatrixGlobalDofs[rowCount];
            if (globalRowDof < this->mNumActiveDofs)
            {
                for (unsigned int colCount = 0; colCount < values[rowCount].size(); colCount++)
                {
                    if (fabs(values[rowCount][colCount])>mToleranceStiffnessEntries)
                    {
                        int globalColumnDof = constraintMatrixGlobalDofs[columns[rowCount][colCount]-constraintMatrix.HasOneBasedIndexing()];
                        if (globalColumnDof < this->mNumActiveDofs)
                        {
                            rMatrixJJ.AddEntry(globalRowDof, globalColumnDof, values[rowCount][colCount]);
                        }
                        else
                        {
                            rMatrixJK.AddEntry(globalRowDof, globalColumnDof - this->mNumActiveDofs, values[rowCount][colCount]);
                        }
                    }
                }
            }
        }
    }
}

//! @brief ... add the contribution of Lagrange multipliers to the global system of equations
//! @param rMatrixJJ ... matrix jj
//! @param rMatrixJK ... matrix jk
//! @param rMatrixKK ... matrix kk
void NuTo::StructureBase::ConstraintBuildGlobalCoefficientSubMatrices0Symmetric(SparseMatrix<double>& rMatrixJJ, SparseMatrix<double>& rMatrixJK, SparseMatrix<double>& rMatrixKK) const
{
    // define variables storing the element contribution outside the loop
    NuTo::SparseMatrixCSRVector2Symmetric<double> constraintMatrix;
    std::vector<int> constraintMatrixGlobalDofs;

    // loop over all elements
    for(boost::ptr_map<int,ConstraintBase>::const_iterator constraintIter = this->mConstraintMap.begin(); constraintIter != this->mConstraintMap.end(); constraintIter++)
    {
        // calculate element contribution
        if (constraintIter->second->IsLinear())
            continue;
        const ConstraintNonlinear* nonlinearPtr(constraintIter->second->AsConstraintNonlinear());
        nonlinearPtr->CalculateCoefficientMatrix_0(constraintMatrix, constraintMatrixGlobalDofs);
        assert(static_cast<unsigned int>(constraintMatrix.GetNumRows()) == constraintMatrixGlobalDofs.size());
        assert(static_cast<unsigned int>(constraintMatrix.GetNumColumns()) == constraintMatrixGlobalDofs.size());

        const std::vector<std::vector<double> >& values(constraintMatrix.GetValues());
        const std::vector<std::vector<int> >& columns(constraintMatrix.GetColumns());

        // write constraint contribution to global matrix
        for (unsigned int rowCount = 0; rowCount < values.size(); rowCount++)
        {
            int globalRowDof = constraintMatrixGlobalDofs[rowCount];
            if (globalRowDof < this->mNumActiveDofs)
            {
                for (unsigned int colCount = 0; colCount < values[rowCount].size(); colCount++)
                {
                    if (fabs(values[rowCount][colCount])>mToleranceStiffnessEntries)
                    {
                        int globalColumnDof = constraintMatrixGlobalDofs[columns[rowCount][colCount]-constraintMatrix.HasOneBasedIndexing()];
                        if (globalColumnDof < this->mNumActiveDofs)
                        {
                            // add upper triangle and diagonal
                            if(globalColumnDof >= globalRowDof)
                            {
                                rMatrixJJ.AddEntry(globalRowDof, globalColumnDof, values[rowCount][colCount]);
                            }
                        }
                        else
                        {
                            rMatrixJK.AddEntry(globalRowDof, globalColumnDof - this->mNumActiveDofs, values[rowCount][colCount]);
                        }
                    }
                }
            }
            else
            {
                for (unsigned int colCount = 0; colCount < values[rowCount].size(); colCount++)
                {
                    if (fabs(values[rowCount][colCount])>mToleranceStiffnessEntries)
                    {
                        int globalColumnDof = constraintMatrixGlobalDofs[columns[rowCount][colCount]-constraintMatrix.HasOneBasedIndexing()];
                        if (globalColumnDof >= this->mNumActiveDofs)
                        {
                            // add upper triangle and diagonal
                            if(globalColumnDof >= globalRowDof)
                            {
                                rMatrixKK.AddEntry(globalRowDof - this->mNumActiveDofs, globalColumnDof - this->mNumActiveDofs, values[rowCount][colCount]);
                            }
                        }
                    }
                }
            }
        }
    }
}

//! @brief ... add the contribution of Lagrange multipliers to the global system of equations
//! @param rActiveDofGradientVector ... gradient of active dofs
//! @param rDependentDofGradientVector ... gradient of dependent dofs
void NuTo::StructureBase::ConstraintBuildGlobalGradientInternalPotentialSubVectors(NuTo::FullMatrix<double>& rActiveDofGradientVector, NuTo::FullMatrix<double>& rDependentDofGradientVector) const
{
    // define variables storing the element contribution outside the loop
    NuTo::FullMatrix<double> constraintVector;
    std::vector<int> constraintVectorGlobalDofs;

    // loop over all constraints
    for (boost::ptr_map<int,ConstraintBase>::const_iterator constraintIter = this->mConstraintMap.begin(); constraintIter != this->mConstraintMap.end(); constraintIter++)
    {
        // calculate element contribution
        if (constraintIter->second->IsLinear())
            continue;
        const ConstraintNonlinear* nonlinearPtr(constraintIter->second->AsConstraintNonlinear());
        nonlinearPtr->CalculateGradientInternalPotential(constraintVector,constraintVectorGlobalDofs);
        assert(static_cast<unsigned int>(constraintVector.GetNumRows()) == constraintVectorGlobalDofs.size());
        assert(static_cast<unsigned int>(constraintVector.GetNumColumns()) == 1);
        //constraintVector.Trans().Info();
        //std::cout << std::endl;

        // write constraint contribution to global vectors
        for (unsigned int rowCount = 0; rowCount < constraintVectorGlobalDofs.size(); rowCount++)
        {
            int globalRowDof = constraintVectorGlobalDofs[rowCount];
            if (globalRowDof < this->mNumActiveDofs)
            {
                rActiveDofGradientVector(globalRowDof,0) += constraintVector(rowCount,0);
            }
            else
            {
                globalRowDof -= this->mNumActiveDofs;
                assert(globalRowDof < this->mNumDofs - this->mNumActiveDofs);
                rDependentDofGradientVector(globalRowDof,0) += constraintVector(rowCount,0);
            }
        }
    }
}

double NuTo::StructureBase::ConstraintTotalGetTotalEnergy()const
{
    double energy(0);
    // loop over all constraints
    for (boost::ptr_map<int,ConstraintBase>::const_iterator constraintIter = this->mConstraintMap.begin(); constraintIter != this->mConstraintMap.end(); constraintIter++)
    {
        // calculate constraint contribution
        if (constraintIter->second->IsLinear()==false)
        {
            const ConstraintNonlinear* constraintPtr (constraintIter->second->AsConstraintNonlinear());

            energy+=constraintPtr->CalculateTotalPotential();
        }
    }
    return energy;
}


//! @brief writes the Lagrange multiplier and Slack variables (inequalities) of a constraint to the prescribed matrix
//! @param ConstraintId constraint id
//! @param rMultiplier Lagrange multiplier (first col Lagrange, evtl. second col Slackvariables)
void NuTo::StructureBase::ConstraintLagrangeGetMultiplier(int ConstraintId, NuTo::FullMatrix<double>& rMultiplier)const
{
    // get iterator
    boost::ptr_map<int,ConstraintBase>::const_iterator it = this->mConstraintMap.find(ConstraintId);
    if(it == this->mConstraintMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::ConstraintLagrangeGetMultiplier] constraint equation does not exist.");
    }
    try
    {
        if (it->second->GetNumLagrangeMultipliers()>0)
        {
            const ConstraintLagrange* constraintPtr (it->second->AsConstraintLagrange());
            constraintPtr->GetLagrangeMultiplier(rMultiplier);
        }
        else
            throw MechanicsException("[NuTo::StructureBase::ConstraintLagrangeGetMultiplier] constraint has no Lagrange multipliers.");
    }
    catch(NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstraintLagrangeGetMultiplier] error getting Lagrange multipliers.");
        throw e;
    }
}

//! @brief sets the penalty stiffness of the augmented Lagragian to the prescribed value
//! @param ConstraintId constraint id
//! @param rPenalty penalty parameter
void NuTo::StructureBase::ConstraintLagrangeSetPenaltyStiffness(int ConstraintId, double rPenalty)
{
    // get iterator
    boost::ptr_map<int,ConstraintBase>::iterator it = this->mConstraintMap.find(ConstraintId);
    if(it == this->mConstraintMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::ConstraintLagrangeSetPenaltyStiffness] constraint equation does not exist.");
    }
    try
    {
        if (it->second->GetNumLagrangeMultipliers()>0)
        {
            ConstraintLagrange* constraintPtr (it->second->AsConstraintLagrange());
            constraintPtr->SetPenaltyStiffness(rPenalty);
        }
        else
            throw MechanicsException("[NuTo::StructureBase::ConstraintLagrangeSetPenaltyStiffness] constraint has no Lagrange multipliers.");
    }
    catch(NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstraintLagrangeSetPenaltyStiffness] error setting penalty stiffness.");
        throw e;
    }
}

//! @brief info about the elements in the Structure
void NuTo::StructureBase::ConstraintInfo(int rVerboseLevel)const
{
    mLogger <<"number of constraints: " << mConstraintMap.size() << "\n";
    if (rVerboseLevel>3)
    {
        for (boost::ptr_map<int,ConstraintBase>::const_iterator it = mConstraintMap.begin(); it!= mConstraintMap.end(); it++)
        {
            if (rVerboseLevel>4)
            {
                mLogger << "\t constraint " << it->first << "\t";
                it->second->Info(rVerboseLevel);
            }
        }
    }
}


//!@brief deletes a constraint equation
//!@param rConstraintEquation id of the constraint equation
//!@param rCrackOpening new crack opening (x,y)
void NuTo::StructureBase::ConstraintDelete(int rConstraintId)
{
    this->mNodeNumberingRequired = true;
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(rConstraintId);
    if (it==mConstraintMap.end())
    {
        throw MechanicsException("[NuTo::StructureBase::ConstraintDelete] Constraint equation does not exist.");
    }
    mConstraintMap.erase(it);
}

//! @brief releases a constraint, (remove from the list but don't delete it)
//!@param rConstraintEquation id of the constraint equation
//! @return ptr to constraint
NuTo::ConstraintBase* NuTo::StructureBase::ConstraintRelease(int rConstraintId)
{
    this->mNodeNumberingRequired = true;
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(rConstraintId);
    if (it==mConstraintMap.end())
    {
        throw MechanicsException("[NuTo::StructureBase::ConstraintRelease] Constraint equation does not exist.");
    }
    boost::ptr_map<int,ConstraintBase>::auto_type ptr = mConstraintMap.release(it);
    return ptr.release();
}

//! @brief adds a constraint to the map
//! @param ConstraintId constraint id
//! @param
void NuTo::StructureBase::ConstraintAdd(int rConstraintId, NuTo::ConstraintBase* rConstraint)
{
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(rConstraintId);
    if (it!=mConstraintMap.end())
    	throw MechanicsException("[NuTo::StructureBase::ConstraintAdd] Id already exists in constraint map.");
    mConstraintMap.insert(rConstraintId, rConstraint);

}
