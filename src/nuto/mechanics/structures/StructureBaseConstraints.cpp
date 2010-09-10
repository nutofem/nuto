// $Id$

#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/constraints/ConstraintDisplacementsPeriodic2D.h"
#include "nuto/mechanics/constraints/ConstraintEquation.h"
#include "nuto/mechanics/constraints/ConstraintNodeDisplacements1D.h"
#include "nuto/mechanics/constraints/ConstraintNodeDisplacements2D.h"
#include "nuto/mechanics/constraints/ConstraintNodeDisplacements3D.h"
#include "nuto/mechanics/constraints/ConstraintNodeGroupDisplacements1D.h"
#include "nuto/mechanics/constraints/ConstraintNodeGroupDisplacements2D.h"
#include "nuto/mechanics/constraints/ConstraintNodeGroupDisplacements3D.h"

//! @brief adds a displacement constraint equation for a node
//! @param rNode pointer to node
//! @param rDirection direction of the constraint (in 2D a point with 2 entries, in 3D 3 entries, in 1D not used)
//! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
//! @return integer id to delete or modify the constraint
int NuTo::StructureBase::ConstraintSetDisplacementNode(NodeBase* rNode, const NuTo::FullMatrix<double>& rDirection, double rValue)
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
        mConstraintMap.insert(id, new NuTo::ConstraintNodeDisplacements1D(rNode,rDirection(0,0),rValue));
        break;
    case 2:
        mConstraintMap.insert(id, new NuTo::ConstraintNodeDisplacements2D(rNode,rDirection,rValue));
        break;
    case 3:
        mConstraintMap.insert(id, new NuTo::ConstraintNodeDisplacements3D(rNode,rDirection,rValue));
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
int  NuTo::StructureBase::ConstraintSetDisplacementNode(int rIdent, const NuTo::FullMatrix<double>& rDirection, double rValue)
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

    return ConstraintSetDisplacementNode(nodePtr,rDirection, rValue);
}

//! @brief adds a displacement constraint equation for a group of node
//! @param rNode pointer to group of nodes
//! @param rDirection direction of the constraint (in 2D a point with 2 entries, in 3D 3 entries, in 1D not used)
//! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
//! @return integer id to delete or modify the constraint
int NuTo::StructureBase::ConstraintSetDisplacementNodeGroup(Group<NodeBase>* rGroup, const NuTo::FullMatrix<double>& rDirection, double rValue)
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
        mConstraintMap.insert(id, new NuTo::ConstraintNodeGroupDisplacements1D(rGroup,rDirection(0,0),rValue));
        break;
    case 2:
        mConstraintMap.insert(id, new NuTo::ConstraintNodeGroupDisplacements2D(rGroup,rDirection,rValue));
        break;
    case 3:
        mConstraintMap.insert(id, new NuTo::ConstraintNodeGroupDisplacements3D(rGroup,rDirection,rValue));
        break;
    default:
        throw MechanicsException("[NuTo::StructureBase::ConstraintSetDisplacementNode] Incorrect dimension of the structure.");
    }
    return id;
}
//! @brief adds a constraint equation for a group of nodes
//! @param rGroupIdent identifier for group of nodes
//! @param rAttribute displacements, rotations, temperatures
//! @param rComponent e.g. the first (count from zero) displacement component
//! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
int NuTo::StructureBase::ConstraintSetDisplacementNodeGroup(int rGroupIdent, const NuTo::FullMatrix<double>& rDirection, double rValue)
{
	this->mNodeNumberingRequired = true;
    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::ConstraintDisplacementNodeGroup] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Nodes)
        throw MechanicsException("[NuTo::Structure::NodeGroupAddSingleConstraint] Group is not a node group.");
    Group<NodeBase> *nodeGroup = dynamic_cast<Group<NodeBase>*>(itGroup->second);
    assert(nodeGroup!=0);

    return ConstraintSetDisplacementNodeGroup(nodeGroup,rDirection, rValue);
}

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::StructureBase::ConstraintGetNumConstraintEquations()const
{
    int numConstraintEquations(0);
    for (boost::ptr_map<int,ConstraintBase>::const_iterator itConstraint = mConstraintMap.begin(); itConstraint != mConstraintMap.end(); itConstraint++ )
    {
        numConstraintEquations+=itConstraint->second-> GetNumConstraintEquations();
    }
    return numConstraintEquations;
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
    int numConstraintEquations = ConstraintGetNumConstraintEquations();
    rConstraintMatrix.Resize(numConstraintEquations,mNumDofs);
    rRHS.Resize(numConstraintEquations,1);
    int curConstraintEquations(0);
    for (boost::ptr_map<int,ConstraintBase>::const_iterator itConstraint = mConstraintMap.begin(); itConstraint != mConstraintMap.end(); itConstraint++ )
    {
        itConstraint->second->AddToConstraintMatrix(curConstraintEquations, rConstraintMatrix, rRHS);
    }

    if (curConstraintEquations!=numConstraintEquations)
    {
        std::cout << "curConstraintEquations " << curConstraintEquations << std::endl;
        std::cout << "numConstraintEquations " << numConstraintEquations << std::endl;
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

//!@brief sets/modifies the strain of a constraint equation (works only for periodic bc)
//!@param rConstraintEquation id of the constraint equation
//!@param rRHS new strain
void NuTo::StructureBase::ConstraintPeriodicSetStrain(int rConstraintEquation, NuTo::FullMatrix<double> rStrain)
{
    this->mNodeNumberingRequired = true;
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(rConstraintEquation);
    if (it==mConstraintMap.end())
    {
        throw MechanicsException("[NuTo::StructureBase::ConstraintSetRHS] Constraint equation does not exist.");
    }
    it->second->SetStrain(rStrain);
}
// create a constraint equation
int NuTo::StructureBase::ConstraintEquationCreate(int rNode, const std::string& rDof, double rCoefficient, double rRHS)
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
    this->ConstraintEquationCreate(id, rNode, rDof, rCoefficient, rRHS);

    // return integer id
    return id;
}

// create a constraint equation
void NuTo::StructureBase::ConstraintEquationCreate(int rConstraint, int rNode, const std::string& rDof, double rCoefficient, double rRHS)
{
	this->mNodeNumberingRequired = true;
    try
    {
        // convert dof string
        NuTo::Node::eAttributes dofType;
        int dofComponent;
        this->ConstraintEquationGetDofInformationFromString(rDof, dofType, dofComponent);

        // create constraint equation
        this->ConstraintEquationCreate(rConstraint, rNode, dofType, dofComponent, rCoefficient, rRHS);
    }
    catch(NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstraintEquationCreate] error creating constraint equation");
        throw e;
    }
}

// create a constraint equation
void NuTo::StructureBase::ConstraintEquationCreate(int rConstraint, int rNode, NuTo::Node::eAttributes rDofType, int rDofComponent, double rCoefficient, double rRHS)
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
        ConstraintBase* constraintPtr = new ConstraintEquation(nodePtr, rDofType, rDofComponent, rCoefficient, rRHS);

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
void NuTo::StructureBase::ConstraintEquationAddTerm(int rConstraint, int rNode, const std::string& rDof, double rCoefficient)
{
	this->mNodeNumberingRequired = true;
    try
    {
        // convert dof string
        NuTo::Node::eAttributes dofType;
        int dofComponent;
        this->ConstraintEquationGetDofInformationFromString(rDof, dofType, dofComponent);

        // add term to constraint equation
        this->ConstraintEquationAddTerm(rConstraint, rNode, dofType, dofComponent, rCoefficient);
    }
    catch(NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstraintEquationAddTerm] error adding term to constraint equation");
        throw e;
    }
}

// add a term to a constraint equation
void NuTo::StructureBase::ConstraintEquationAddTerm(int rConstraint, int rNode, NuTo::Node::eAttributes rDofType, int rDofComponent, double rCoefficient)
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
        ConstraintEquation* equationPtr = dynamic_cast<ConstraintEquation*>(it->second);
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
int NuTo::StructureBase::ConstraintDisplacementsSetPeriodic2D(double rAngle, NuTo::FullMatrix<double> rStrain,
        int rNodeGroupUpperId, int rNodeGroupLowerId, int rNodeGroupLeftId, int rNodeGroupRightId)
{
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
        ConstraintBase* constraintPtr = new NuTo::ConstraintDisplacementsPeriodic2D(this, rAngle, rStrain, nodeGroupUpperPtr, nodeGroupLowerPtr, nodeGroupLeftPtr, nodeGroupRightPtr);

        // insert constraint equation into map
        this->mConstraintMap.insert(id, constraintPtr);
    }
    catch(NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstraintSetPeriodicBoundaryConditions2D] error creating periodic boundary conditions in 2D");
        throw e;
    }

    // return integer id
    return id;

}

