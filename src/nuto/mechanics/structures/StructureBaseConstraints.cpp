// $Id$

#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/constraints/ConstraintEquation.h"
#include "nuto/mechanics/constraints/ConstraintNodeDisplacements1D.h"
#include "nuto/mechanics/constraints/ConstraintNodeDisplacements3D.h"
#include "nuto/mechanics/constraints/ConstraintNodeGroupDisplacements1D.h"
#include "nuto/mechanics/constraints/ConstraintNodeGroupDisplacements3D.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"

//! @brief adds a displacement constraint equation for a node
//! @param rNode pointer to node
//! @param rDirection direction of the constraint (in 2D a point with 2 entries, in 3D 3 entries, in 1D not used)
//! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
//! @return integer id to delete or modify the constraint
int NuTo::StructureBase::ConstraintSetDisplacementNode(NodeBase* rNode, const NuTo::FullMatrix<double>& rDirection, double rValue)
{
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
        throw MechanicsException("[NuTo::StructureBase::ConstraintDisplacementNode] To be implemented.");
        break;
    case 3:
        mConstraintMap.insert(id, new NuTo::ConstraintNodeDisplacements3D(rNode,rDirection,rValue));
        break;
    default:
        throw MechanicsException("[NuTo::StructureBase::ConstraintDisplacementNode] Incorrect dimension of the structure.");
    }
    return id;
}

//! @brief adds a displacement constraint equation for a node
//! @param rNode identifier for node
//! @param rComponent e.g. the first (count from zero) displacement component
//! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
int  NuTo::StructureBase::ConstraintSetDisplacementNode(int rIdent, const NuTo::FullMatrix<double>& rDirection, double rValue)
{
    NodeBase* nodePtr;
    try
    {
        nodePtr = NodeGetNodePtr(rIdent);
    }
    catch (NuTo::MechanicsException &e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstraintDisplacementNode] Node with the given identifier could not be found.");
        throw e;
    }
    catch (...)
    {
        throw MechanicsException("[NuTo::StructureBase::ConstraintDisplacementNode] Node with the given identifier could not be found.");
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
        throw MechanicsException("[NuTo::StructureBase::ConstraintDisplacementNode] To be implemented.");
        break;
    case 3:
        mConstraintMap.insert(id, new NuTo::ConstraintNodeGroupDisplacements3D(rGroup,rDirection,rValue));
        break;
    default:
        throw MechanicsException("[NuTo::StructureBase::ConstraintDisplacementNode] Incorrect dimension of the structure.");
    }
    return id;
}
//! @brief adds a constraint equation for a group of nodes
//! @param rGroupIdent identifier for group of nodes
//! @param rAttribute displacements, rotations, temperatures
//! @param rComponent e.g. the first (count from zero) displacement component
//! @param rValue prescribed value (e.g. zero to fix a displacement to zero)
int NuTo::StructureBase::ConstraintSetDisplacementNodeGroup(std::string rGroupIdent, const NuTo::FullMatrix<double>& rDirection, double rValue)
{
    boost::ptr_map<std::string,GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::ConstraintDisplacementNodeGroup] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::GroupBase::Nodes)
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
        throw MechanicsException("[NuTo::StructureBase::ConstraintGetConstraintMatrix] Internal error, there is something wrong with the constraint equations.");
}

//!@brief sets/modifies the right hand side of the constraint equations
//!@param rRHS new right hand side
//!@param rRHS new right hand side
void NuTo::StructureBase::ConstraintSetRHS(int rConstraintEquation, double rRHS)
{
    //find unused integer id
    //int id(0);
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(rConstraintEquation);
    if (it==mConstraintMap.end())
    {
    	throw MechanicsException("[NuTo::StructureBase::ConstraintSetRHS] Constraint equation does not exist.");
    }
    it->second->SetRHS(rRHS);
}

// create a constraint equation
int NuTo::StructureBase::ConstraintEquationCreate(int rNode, const std::string& rDof, double rCoefficient, double rRHS)
{
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
    try
    {
        // convert dof string
        NuTo::NodeBase::eAttributes dofType;
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
void NuTo::StructureBase::ConstraintEquationCreate(int rConstraint, int rNode, NuTo::NodeBase::eAttributes rDofType, int rDofComponent, double rCoefficient, double rRHS)
{
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
    try
    {
        // convert dof string
        NuTo::NodeBase::eAttributes dofType;
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
void NuTo::StructureBase::ConstraintEquationAddTerm(int rConstraint, int rNode, NuTo::NodeBase::eAttributes rDofType, int rDofComponent, double rCoefficient)
{
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

void NuTo::StructureBase::ConstraintEquationGetDofInformationFromString(const std::string& rDof, NuTo::NodeBase::eAttributes& rDofType, int& rDofComponent)
{
    // convert string to upper-case
    std::string dofString;
    std::transform(rDof.begin(), rDof.end(), std::back_inserter(dofString), (int(*)(int)) toupper);
    if(dofString == "X_DISPLACEMENT")
    {
        rDofType = NuTo::NodeBase::DISPLACEMENTS;
        rDofComponent = 0;
    }
    else if(dofString == "Y_DISPLACEMENT")
    {
        if(this->mDimension < 2)
        {
            throw NuTo::MechanicsException("[NuTo::StructureBase::ConstraintEquationGetDofInformationFromString] y-displacement dofs are only available in 2D or 3D structures.");
        }
        rDofType = NuTo::NodeBase::DISPLACEMENTS;
        rDofComponent = 1;
    }
    else if(dofString == "Z_DISPLACEMENT")
    {
        if(this->mDimension < 3)
        {
            throw NuTo::MechanicsException("[NuTo::StructureBase::ConstraintEquationGetDofInformationFromString] z-displacement dofs are only available in 3D structures.");
        }
        rDofType = NuTo::NodeBase::DISPLACEMENTS;
        rDofComponent = 2;
    }
    else if(dofString == "X_ROTATION")
    {
        if(this->mDimension < 3)
        {
            throw NuTo::MechanicsException("[NuTo::StructureBase::ConstraintEquationGetDofInformationFromString] x-rotation dofs are only available in 3D structures.");
        }
        rDofType = NuTo::NodeBase::ROTATIONS;
        rDofComponent = 0;
    }
    else if(dofString == "Y_ROTATION")
    {
        if(this->mDimension < 3)
        {
            throw NuTo::MechanicsException("[NuTo::StructureBase::ConstraintEquationGetDofInformationFromString] y-rotation dofs are only available in 3D structures.");
        }
        rDofType = NuTo::NodeBase::ROTATIONS;
        rDofComponent = 1;
    }
    else if(dofString == "Z_ROTATION")
    {
        rDofType = NuTo::NodeBase::ROTATIONS;
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
        rDofType = NuTo::NodeBase::TEMPERATURES;
        rDofComponent = 0;
    }
    else
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::ConstraintEquationGetDofInformationFromString] invalid dof string.");
    }
}
