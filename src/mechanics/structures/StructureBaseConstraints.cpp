// $Id$
#include <eigen3/Eigen/Dense>
#include "base/Timer.h"

#include <boost/foreach.hpp>
#include "math/SparseMatrixCSRVector2.h"

#include "mechanics/structures/StructureBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/elements/ContinuumElement.h"
#include "mechanics/constraints/ConstraintLinearDerivativeNonlocalTotalStrain1D.h"
#include "mechanics/constraints/ConstraintLinearDisplacementsPeriodic2D.h"
#include "mechanics/constraints/ConstraintLinearEquation.h"
#include "mechanics/constraints/ConstraintLinearNodeDisplacements1D.h"
#include "mechanics/constraints/ConstraintLinearNodeDisplacements2D.h"
#include "mechanics/constraints/ConstraintLinearNodeDisplacements3D.h"
#include "mechanics/constraints/ConstraintLinearNodeGroupDisplacements1D.h"
#include "mechanics/constraints/ConstraintLinearNodeGroupDisplacements2D.h"
#include "mechanics/constraints/ConstraintLinearNodeGroupDisplacements3D.h"
#include "mechanics/constraints/ConstraintLinearNodeGroupRotations2D.h"
#include "mechanics/constraints/ConstraintLinearNodeGroupTemperature.h"
#include "mechanics/constraints/ConstraintLinearNodeRelativeHumidity.h"
#include "mechanics/constraints/ConstraintLinearNodeRotations2D.h"
#include "mechanics/constraints/ConstraintLinearNodeTemperature.h"
#include "mechanics/constraints/ConstraintLinearNodeWaterVolumeFraction.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/interpolationtypes/InterpolationBase.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"


int NuTo::StructureBase::ConstraintLinearSetDisplacementNode(NodeBase* rNode, const Eigen::VectorXd& rDirection, double rValue)
{
	this->mNodeNumberingRequired = true;

    int id = GetUnusedId(mConstraintMap);

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
        throw MechanicsException(__PRETTY_FUNCTION__,"Incorrect dimension of the structure.");
    }
    return id;
}


int NuTo::StructureBase::ConstraintLinearSetRotationNode(NodeBase* rNode, double rValue)
{
	this->mNodeNumberingRequired = true;

    int id = GetUnusedId(mConstraintMap);

    switch (mDimension)
    {
    case 1:
        throw MechanicsException(__PRETTY_FUNCTION__,"not implemented for 1D.");
        break;
    case 2:
        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeRotations2D(rNode,rValue));
        break;
    case 3:
        throw MechanicsException(__PRETTY_FUNCTION__,"not implemented for 3D.");
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__,"Incorrect dimension of the structure.");
    }
    return id;
}

int  NuTo::StructureBase::ConstraintLinearSetDisplacementNode(int rIdent, const Eigen::VectorXd& rDirection, double rValue)
{
	this->mNodeNumberingRequired = true;
    NodeBase* nodePtr;
    try
    {
        nodePtr = NodeGetNodePtr(rIdent);
    }
    catch (NuTo::MechanicsException &e)
    {
        e.AddMessage(__PRETTY_FUNCTION__,"Node with the given identifier could not be found.");
        throw;
    }
    catch (...)
    {
        throw MechanicsException(__PRETTY_FUNCTION__,"Node with the given identifier could not be found.");
    }

    return ConstraintLinearSetDisplacementNode(nodePtr,rDirection, rValue);
}

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
        e.AddMessage(__PRETTY_FUNCTION__,"Node with the given identifier could not be found.");
        throw;
    }
    catch (...)
    {
        throw MechanicsException(__PRETTY_FUNCTION__,"Node with the given identifier could not be found.");
    }

    return ConstraintLinearSetRotationNode(nodePtr, rValue);
}

int NuTo::StructureBase::ConstraintLinearSetRelativeHumidityNode(NodeBase* rNode, double rValue)
{
    this->mNodeNumberingRequired = true;
    int id = GetUnusedId(mConstraintMap);

    switch (mDimension)
    {
    case 1:
    case 2:
    case 3:
        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeRelativeHumidity(rNode,rValue));
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__,"Incorrect dimension of the structure.");
    }
    return id;
}

int NuTo::StructureBase::ConstraintLinearSetRelativeHumidityNode(int rIdent, double rValue)
{
    this->mNodeNumberingRequired = true;
    NodeBase* nodePtr;
    try
    {
        nodePtr = NodeGetNodePtr(rIdent);
    }
    catch (NuTo::MechanicsException &e)
    {
        e.AddMessage(__PRETTY_FUNCTION__,"Node with the given identifier could not be found.");
        throw;
    }
    catch (...)
    {
        throw MechanicsException(__PRETTY_FUNCTION__,"Node with the given identifier could not be found.");
    }

    return ConstraintLinearSetRelativeHumidityNode(nodePtr, rValue);
}


int NuTo::StructureBase::ConstraintLinearSetTemperatureNode(NodeBase* rNode, double rValue)
{
    this->mNodeNumberingRequired = true;

    int id = GetUnusedId(mConstraintMap);

    switch (mDimension)
    {
    case 1:
    case 2:
    case 3:
        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeTemperature(rNode,rValue));
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__,"Incorrect dimension of the structure.");
    }
    return id;
}

int NuTo::StructureBase::ConstraintLinearSetTemperatureNode(int rIdent, double rValue)
{
    this->mNodeNumberingRequired = true;
    NodeBase* nodePtr;
    try
    {
        nodePtr = NodeGetNodePtr(rIdent);
    }
    catch (NuTo::MechanicsException &e)
    {
        e.AddMessage(__PRETTY_FUNCTION__,"Node with the given identifier could not be found.");
        throw;
    }
    catch (...)
    {
        throw MechanicsException(__PRETTY_FUNCTION__,"Node with the given identifier could not be found.");
    }

    return ConstraintLinearSetTemperatureNode(nodePtr, rValue);
}


int NuTo::StructureBase::ConstraintLinearSetWaterVolumeFractionNode(NodeBase* rNode, double rValue)
{
    this->mNodeNumberingRequired = true;

    int id = GetUnusedId(mConstraintMap);

    switch (mDimension)
    {
    case 1:
        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeWaterVolumeFraction(rNode,rValue));
        break;
    case 2:
        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeWaterVolumeFraction(rNode,rValue));
        break;
    case 3:
        mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeWaterVolumeFraction(rNode,rValue));
        break;
    default:
        throw MechanicsException("[NuTo::StructureBase::ConstraintLinearSetWaterVolumeFractionNode] Incorrect dimension of the structure.");
    }
    return id;
}

int NuTo::StructureBase::ConstraintLinearSetWaterVolumeFractionNode(int rIdent, double rValue)
{
    this->mNodeNumberingRequired = true;
    NodeBase* nodePtr;
    try
    {
        nodePtr = NodeGetNodePtr(rIdent);
    }
    catch (NuTo::MechanicsException &e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstraintLinearSetWaterVolumeFractionNode] Node with the given identifier could not be found.");
        throw;
    }
    catch (...)
    {
        throw MechanicsException("[NuTo::StructureBase::ConstraintLinearSetWaterVolumeFractionNode] Node with the given identifier could not be found.");
    }

    return ConstraintLinearSetWaterVolumeFractionNode(nodePtr, rValue);
}


int NuTo::StructureBase::ConstraintLinearSetDisplacementNodeGroup(Group<NodeBase>* rGroup, const Eigen::VectorXd& rDirection, double rValue)
{
	this->mNodeNumberingRequired = true;
    int id = GetUnusedId(mConstraintMap);

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

int NuTo::StructureBase::ConstraintLinearSetRotationNodeGroup(Group<NodeBase>* rGroup, double rValue)
{
	this->mNodeNumberingRequired = true;
    int id = GetUnusedId(mConstraintMap);

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


int NuTo::StructureBase::ConstraintLinearSetDisplacementNodeGroup(int rGroupIdent, const Eigen::VectorXd& rDirection, double rValue)
{
	this->mNodeNumberingRequired = true;
    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::ConstraintSetDisplacementNodeGroup] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::eGroupId::Nodes)
        throw MechanicsException("[NuTo::Structure::ConstraintSetDisplacementNodeGroup] Group is not a node group.");
    Group<NodeBase> *nodeGroup = dynamic_cast<Group<NodeBase>*>(itGroup->second);
    assert(nodeGroup!=0);

    return ConstraintLinearSetDisplacementNodeGroup(nodeGroup,rDirection, rValue);
}

int NuTo::StructureBase::ConstraintLinearSetRotationNodeGroup(int rGroupIdent, double rValue)
{
	this->mNodeNumberingRequired = true;
    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::ConstraintLinearSetRotationNodeGroup] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::eGroupId::Nodes)
        throw MechanicsException("[NuTo::Structure::ConstraintLinearSetRotationNodeGroup] Group is not a node group.");
    Group<NodeBase> *nodeGroup = dynamic_cast<Group<NodeBase>*>(itGroup->second);
    assert(nodeGroup!=0);

    return ConstraintLinearSetRotationNodeGroup(nodeGroup, rValue);
}

int NuTo::StructureBase::ConstraintLinearSetTemperatureNodeGroup(int rGroupIdent, double rValue)
{
	this->mNodeNumberingRequired = true;
    boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupIdent);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::ConstraintLinearSetTemperatureNodeGroup] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::eGroupId::Nodes)
        throw MechanicsException("[NuTo::Structure::ConstraintLinearSetTemperatureNodeGroup] Group is not a node group.");
    Group<NodeBase> *nodeGroup = dynamic_cast<Group<NodeBase>*>(itGroup->second);
    assert(nodeGroup!=0);

    return ConstraintLinearSetTemperatureNodeGroup(nodeGroup, rValue);
}


int NuTo::StructureBase::ConstraintLinearSetTemperatureNodeGroup(Group<NodeBase>* rGroup, double rValue)
{
	this->mNodeNumberingRequired = true;
    int id = GetUnusedId(mConstraintMap);

    mConstraintMap.insert(id, new NuTo::ConstraintLinearNodeGroupTemperature(rGroup,rValue));
    return id;
}

int NuTo::StructureBase::ConstraintGetNumLinearConstraints(Node::eDof rDof) const
{
    int numLinearConstraints = 0;
    BOOST_FOREACH(auto itConstraint, mConstraintMap)
    {
        const auto& constraint = itConstraint.second;
        if (constraint->GetDofType() == rDof)
            numLinearConstraints += constraint->GetNumLinearConstraints();
    }
    return numLinearConstraints;
}


NuTo::BlockSparseMatrix NuTo::StructureBase::ConstraintGetConstraintMatrixBeforeGaussElimination() const
{
    if (mNodeNumberingRequired)
    {
        throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] build global numbering first");
    }

    BlockSparseMatrix constraintMatrix(GetDofStatus(), false);

    for (auto dofType : DofTypesGet())
    {
        int numLinearConstraints = ConstraintGetNumLinearConstraints(dofType);
        int curConstraintEquations = 0;

        constraintMatrix(dofType, dofType).Resize(numLinearConstraints,GetNumDofs(dofType));

        BOOST_FOREACH(auto itConstraint, mConstraintMap)
        {
            if (itConstraint->second->GetNumLinearConstraints()>0)
            {
                try
                {
                    if (itConstraint->second->GetDofType() == dofType)
                        itConstraint->second->AsConstraintLinear()->AddToConstraintMatrix(curConstraintEquations, constraintMatrix(dofType, dofType));
                }
                catch (MechanicsException& e)
                {
                    e.AddMessage(std::string("[") + __PRETTY_FUNCTION__ + "] mechanics exception while building constraint matrix for constraint with nonzero number of linear components.");
                    throw;
                }
                catch (...)
                {
                    throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] error building constraint matrix for constraint with nonzero number of linear components.");
                }
            }

        }

        if (curConstraintEquations != numLinearConstraints)
        {
            std::cout << "curConstraintEquations " << curConstraintEquations << std::endl;
            std::cout << "numConstraintEquations " << numLinearConstraints << std::endl;
            throw MechanicsException("[NuTo::StructureBase::ConstraintGetConstraintMatrixBeforeGaussEliminationDof] Internal error, there is something wrong with the constraint equations.");
        }
    }
    return constraintMatrix;
}

const NuTo::BlockFullVector<double>& NuTo::StructureBase::ConstraintGetRHSAfterGaussElimination() const
{
    if (mNodeNumberingRequired)
    {
        throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] build global numbering first");
    }
    return mConstraintRHS;
}


NuTo::BlockFullVector<double> NuTo::StructureBase::ConstraintGetRHSBeforeGaussElimination()
{
    if (mNodeNumberingRequired)
    {
        throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] build global numbering first");
    }

    NuTo::BlockFullVector<double> rhsBeforeGaussElimination(GetDofStatus());

    for (auto dof : DofTypesGet())
    {

        int numLinearConstraints = ConstraintGetNumLinearConstraints(dof);

        rhsBeforeGaussElimination[dof].resize(numLinearConstraints);

        //calculate the rhs vector of the constraint equations before the Gauss elimination
        int curConstraintEquations = 0;
        BOOST_FOREACH(auto itConstraint, mConstraintMap)
        {
            if (itConstraint.second->GetNumLinearConstraints()>0)
            {
                try
                {
                    if (itConstraint.second->GetDofType() == dof)
                        itConstraint.second->AsConstraintLinear()->GetRHS(curConstraintEquations, rhsBeforeGaussElimination[dof]);
                }
                catch (MechanicsException& e)
                {
                    e.AddMessage(std::string("[") + __PRETTY_FUNCTION__ + "] mechanics exception while building rhs vector after gauss elimination.");
                    throw;
                }
                catch (...)
                {
                    throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] mechanics exception while building rhs vector after gauss elimination.");
                }
            }
        }

        if (curConstraintEquations!=numLinearConstraints)
        {
            std::cout << "curConstraintEquations " << curConstraintEquations << std::endl;
            std::cout << "numConstraintEquations " << numLinearConstraints << std::endl;
            throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Internal error, there is something wrong with the constraint equations.");
        }
    }
    return rhsBeforeGaussElimination;
}

void NuTo::StructureBase::ConstraintUpdateRHSAfterGaussElimination()
{
    if (mNodeNumberingRequired)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "build global numbering first");
    }

    BlockFullVector<double> rhsBeforeGaussElimination = ConstraintGetRHSBeforeGaussElimination();

    if (mConstraintMappingRHS.GetNumColumns()!=rhsBeforeGaussElimination.GetNumRows())
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "here is something wrong in the implementation.");
    }

    //calculate the rhs vector of the constraint equations after the Gauss elimination using the mapping matrix
    mConstraintRHS = mConstraintMappingRHS * rhsBeforeGaussElimination;
}


void NuTo::StructureBase::ConstraintSetRHS(int rConstraintEquation, double rRHS)
{
    auto  it = mConstraintMap.find(rConstraintEquation);

    if (it == mConstraintMap.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Constraint equation does not exist.");

    it->second->SetRHS(rRHS);

    //since the rhs before Gauss elimination has changed, update the rhs after Gauss elimination using the mapping matrix
    ConstraintUpdateRHSAfterGaussElimination();
}

double NuTo::StructureBase::ConstraintGetRHS(int rConstraintEquation)const
{
    boost::ptr_map<int,ConstraintBase>::const_iterator it = mConstraintMap.find(rConstraintEquation);
    if (it==mConstraintMap.end())
    {
    	throw MechanicsException(__PRETTY_FUNCTION__, "Constraint equation does not exist.");
    }
    return it->second->GetRHS();
}

// create a constraint equation
int NuTo::StructureBase::ConstraintLinearEquationCreate(int rNode, const std::string& rDof, double rCoefficient, double rRHS)
{
	this->mNodeNumberingRequired = true;
    int id = GetUnusedId(mConstraintMap);

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
        NuTo::Node::eDof dofType;
        int dofComponent;
        this->ConstraintEquationGetDofInformationFromString(rDof, dofType, dofComponent);

        // create constraint equation
        this->ConstraintLinearEquationCreate(rConstraint, rNode, dofType, dofComponent, rCoefficient, rRHS);
    }
    catch(NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstraintEquationCreate] error creating constraint equation");
        throw;
    }
}

// create a constraint equation
void NuTo::StructureBase::ConstraintLinearEquationCreate(int rConstraint, int rNode, NuTo::Node::eDof rDofType, int rDofComponent, double rCoefficient, double rRHS)
{
	this->mNodeNumberingRequired = true;
    // check if constraint equation already exists
    if(this->mConstraintMap.find(rConstraint) != this->mConstraintMap.end())
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "constraint equation already exist.");
    }

    try
    {
        // get node pointer
        NodeBase* nodePtr = this->NodeGetNodePtr(rNode);

        // create new constraint equation term
        ConstraintBase* constraintPtr = new ConstraintLinearEquation(nodePtr, rDofType, rDofComponent, rCoefficient, rRHS);

        // insert constraint equation into map
        this->ConstraintAdd(rConstraint, constraintPtr);
    }
    catch(NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstraintEquationCreate] error creating constraint equation");
        throw;
    }
}

void NuTo::StructureBase::ConstraintLinearEquationNodeToElementCreate(int rNode, int rElementGroup,
                                                                      NuTo::Node::eDof rDofType,
                                                                      const double rTolerance,
                                                                      Eigen::Vector3d rNodeCoordOffset)
{
    this->mNodeNumberingRequired = true;

    const int dim = GetDimension();

    Eigen::VectorXd queryNodeCoords = NodeGetNodePtr(rNode)->Get(Node::eDof::COORDINATES);
    queryNodeCoords = queryNodeCoords +  rNodeCoordOffset.head(dim);


    std::vector<int> elementGroupIds = GroupGetMemberIds(rElementGroup);

    ElementBase* elementPtr = nullptr;
    Eigen::VectorXd elementNaturalNodeCoords;
    bool nodeInElement = false;
    for (auto const& eleId : elementGroupIds) {
        elementPtr = ElementGetElementPtr(eleId);

        // Coordinate interpolation must be linear so the shape function derivatives are constant!
        assert(elementPtr->GetInterpolationType().Get(Node::eDof::COORDINATES).GetTypeOrder() ==
               Interpolation::eTypeOrder::EQUIDISTANT1);
        const Eigen::MatrixXd &derivativeShapeFunctionsGeometryNatural = elementPtr->GetInterpolationType().Get(
                Node::eDof::COORDINATES).GetDerivativeShapeFunctionsNatural(0);

        // real coordinates of every node in rElement
        Eigen::VectorXd elementNodeCoords = elementPtr->ExtractNodeValues(NuTo::Node::eDof::COORDINATES);

        switch (mDimension) {
            case 2: {
                Eigen::Matrix2d invJacobian = elementPtr->AsContinuumElement2D().CalculateJacobian(
                        derivativeShapeFunctionsGeometryNatural, elementNodeCoords).inverse();

                elementNaturalNodeCoords = invJacobian * (queryNodeCoords - elementNodeCoords.head(2));
            }
                break;
            case 3: {
                Eigen::Matrix3d invJacobian = elementPtr->AsContinuumElement3D().CalculateJacobian(
                        derivativeShapeFunctionsGeometryNatural, elementNodeCoords).inverse();

                elementNaturalNodeCoords = invJacobian * (queryNodeCoords - elementNodeCoords.head(3));


            }
                break;

            default:
                throw NuTo::MechanicsException(
                        std::string(__PRETTY_FUNCTION__) + ": \t Only implemented for 2D and 3D");

        }


        if ( (elementNaturalNodeCoords.array() > -rTolerance).all() and elementNaturalNodeCoords.sum() <= 1. + rTolerance)
        {
            nodeInElement = true;
            break;
        }

    }

    if (not nodeInElement)
    {
        GetLogger() << "Natural node coordinates: \n" << elementNaturalNodeCoords << "\n";
        throw MechanicsException(__PRETTY_FUNCTION__, "Node is not inside any element.");
    }

    auto shapeFunctions = elementPtr->GetInterpolationType().Get(Node::eDof::DISPLACEMENTS).CalculateShapeFunctions(elementNaturalNodeCoords);

    //find unused integer id
    std::vector<int> unusedId(dim);
    for (int iDim = 0; iDim < dim; ++iDim)
    {
        unusedId[iDim] = GetUnusedId(mConstraintMap);
        ConstraintLinearEquationCreate(unusedId[iDim], rNode, NuTo::Node::eDof::DISPLACEMENTS, iDim, 1.0, 0.0);
    }


    for (int iNode = 0; iNode < shapeFunctions.rows(); ++iNode)
    {
        int localNodeId = elementPtr->GetInterpolationType().Get(Node::eDof::DISPLACEMENTS).GetNodeIndex(iNode);
        int globalNodeId = NodeGetId(elementPtr->GetNode(localNodeId, Node::eDof::DISPLACEMENTS));
//        std::cout << "globalNodeId \t" << globalNodeId << std::endl;
        double coefficient = -shapeFunctions(iNode, 0);

        for (int iDim = 0; iDim < dim; ++iDim)
            ConstraintLinearEquationAddTerm(unusedId[iDim], globalNodeId, Node::eDof::DISPLACEMENTS, iDim, coefficient);

    }


}

// add a term to a constraint equation
void NuTo::StructureBase::ConstraintLinearEquationAddTerm(int rConstraint, int rNode, const std::string& rDof, double rCoefficient)
{
    try
    {
        // convert dof string
        NuTo::Node::eDof dofType;
        int dofComponent;
        this->ConstraintEquationGetDofInformationFromString(rDof, dofType, dofComponent);

        // add term to constraint equation
        this->ConstraintLinearEquationAddTerm(rConstraint, rNode, dofType, dofComponent, rCoefficient);
    }
    catch(NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstraintEquationAddTerm] error adding term to constraint equation");
        throw;
    }
}

// add a term to a constraint equation
void NuTo::StructureBase::ConstraintLinearEquationAddTerm(int rConstraint, int rNode, NuTo::Node::eDof rDofType, int rDofComponent, double rCoefficient)
{
	this->mNodeNumberingRequired = true;
    // get iterator
    boost::ptr_map<int,ConstraintBase>::iterator it = this->mConstraintMap.find(rConstraint);
    if(it == this->mConstraintMap.end())
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "constraint equation does not exist.");
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
        throw;
    }
}

void NuTo::StructureBase::ConstraintEquationGetDofInformationFromString(const std::string& rDof, NuTo::Node::eDof& rDofType, int& rDofComponent)
{
    // convert string to upper-case
    std::string dofString;
    std::transform(rDof.begin(), rDof.end(), std::back_inserter(dofString), (int(*)(int)) toupper);
    if(dofString == "X_DISPLACEMENT")
    {
        rDofType = NuTo::Node::eDof::DISPLACEMENTS;
        rDofComponent = 0;
    }
    else if(dofString == "Y_DISPLACEMENT")
    {
        if(this->mDimension < 2)
        {
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "y-displacement dofs are only available in 2D or 3D structures.");
        }
        rDofType = NuTo::Node::eDof::DISPLACEMENTS;
        rDofComponent = 1;
    }
    else if(dofString == "Z_DISPLACEMENT")
    {
        if(this->mDimension < 3)
        {
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "z-displacement dofs are only available in 3D structures.");
        }
        rDofType = NuTo::Node::eDof::DISPLACEMENTS;
        rDofComponent = 2;
    }
    else if(dofString == "X_ROTATION")
    {
        if(this->mDimension < 3)
        {
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "x-rotation dofs are only available in 3D structures.");
        }
        rDofType = NuTo::Node::eDof::ROTATIONS;
        rDofComponent = 0;
    }
    else if(dofString == "Y_ROTATION")
    {
        if(this->mDimension < 3)
        {
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "y-rotation dofs are only available in 3D structures.");
        }
        rDofType = NuTo::Node::eDof::ROTATIONS;
        rDofComponent = 1;
    }
    else if(dofString == "Z_ROTATION")
    {
        rDofType = NuTo::Node::eDof::ROTATIONS;
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
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "z-rotation dofs are only available in 2D and 3D structures.");
        }
    }
    else if(dofString == "TEMPERATURE")
    {
        rDofType = NuTo::Node::eDof::TEMPERATURE;
        rDofComponent = 0;
    }
    else if(dofString == "NONLOCALEQSTRAIN")
    {
        rDofType = NuTo::Node::eDof::NONLOCALEQSTRAIN;
        rDofComponent = 0;
    }
    else
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "invalid dof string.");
    }
}

int NuTo::StructureBase::ConstraintLinearDisplacementsSetPeriodic2D(double rAngle, Eigen::MatrixXd rStrain,
        double rRadiusToCrackWithoutConstraints,
        int rNodeGroupUpperId, int rNodeGroupLowerId, int rNodeGroupLeftId, int rNodeGroupRightId)
{
    NuTo::Timer timer(__FUNCTION__, GetShowTime(), GetLogger());
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
    if (itGroup->second->GetType()!=NuTo::eGroupId::Nodes)
        throw MechanicsException("[NuTo::Structure::ConstraintSetPeriodicBoundaryConditions2D] Group of upper nodes is not a node group.");
    Group<NodeBase> *nodeGroupUpperPtr = dynamic_cast<Group<NodeBase>*>(itGroup->second);
    assert(nodeGroupUpperPtr!=0);

     itGroup = mGroupMap.find(rNodeGroupLowerId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::ConstraintSetPeriodicBoundaryConditions2D] Group of lower nodes with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::eGroupId::Nodes)
        throw MechanicsException("[NuTo::Structure::ConstraintSetPeriodicBoundaryConditions2D] Group of lower nodes is not a node group.");
    Group<NodeBase> *nodeGroupLowerPtr = dynamic_cast<Group<NodeBase>*>(itGroup->second);
    assert(nodeGroupLowerPtr!=0);

    itGroup = mGroupMap.find(rNodeGroupLeftId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::ConstraintSetPeriodicBoundaryConditions2D] Group of left nodes with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::eGroupId::Nodes)
        throw MechanicsException("[NuTo::Structure::ConstraintSetPeriodicBoundaryConditions2D] Group of left nodes is not a node group.");
    Group<NodeBase> *nodeGroupLeftPtr = dynamic_cast<Group<NodeBase>*>(itGroup->second);
    assert(nodeGroupLeftPtr!=0);

    itGroup = mGroupMap.find(rNodeGroupRightId);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::ConstraintSetPeriodicBoundaryConditions2D] Group of right nodes with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::eGroupId::Nodes)
        throw MechanicsException("[NuTo::Structure::ConstraintSetPeriodicBoundaryConditions2D] Group of right nodes is not a node group.");
    Group<NodeBase> *nodeGroupRightPtr = dynamic_cast<Group<NodeBase>*>(itGroup->second);
    assert(nodeGroupRightPtr!=0);

    try
    {
        // create new constraint equation term
        if (rStrain.rows()!=3 || rStrain.cols()!=1)
            throw MechanicsException("[NuTo::ConstraintNodeDisplacementsPeriodic2D::ConstraintNodeDisplacementsPeriodic2D] the strain is matrix (3,1) with (e_xx, e_yy, gamma_xy)");

        EngineeringStrain<2> engineeringStrain;
        engineeringStrain.AsVector() = rStrain;
        ConstraintBase* constraintPtr = new NuTo::ConstraintLinearDisplacementsPeriodic2D(this, rAngle, engineeringStrain, rRadiusToCrackWithoutConstraints,
                   nodeGroupUpperPtr, nodeGroupLowerPtr, nodeGroupLeftPtr, nodeGroupRightPtr);

        // insert constraint equation into map
        this->ConstraintAdd(id, constraintPtr);
    }
    catch(NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::ConstraintSetPeriodicBoundaryConditions2D] error creating periodic boundary conditions in 2D");
        throw;
    }
    return id;
}


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

void NuTo::StructureBase::ConstraintDelete(int rConstraintId)
{
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(rConstraintId);
    if (it==mConstraintMap.end())
    {
        throw MechanicsException("[NuTo::StructureBase::ConstraintDelete] Constraint equation does not exist.");
    }
    mConstraintMap.erase(it);
    this->mNodeNumberingRequired = true;
}


int NuTo::StructureBase::ConstraintGetNumLinearConstraints(std::string rDof) const
{
    return ConstraintGetNumLinearConstraints(Node::DofToEnum(rDof));
}

NuTo::ConstraintBase* NuTo::StructureBase::ConstraintRelease(int rConstraintId)
{
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(rConstraintId);
    if (it==mConstraintMap.end())
    {
        throw MechanicsException("[NuTo::StructureBase::ConstraintRelease] Constraint equation does not exist.");
    }
    boost::ptr_map<int,ConstraintBase>::auto_type ptr = mConstraintMap.release(it);
    this->mNodeNumberingRequired = true;
    return ptr.release();
}


int NuTo::StructureBase::ConstraintLinearSetNode(NuTo::Node::eDof rDOFType, NuTo::NodeBase *rNode, double rValue)
{
    Eigen::MatrixXd direction(1,1);
    return ConstraintLinearSetNode(rDOFType,rNode,direction,rValue);
}

int NuTo::StructureBase::ConstraintLinearSetNode(NuTo::Node::eDof rDOFType,
                                                 NuTo::NodeBase *rNode,
                                                 const Eigen::VectorXd& rDirection,
                                                 double rValue)
{
    switch(rDOFType)
    {
    case Node::eDof::DISPLACEMENTS:
        return ConstraintLinearSetDisplacementNode(rNode,rDirection,rValue);

    case Node::eDof::RELATIVEHUMIDITY:
        return ConstraintLinearSetRelativeHumidityNode(rNode,rValue);

    case Node::eDof::TEMPERATURE:
        return ConstraintLinearSetTemperatureNode(rNode, rValue);

    default:
        throw MechanicsException(__PRETTY_FUNCTION__,std::string("not implemented for dof ")+Node::DofToString(rDOFType));
    }
}

void NuTo::StructureBase::ConstraintAdd(int rConstraintId, NuTo::ConstraintBase* rConstraint)
{
    boost::ptr_map<int,ConstraintBase>::iterator it = mConstraintMap.find(rConstraintId);
    if (it!=mConstraintMap.end())
    	throw MechanicsException("[NuTo::StructureBase::ConstraintAdd] Id already exists in constraint map.");
    mConstraintMap.insert(rConstraintId, rConstraint);
    mNodeNumberingRequired = true;
}
