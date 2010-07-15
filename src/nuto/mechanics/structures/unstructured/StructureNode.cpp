#include <boost/tokenizer.hpp>

#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/nodes/NodeCoordinates1D.h"
#include "nuto/mechanics/nodes/NodeCoordinates2D.h"
#include "nuto/mechanics/nodes/NodeCoordinates3D.h"
#include "nuto/mechanics/nodes/NodeCoordinatesDisplacements1D.h"
#include "nuto/mechanics/nodes/NodeCoordinatesDisplacements2D.h"
#include "nuto/mechanics/nodes/NodeCoordinatesDisplacements3D.h"
#include "nuto/mechanics/nodes/NodeCoordinatesDisplacementsRotations2D.h"
#include "nuto/mechanics/nodes/NodeCoordinatesDisplacementsRotations3D.h"
#include "nuto/mechanics/nodes/NodeCoordinatesTemperature1D.h"
#include "nuto/mechanics/nodes/NodeCoordinatesTemperature2D.h"
#include "nuto/mechanics/nodes/NodeCoordinatesTemperature3D.h"
#include "nuto/mechanics/nodes/NodeCoordinatesDisplacementsNonlocalData2D.h"
#include "nuto/mechanics/nodes/NodeCoordinatesDisplacementsNonlocalData3D.h"

//! @brief returns the number of nodes
//! @return number of nodes
int NuTo::Structure::GetNumNodes() const
{
    return mNodeMap.size();
}

//! @brief a reference to a node
//! @param identifier
//! @return reference to a node
NuTo::NodeBase* NuTo::Structure::NodeGetNodePtr(int rIdent)
{
    boost::ptr_map<int,NodeBase>::iterator it = mNodeMap.find(rIdent);
    if (it!=mNodeMap.end())
        return it->second;
    else
    {
    	std::stringstream out;
    	out << rIdent;
    	throw MechanicsException("[NuTo::Structure::NodeGetNodePtr] Node with id " + out.str() +" does not exist.");
    }
}

//! @brief a reference to a node
//! @param identifier
//! @return reference to a node
const NuTo::NodeBase* NuTo::Structure::NodeGetNodePtr(int rIdent)const
{
    boost::ptr_map<int,NodeBase>::const_iterator it = mNodeMap.find(rIdent);
    if (it!=mNodeMap.end())
        return it->second;
    else
    {
    	std::stringstream out;
    	out << rIdent;
    	throw MechanicsException("[NuTo::Structure::NodeGetNodePtr] Node with id " + out.str() +" does not exist.");
    }
}

//! @brief gives the identifier of a node
//! @param reference to a node
//! @return identifier
int NuTo::Structure::NodeGetId(const NodeBase* rNode)const
{
    for (boost::ptr_map<int,NodeBase>::const_iterator
            it = mNodeMap.begin(); it!= mNodeMap.end(); it++)
    {
        if (it->second==rNode)
            return it->first;
    }
    throw MechanicsException("[NuTo::Structure::GetNodeId] Node does not exist.");
}

//! @brief info about the elements in the Structure
void NuTo::Structure::NodeInfo(int rVerboseLevel)const
{
    std::cout<<"number of nodes   : " << mNodeMap.size() <<std::endl;
    if (rVerboseLevel>3)
    {
    	std::cout << "\t\tnodes :" <<std::endl;
    	for (boost::ptr_map<int,NodeBase>::const_iterator it = mNodeMap.begin(); it!= mNodeMap.end(); it++)
    	{
			std::cout << "\t\t" << it->first;
			if (rVerboseLevel>4)
			{
				std::cout << "\t:";
				for(unsigned short iDof=0; iDof<it->second->GetNumCoordinates(); ++iDof)
					std::cout << "\t" << it->second->GetCoordinate(iDof);
			}
			std::cout << std::endl;
    	}

    }
}

// creates a node
int NuTo::Structure::NodeCreate(std::string rDOFs, NuTo::FullMatrix<double>& rCoordinates)
{

    //find unused integer id
    int id(mNodeMap.size());
    boost::ptr_map<int,NodeBase>::iterator it = mNodeMap.find(id);
    while (it!=mNodeMap.end())
    {
        id++;
        it = mNodeMap.find(id);
    }
    this->NodeCreate(id, rDOFs, rCoordinates);

    //return int identifier of the new node
    return id;
}

void NuTo::Structure::NodeCreate(int rNodeNumber, std::string rDOFs, NuTo::FullMatrix<double>& rCoordinates)
{
	// check node number
	boost::ptr_map<int,NodeBase>::iterator it = this->mNodeMap.find(rNodeNumber);
	if(it != this->mNodeMap.end())
	{
		throw MechanicsException("[NuTo::Structure::NodeCreate] Node already exists.");
	}

    // transform string to uppercase
    std::transform(rDOFs.begin(), rDOFs.end(), rDOFs.begin(), toupper);

    // check all values
    int attributes(1 << Node::COORDINATES);
    //! bit 0 : coordinates
    //! bit 1 : displacements
    //! bit 2 : rotations
    //! bit 3 : temperatures
    //! bit 4 : nonlocal data
    boost::tokenizer<> tok(rDOFs);
    for (boost::tokenizer<>::iterator beg=tok.begin(); beg!=tok.end(); ++beg)
    {
        if (*beg=="DISPLACEMENTS")
            attributes = attributes | 1 << Node::DISPLACEMENTS;
        if (*beg=="ROTATIONS")
            attributes = attributes | 1 << Node::ROTATIONS;
        if (*beg=="TEMPERATURES")
            attributes = attributes | 1 << Node::TEMPERATURES;
        if (*beg=="NONLOCALDATA")
            attributes = attributes | 1 << Node::NONLOCALDATA;
    }
    if (rCoordinates.GetNumRows()!=mDimension || rCoordinates.GetNumColumns()!=1)
        throw MechanicsException("[NuTo::Structure::NodeCreate]\
 Dimension of the coordinate vector does not fit the dimension of the structure");

    NodeBase* nodePtr;
    switch (attributes)
    {
        // the << shifts the 1 bitwise to the left, so 1<<n = 2^n
        // it actually sets the n-th bit (from the right) to 1, and all the other to zero
    case (1 << Node::COORDINATES):
        // reference node having only coordinates
        switch (mDimension)
        {
        case 1:
        	nodePtr = new NuTo::NodeCoordinates1D();
            break;
        case 2:
        	nodePtr = new NuTo::NodeCoordinates2D();
            break;
        case 3:
        	nodePtr = new NuTo::NodeCoordinates3D();
            break;
        default:
            throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
        }
        break;
    case (1 << Node::COORDINATES) | (1 << Node::DISPLACEMENTS):
        // coordinates and displacements
        switch (mDimension)
        {
        case 1:
        	nodePtr = new NuTo::NodeCoordinatesDisplacements1D();
            break;
        case 2:
        	nodePtr = new NuTo::NodeCoordinatesDisplacements2D();
            break;
        case 3:
        	nodePtr = new NuTo::NodeCoordinatesDisplacements3D();
            break;
        default:
            throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
        }
        break;
    case (1 << Node::COORDINATES) | (1 << Node::DISPLACEMENTS) | (1 << Node::ROTATIONS):
        // coordinates and displacements
        switch (mDimension)
        {
        case 1:
            throw MechanicsException("[NuTo::Structure::NodeCreate] A 1D node with a rotational DOF is not allowed.");
            break;
        case 2:
        	nodePtr = new NuTo::NodeCoordinatesDisplacementsRotations2D();
            break;
        case 3:
        	nodePtr = new NuTo::NodeCoordinatesDisplacementsRotations3D();
            break;
        default:
            throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
        }
        break;
    case (1 << Node::COORDINATES) | (1 << Node::TEMPERATURES):
        // coordinates and temperatures
        switch (mDimension)
        {
        case 1:
        	nodePtr = new NuTo::NodeCoordinatesTemperature1D();
            break;
        case 2:
        	nodePtr = new NuTo::NodeCoordinatesTemperature2D();
            break;
        case 3:
        	nodePtr = new NuTo::NodeCoordinatesTemperature3D();
            break;
        default:
            throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
        }
        break;
	case (1 << Node::COORDINATES) | (1 << Node::DISPLACEMENTS) | (1 << Node::NONLOCALDATA):
		// coordinates, displacements and nonlocal data
		switch (mDimension)
		{
		case 1:
			throw MechanicsException("[NuTo::Structure::NodeCreate] A 1D node with nonlocal data is not implemented");
			break;
		case 2:
			nodePtr = new NuTo::NodeCoordinatesDisplacementsNonlocalData2D();
			break;
		case 3:
			nodePtr = new NuTo::NodeCoordinatesDisplacementsNonlocalData3D();
			break;
		default:
			throw MechanicsException("[NuTo::Structure::NodeCreate] Dimension of the structure is not valid.");
		}
		break;
    default:
        throw MechanicsException("[NuTo::Structure::NodeCreate] This combination of attributes is not implemented (just add in the source file the relevant combination).");
    }

    //add coordinates
    try
    {
        switch (mDimension)
        {
        case 1:
            nodePtr->SetCoordinates1D(rCoordinates.mEigenMatrix.data());
            break;
        case 2:
            nodePtr->SetCoordinates2D(rCoordinates.mEigenMatrix.data());
            break;
        case 3:
            nodePtr->SetCoordinates3D(rCoordinates.mEigenMatrix.data());
            break;
        case 0:
            throw MechanicsException("[NuTo::StructureBase::NodeCreate] Node has no displacements.");
            break;
        }
    }
    catch (std::bad_cast & b)
    {
        throw MechanicsException("[NuTo::StructureBase::NodeCreate] Node has no coordinates or its dimension is not equivalent to the dimension of the input matrix.");
    }
    catch (NuTo::MechanicsException & b)
    {
        b.AddMessage("[NuTo::StructureBase::NodeCreate] Error setting coordinates.");
        throw b;
    }
    catch (...)
    {
        throw MechanicsException("[NuTo::StructureBase::NodeCreate] Error setting coordinates of node (unspecified exception).");
    }
    // add node to map
    this->mNodeMap.insert(rNodeNumber, nodePtr);

    //renumbering of dofs for global matrices required
    this->mNodeNumberingRequired  = true;
}

//! @brief create multiple nodes
//! @param reference to a FullMatrix containing the node coordinates (row->coordinate; col->nodes)
//! @return a FullMatrix<int> with the identifiers
NuTo::FullMatrix<int> NuTo::Structure::NodesCreate(std::string rDOFs, NuTo::FullMatrix<double>& rCoordinates)
{
	std::vector<int> idVec;
	/// go through the nodes
	for(size_t i=0 ; i<(size_t)rCoordinates.GetNumColumns(); ++i)
	{
		NuTo::FullMatrix<double> coordinate(rCoordinates.GetColumn(i));
		idVec.push_back(this->NodeCreate(rDOFs, coordinate));
	}

    //return int identifiers of the new nodes as FullMatrix
	NuTo::FullMatrix<int> ids(idVec);

    return ids;
}

//! @brief Deletes a nodes
//! @param rNodeNumber identifier for the node
void NuTo::Structure::NodeDelete(const int rNodeNumber)
{
    boost::ptr_map<int,NodeBase>::iterator itNode = mNodeMap.find(rNodeNumber);
    if (itNode==this->mNodeMap.end())
    {
        throw MechanicsException("[NuTo::Structure::NodeDelete] Node with the given id does not exist.");
    }
    else
    {
    	/// Search for node in elements: using a loop over all elements
        for(boost::ptr_map<int,ElementBase>::const_iterator elemIt=mElementMap.begin(); elemIt!=mElementMap.end(); ++elemIt){
        	// loop over all element nodes
        	for(unsigned short iNode=0; iNode<elemIt->second->GetNumNodes(); ++iNode){
            	// if the id of the node to be deleted is found, throw an error
        		if(rNodeNumber == this->NodeGetId(elemIt->second->GetNode(iNode)) ){
        	        throw MechanicsException( "[NuTo::Structure::NodeDelete] Node with given id is already in use.");
        		}
        	}
        }

    	//! Search for node in groups
        //! using a loop over all groups
        //! remove the entry from all groups
        for(boost::ptr_map<int,GroupBase>::iterator groupIt=mGroupMap.begin();groupIt!=mGroupMap.end(); ++groupIt){
        	if(groupIt->second->GetType()==NuTo::Groups::Nodes){
        		if(groupIt->second->Contain(itNode->second)){
        			groupIt->second->RemoveMember(itNode->second);
        		}
        	}
        }

        // delete element from map
        this->mNodeMap.erase(itNode);

        this->mNodeNumberingRequired = true;
    }

}

//! @brief number the dofs in the structure
void NuTo::Structure::NodeBuildGlobalDofs()
{
    // build initial node numbering
    this->mNumDofs = 0;
    for (boost::ptr_map<int,NodeBase>::iterator it = mNodeMap.begin(); it!= mNodeMap.end(); it++)
    {
        it->second->SetGlobalDofs(this->mNumDofs);
    }

    // build constraint matrix
    this->mNodeNumberingRequired = false;
    this->ConstraintGetConstraintMatrix(this->mConstraintMatrix, this->mConstraintRHS);
    this->mNodeNumberingRequired = true;

    // perform gauss algorithm
    std::vector<int> mappingInitialToNewOrdering;
    std::vector<int> mappingNewToInitialOrdering;
    this->mConstraintMatrix.Gauss(this->mConstraintRHS, mappingNewToInitialOrdering, mappingInitialToNewOrdering);

    // move dependent dofs at the end
    // Warning!!! after this loop mappingNewToInitialOrdering is no longer valid !!!
    unsigned int numDependentDofs = this->mConstraintMatrix.GetNumRows();
    this->mNumActiveDofs = this->mNumDofs - numDependentDofs;
    std::vector<int> tmpMapping;
    for (unsigned int dependentDofCount = 0; dependentDofCount < numDependentDofs; dependentDofCount++)
    {
        tmpMapping.push_back(this->mNumActiveDofs + dependentDofCount);
        mappingInitialToNewOrdering[mappingNewToInitialOrdering[dependentDofCount]] += this->mNumActiveDofs;
    }
    for (int activeDofCount = numDependentDofs; activeDofCount < this->mNumDofs; activeDofCount++)
    {
        tmpMapping.push_back(activeDofCount - numDependentDofs);
        mappingInitialToNewOrdering[mappingNewToInitialOrdering[activeDofCount]] -= numDependentDofs;
    }
    mappingNewToInitialOrdering.clear();

    // reorder columns
    this->mConstraintMatrix.ReorderColumns(tmpMapping);

    // remove columns of dependent dofs
    // check if the submatrix which is removed is a diagonal matrix
    const std::vector<int>& constraintMatrixRowIndex = this->mConstraintMatrix.GetRowIndex();
    const std::vector<int>& constraintMatrixColumns = this->mConstraintMatrix.GetColumns();
    for (int row = 0; row < this->mConstraintMatrix.GetNumRows(); row++)
    {
        for (int pos = constraintMatrixRowIndex[row]; pos < constraintMatrixRowIndex[row + 1]; pos++)
        {
            if ((constraintMatrixColumns[pos] > this->mNumActiveDofs) && (constraintMatrixColumns[pos] != row + this->mNumActiveDofs))
            {
                throw MechanicsException("[NuTo::Structure::NodeBuildGlobalDofs] invalid matrix structure.");
            }
        }
    }
    this->mConstraintMatrix.RemoveLastColumns(numDependentDofs);

    // renumber dofs
    for (boost::ptr_map<int,NodeBase>::iterator it = mNodeMap.begin(); it!= mNodeMap.end(); it++)
    {
        it->second->RenumberGlobalDofs(mappingInitialToNewOrdering);
    }

    mNodeNumberingRequired = false;
}

// merge dof values
void NuTo::Structure::NodeMergeActiveDofValues(const FullMatrix<double>& rActiveDofValues)
{
    if (this->mNodeNumberingRequired)
    {
        throw MechanicsException("[NuTo::Structure::NodeMergeActiceDofValues] a valid dof numbering was not found (build dof numbering using NodeBuildGlobalDofs).");
    }
    if ((rActiveDofValues.GetNumRows() != this->mNumActiveDofs) || rActiveDofValues.GetNumColumns() != 1)
    {
        throw MechanicsException("[NuTo::Structure::NodeMergeActiceDofValues] invalid dimension of input object (number of active dofs,1).");
    }
	this->mUpdateTmpStaticDataRequired=true;

    // calculate dependent dof values
    FullMatrix<double> dependentDofValues = this->mConstraintRHS - this->mConstraintMatrix * rActiveDofValues;

    // write dof values to the nodes
    for (boost::ptr_map<int,NodeBase>::iterator it = mNodeMap.begin(); it!= mNodeMap.end(); it++)
    {
        it->second->SetGlobalDofValues(rActiveDofValues, dependentDofValues);
    }
}

// extract dof values
void NuTo::Structure::NodeExtractDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
    if (this->mNodeNumberingRequired)
    {
        throw MechanicsException("[NuTo::Structure::NodeExtractDofValues] a valid dof numbering was not found (build dof numbering using NodeBuildGlobalDofs).");
    }
    rActiveDofValues.Resize(this->mNumActiveDofs,1);
    rDependentDofValues.Resize(this->mNumDofs - this->mNumActiveDofs,1);

    // extract dof values from nodes
    for (boost::ptr_map<int,NodeBase>::const_iterator it = mNodeMap.begin(); it!= mNodeMap.end(); it++)
    {
        it->second->GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }
}

void NuTo::Structure::NodeGetInternalForce(const NodeBase* rNode, NuTo::FullMatrix<double>& rNodeForce)const
{
    // define variables storing the element contribution outside the loop
    NuTo::FullMatrix<double> elementVector;
    std::vector<int> elementVectorGlobalDofs;

    //iterate over the elements, check if the node is part of the element
	for (boost::ptr_map<int, ElementBase>::const_iterator itElement = mElementMap.begin(); itElement!=mElementMap.end(); itElement++)
	{
		const ElementBase* elem_ptr(itElement->second);
		for (int local_node=0; local_node<elem_ptr->GetNumNodes(); local_node++)
		{
			if (elem_ptr->GetNode(local_node)==rNode)
			{
				//node has been found, build internal force vector of element and add relevant parts to the global vector
				elem_ptr->CalculateGradientInternalPotential(elementVector, elementVectorGlobalDofs);
				if (rNode->GetNumDisplacements()==0)
				{
				    std::stringstream message;
				    message << "[NuTo::Structure::NodeGetInternalForce] Node" << NodeGetId(rNode) << " has no displacement degrees of freedom.";
					throw MechanicsException(message.str());
				}

				for (int theDisp=0; theDisp<rNode->GetNumDisplacements(); theDisp++)
				{
					// find global dof number
					int theDOF = rNode->GetDofDisplacement(theDisp);

					for (unsigned int local_dof=0; local_dof<elementVectorGlobalDofs.size(); local_dof++)
					{
						if (elementVectorGlobalDofs[local_dof]==theDOF)
						{
							// add node force
							rNodeForce(theDisp,0)+=elementVector(local_dof,0);
						}
					}
				}
			}
		}
	}
}

//! @brief calculates the internal force vector for a given node
//! @param rNodeId node id
//! @param rNodeForce return value
void NuTo::Structure::NodeGetInternalForce(int rNodeId, NuTo::FullMatrix<double>& rNodeForce)const
{
	const NodeBase* node_ptr(NodeGetNodePtr(rNodeId));

	// initialize matrix
	rNodeForce.Resize(node_ptr->GetNumDisplacements(),1);

	// add internal force
	NodeGetInternalForce(node_ptr, rNodeForce);
}

//! @brief calculates the internal force vector for a given node group
//! @param rGroupIdent group identifier
//! @param rNodeForce return value
void NuTo::Structure::NodeGroupGetInternalForce(int rIdentGroup, NuTo::FullMatrix<double>& rNodeForce)const
{
    boost::ptr_map<int,GroupBase>::const_iterator itGroup = mGroupMap.find(rIdentGroup);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::Structure::NodeGroupGetInternalForce] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=Groups::Nodes)
        throw MechanicsException("[NuTo::Structure::NodeGroupGetInternalForce] Group is not a node group.");

    const Group<NodeBase>* nodeGroup = dynamic_cast<const Group<NodeBase>* >(itGroup->second);

	// initialize matrix with first node
	rNodeForce.Resize((*(nodeGroup->begin()))->GetNumDisplacements(),1);
    for (Group<NodeBase>::iterator itNode=nodeGroup->begin(); itNode!=nodeGroup->end();itNode++)
    {
        try
        {
        	NodeGetInternalForce(*itNode, rNodeForce);
        }
        catch(NuTo::MechanicsException e)
        {
            std::stringstream ss;
            ss << NodeGetId(*itNode);
            e.AddMessage("[NuTo::Structure::NodeGroupGetInternalForce] Error calculating internal force for node "
            	+ ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            ss << NodeGetId(*itNode);
       	    throw NuTo::MechanicsException
        	   ("[NuTo::Structure::NodeGroupGetInternalForce] Error calculating internal force for node " + ss.str() + ".");
        }
    }

}

