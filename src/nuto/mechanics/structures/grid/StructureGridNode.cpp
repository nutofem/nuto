#include <boost/tokenizer.hpp>
#include <sstream>

#include "nuto/mechanics/structures/grid/StructureGrid.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/nodes/NodeGrid1D.h"
#include "nuto/mechanics/nodes/NodeGrid2D.h"
#include "nuto/mechanics/nodes/NodeGrid3D.h"
#include "nuto/mechanics/nodes/NodeGridDisplacements1D.h"
#include "nuto/mechanics/nodes/NodeGridDisplacements2D.h"
#include "nuto/mechanics/nodes/NodeGridDisplacements3D.h"

//! @brief returns the number of nodes
//! @return number of nodes
int NuTo::StructureGrid::GetNumNodes() const
{
    return mNodeVec.size();
}

//! @brief a reference to a node
//! @param identifier
//! @return reference to a node
NuTo::NodeBase* NuTo::StructureGrid::NodeGetNodePtr(int rNodeNumber)
{
    if (rNodeNumber<0 || rNodeNumber>=GetNumNodes())
         throw MechanicsException("[NuTo::StructureGrid::NodeGetNodePtr] Conversion from string to int did not yield valid node number.");
     return &mNodeVec[rNodeNumber];
}

//! @brief a reference to a node
//! @param identifier
//! @return reference to a node
const NuTo::NodeBase* NuTo::StructureGrid::NodeGetNodePtr(int rNodeNumber) const
{
    if (rNodeNumber<0 || rNodeNumber>=GetNumNodes())
        throw MechanicsException("[NuTo::StructureGrid::NodeGetNodePtr] Conversion from string to int did not yield valid node number.");
    return &mNodeVec[rNodeNumber];
 }

//! @brief a reference to a node
//! @param node ID
//! @return reference to a node
NuTo::NodeBase* NuTo::StructureGrid::NodeGetNodePtrFromId(int rNodeId)
{
    unsigned int numGridNodes=(mGridDimension[0]+1)*(mGridDimension[1]+1)*(mGridDimension[2]+1);//all nodes of the grid
    if (rNodeId<0 || rNodeId>=(int) numGridNodes)
         throw MechanicsException("[NuTo::StructureGrid::NodeGetNodePtrFromId] Grid Node number is not valid.");
    int nodeNumber(0);
    boost::ptr_vector<NodeBase>::const_iterator it;
    for (it = mNodeVec.begin(); it!= mNodeVec.end(); it++,nodeNumber++)
    {
        if (it->GetNodeId()==rNodeId)
        {
            break;
        }
    }
    if (it== mNodeVec.end())
        throw MechanicsException("[NuTo::StructureGrid::NodeGetNodePtrFromId] Node with this id does not exist.");
    return &mNodeVec[nodeNumber];
}

//! @brief a reference to a node
//! @param node ID
//! @return reference to a node
const NuTo::NodeBase* NuTo::StructureGrid::NodeGetNodePtrFromId(int rNodeId) const
{
    unsigned int numGridNodes=(mGridDimension[0]+1)*(mGridDimension[1]+1)*(mGridDimension[2]+1);//all nodes of the grid
    if (rNodeId<0 || rNodeId>=(int) numGridNodes)
         throw MechanicsException("[NuTo::StructureGrid::NodeGetNodePtrFromId] Grid Node number is not valid.");
    int nodeNumber(0);
     boost::ptr_vector<NodeBase>::const_iterator it;
     for (it = mNodeVec.begin(); it!= mNodeVec.end(); it++,nodeNumber++)
     {
         if (it->GetNodeId()==rNodeId)
              break;
     }
     if (it== mNodeVec.end())
           throw MechanicsException("[NuTo::StructureGrid::NodeGetNodePtrFromId] Node with this node id  does not exist.");
    return &mNodeVec[nodeNumber];
}

//! @brief a reference to a node
//! @param node ID
//! @return reference to a node
int NuTo::StructureGrid::NodeGetNodeNumberFromId(int rNodeId)
{
    unsigned int numGridNodes=(mGridDimension[0]+1)*(mGridDimension[1]+1)*(mGridDimension[2]+1);//all nodes of the grid
    if (rNodeId<0 || rNodeId>=(int) numGridNodes)
         throw MechanicsException("[NuTo::StructureGrid::NodeGetNodePtrFromId] Grid Node number is not valid.");
    int nodeNumber(0);
     boost::ptr_vector<NodeBase>::const_iterator it;
     for (it = mNodeVec.begin(); it!= mNodeVec.end(); it++,nodeNumber++)
     {
         if (it->GetNodeId()==rNodeId)
             break;
      }
     if (it== mNodeVec.end())
           throw MechanicsException("[NuTo::StructureGrid::NodeGetNodeNumberFromId] Node with this node id does not exist.");
    return nodeNumber;
}

//! @brief a reference to a node
//! @param node ID
//! @return reference to a node
const int NuTo::StructureGrid::NodeGetNodeNumberFromId(int rNodeId) const
{
    unsigned int numGridNodes=(mGridDimension[0]+1)*(mGridDimension[1]+1)*(mGridDimension[2]+1);//all nodes of the grid
    if (rNodeId<0 || rNodeId>=(int) numGridNodes)
         throw MechanicsException("[NuTo::StructureGrid::NodeGetNodePtrFromId] Grid Node number is not valid.");
    int nodeNumber(0);
     boost::ptr_vector<NodeBase>::const_iterator it;
     for (it = mNodeVec.begin(); it!= mNodeVec.end(); it++,nodeNumber++)
     {
         if (it->GetNodeId()==rNodeId)
             break;
     }
     if (it== mNodeVec.end())
           throw MechanicsException("[NuTo::StructureGrid::NodeGetNodeNumberFromId] Node with this node id does not exist.");
     return nodeNumber;
}


//! @brief gives the identifier of a node
//! @param reference to a node
//! @return identifier
int NuTo::StructureGrid::NodeGetId(const NodeBase* rNode)const
{
    return  rNode->GetNodeId();
}

//! @brief gives the identifier of a node
//! @param reference to a node
//! @return identifier
int NuTo::StructureGrid::NodeGetId(int rNodeNumber)const
{
    if (rNodeNumber<0 || rNodeNumber>=GetNumNodes())
        throw MechanicsException("[NuTo::StructureGrid::NodeGetID] Node number is not valid.");
    const NodeBase* myNode = &mNodeVec[rNodeNumber];
    return myNode->GetNodeId();
}
//! @brief info about the nodes in the grid structure
void NuTo::StructureGrid::NodeInfo(int mVerboseLevel) const
{
    std::cout<<"number of nodes   : " << mNodeVec.size() <<std::endl;
}

/*void NuTo::StructureGrid::GetNodeCoordinates(int rNodeNumber)
{

}*/
/*void NuTo::StructureGrid::CreateNodeGrid(std::string rDOFs)
{
    int numNodes=mGridDimension[0]+1*mGridDimension[1]+1*mGridDimension[2]+1;
    for (int nodeID =0; nodeID<numNodes;nodeID++)
    {
        NuTo::StructureGrid::NodeCreate(nodeID,rDOFs);
    }
}*/

//! @brief create node grid without data free nodes
void NuTo::StructureGrid::CreateNodeGrid(std::string rDOFs)
{
    unsigned int numGridNodes=(mGridDimension[0]+1)*(mGridDimension[1]+1)*(mGridDimension[2]+1);//all nodes of the grid
    unsigned int numMatNodes=0;//all existing nodes (with material)
    int numVoxel=mGridDimension[0]*mGridDimension[1]*mGridDimension[2];
    NuTo::FullMatrix<int> imageValues (numVoxel,1);
    imageValues.FullMatrix<int>::ImportFromVtkASCIIFile(mImageDataFile);
    for (unsigned int countNodes =0; countNodes<numGridNodes;countNodes++)//countNodes correspond to nodeID
        {
         //get coincident voxels for each node, check if one voxel has material, then create node
         TCoincidentVoxelList coincidentVoxels=GetCoincidenceVoxelIDs(countNodes);
         int flag=0;
         for (int count =0; count<8; count++)
         {
             // voxel exist (for boundary nodes)
             if (coincidentVoxels[count]>-1)
             {
                   // voxel has material
                 if(imageValues(coincidentVoxels[count],0)>130)//@TODO replace with variable for material boundary values
                 {
                     flag=1;
                     count=8;
                 }
             }
         }
         if (flag) //flag!=0
         {
             NuTo::StructureGrid::NodeCreate(numMatNodes,countNodes,rDOFs);//node number, node id, attr.
             numMatNodes++;
         }
     }
}

//! @brief get coincident Voxels from one node in following order: against the clock, first bottom voxels, first voxel BSW = BottomSouthWest
//! @param Node ID
//! @return vector of Voxel IDSs coincident to this node
/* voxels
 * up:
 *      7       6
 *
 *      4       5
 * bottom:
 *      3       2
 *
 *      0       1
 */
NuTo::StructureGrid::TCoincidentVoxelList  NuTo::StructureGrid::GetCoincidenceVoxelIDs(unsigned int rNodeID)
{
    TCoincidentVoxelList coincidentVoxels(8);
    unsigned int numDim2;
    unsigned int numDim3;
    numDim3=rNodeID/((mGridDimension[0]+1)*(mGridDimension[1]+1));
    numDim2=(rNodeID-(mGridDimension[0]+1)*(mGridDimension[1]+1)*numDim3)/(mGridDimension[0]+1);

    // for all nodes
    coincidentVoxels[0]=(rNodeID-numDim2-numDim3 *( mGridDimension[1] + mGridDimension[0]+2)-mGridDimension[0]-1 -mGridDimension[0]*mGridDimension[1]);
    coincidentVoxels[1]=(rNodeID-numDim2-numDim3 *( mGridDimension[1] + mGridDimension[0]+2)-mGridDimension[0] -mGridDimension[0]*mGridDimension[1]);
    coincidentVoxels[2]=(rNodeID-numDim2-numDim3 *( mGridDimension[1] + mGridDimension[0]+2)-mGridDimension[0]*mGridDimension[1]);
    coincidentVoxels[3]=(rNodeID-numDim2-numDim3 *( mGridDimension[1] + mGridDimension[0]+2) -1 -mGridDimension[0]*mGridDimension[1]);

    coincidentVoxels[4]=(rNodeID-numDim2-numDim3 *( mGridDimension[1] + mGridDimension[0]+2)-mGridDimension[0]-1);
    coincidentVoxels[5]=(rNodeID-numDim2-numDim3 *( mGridDimension[1] + mGridDimension[0]+2)-mGridDimension[0]);
    coincidentVoxels[6]=(rNodeID-numDim2-numDim3 *( mGridDimension[1] + mGridDimension[0]+2));
    coincidentVoxels[7]=(rNodeID-numDim2-numDim3 *( mGridDimension[1] + mGridDimension[0]+2) -1);

    // for nodes in first level related to z
    if (numDim3==0)
    {
        for (int count =0;count<4;count++)
            coincidentVoxels[count]=-1;
    }
    // for nodes in first last related to z
    else if (numDim3==mGridDimension[2])
    {
        for (int count =4;count<8;count++)
            coincidentVoxels[count]=-1;
    }
    if(rNodeID == 0)
        coincidentVoxels[7]=-1;

    // for nodes with dim0=mGridDimension[0]!!!
    else if ((rNodeID-(mGridDimension[0]+1)*(mGridDimension[1]+1)*numDim3) % (mGridDimension[0]+1)==0 )
    {
         coincidentVoxels[1]=-1;
         coincidentVoxels[2]=-1;
         coincidentVoxels[5]=-1;
         coincidentVoxels[6]=-1;
    }
    // for nodes with dim0=0
    else if ((rNodeID-(mGridDimension[0]+1)*(mGridDimension[1]+1)*numDim3) % (mGridDimension[0]+1)==1 )
    {
         coincidentVoxels[0]=-1;
         coincidentVoxels[3]=-1;
         coincidentVoxels[4]=-1;
         coincidentVoxels[7]=-1;
    }
    // for node with dim1=0
    if (numDim2==0)
    {
        coincidentVoxels[0]=-1;
        coincidentVoxels[1]=-1;
        coincidentVoxels[4]=-1;
        coincidentVoxels[5]=-1;
    }
    // for node with dim1=MGridDimension[1]
    else if (numDim2==mGridDimension[1])
    {
        coincidentVoxels[2]=-1;
        coincidentVoxels[3]=-1;
        coincidentVoxels[6]=-1;
        coincidentVoxels[7]=-1;
    }
/*
    std::cout<<"Knoten"<<rNodeID<< "Voxels:"<<std::endl;
    for (int count=0;count<8;count++)
        std::cout<< coincidentVoxels[count]<<" ";
    std::cout<<" "<<std::endl;
 */
    return coincidentVoxels;
}
//! @brief creates a node
//! @param rNodeID - Number in grid; rDoFs - kind of degrees of freedom

void NuTo::StructureGrid::NodeCreate(unsigned int rNodeNumber,unsigned int rNodeID,std::string rDOFs)
{
    // check node number
    unsigned int numNodes=(mGridDimension[0]+1)*(mGridDimension[1]+1)*(mGridDimension[2]+1);
    if (rNodeID >numNodes)
    {
        throw MechanicsException("[NuTo::StructureGrid::NodeCreate] Node number is outside the grid.");
    }
    if (rNodeNumber >  mNodeVec.size())
    {
        throw MechanicsException("[NuTo::StructureGrid::NodeCreate] Node exist already.");

    }

    // transform string to uppercase
    std::transform(rDOFs.begin(), rDOFs.end(), rDOFs.begin(), toupper);

    // check all values
    int attributes(1 << Node::COORDINATES);
    //! bit 0 : coordinates
    //! bit 1 : displacements
    //! bit 2 : rotations
    //! bit 3 : temperatures
    boost::tokenizer<> tok(rDOFs);
    for (boost::tokenizer<>::iterator beg=tok.begin(); beg!=tok.end(); ++beg)
    {
        if (*beg=="DISPLACEMENTS")
            attributes = attributes |  1 << Node::DISPLACEMENTS;
        if (*beg=="ROTATIONS")
            attributes = attributes |  1 << Node::ROTATIONS;
        if (*beg=="TEMPERATURES")
            attributes = attributes |  1 << Node::TEMPERATURES;
    }

    NodeBase* nodePtr;
    switch (attributes)
    {
        // the << shifts the 1 bitwise to the left, so 1<<n = 2^n
        // it actually sets the n-th bit (from the right) to 1, and all the other to zero
    case (1 << Node::COORDINATES):
        // for grid nodes COORDINATES are replaced with NodeID
        switch (mDimension)
        {
        case 1:
            nodePtr = new NuTo::NodeGrid1D(rNodeID);
            break;
        case 2:
            nodePtr = new NuTo::NodeGrid2D(rNodeID);
            break;
        case 3:
            nodePtr = new NuTo::NodeGrid3D(rNodeID);
            break;
        default:
            throw MechanicsException("[NuTo::StructureGrid::NodeCreate] Dimension of the structure is not valid.");
        }
       break;
    case (1 << Node::COORDINATES) | (1 << Node::DISPLACEMENTS):
        switch (mDimension)
         {
         case 1:
             nodePtr = new NuTo::NodeGridDisplacements1D(rNodeID);
             break;
         case 2:
             nodePtr = new NuTo::NodeGridDisplacements2D(rNodeID);
             break;
         case 3:
             nodePtr = new NuTo::NodeGridDisplacements3D(rNodeID);
             break;
         default:
             throw MechanicsException("[NuTo::StructureGrid::NodeCreate] Dimension of the structure is not valid.");
         }
        break;
//! @Todo: add NodeGridCoordinatesDisplacementsRotations.h
    case (1 << Node::COORDINATES) | (1 << Node::DISPLACEMENTS) | (1 << Node::ROTATIONS):
         throw MechanicsException("[NuTo::StructureGrid::NodeCreate] Rotation not implemented.");
         break;
   }
    //add grid node number to node
   //add node to vector, richtige id ????????????????????????
    this->mNodeVec.push_back(nodePtr);

    //renumbering of dofs for global matrices not required, why not???
    this->mNodeNumberingRequired  = true;


}


//! @brief Deletes a node
//! @param rElementIdent identifier for the node
void NuTo::StructureGrid::NodeDelete(const int rNodeNumber)
{

	// @TODO [NuTo::StructureGrid::NodeDelete] has to be implemented
    throw MechanicsException("[NuTo::StructureGrid::NodeDelete] Not implemented yet!!!");

}


//! @brief number the dofs in the structure
void NuTo::StructureGrid::NodeBuildGlobalDofs()
{
    // build initial node numbering
    this->mNumDofs = 0;
    for (boost::ptr_vector<NodeBase>::iterator it = mNodeVec.begin(); it!= mNodeVec.end(); it++)
    {
        it->SetGlobalDofs(this->mNumDofs);
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
                throw MechanicsException("[NuTo::StructureGrid::NodeBuildGlobalDofs] invalid matrix structure.");
            }
        }
    }
    this->mConstraintMatrix.RemoveLastColumns(numDependentDofs);

    // renumber dofs
    for (boost::ptr_vector<NodeBase>::iterator it = mNodeVec.begin(); it!= mNodeVec.end(); it++)
    {
        it->RenumberGlobalDofs(mappingInitialToNewOrdering);
    }

    mNodeNumberingRequired = false;
}

// merge dof values
void NuTo::StructureGrid::NodeMergeActiveDofValues(const NuTo::FullMatrix<double>& rActiveDofValues)
{
    if (this->mNodeNumberingRequired)
    {
        throw MechanicsException("[NuTo::StructureGrid::NodeMergeActiceDofValues] a valid dof numbering was not found (build dof numbering using NodeBuildGlobalDofs).");
    }
    if ((rActiveDofValues.GetNumRows() != this->mNumActiveDofs) || rActiveDofValues.GetNumColumns() != 1)
    {
        throw MechanicsException("[NuTo::StructureGrid::NodeMergeActiceDofValues] invalid dimension of input object (number of active dofs,1).");
    }

    // calculate dependent dof values
    NuTo::FullMatrix<double> dependentDofValues = this->mConstraintRHS - this->mConstraintMatrix * rActiveDofValues;

    // write dof values to the nodes
    for (boost::ptr_vector<NodeBase>::iterator it = mNodeVec.begin(); it!= mNodeVec.end(); it++)
    {
        it->SetGlobalDofValues(rActiveDofValues, dependentDofValues);
    }
}

//! @brief extract dof values (e.g. displacements, temperatures to the nodes)
void NuTo::StructureGrid::NodeExtractDofValues(NuTo::FullMatrix<double>& rActiveDofValues, NuTo::FullMatrix<double>& rDependentDofValues) const
{

    try
    {
        if (this->mNodeNumberingRequired)
        {
           throw MechanicsException("[NuTo::StructureGrid::NodeExtractDofValues] a valid dof numbering was not found (build dof numbering using NodeBuildGlobalDofs).");
        }
    }
    catch(MechanicsException& e)
    {
        throw;
    }

    rActiveDofValues.Resize(this->mNumActiveDofs,1);
    rDependentDofValues.Resize(this->mNumDofs - this->mNumActiveDofs,1);

    // extract dof values from nodes
    for (boost::ptr_vector<NodeBase>::const_iterator it = this->mNodeVec.begin(); it!= this->mNodeVec.end(); it++)
    {
        it->GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }
}
/*
void NuTo::StructureGrid::NodeGetInternalForce(const NodeBase* rNode, NuTo::FullMatrix<double>& rNodeForce)const
{
   // define variables storing the element contribution outside the loop
    NuTo::FullMatrix<double> elementVector;
    std::vector<int> elementVectorGlobalDofs;

    //iterate over the elements, check if the node is part of the element
    for (boost::ptr_vector<ElementBase>::const_iterator itElement = mElementVec.begin(); itElement!=mElementVec.end(); itElement++)
    {
        const ElementBase* elem_ptr(itElement);
        for (int local_node=0; local_node<elem_ptr->GetNumNodes(); local_node++)
        {
            if (elem_ptr->GetNode(local_node)==rNode)
            {
                //node has been found, build internal force vector of element and add relevant parts to the global vector
                elem_ptr->CalculateGradientInternalPotential(elementVector, elementVectorGlobalDofs);
                const NodeDisplacementsBase *theDispNode(dynamic_cast<const NodeDisplacementsBase* > (rNode));
                if (theDispNode==0)
                {
                    std::stringstream message;
                    message << "[NuTo::StructureGrid::NodeGetInternalForce] Node" << NodeGetId(rNode) << " has no displacement degrees of freedom.";
                    throw MechanicsException(message.str());
                }

                for (int theDisp=0; theDisp<theDispNode->GetNumDisplacements(); theDisp++)
                {
                    // find global dof number
                    int theDOF = theDispNode->GetDofDisplacement(theDisp);

                    for (int local_dof=0; local_dof<elementVectorGlobalDofs.size(); local_dof++)
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
void NuTo::StructureGrid::NodeGetInternalForce(int rNodeId, FullMatrix<double>& rNodeForce)const
{
    const NodeBase* node_ptr(NodeGetNodePtr(rNodeId));

    // initialize matrix
    rNodeForce.Resize(node_ptr->GetNumDisplacements(),1);

    // add internal force
    NodeGetInternalForce(node_ptr, rNodeForce);
}
*/

