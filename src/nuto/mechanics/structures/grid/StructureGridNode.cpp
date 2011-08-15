// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/vector.hpp>
#endif //ENABLE_SERIALIZATION
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
NuTo::NodeGrid3D* NuTo::StructureGrid::NodeGridGetNodePtr(int rIdent)
{
	if (rIdent<0 || rIdent>=GetNumNodes())
	{
		std::cout<<"node number "<<rIdent<<std::endl;
		throw MechanicsException("[NuTo::StructureGrid::NodeGridGetNodePtr] (sec) Conversion from string to int did not yield valid node number.");
	}
	return mNodeVec[rIdent];
}
//! @brief a reference to a node
//! @param identifier
//! @return reference to a node
const NuTo::NodeGrid3D* NuTo::StructureGrid::NodeGridGetNodePtr(int rIdent) const
{
    if (rIdent<0 || rIdent>=GetNumNodes())
     	throw MechanicsException("[NuTo::StructureGrid::NodeGetNodePtr] (third) Conversion from string to int did not yield valid node number.");
    return mNodeVec[rIdent];
 }

//! @brief info about the nodes in the grid structure
void NuTo::StructureGrid::NodeInfo(int mVerboseLevel) const
{
    std::cout<<"number of nodes   : " << mNodeVec.size() <<std::endl;
}


//! @brief create node grid without data free nodes
void NuTo::StructureGrid::CreateNodeGrid(std::string rDOFs,int rThresholdMaterialValue)
{
    int numGridNodes=(mGridDimension[0]+1)*(mGridDimension[1]+1)*(mGridDimension[2]+1);//all nodes of the grid
    //int numMatNodes=0;//all existing nodes (with material)
    int numVoxel=mGridDimension[0]*mGridDimension[1]*mGridDimension[2];
    int numGlobDofs=3*numGridNodes;
    mDofIsNotConstraint = new bool [numGlobDofs];
    memset(mDofIsNotConstraint,true,sizeof (mDofIsNotConstraint[0])*numGlobDofs );
    NuTo::FullMatrix<int> imageValues (numVoxel,1);
    imageValues.FullMatrix<int>::ImportFromVtkASCIIFile(mImageDataFile);
    int * coincidentVoxels=new int[8];
    for (int countNodes =0; countNodes<numGridNodes;countNodes++)//countNodes correspond to nodeID
    {
         //get coincident voxels for each node, check if one voxel has material, then create node
         coincidentVoxels=GetCoincidenceVoxelIDs(countNodes);
         bool flag=0;
         for (int count =0; count<8; count++)
         {
        	 // voxel exist (for boundary nodes)
             if (coincidentVoxels[count]>-1)
             {
                 // voxel has material
            	 // color value 0 is material, 255 is air
            	 // material value smaller than thresholdvalue
                 if(imageValues(coincidentVoxels[count],0)<rThresholdMaterialValue)
                 {
                     flag=true;
                     count=8;
					 //numMatNodes++;
                 }
             }
         }
         NuTo::StructureGrid::NodeCreate(flag,countNodes,rDOFs);//flag,gridNode, attr.
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
int*  NuTo::StructureGrid::GetCoincidenceVoxelIDs(int rNodeID)
{
    int* coincidentVoxels=new int[8];
    //get the number of nodes in the actual x-y-dim
    int numDimxy=rNodeID/((mGridDimension[0]+1)*(mGridDimension[1]+1));
    int numDimx=0;
    int residual1=rNodeID%((mGridDimension[0]+1)*(mGridDimension[1]+1));
    int residual2=0;
    numDimx=residual1/(mGridDimension[0]+1);
    residual2=residual1%(mGridDimension[0]+1);

    //get
    //numDimx=(rNodeID-(mGridDimension[0]+1)*(mGridDimension[1]+1)*numDimxy)/(mGridDimension[0]+1);
    if (mVerboseLevel>3)
    	std::cout<<__FILE__<<" "<<__LINE__<<" res1 "<< residual1 <<" numDimx "<<numDimx<<" numDimxy "<<numDimxy<<" res2 "<<residual2 << std::endl;
    //std::cout<<__FILE__<<" "<<__LINE__<<" Rest "<<(rNodeID-(mGridDimension[0]+1)*(mGridDimension[1]+1)*numDimxy) % (mGridDimension[0]+1)<<std::endl;
    //residual2 = (rNodeID-(mGridDimension[0]+1)*(mGridDimension[1]+1)*numDimxy) % (mGridDimension[0]+1);
    // for all nodes
/*    coincidentVoxels[0]=(rNodeID-numDimx  -numDimxy *( mGridDimension[1] * mGridDimension[0]+2)-mGridDimension[0]-1 -mGridDimension[0]*mGridDimension[1]);
    coincidentVoxels[1]=(rNodeID-numDimx  -numDimxy *( mGridDimension[1] * mGridDimension[0]+2)-mGridDimension[0] -mGridDimension[0]*mGridDimension[1]);
    coincidentVoxels[2]=(rNodeID-numDimx  -numDimxy *( mGridDimension[1] * mGridDimension[0]+2)-mGridDimension[0]*mGridDimension[1]);
    coincidentVoxels[3]=(rNodeID-numDimx  -numDimxy *( mGridDimension[1] * mGridDimension[0]+2) -1 -mGridDimension[0]*mGridDimension[1]);

    coincidentVoxels[4]=(rNodeID-numDimx  -numDimxy *( mGridDimension[1] * mGridDimension[0]+2)-mGridDimension[0]-1);
    coincidentVoxels[5]=(rNodeID-numDimx  -numDimxy *( mGridDimension[1] * mGridDimension[0]+2)-mGridDimension[0]);
    //coincidentVoxels[6]=(rNodeID-numDimx  -numDimxy *( mGridDimension[1] * mGridDimension[0]+2));
    coincidentVoxels[7]=(rNodeID-numDimx  -numDimxy *( mGridDimension[1] * mGridDimension[0]+2) -1);
*/

    coincidentVoxels[0]=(numDimx * mGridDimension[0] + (numDimxy - 1) *( mGridDimension[1] * mGridDimension[0]) + residual2 ) - mGridDimension[0]-1;
    coincidentVoxels[1]=(numDimx * mGridDimension[0] + (numDimxy - 1) *( mGridDimension[1] * mGridDimension[0]) + residual2 ) - mGridDimension[0];
    coincidentVoxels[2]=(numDimx * mGridDimension[0] + (numDimxy - 1) *( mGridDimension[1] * mGridDimension[0]) + residual2 );
    coincidentVoxels[3]=(numDimx * mGridDimension[0] + (numDimxy - 1) *( mGridDimension[1] * mGridDimension[0]) + residual2 ) -1;

    coincidentVoxels[4]=(numDimx * mGridDimension[0] + numDimxy *( mGridDimension[1] * mGridDimension[0]) + residual2 ) - mGridDimension[0]-1;
    coincidentVoxels[5]=(numDimx * mGridDimension[0] + numDimxy *( mGridDimension[1] * mGridDimension[0]) + residual2 ) - mGridDimension[0];
    coincidentVoxels[6]=(numDimx * mGridDimension[0] + numDimxy *( mGridDimension[1] * mGridDimension[0]) + residual2 );
    coincidentVoxels[7]=(numDimx * mGridDimension[0] + numDimxy *( mGridDimension[1] * mGridDimension[0]) + residual2 ) -1;

/*
    std::cout<<__FILE__ <<" "<<__LINE__<<" Knoten "<<rNodeID<< " Voxels: ";
	for (int count=0;count<8;count++)
		std::cout<< coincidentVoxels[count]<<" ";
	std::cout<<" "<<std::endl;
*/
	// for nodes in first level related to z
    if (numDimxy==0)
    {
        for (int count =0;count<4;count++)
            coincidentVoxels[count]=-1;
    }
    // for nodes in first last related to z
    else if (numDimxy==mGridDimension[2])
    {
        for (int count =4;count<8;count++)
            coincidentVoxels[count]=-1;
    }
    // for nodes with dim0=0!!!
    if ((rNodeID-(mGridDimension[0]+1)*(mGridDimension[1]+1)*numDimxy) % (mGridDimension[0]+1)==0 )
    {
         coincidentVoxels[0]=-1;
         coincidentVoxels[3]=-1;
         coincidentVoxels[4]=-1;
         coincidentVoxels[7]=-1;
    }
    // for nodes with dim0=dimension[0]
    else if ((rNodeID-(mGridDimension[0]+1)*(mGridDimension[1]+1)*numDimxy) % (mGridDimension[0]+1)==mGridDimension[0] )
    {
         coincidentVoxels[1]=-1;
         coincidentVoxels[2]=-1;
         coincidentVoxels[5]=-1;
         coincidentVoxels[6]=-1;
    }
    // for node with dim1=0
    if (numDimx==0)
    {
        coincidentVoxels[0]=-1;
        coincidentVoxels[1]=-1;
        coincidentVoxels[4]=-1;
        coincidentVoxels[5]=-1;
    }
    // for node with dim1=MGridDimension[1]
    else if (numDimx==mGridDimension[1])
    {
        coincidentVoxels[2]=-1;
        coincidentVoxels[3]=-1;
        coincidentVoxels[6]=-1;
        coincidentVoxels[7]=-1;
    }
    if (mVerboseLevel>3)
    {
    	std::cout<<__FILE__ <<" "<<__LINE__<<" Knoten "<<rNodeID<< " Voxels: ";
    	for (int count=0;count<8;count++)
    		std::cout<< coincidentVoxels[count]<<" ";
    	std::cout<<" "<<std::endl;
    }
    return coincidentVoxels;
}
//! @brief creates a node
//! @param flag - node exists really; rDoFs - kind of degrees of freedom

void NuTo::StructureGrid::NodeCreate(bool flag,int rNodeGridNum,std::string rDOFs)
{
   // check node number
    int numNodes=(mGridDimension[0]+1)*(mGridDimension[1]+1)*(mGridDimension[2]+1);
    if (rNodeGridNum != (int) mNodeVec.size())
        throw MechanicsException("[NuTo::StructureGrid::NodeCreate] Wrong next node number.");

    if (rNodeGridNum >numNodes)
        throw MechanicsException("[NuTo::StructureGrid::NodeCreate] Node number is outside the grid.");

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

    NodeGrid3D* nodePtr(0);
    if(flag)
    {
		switch (attributes)
		{
			// the << shifts the 1 bitwise to the left, so 1<<n = 2^n
			// it actually sets the n-th bit (from the right) to 1, and all the other to zero
		case (1 << Node::COORDINATES):
			// for grid nodes COORDINATES are replaced with NodeID
			switch (mDimension)
			{
			case 1:
				//nodePtr = new NuTo::NodeGrid1D(rNodeGridNum);
				break;
			case 2:
				//nodePtr = new NuTo::NodeGrid2D(rNodeGridNum);
				break;
			case 3:
				nodePtr = new NuTo::NodeGrid3D();
				break;
			default:
				throw MechanicsException("[NuTo::StructureGrid::NodeCreate] Dimension of the structure is not valid.");
			}
		   break;
		case (1 << Node::COORDINATES) | (1 << Node::DISPLACEMENTS):
			switch (mDimension)
			 {
			 case 1:
				 //nodePtr = new NuTo::NodeGridDisplacements1D(rNodeGridNum);
				 break;
			 case 2:
				 //nodePtr = new NuTo::NodeGridDisplacements2D(rNodeGridNum);
				 break;
			 case 3:
				 nodePtr = new NuTo::NodeGridDisplacements3D();
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
    }
	//add node to vector
	this->mNodeVec.push_back(nodePtr);
}
//! @brief NodeGetConstraintSwitch
//! @param rGlobalDof
//! @return switch for constraint
bool NuTo::StructureGrid::NodeGetConstraintSwitch(int rGlobalDof)
{
	return mDofIsNotConstraint[rGlobalDof];
}

//! @brief NodeSetConstraintSwitch
//! @brief rDirection goes from 0 to 2
//! @param rGridNodeNum, rConstraint
void NuTo::StructureGrid::NodeSetConstraintSwitch(int rGridNodeNum, int rDirection, bool rConstraint)
{
	assert(rDirection<mDimension);
	assert(rGridNodeNum<(mGridDimension[0]+1)*(mGridDimension[1]+1)*(mGridDimension[2]+1));
	mDofIsNotConstraint[3*rGridNodeNum+rDirection]=rConstraint;
	if (mVerboseLevel>3)
		std::cout<<__FILE__<<__LINE__<<"dof is constraint: GridNode "<<rGridNodeNum<<" direction "<<rDirection<<" place in mDofIsNotConstraint "<<3*rGridNodeNum+rDirection<<std::endl;
}


//! @brief Deletes a node
//! @param rElementIdent identifier for the node
void NuTo::StructureGrid::NodeDelete(const int rIdent)
{
	   if (rIdent<0 || rIdent>=GetNumNodes())
	     	throw MechanicsException("[NuTo::StructureGrid::NodeGetNodePtr] (third) Conversion from string to int did not yield valid node number.");
	   mNodeVec[rIdent]=NULL;

	// @TODO [NuTo::StructureGrid::NodeDelete] has to be implemented
    throw MechanicsException("[NuTo::StructureGrid::NodeDelete] Not implemented yet!!!");

}

//! @brief gets the displacements of a node
//! @param rIdent node identifier
//! @param rDisplacements matrix (one column) with the displacements
void NuTo::StructureGrid::NodeGetDisplacements(int rNode, FullMatrix<double>& rDisplacements)const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	const NodeBase* nodePtr = NodeGridGetNodePtr(rNode);
	if (nodePtr)
	{
		rDisplacements.Resize(3,1);
		nodePtr->GetDisplacements3D(rDisplacements.mEigenMatrix.data());
	}
	//else: displacements for non existing node remain zero
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime && mVerboseLevel>3)
        std::cout<<"[NuTo::StructureBase::NodeGetDisplacements] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}
//! @brief sets the displacements of a node
//! @param rIdent node identifier
//! @param rDisplacements matrix (one column) with the displacements
void NuTo::StructureGrid::NodeSetDisplacements(int rNode, const FullMatrix<double>& rDisplacements)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	NodeBase* nodePtr=NodeGridGetNodePtr(rNode);

	if (rDisplacements.GetNumColumns()!=1)
	throw MechanicsException("[NuTo::StructureGrid::NodeSetDisplacements] Displacement matrix has to have a single column.");
	try
	{
		switch (rDisplacements.GetNumRows())
		{
		case 3:
			nodePtr->SetDisplacements3D(rDisplacements.mEigenMatrix.data());
		break;
		default:
			throw MechanicsException("[NuTo::StructureGrid::NodeSetDisplacements] The number of displacement components is either 1, 2 or 3.");
		}
	}
    catch(NuTo::MechanicsException & b)
	{
    	b.AddMessage("[NuTo::StructureGrid::NodeSetDisplacements] Error setting displacements.");
    	throw b;
	}
    catch(...)
	{
	    throw MechanicsException("[NuTo::StructureGrid::NodeSetDisplacements] Error setting displacements of node (unspecified exception).");
	}
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime && mVerboseLevel>3)
        std::cout<<"[NuTo::StructureGrid::NodeSetDisplacements] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

/*
//! @brief number the dofs in the structure
void NuTo::StructureGrid::NodeBuildGlobalDofs()
{
    // build initial node numbering
    this->mNumDofs = 0;
     for (std::vector<NodeGrid3D*>::iterator it = mNodeVec.begin(); it!= mNodeVec.end(); it++)
    {
        it->SetGlobalDofs(this->mNumDofs);
    }
    // build constraint matrix
    this->mNodeNumberingRequired = false;
    //this->ConstraintGetConstraintMatrix(this->mConstraintMatrix, this->mConstraintRHS);
 //   std::cout<<__FILE__<<" "<<__LINE__<<" constraint matrix \n"<<(NuTo::FullMatrix<double> (this->mConstraintMatrix))<<"\n rhs \n"<<this->mConstraintRHS<<std::endl;
   //get constraint dofs
    // int numDependentDofs = this->mConstraintMatrix.GetNumRows();

    int numDependentDofs=0;
    int numGridNodes = (this->mGridDimension[0]+1)*(this->mGridDimension[1]+1)*(this->mGridDimension[2]+1);
    for (int count = 0;count<3*numGridNodes;++count)
    {
    	if (!mDofIsNotConstraint[count])
    		numDependentDofs++;
    }

    this->mNumActiveDofs = this->mNumDofs - numDependentDofs;
    //this->mConstraintMatrix=0;
    //this->mNodeNumberingRequired = true;


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
    for (std::vector<NodeGrid3D*>::iterator it = mNodeVec.begin(); it!= mNodeVec.end(); it++)
    {
        it->SetGlobalDofValues(rActiveDofValues, dependentDofValues);
    }
}
*/
/*
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
    for (std::vector<NodeGrid3D>::const_iterator it = this->mNodeVec.begin(); it!= this->mNodeVec.end(); it++)
    {
        it->GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }
}
*/
/*
void NuTo::StructureGrid::NodeGetInternalForce(const NodeGrid3D* rNode, NuTo::FullMatrix<double>& rNodeForce)const
{
   // define variables storing the element contribution outside the loop
    NuTo::FullMatrix<double> elementVector;
    std::vector<int> elementVectorGlobalDofs;

    //iterate over the elements, check if the node is part of the element
    for (std::vector<ElementBase>::const_iterator itElement = mElementVec.begin(); itElement!=mElementVec.end(); itElement++)
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
    const NodeGrid3D* node_ptr(NodeGetNodePtr(rNodeId));

    // initialize matrix
    rNodeForce.Resize(node_ptr->GetNumDisplacements(),1);

    // add internal force
    NodeGetInternalForce(node_ptr, rNodeForce);
}
*/

