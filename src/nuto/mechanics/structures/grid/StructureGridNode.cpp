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
    return (int) mEdgeId.size();
}

//! @brief a reference to a node
//! @param identifier
//! @return reference to a node
NuTo::NodeBase* NuTo::StructureGrid::NodeGetNodePtr(int rIdent)
{
    throw MechanicsException("[NuTo::StructureGrid::NodeGetNodePtr] Not implemented.");
//    boost::ptr_map<int,NodeBase>::iterator it = mNodeMap.find(rIdent);
//    if (it!=mNodeMap.end())
//        return it->second;
//    else
//    {
//    	std::stringstream out;
//    	out << rIdent;
//    	throw MechanicsException("[NuTo::StructureGrid::NodeGetNodePtr] Node with id " + out.str() +" does not exist.");
//    }
}

//! @brief a reference to a node
//! @param identifier
//! @return reference to a node
const NuTo::NodeBase* NuTo::StructureGrid::NodeGetNodePtr(int rIdent)const
{
    throw MechanicsException("[NuTo::StructureGrid::NodeGetNodePtr] Node does not exist.");
//    boost::ptr_map<int,NodeBase>::const_iterator it = mNodeMap.find(rIdent);
//    if (it!=mNodeMap.end())
//        return it->second;
//    else
//    {
//    	std::stringstream out;
//    	out << rIdent;
//    	throw MechanicsException("[NuTo::StructureGrid::NodeGetNodePtr] Node with id " + out.str() +" does not exist.");
//    }
}

//! @brief gives the identifier of a node
//! @param reference to a node
//! @return identifier
int NuTo::StructureGrid::NodeGetId(const NodeBase* rNode)const
{
    throw MechanicsException("[NuTo::StructureGrid::NodeGetId] Node does not exist.");
//    for (boost::ptr_map<int,NodeBase>::const_iterator
//            it = mNodeMap.begin(); it!= mNodeMap.end(); it++)
//    {
//        if (it->second==rNode)
//            return it->first;
//    }
//    throw MechanicsException("[NuTo::StructureGrid::GetNodeId] Node does not exist.");
}

//! @brief info about the nodes in the grid StructureGrid
void NuTo::StructureGrid::NodeInfo(int mVerboseLevel) const
{
    std::cout<<"number of nodes   : " << mEdgeId.size() <<std::endl;
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
    if (StructureBase::mVerboseLevel>3)
    	std::cout<<__FILE__<<" "<<__LINE__<<" res1 "<< residual1 <<" numDimx "<<numDimx<<" numDimxy "<<numDimxy<<" res2 "<<residual2 << std::endl;
    coincidentVoxels[0]=(numDimx * mGridDimension[0] + (numDimxy - 1) *( mGridDimension[1] * mGridDimension[0]) + residual2 ) - mGridDimension[0]-1;
    coincidentVoxels[1]=(numDimx * mGridDimension[0] + (numDimxy - 1) *( mGridDimension[1] * mGridDimension[0]) + residual2 ) - mGridDimension[0];
    coincidentVoxels[2]=(numDimx * mGridDimension[0] + (numDimxy - 1) *( mGridDimension[1] * mGridDimension[0]) + residual2 );
    coincidentVoxels[3]=(numDimx * mGridDimension[0] + (numDimxy - 1) *( mGridDimension[1] * mGridDimension[0]) + residual2 ) -1;

    coincidentVoxels[4]=(numDimx * mGridDimension[0] + numDimxy *( mGridDimension[1] * mGridDimension[0]) + residual2 ) - mGridDimension[0]-1;
    coincidentVoxels[5]=(numDimx * mGridDimension[0] + numDimxy *( mGridDimension[1] * mGridDimension[0]) + residual2 ) - mGridDimension[0];
    coincidentVoxels[6]=(numDimx * mGridDimension[0] + numDimxy *( mGridDimension[1] * mGridDimension[0]) + residual2 );
    coincidentVoxels[7]=(numDimx * mGridDimension[0] + numDimxy *( mGridDimension[1] * mGridDimension[0]) + residual2 ) -1;

    return coincidentVoxels;
}
//! @brief creates a node
//! @param flag - node exists really; rDoFs - kind of degrees of freedom

void NuTo::StructureGrid::NodeCreate(int rNodeNum,std::string rDOFs)
{
    throw MechanicsException("[NuTo::StructureGrid::NodeCreate] Node does not exist.");   // check node number
//	boost::ptr_map<int,NodeBase>::iterator it = this->mNodeMap.find(rNodeNum);
//	if(it != this->mNodeMap.end())
//	{
//		throw MechanicsException("[NuTo::StructureGrid::NodeCreate] Node already exists.");
//	}
//    int numGridNodes=(mGridDimension[0]+1)*(mGridDimension[1]+1)*(mGridDimension[2]+1);
//
//    if (rNodeNum >numGridNodes)
//        throw MechanicsException("[NuTo::StructureGrid::NodeCreate] Node number is outside the grid.");
//
//    // transform string to uppercase
//    std::transform(rDOFs.begin(), rDOFs.end(), rDOFs.begin(), toupper);
//
//    // check all values
//    int attributes(1 << Node::COORDINATES);
//    //! bit 0 : coordinates
//    //! bit 1 : displacements
//    //! bit 2 : rotations
//    //! bit 3 : temperatures
//    boost::tokenizer<> tok(rDOFs);
//    for (boost::tokenizer<>::iterator beg=tok.begin(); beg!=tok.end(); ++beg)
//    {
//        if (*beg=="DISPLACEMENTS")
//            attributes = attributes |  1 << Node::DISPLACEMENTS;
//        if (*beg=="ROTATIONS")
//            attributes = attributes |  1 << Node::ROTATIONS;
//        if (*beg=="TEMPERATURES")
//            attributes = attributes |  1 << Node::TEMPERATURES;
//    }
//
//    NodeGrid3D* nodePtr(0);
//	switch (attributes)
//	{
//		// the << shifts the 1 bitwise to the left, so 1<<n = 2^n
//		// it actually sets the n-th bit (from the right) to 1, and all the other to zero
//	case (1 << Node::COORDINATES):
//		// for grid nodes COORDINATES are replaced with NodeID
//		throw MechanicsException("[NuTo::StructureGrid::NodeCreate] Grid nodes do not have coordinates.");
//		break;
//	case (1 << Node::COORDINATES) | (1 << Node::DISPLACEMENTS):
//		switch (mDimension)
//		 {
//		 case 1:
//			 //nodePtr = new NuTo::NodeGridDisplacements1D(rNodeGridNum);
//			 break;
//		 case 2:
//			 //nodePtr = new NuTo::NodeGridDisplacements2D(rNodeGridNum);
//			 break;
//		 case 3:
//			 nodePtr = new NuTo::NodeGridDisplacements3D();
//			 break;
//		 default:
//			 throw MechanicsException("[NuTo::StructureGrid::NodeCreate] Dimension of the StructureGrid is not valid.");
//		 }
//		break;
////! @Todo: add NodeGridCoordinatesDisplacementsRotations.h
//	case (1 << Node::COORDINATES) | (1 << Node::DISPLACEMENTS) | (1 << Node::ROTATIONS):
//		 throw MechanicsException("[NuTo::StructureGrid::NodeCreate] Rotation not implemented.");
//		 break;
//   }
//    // add node to map
//    this->mNodeMap.insert(rNodeNum, nodePtr);
}

//! @brief NodeGetConstraintSwitch true is constraint
//! @param rGlobalDof
//! @return bool ... switch for constraint
bool NuTo::StructureGrid::NodeGetConstraintSwitch(int rGlobalDof)
{
	return mDofIsConstraint[rGlobalDof];
}

//! @brief Deletes a node
//! @param rElementIdent identifier for the node
void NuTo::StructureGrid::NodeDelete(const int rIdent)
{
    throw MechanicsException("[NuTo::StructureGrid::NodeDelete] Not implemented yet!!!");
}

//! @brief gets the displacements of a node
//! @param rIdent node identifier
//! @param rDisplacements matrix (one column) with the displacements
void NuTo::StructureGrid::NodeGetDisplacements(int rNode, FullVector<double,Eigen::Dynamic>& rDisplacements)const
{
    throw MechanicsException("[NuTo::StructureGrid::NodeGetDisplacements] Not implemented.");
//#ifdef SHOW_TIME
//    std::clock_t start,end;
//    start=clock();
//#endif
//	const NodeBase* nodePtr = NodeGetNodePtr(rNode);
//	if (nodePtr)
//	{
//		rDisplacements.Resize(3,1);
//		nodePtr->GetDisplacements3D(rDisplacements.mEigenMatrix.data());
//	}
//	//else: displacements for non existing node remain zero
//#ifdef SHOW_TIME
//    end=clock();
//    if (mShowTime && mVerboseLevel>3)
//        std::cout<<"[NuTo::StructureBase::NodeGetDisplacements] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
//#endif
}
//! @brief sets the displacements of a node
//! @param rIdent node identifier
//! @param rDisplacements matrix (one column) with the displacements
void NuTo::StructureGrid::NodeSetDisplacements(int rNode, const FullVector<double,Eigen::Dynamic>& rDisplacements)
{
	   throw MechanicsException("[NuTo::StructureGrid::NodeSetDisplacements] Not implemented.");
//#ifdef SHOW_TIME
//    std::clock_t start,end;
//    start=clock();
//#endif
//	NodeBase* nodePtr=NodeGetNodePtr(rNode);
//
//	if (rDisplacements.GetNumColumns()!=1)
//	throw MechanicsException("[NuTo::StructureGrid::NodeSetDisplacements] Displacement matrix has to have a single column.");
//	try
//	{
//		switch (rDisplacements.GetNumRows())
//		{
//		case 3:
//			nodePtr->SetDisplacements3D(rDisplacements.mEigenMatrix.data());
//		break;
//		default:
//			throw MechanicsException("[NuTo::StructureGrid::NodeSetDisplacements] The number of displacement components is either 1, 2 or 3.");
//		}
//	}
//    catch(NuTo::MechanicsException & b)
//	{
//    	b.AddMessage("[NuTo::StructureGrid::NodeSetDisplacements] Error setting displacements.");
//    	throw b;
//	}
//    catch(...)
//	{
//	    throw MechanicsException("[NuTo::StructureGrid::NodeSetDisplacements] Error setting displacements of node (unspecified exception).");
//	}
//#ifdef SHOW_TIME
//    end=clock();
//    if (mShowTime && mVerboseLevel>3)
//        std::cout<<"[NuTo::StructureGrid::NodeSetDisplacements] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
//#endif

}

