// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/vector.hpp>
#endif //ENABLE_SERIALIZATION
#include <boost/tokenizer.hpp>
#include <sstream>

#include "nuto/mechanics/structures/grid/StructureGrid.h"
#include "nuto/math/FullMatrix.h"


//! @brief returns the number of nodes
//! @return number of nodes
size_t NuTo::StructureGrid::GetNumNodes() const
{
    return mEdgeId.size();
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
    if (mVerboseLevel>3)
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

//! @brief NodeGetConstraintSwitch true is constraint
//! @param rGlobalDof
//! @return bool ... switch for constraint
bool NuTo::StructureGrid::NodeGetConstraintSwitch(int rGlobalDof)
{
	return mDofIsConstraint[rGlobalDof];
}
