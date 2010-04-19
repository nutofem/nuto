#include <boost/tokenizer.hpp>

#include "nuto/mechanics/structures/grid/StructureGrid.h"
#include "nuto/math/FullMatrix.h"
#include <sstream>

#include "nuto/mechanics/nodes/NodeGridCoordinates.h"
#include "nuto/mechanics/nodes/NodeGridCoordinatesDisplacements.h"
#include "nuto/mechanics/nodes/NodeDisplacements.h"
#include "nuto/mechanics/nodes/NodeTemperatures.h"

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

//! @brief gives the identifier of a node
//! @param reference to a node
//! @return identifier
int NuTo::StructureGrid::NodeGetId(const NodeBase* rNode)const
{
    int nodeNumber(0);
    boost::ptr_vector<NodeGridCoordinates>::const_iterator it;
    for (it = mNodeVec.begin(); it!= mNodeVec.end(); it++,nodeNumber++)
    {
        if (&(*it)==rNode)
        {
            break;
        }
    }
    if (it== mNodeVec.end())
        throw MechanicsException("[NuTo::StructureGrid::GetNodeId] Node does not exist.");
    const NuTo::NodeGridCoordinates* myNode = dynamic_cast<const NuTo::NodeGridCoordinates*>(rNode);
    return  myNode->GetNodeID();
    //   return ((const NuTo::NodeGridCoordinates*)rNode)->GetNodeID();
   //return rNode->GetNumCoordinates();
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
    std::cout<<"numVoxels: "<<numVoxel <<std::endl;
   NuTo::FullMatrix<int> imageValues (numVoxel,1);
    imageValues.FullMatrix<int>::ImportFromVtkASCIIFile(mImageDataFile);
    for (unsigned int countNodes =0; countNodes<numGridNodes;countNodes++)//countNodes correpond to nodeID
    {
         //get coincident voxels for each node, check if one voxel has material, then create node
         std::cout<<"geht in coincident routine!"<<std::endl;
         TCoincidentVoxelList coincidentVoxels=GetCoincidenceVoxelIDs(countNodes);
         int flag=0;
         for (int count =0; count<8; count++)
         {
             // voxel exist (for boundary nodes)
             if (coincidentVoxels[count]>-1)
             {
                   // voxel has material
                 std::cout<<"VoxelValue: "<<imageValues(coincidentVoxels[count],0)<<std::endl;

                 if(imageValues(coincidentVoxels[count],0)>180)//@TODO replace with variable for material boundary values
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
             std::cout<<"Knoten erstellt Nr.:"<<numMatNodes-1<<std::endl;
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
    std::cout<<"numDim3 "<<numDim3<< "; numDim2 "<<numDim2<< std::endl;


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
    {
        std::cout<<"Schleife rNodeId =0"<<std::endl;
        coincidentVoxels[7]=-1;
    }
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
    std::cout<<"Knoten"<<rNodeID<< "Voxels:"<<std::endl;
    for (int count=0;count<8;count++)
        std::cout<< coincidentVoxels[count]<<" ";
    std::cout<<" "<<std::endl;
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
    int attributes(1 << NodeBase::COORDINATES);
    //! bit 0 : coordinates
    //! bit 1 : displacements
    //! bit 2 : rotations
    //! bit 3 : temperatures
    boost::tokenizer<> tok(rDOFs);
    for (boost::tokenizer<>::iterator beg=tok.begin(); beg!=tok.end(); ++beg)
    {
        if (*beg=="DISPLACEMENTS")
            attributes = attributes |  1 << NodeBase::DISPLACEMENTS;
        if (*beg=="ROTATIONS")
            attributes = attributes |  1 << NodeBase::ROTATIONS;
        if (*beg=="TEMPERATURES")
            attributes = attributes |  1 << NodeBase::TEMPERATURES;
    }

    NodeBase* nodePtr;
    switch (attributes)
    {
        // the << shifts the 1 bitwise to the left, so 1<<n = 2^n
        // it actually sets the n-th bit (from the right) to 1, and all the other to zero
    case (1 << NodeBase::COORDINATES):
        // for grid nodes COORDINATES are replaced with NodeID
        switch (mDimension)
        {
        case 1:
            nodePtr = new NuTo::NodeGridCoordinates();
            break;
        case 2:
            nodePtr = new NuTo::NodeGridCoordinates();
            break;
        case 3:
            nodePtr = new NuTo::NodeGridCoordinates();
            break;
        default:
            throw MechanicsException("[NuTo::StructureGrid::NodeCreate] Dimension of the structure is not valid.");
        }
       break;
    case (1 << NodeBase::COORDINATES) | (1 << NodeBase::DISPLACEMENTS):
        switch (mDimension)
         {
         case 1:
             nodePtr = new NuTo::NodeGridCoordinatesDisplacements<1>();
             break;
         case 2:
             nodePtr = new NuTo::NodeGridCoordinatesDisplacements<2>();
             break;
         case 3:
             nodePtr = new NuTo::NodeGridCoordinatesDisplacements<3>();
             break;
         default:
             throw MechanicsException("[NuTo::StructureGrid::NodeCreate] Dimension of the structure is not valid.");
         }
        break;
//! @Todo: add NodeGridCoordinatesDisplacementsRotations.h
    case (1 << NodeBase::COORDINATES) | (1 << NodeBase::DISPLACEMENTS) | (1 << NodeBase::ROTATIONS):
         throw MechanicsException("[NuTo::StructureGrid::NodeCreate] Rotation not implemented.");
         break;
   }
    //add grid node number to node
   //add node to vector, richtige id ????????????????????????
    this->mNodeVec.push_back(nodePtr);

    //renumbering of dofs for global matrices not required
    this->mNodeNumberingRequired  = false;


}

// extract dof values (e.g. displacements, temperatures to the nodes)
void NuTo::StructureGrid::NodeExtractDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{

    try
    {
        if (this->mNodeNumberingRequired)
        {
           throw MechanicsException("[NuTo::GridStructure::NodeExtractDofValues] a valid dof numbering was not found (build dof numbering using NodeBuildGlobalDofs).");
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
void NuTo::StructureGrid::NodeBuildGlobalDofs()
{
    throw MechanicsException("[NuTo::GridStructure::NodeBuildGlobalDofs] routine is not implemented");
}



