// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/vector.hpp>
#endif //ENABLE_SERIALIZATION

#include "nuto/mechanics/structures/grid/StructureGrid.h"
#include "nuto/mechanics/elements/ElementBase.h"
//#include "nuto/mechanics/elements/Voxel8N.h"
#include <sstream>

//! @brief returns the number of nodes
//! @return number of nodes
int NuTo::StructureGrid::GetNumElements() const
{
    return mVoxelId.size();
}



//! @brief info about the elements in the StructureGrid
void NuTo::StructureGrid::ElementInfo(int mVerboseLevel)const
{
//    std::cout<<"number of elements: " << mElementMap.size() <<std::endl;
    std::cout<<"number of elements: " << mVoxelId.size() <<std::endl;
}

//! @brief create element grid without data free elements
//! @param reference to a base coefficient matrix, to a ColorToMaterialMatrix and to an element type
void NuTo::StructureGrid::CreateElementGrid( NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rBaseCoefficientMatrix0,
const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rColorToMaterialData,const std::string& rElementType)
{
    throw MechanicsException("[NuTo::StructureGrid::CreateElementGrid] Not implemented.");
}
