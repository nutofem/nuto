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

//! @brief a reference to an element
//! @param identifier
//! @return reference to an element
NuTo::ElementBase* NuTo::StructureGrid::ElementGetElementPtr(int rIdent)
{
    throw MechanicsException("[NuTo::StructureGrid::ElementGetElementPtr] not implemented.");
//    boost::ptr_map<int,ElementBase>::iterator it = mElementMap.find(rIdent);
//    if (it!=mElementMap.end())
//        return it->second;
//    else
//    {
//    	std::stringstream message;
//    	message << "[NuTo::StructureGrid::ElementGetElementPtr] Element with identifier " << rIdent << " does not exist." << std::endl;
//    	throw MechanicsException(message.str());
//    }
}

//! @brief a reference to an element
//! @param identifier
//! @return reference to an element
const NuTo::ElementBase* NuTo::StructureGrid::ElementGetElementPtr(const int rIdent)const
{
    throw MechanicsException("[NuTo::StructureGrid::ElementGetElementPtr] not implemented.");
//    boost::ptr_map<int,ElementBase>::const_iterator it = mElementMap.find(rIdent);
//    if (it!=mElementMap.end())
//        return it->second;
//    else
//    {
//    	std::stringstream message;
//    	message << "[NuTo::StructureGrid::ElementGetElementPtr] Element with identifier " << rIdent << " does not exist." << std::endl;
//    	throw MechanicsException(message.str());
//    }
}

//! @brief gives the identifier of an element
//! @param reference to an element
//! @return element number
int NuTo::StructureGrid::ElementGetId(const ElementBase* rElement)const
{
//    for (boost::ptr_map<int,ElementBase>::const_iterator
//            it = mElementMap.begin(); it!= mElementMap.end(); it++)
//    {
//        if (it->second==rElement)
//            return it->first;
//    }
//	    throw MechanicsException("[NuTo::StructureGrid::ElementGetId] Element does not exist.");
    throw MechanicsException("[NuTo::StructureGrid::ElementGetId] not implemented.");
}

//! @brief info about the elements in the StructureGrid
void NuTo::StructureGrid::ElementInfo(int mVerboseLevel)const
{
//    std::cout<<"number of elements: " << mElementMap.size() <<std::endl;
    std::cout<<"number of elements: " << mVoxelId.size() <<std::endl;
}

//! @brief create element grid without data free elements
//! @param reference to a base coefficient matrix, to a ColorToMaterialMatrix and to an element type
void NuTo::StructureGrid::CreateElementGrid( NuTo::FullMatrix<double>& rBaseCoefficientMatrix0,
const NuTo::FullMatrix<double>& rColorToMaterialData,const std::string& rElementType)
{
    throw MechanicsException("[NuTo::StructureGrid::CreateElementGrid] Not implemented.");
}

//! @brief Creates an element
//! @param rElementIdent identifier for the element
//! @param rElementType element type
//! @param rNodeIdents Identifier for the corresponding nodes
int NuTo::StructureGrid::ElementCreate (int materialNumber)
{
   	throw MechanicsException("[NuTo::StructureGrid::ElementCreate] Not implemented.");
//	//find unused integer id
//	int elementNumber(mElementMap.size());
//	boost::ptr_map<int,ElementBase>::iterator it = mElementMap.find(elementNumber);
//	while (it!=mElementMap.end())
//	{
//		elementNumber++;
//		it = mElementMap.find(elementNumber);
//	}
//
//	// create element
//	this->ElementCreate(elementNumber, materialNumber);
//
//	// return element number
//	return elementNumber;
}
void NuTo::StructureGrid::ElementCreate (int rElementNumber, int rMaterialNumber)
{
   	throw MechanicsException("[NuTo::StructureGrid::ElementCreate] Not implemented.");
//	// check node number
//	boost::ptr_map<int,ElementBase>::iterator it = mElementMap.find(rElementNumber);
//	if(it != this->mElementMap.end())
//	{
//    	throw MechanicsException("[NuTo::Structure::ElementCreate] Element already exists.");
//	}
//	if(mDimension!=3)
//	   	throw MechanicsException("[NuTo::StructureGrid::ElementCreate] Only 3D implemented.");
//
//  	ElementBase* ptrElement=0;
//  	ptrElement= new NuTo::Voxel8N(this,  rMaterialNumber);
// 	 mElementMap.insert(rElementNumber, ptrElement);
}

//! @brief Deletes an element
//! @param rElementIdent identifier for the element
void NuTo::StructureGrid::ElementDelete(int rElementNumber)
{

	// @TODO [NuTo::StructureGrid::ElementDelete] has to be implemented
    throw MechanicsException("[NuTo::StructureGrid::ElementDelete] Not implemented yet!!!");

}



