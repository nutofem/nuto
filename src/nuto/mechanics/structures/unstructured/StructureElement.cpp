// $Id$

#include <assert.h>
#include <boost/tokenizer.hpp>
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/elements/Truss1D2N.h"
#include "nuto/mechanics/elements/Brick8N.h"
#include "nuto/mechanics/elements/Plane2D3N.h"
#include "nuto/mechanics/elements/Plane2D6N.h"
#include "nuto/mechanics/elements/Plane2D4N.h"
#include "nuto/mechanics/elements/Truss1D3N.h"
#include "nuto/mechanics/elements/Tetrahedron10N.h"
#include "nuto/mechanics/groups/Group.h"
//! @brief returns the number of nodes
//! @return number of nodes
int NuTo::Structure::GetNumElements() const
{
    return mElementMap.size();
}

//! @brief a reference to an element
//! @param identifier
//! @return reference to an element
NuTo::ElementBase* NuTo::Structure::ElementGetElementPtr(int rIdent)
{
    boost::ptr_map<int,ElementBase>::iterator it = mElementMap.find(rIdent);
    if (it!=mElementMap.end())
        return it->second;
    else
    {
    	std::stringstream message;
    	message << "[NuTo::Structure::ElementGetElementPtr] Element with identifier " << rIdent << " does not exist." << std::endl;
    	throw MechanicsException(message.str());
    }
}

//! @brief a reference to an element
//! @param identifier
//! @return reference to an element
const NuTo::ElementBase* NuTo::Structure::ElementGetElementPtr(int rIdent)const
{
    boost::ptr_map<int,ElementBase>::const_iterator it = mElementMap.find(rIdent);
    if (it!=mElementMap.end())
        return it->second;
    else
    {
    	std::stringstream message;
    	message << "[NuTo::Structure::ElementGetElementPtr] Element with identifier " << rIdent << " does not exist." << std::endl;
    	throw MechanicsException(message.str());
    }
}

//! @brief gives the identifier of an element
//! @param reference to an element
//! @return element number
int NuTo::Structure::ElementGetId(const ElementBase* rElement)const
{
    for (boost::ptr_map<int,ElementBase>::const_iterator
            it = mElementMap.begin(); it!= mElementMap.end(); it++)
    {
        if (it->second==rElement)
            return it->first;
    }
    throw MechanicsException("[NuTo::Structure::GetElementId] Element does not exist.");
}

//! @brief info about one single element
//! @param rElement (Input) ... pointer to the element
//! @param rVerboseLevel (Input) ... level of verbosity
void NuTo::Structure::ElementInfo(const ElementBase* rElement, int rVerboseLevel)const
{
	std::cout << "element : " << rElement->ElementGetId() << std::endl;
	if (rVerboseLevel>2)
	{
		std::cout << "\tenum::type=" << rElement->GetEnumType() << std::endl;
		if (rVerboseLevel>3)
		{
			std::cout << "\tNodes:";
			for(unsigned short iNode=0; iNode<rElement->GetNumNodes(); ++iNode)
				std::cout << "\t" << this->NodeGetId(rElement->GetNode(iNode));
			std::cout << std::endl;
			if (rVerboseLevel>4)
			{
				std::cout << "\tintegration points :" << std::endl;
				for(int iIp=0; iIp<rElement->GetNumIntegrationPoints(); ++iIp)
				{
					double coor[3];
					rElement->GetGlobalIntegrationPointCoordinates(iIp,coor);
					std::cout << "\t\t" << iIp << ": [" << coor[0] << ";" << coor[1] << ";" << coor[2] << "]" << std::endl;
				}
			}
		}
	}
}

//! @brief info about the elements in the Structure
void NuTo::Structure::ElementInfo(int rVerboseLevel)const
{
    std::cout<<"number of elements: " << mElementMap.size() <<std::endl;
    if (rVerboseLevel>3)
    {
    	std::cout << "\t\telements :" <<std::endl;
    	for (boost::ptr_map<int,ElementBase>::const_iterator it = mElementMap.begin(); it!= mElementMap.end(); it++)
    	{
    		std::cout << "\t\t" << it->first;
    	    if (rVerboseLevel>4)
    	    {
				std::cout << "\t:";
				for(unsigned short iNode=0; iNode<it->second->GetNumNodes(); ++iNode)
					std::cout << "\t" << this->NodeGetId(it->second->GetNode(iNode));
    	    }
        	std::cout << std::endl;
    	}
    }
}

//! @brief Creates an element
//! @param rElementIdent identifier for the element
//! @param rElementType element type
//! @param rNodeIdents Identifier for the corresponding nodes
int NuTo::Structure::ElementCreate (const std::string& rElementType,
        const NuTo::FullVector<int,Eigen::Dynamic>& rNodeNumbers)
{
	return ElementCreate(rElementType,rNodeNumbers,std::string("CONSTITUTIVELAWIP"),std::string("NOIPDATA") );
}

//! @brief Creates an element
//! @param rElementIdent identifier for the element
//! @param rElementType element type
//! @param rNodeIdents Identifier for the corresponding nodes
int NuTo::Structure::ElementCreate (const std::string& rElementType,
        const NuTo::FullVector<int,Eigen::Dynamic>& rNodeNumbers, const std::string& rElementDataType, const std::string& rIpDataType)
{
	//find unused integer id
	int elementNumber(mElementMap.size());
	boost::ptr_map<int,ElementBase>::iterator it = mElementMap.find(elementNumber);
	while (it!=mElementMap.end())
	{
		elementNumber++;
		it = mElementMap.find(elementNumber);
	}

	// create element
	this->ElementCreate(elementNumber, rElementType, rNodeNumbers,rElementDataType,rIpDataType);

	// return element number
	return elementNumber;
}
void NuTo::Structure::ElementCreate (int rElementNumber, const std::string& rElementType,
        const NuTo::FullVector<int,Eigen::Dynamic> &rNodeNumbers)
{
	ElementCreate(rElementNumber,rElementType,rNodeNumbers,std::string("CONSTITUTIVELAWIP"),std::string("NOIPDATA"));
}

void NuTo::Structure::ElementCreate (int rElementNumber, const std::string& rElementType,
        const NuTo::FullVector<int,Eigen::Dynamic> &rNodeNumbers, const std::string& rElementDataType, const std::string& rIpDataType)
{
	// check node number
	boost::ptr_map<int,ElementBase>::iterator it = mElementMap.find(rElementNumber);
	if(it != this->mElementMap.end())
	{
    	throw MechanicsException("[NuTo::Structure::ElementCreate] Element already exists.");
	}

	// convert node numbers to pointer
	if (rNodeNumbers.GetNumColumns()!=1)
    	throw MechanicsException("[NuTo::Structure::ElementCreate] Matrix with node numbers should have a single column.");
	std::vector<NodeBase*> nodeVector;
    for (int count=0; count<rNodeNumbers.GetNumRows(); count++)
    {
        nodeVector.push_back(NodeGetNodePtr(rNodeNumbers(count,0)));
    }

    // get element type
    std::string upperCaseElementType;
    std::transform(rElementType.begin(), rElementType.end(), std::back_inserter(upperCaseElementType), (int(*)(int)) toupper);

    Element::eElementType elementType;
    if (upperCaseElementType=="TRUSS1D2N")
    {
    	elementType = NuTo::Element::TRUSS1D2N;
    }
    else if (upperCaseElementType=="TRUSS1D3N")
	{
    	elementType = NuTo::Element::TRUSS1D3N;
	}
    else if (upperCaseElementType=="BRICK8N")
    {
    	elementType = NuTo::Element::BRICK8N;
    }
    else if (upperCaseElementType=="PLANE2D3N")
    {
    	elementType = NuTo::Element::PLANE2D3N;
    }
    else if (upperCaseElementType=="PLANE2D4N")
    {
    	elementType = NuTo::Element::PLANE2D4N;
    }
    else if (upperCaseElementType=="PLANE2D6N")
    {
    	elementType = NuTo::Element::PLANE2D6N;
    }
    else if (upperCaseElementType=="TETRAHEDRON10N")
    {
    	elementType = NuTo::Element::TETRAHEDRON10N;
    }
    else
    {
    	throw MechanicsException("[NuTo::Structure::ElementCreate] Element type "+upperCaseElementType +" does not exist.");
    }

    // check element data
    std::string upperCaseElementDataType;
    std::transform(rElementDataType.begin(), rElementDataType.end(), std::back_inserter(upperCaseElementDataType), (int(*)(int)) toupper);

    ElementData::eElementDataType elementDataType;
    if (upperCaseElementDataType=="CONSTITUTIVELAWIP")
    {
    	elementDataType = NuTo::ElementData::CONSTITUTIVELAWIP;
    }
    else if (upperCaseElementDataType=="CONSTITUTIVELAWIPCRACK")
   	{
		elementDataType = NuTo::ElementData::CONSTITUTIVELAWIPCRACK;
   	}
    else if (upperCaseElementDataType=="CONSTITUTIVELAWIPNONLOCAL")
   	{
		elementDataType = NuTo::ElementData::CONSTITUTIVELAWIPNONLOCAL;
   	}
    else
	{
		throw MechanicsException("[NuTo::Structure::ElementCreate] Element data type "+upperCaseElementDataType +" does not exist.");
	}

    // check ip data
    std::string upperCaseIpDataType;
    std::transform(rIpDataType.begin(), rIpDataType.end(), std::back_inserter(upperCaseIpDataType), (int(*)(int)) toupper);

    IpData::eIpDataType ipDataType;
    if (upperCaseIpDataType=="NOIPDATA")
    {
    	ipDataType = NuTo::IpData::NOIPDATA;
    }
    else if (upperCaseIpDataType=="STATICDATA")
   	{
    	ipDataType = NuTo::IpData::STATICDATA;
   	}
    else if (upperCaseIpDataType=="STATICDATANONLOCAL")
   	{
    	ipDataType = NuTo::IpData::STATICDATANONLOCAL;
   	}
	else
	{
		throw MechanicsException("[NuTo::Structure::ElementCreate] Ip data type "+upperCaseIpDataType +" does not exist.");
	}

    // create element
    this->ElementCreate(rElementNumber, elementType, nodeVector, elementDataType, ipDataType);

}

//! @brief Creates an element
//! @param rElementNumber element number
//! @param rElementType element type
//! @param rNodeIdents pointers to the corresponding nodes
//! @return int rElementNumber
int NuTo::Structure::ElementCreate(Element::eElementType rType,
        std::vector<NodeBase*> rNodeVector, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType)
{
	//find unused integer id
	int elementNumber(mElementMap.size());
	boost::ptr_map<int,ElementBase>::iterator it = mElementMap.find(elementNumber);
	while (it!=mElementMap.end())
	{
		elementNumber++;
		it = mElementMap.find(elementNumber);
	}

	// create element
	this->ElementCreate(elementNumber,rType,rNodeVector,rElementDataType,rIpDataType);

	// return element number
	return elementNumber;
}
//! @brief Creates an element
//! @param rElementIdent identifier for the element
//! @param rElementType element type
//! @param rNodeIdents pointers to the corresponding nodes
//! @return element number
void NuTo::Structure::ElementCreate(int rElementNumber, Element::eElementType rType,
        std::vector<NodeBase*> rNodeVector, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType)
{
	ElementBase* ptrElement(0);
    switch (rType)
	{
    case NuTo::Element::TRUSS1D2N:
		if (1!=mDimension)
			throw MechanicsException("[NuTo::Structure::ElementCreate] TRUSS1D2N is only a 1D element, either change the dimension of the structure to one or use TRUSS3D2N.");
		ptrElement = new NuTo::Truss1D2N(this, rNodeVector, rElementDataType, rIpDataType);
		break;
    case NuTo::Element::TRUSS1D3N:
		if (1!=mDimension)
			throw MechanicsException("[NuTo::Structure::ElementCreate] TRUSS1D3N is only a 1D element, either change the dimension of the structure to one or use TRUSS3D3N.");
		ptrElement = new NuTo::Truss1D3N(this, rNodeVector, rElementDataType, rIpDataType);
		break;
    case NuTo::Element::BRICK8N:
        if (this->mDimension != 3)
        {
            throw MechanicsException("[NuTo::Structure::ElementCreate] Brick8N is a 3D element.");
        }
        ptrElement = new NuTo::Brick8N(this, rNodeVector, rElementDataType, rIpDataType);
        break;
    case NuTo::Element::PLANE2D3N:
        if (this->mDimension != 2)
        {
            throw MechanicsException("[NuTo::Structure::ElementCreate] PLANE2D3N is a 2D element.");
        }
        ptrElement = new NuTo::Plane2D3N(this, rNodeVector, rElementDataType, rIpDataType);
        break;
    case NuTo::Element::PLANE2D4N:
        if (this->mDimension != 2)
        {
            throw MechanicsException("[NuTo::Structure::ElementCreate] PLANE2D4N is a 2D element.");
        }
        ptrElement = new NuTo::Plane2D4N(this, rNodeVector, rElementDataType, rIpDataType);
        break;
    case NuTo::Element::PLANE2D6N:
        if (this->mDimension != 2)
        {
            throw MechanicsException("[NuTo::Structure::ElementCreate] PLANE2D6N is a 2D element.");
        }
        ptrElement = new NuTo::Plane2D6N(this, rNodeVector, rElementDataType, rIpDataType);
        break;
    case NuTo::Element::TETRAHEDRON10N:
        if (this->mDimension != 3)
        {
            throw MechanicsException("[NuTo::Structure::ElementCreate] Tetrahedron10N is a 3D element.");
        }
        ptrElement = new NuTo::Tetrahedron10N(this, rNodeVector, rElementDataType, rIpDataType);
        break;
	default:
		throw NuTo::MechanicsException("[NuTo::Structure::ElementCreate] Invalid element type.");
	}
    mElementMap.insert(rElementNumber, ptrElement);
}

//! @brief creates multiple elements
//! @param rElementType element type
//! @param rNodeIdents Identifier for the corresponding nodes (Incidences have to be stored column-wise)
//! @return a NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> containing the element numbers
NuTo::FullVector<int,Eigen::Dynamic> NuTo::Structure::ElementsCreate (const std::string& rElementType, NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> & rNodeNumbers)
{
	std::vector<int> idVec;
	/// go through the elements
	for(size_t i=0 ; i<(size_t)rNodeNumbers.GetNumColumns(); ++i)
	{
		const NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> incidence(rNodeNumbers.GetColumn(i));
		idVec.push_back(this-> ElementCreate(rElementType,incidence ));
	}

    //return int identifiers of the new elements as FullMatrix
	NuTo::FullVector<int,Eigen::Dynamic> ids(idVec);
    return ids;
}

//! @brief creates multiple elements
//! @param rElementType element type
//! @param rNodeIdents Identifier for the corresponding nodes (Incidences have to be stored column-wise)
//! @param rElementDataType Element data for the elements
//! @param rIpDataType Integration point data for the elements
//! @return a NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> containing the element numbers
NuTo::FullVector<int,Eigen::Dynamic> NuTo::Structure::ElementsCreate (const std::string& rElementType, NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> & rNodeNumbers, const std::string& rElementDataType, const std::string& rIpDataType)
{
	std::vector<int> idVec;
	/// go through the elements
	for(size_t i=0 ; i<(size_t)rNodeNumbers.GetNumColumns(); ++i)
	{
		const NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> incidence(rNodeNumbers.GetColumn(i));
		idVec.push_back(this-> ElementCreate(rElementType, incidence, rElementDataType, rIpDataType));
	}

    //return int identifiers of the new elements as FullMatrix
	NuTo::FullVector<int,Eigen::Dynamic> ids(idVec);
    return ids;
}

//! @brief Deletes an element
//! @param rElementIdent identifier for the element
void NuTo::Structure::ElementDelete(int rElementNumber)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
   	ElementDeleteInternal(rElementNumber);
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::Structure::ElementDelete] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief Deletes a group of elements element
//! @param rGroupNumber group number
void NuTo::Structure::ElementGroupDelete (int rGroupNumber, bool deleteNodes)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	boost::ptr_map<int,GroupBase>::iterator itGroup = mGroupMap.find(rGroupNumber);
    if (itGroup==mGroupMap.end())
        throw MechanicsException("[NuTo::StructureBase::ElementGroupDelete] Group with the given identifier does not exist.");
    if (itGroup->second->GetType()!=NuTo::Groups::Elements)
    	throw MechanicsException("[NuTo::StructureBase::ElementGroupDelete] Group is not an element group.");

    //the group has to be copied, since the elements are removed from this group, which invalidates the iterators
    Group<ElementBase> copyOfElementGroup = *(dynamic_cast<Group<ElementBase>*>(itGroup->second));

    std::set<NodeBase*> potentialNodesToBeRemoved;
    for (Group<ElementBase>::iterator itElement=copyOfElementGroup.begin(); itElement!=copyOfElementGroup.end();itElement++)
    {
        try
        {
        	//save the nodes, which are eventually to be removed
        	if (deleteNodes)
        	{
        		for (int countNode=0; countNode<itElement->second->GetNumNodes(); countNode++)
        		{
        			NodeBase* nodePtr = itElement->second->GetNode(countNode);
        			potentialNodesToBeRemoved.insert(nodePtr);

        		}
        	}
        	ElementDeleteInternal(itElement->second->ElementGetId());
        }
        catch(NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
            e.AddMessage("[NuTo::StructureBase::ElementGroupDelete] Error deleting element "
            	+ ss.str() + ".");
            throw e;
        }
        catch(...)
        {
            std::stringstream ss;
            assert(ElementGetId(itElement->second)==itElement->first);
            ss << itElement->first;
        	throw NuTo::MechanicsException
        	   ("[NuTo::StructureBase::ElementGroupDelete] Error deleting element " + ss.str() + ".");
        }
    }

    //check all the other elements and see, if they have one of the potential Nodes To Be Removed as valid node
    for (boost::ptr_map<int,ElementBase>::iterator itElement = mElementMap.begin(); itElement != mElementMap.end(); itElement++)
    {
    	for (int countNode=0; countNode<(itElement->second)->GetNumNodes(); countNode++)
		{
			NodeBase* nodePtr = (itElement->second)->GetNode(countNode);
			//int numRemoved = potentialNodesToBeRemoved.erase(nodePtr);
			potentialNodesToBeRemoved.erase(nodePtr);
		}
    }

    for (std::set<NodeBase*>::iterator itNode = potentialNodesToBeRemoved.begin(); itNode != potentialNodesToBeRemoved.end(); itNode++)
    {
    	int nodeId(NodeGetId(*itNode));
    	NodeDelete(nodeId,false);
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureBase::ElementGroupDelete] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

//! @brief Deletes an element
//! @param rItElement iterator of the map
void NuTo::Structure::ElementDeleteInternal(int rElementId)
{
	// find element
	boost::ptr_map<int,ElementBase>::iterator itElement = mElementMap.find(rElementId);
    if (itElement == this->mElementMap.end())
    {
        throw MechanicsException("[NuTo::Structure::ElementDeleteInternal] Element does not exist.");
    }
    else
    {
		// Search for elements in groups: using a loop over all groups
		for(boost::ptr_map<int,GroupBase>::iterator groupIt=mGroupMap.begin();groupIt!=mGroupMap.end(); ++groupIt)
		{
			if(groupIt->second->GetType()==NuTo::Groups::Elements)
			{
				if(groupIt->second->Contain(rElementId))
				{
					groupIt->second->RemoveMember(rElementId);
				}
			}
		}

		// delete element from map
		this->mElementMap.erase(itElement);
    }
}

// store all elements of a structure in a vector
void NuTo::Structure::GetElementsTotal(std::vector<const ElementBase*>& rElements) const
{
    rElements.reserve(mElementMap.size());
    rElements.resize(0);
	boost::ptr_map<int,ElementBase>::const_iterator ElementIter = this->mElementMap.begin();
    while (ElementIter != this->mElementMap.end())
    {
        rElements.push_back(ElementIter->second);
        ElementIter++;
    }
}

// store all elements of a structure in a vector
void NuTo::Structure::GetElementsTotal(std::vector<std::pair<int, const ElementBase*> >& rElements) const
{
	rElements.reserve(mElementMap.size());
	rElements.resize(0);
	boost::ptr_map<int,ElementBase>::const_iterator ElementIter = this->mElementMap.begin();
    while (ElementIter != this->mElementMap.end())
    {
    	rElements.push_back(std::pair<int, const ElementBase*>(ElementIter->first,ElementIter->second));
        ElementIter++;
    }
}

// store all elements of a structure in a vector
void NuTo::Structure::GetElementsTotal(std::vector<ElementBase*>& rElements)
{
    rElements.reserve(mElementMap.size());
    rElements.resize(0);
    boost::ptr_map<int,ElementBase>::iterator ElementIter = this->mElementMap.begin();
    while (ElementIter != this->mElementMap.end())
    {
        rElements.push_back(ElementIter->second);
        ElementIter++;
    }
}

// store all elements of a structure in a vector
void NuTo::Structure::GetElementsTotal(std::vector<std::pair<int, ElementBase*> >& rElements)
{
	rElements.reserve(mElementMap.size());
	rElements.resize(0);
	boost::ptr_map<int,ElementBase>::iterator ElementIter = this->mElementMap.begin();
    while (ElementIter != this->mElementMap.end())
    {
    	rElements.push_back(std::pair<int, ElementBase*>(ElementIter->first,ElementIter->second));
        ElementIter++;
    }
}

