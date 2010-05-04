#include <assert.h>
#include <boost/tokenizer.hpp>
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/elements/Truss1D2N.h"
#include "nuto/mechanics/elements/Brick8N.h"
#include "nuto/mechanics/elements/Plane2D4N.h"
#include "nuto/mechanics/elements/Truss1D3N.h"
#include "nuto/mechanics/elements/Tetrahedron10N.h"
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

//! @brief info about the elements in the Structure
void NuTo::Structure::ElementInfo(int mVerboseLevel)const
{
    std::cout<<"number of elements: " << mElementMap.size() <<std::endl;
}

//! @brief Creates an element
//! @param rElementIdent identifier for the element
//! @param rElementType element type
//! @param rNodeIdents Identifier for the corresponding nodes
int NuTo::Structure::ElementCreate (const std::string& rElementType,
        const NuTo::FullMatrix<int>& rNodeNumbers)
{
	return ElementCreate(rElementType,rNodeNumbers,std::string("CONSTITUTIVELAWIP"),std::string("NOIPDATA") );
}

//! @brief Creates an element
//! @param rElementIdent identifier for the element
//! @param rElementType element type
//! @param rNodeIdents Identifier for the corresponding nodes
int NuTo::Structure::ElementCreate (const std::string& rElementType,
        const NuTo::FullMatrix<int>& rNodeNumbers, const std::string& rElementDataType, const std::string& rIpDataType)
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
        const NuTo::FullMatrix<int> &rNodeNumbers)
{
	ElementCreate(rElementNumber,rElementType,rNodeNumbers,std::string("CONSTITUTIVELAWIP"),std::string("NOIPDATA"));
}

void NuTo::Structure::ElementCreate (int rElementNumber, const std::string& rElementType,
        const NuTo::FullMatrix<int> &rNodeNumbers, const std::string& rElementDataType, const std::string& rIpDataType)
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
    else if (upperCaseElementDataType=="CONSTITUTIVELAWIP")
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
		throw MechanicsException("[NuTo::Structure::ElementCreate] Element data type "+upperCaseIpDataType +" does not exist.");
	}

    // create element
    this->ElementCreate(rElementNumber, elementType, nodeVector, elementDataType, ipDataType);

}

//! @brief Creates an element
//! @param rElementIdent identifier for the element
//! @param rElementType element type
//! @param rNodeIdents pointers to the corresponding nodes
//! @return element number
void NuTo::Structure::ElementCreate(int rElementNumber, Element::eElementType rType,
        std::vector<NodeBase*> rNodeVector, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType)
{

	//const IntegrationTypeBase *ptrIntegrationType;
	ElementBase* ptrElement;
    switch (rType)
	{
    case NuTo::Element::TRUSS1D2N:
		// get the integration type pointer, if not existent, create the integration type
		//ptrIntegrationType = GetPtrIntegrationType(NuTo::Truss1D2N::GetStandardIntegrationType());
		if (1!=mDimension)
			throw MechanicsException("[NuTo::Structure::ElementCreate] TRUSS1D2N is only a 1D element, either change the dimension of the structure to one or use TRUSS3D2N.");
		ptrElement = new NuTo::Truss1D2N(this, rNodeVector, rElementDataType, rIpDataType);
		break;
    case NuTo::Element::TRUSS1D3N:
		// get the integration type pointer, if not existent, create the integration type
		//ptrIntegrationType = GetPtrIntegrationType(NuTo::Truss1D3N::GetStandardIntegrationType());
		if (1!=mDimension)
			throw MechanicsException("[NuTo::Structure::ElementCreate] TRUSS1D3N is only a 1D element, either change the dimension of the structure to one or use TRUSS3D3N.");
		ptrElement = new NuTo::Truss1D3N(this, rNodeVector, rElementDataType, rIpDataType);
		break;
    case NuTo::Element::BRICK8N:
        // get the integration type pointer, if not existent, create the integration type
        //ptrIntegrationType = GetPtrIntegrationType(NuTo::Brick8N::GetStandardIntegrationType());
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
        //ptrElement = new NuTo::Plane2D3N(this, rNodeVector, rElementDataType, rIpDataType);
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
        //ptrElement = new NuTo::Plane2D6N(this, rNodeVector, rElementDataType, rIpDataType);
        break;
    case NuTo::Element::TETRAHEDRON10N:
        // get the integration type pointer, if not existent, create the integration type
        //ptrIntegrationType = GetPtrIntegrationType(NuTo::Tetrahedron10N::GetStandardIntegrationType());
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
//! @return a NuTo::FullMatrix<int> containing the element numbers
NuTo::FullMatrix<int> NuTo::Structure::ElementsCreate (const std::string& rElementType, NuTo::FullMatrix<int> & rNodeNumbers)
{
	std::vector<int> idVec;
	/// go through the elements
	for(size_t i=0 ; i<(size_t)rNodeNumbers.GetNumColumns(); ++i)
	{
		const NuTo::FullMatrix<int> incidence(rNodeNumbers.GetColumn(i));
		idVec.push_back(this-> ElementCreate(rElementType,incidence ));
	}

    //return int identifiers of the new elements as FullMatrix
	NuTo::FullMatrix<int> ids(idVec);
    return ids;
}
