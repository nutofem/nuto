#include <assert.h>
#include <boost/tokenizer.hpp>
#include "nuto/mechanics/elements/Truss1D2N.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/elements/Brick8N.h"
#include "nuto/mechanics/elements/Plane2D4N.h"
#include "nuto/mechanics/elements/Truss1D3N.h"
#include "nuto/mechanics/elements/Tetrahedron10N.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

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
	return ElementCreate(rElementType,rNodeNumbers,std::string("CONSTITUTIVELAWELEMENT_NOSTATICDATA"));
}

//! @brief Creates an element
//! @param rElementIdent identifier for the element
//! @param rElementType element type
//! @param rNodeIdents Identifier for the corresponding nodes
int NuTo::Structure::ElementCreate (const std::string& rElementType,
        const NuTo::FullMatrix<int>& rNodeNumbers, const std::string& rElementDataType)
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
	this->ElementCreate(elementNumber, rElementType, rNodeNumbers,rElementDataType);

	// return element number
	return elementNumber;
}
void NuTo::Structure::ElementCreate (int rElementNumber, const std::string& rElementType,
        const NuTo::FullMatrix<int> &rNodeNumbers)
{
    std::cout<<__FILE__<<" "<<" test 6"<<std::endl;
	ElementCreate(rElementNumber,rElementType,rNodeNumbers,std::string("CONSTITUTIVELAWELEMENT_NOSTATICDATA"));
}

void NuTo::Structure::ElementCreate (int rElementNumber, const std::string& rElementType,
        const NuTo::FullMatrix<int> &rNodeNumbers, const std::string& rElementDataType)
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

    ElementBase::eElementType elementType;
    if (upperCaseElementType=="TRUSS1D2N")
    {
    	elementType = NuTo::ElementBase::TRUSS1D2N;
    }
    else if (upperCaseElementType=="TRUSS1D3N")
	{
    	elementType = NuTo::ElementBase::TRUSS1D3N;
	}
    else if (upperCaseElementType=="BRICK8N")
    {
    	elementType = NuTo::ElementBase::BRICK8N;
    }
    else if (upperCaseElementType=="PLANE2D3N")
    {
    	elementType = NuTo::ElementBase::PLANE2D3N;
    }
    else if (upperCaseElementType=="PLANE2D4N")
    {
    	elementType = NuTo::ElementBase::PLANE2D4N;
    }
    else if (upperCaseElementType=="PLANE2D6N")
    {
    	elementType = NuTo::ElementBase::PLANE2D6N;
    }
    else if (upperCaseElementType=="TETRAHEDRON10N")
    {
    	elementType = NuTo::ElementBase::TETRAHEDRON10N;
    }
    else
    {
    	throw MechanicsException("[NuTo::Structure::ElementCreate] Element type "+upperCaseElementType +" does not exist.");
    }

    // check element data
    std::string upperCaseElementDataType;
    std::transform(rElementDataType.begin(), rElementDataType.end(), std::back_inserter(upperCaseElementDataType), (int(*)(int)) toupper);

    ElementDataBase::eElementDataType elementDataType;
    if (upperCaseElementDataType=="CONSTITUTIVELAWELEMENT_NOSTATICDATA")
    {
    	elementDataType = NuTo::ElementDataBase::CONSTITUTIVELAWELEMENT_NOSTATICDATA;
    }
    else if (upperCaseElementDataType=="CONSTITUTIVELAWELEMENT_STATICDATA")
   	{
		elementDataType = NuTo::ElementDataBase::CONSTITUTIVELAWELEMENT_STATICDATA;
   	}
	else if (upperCaseElementType=="CONSTITUTIVELAWIP_NOSTATICDATA")
	{
		elementDataType = NuTo::ElementDataBase::CONSTITUTIVELAWIP_NOSTATICDATA;
	}
	else if (upperCaseElementType=="CONSTITUTIVELAWIP_STATICDATA")
	{
		elementDataType = NuTo::ElementDataBase::CONSTITUTIVELAWIP_STATICDATA;
	}
	else
	{
		throw MechanicsException("[NuTo::Structure::ElementCreate] Element data type "+upperCaseElementDataType +" does not exist.");
	}

    // create element
    this->ElementCreate(rElementNumber, elementType, nodeVector, elementDataType);

}

//! @brief Creates an element
//! @param rElementIdent identifier for the element
//! @param rElementType element type
//! @param rNodeIdents pointers to the corresponding nodes
//! @return element number
void NuTo::Structure::ElementCreate(int rElementNumber, ElementBase::eElementType rType,
        std::vector<NodeBase*> rNodeVector, ElementDataBase::eElementDataType rElementDataType)
{

	//const IntegrationTypeBase *ptrIntegrationType;
	ElementBase* ptrElement;
    switch (rType)
	{
    case NuTo::ElementBase::TRUSS1D2N:
		// get the integration type pointer, if not existent, create the integration type
		//ptrIntegrationType = GetPtrIntegrationType(NuTo::Truss1D2N::GetStandardIntegrationType());
		if (1!=mDimension)
			throw MechanicsException("[NuTo::Structure::ElementCreate] TRUSS1D2N is only a 1D element, either change the dimension of the structure to one or use TRUSS3D2N.");
		ptrElement = new NuTo::Truss1D2N(this, rNodeVector, rElementDataType);
		break;
    case NuTo::ElementBase::TRUSS1D3N:
		// get the integration type pointer, if not existent, create the integration type
		//ptrIntegrationType = GetPtrIntegrationType(NuTo::Truss1D3N::GetStandardIntegrationType());
		if (1!=mDimension)
			throw MechanicsException("[NuTo::Structure::ElementCreate] TRUSS1D3N is only a 1D element, either change the dimension of the structure to one or use TRUSS3D3N.");
		ptrElement = new NuTo::Truss1D3N(this, rNodeVector, rElementDataType);
		break;
    case NuTo::ElementBase::BRICK8N:
        // get the integration type pointer, if not existent, create the integration type
        //ptrIntegrationType = GetPtrIntegrationType(NuTo::Brick8N::GetStandardIntegrationType());
        if (this->mDimension != 3)
        {
            throw MechanicsException("[NuTo::Structure::ElementCreate] Brick8N is a 3D element.");
        }
        ptrElement = new NuTo::Brick8N(this, rNodeVector, rElementDataType);
        break;
    case NuTo::ElementBase::PLANE2D3N:
        if (this->mDimension != 2)
        {
            throw MechanicsException("[NuTo::Structure::ElementCreate] PLANE2D3N is a 2D element.");
        }
        //ptrElement = new NuTo::Plane2D3N(this, rNodeVector, rElementDataType);
        break;
    case NuTo::ElementBase::PLANE2D4N:
        if (this->mDimension != 2)
        {
            throw MechanicsException("[NuTo::Structure::ElementCreate] PLANE2D4N is a 2D element.");
        }
        ptrElement = new NuTo::Plane2D4N(this, rNodeVector, rElementDataType);
        break;
    case NuTo::ElementBase::PLANE2D6N:
        if (this->mDimension != 2)
        {
            throw MechanicsException("[NuTo::Structure::ElementCreate] PLANE2D6N is a 2D element.");
        }
        //ptrElement = new NuTo::Plane2D6N(this, rNodeVector, rElementDataType);
        break;
    case NuTo::ElementBase::TETRAHEDRON10N:
        // get the integration type pointer, if not existent, create the integration type
        //ptrIntegrationType = GetPtrIntegrationType(NuTo::Tetrahedron10N::GetStandardIntegrationType());
        if (this->mDimension != 3)
        {
            throw MechanicsException("[NuTo::Structure::ElementCreate] Tetrahedron10N is a 3D element.");
        }
        ptrElement = new NuTo::Tetrahedron10N(this, rNodeVector, rElementDataType);
        break;
	default:
		throw NuTo::MechanicsException("[NuTo::Structure::ElementCreate] Invalid element type.");
	}
    mElementMap.insert(rElementNumber, ptrElement);
}
