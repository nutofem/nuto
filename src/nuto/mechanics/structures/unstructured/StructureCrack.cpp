// $Id$

#include <sstream>

#include <boost/foreach.hpp>

#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/cracks/CrackExplicit2D.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/MechanicsException.h"

#include "nuto/mechanics/elements/ElementDataConstitutiveIpCrack.h"

//! @brief returns the number of cracks
//! @return number of cracks
unsigned int NuTo::Structure::GetNumCracks() const
{
    return mCrackMap.size();
}

//! @brief returns a reference to a crack
//! @param identifier
//! @return reference to a crack
NuTo::CrackBase* NuTo::Structure::CrackGetCrackPtr(int rIdent)
{
    crackMap_t::iterator it = mCrackMap.find(rIdent);
    if (it!=mCrackMap.end())
        return it->second;
    else
    {
    	std::stringstream message;
    	message << "[NuTo::Structure::CrackGetCrackPtr] Crack with identifier " << rIdent << " does not exist." << std::endl;
    	throw MechanicsException(message.str());
    }
}

//! @brief returns a const reference to a crack
//! @param identifier
//! @return const reference to a crack
const NuTo::CrackBase* NuTo::Structure::CrackGetCrackPtr(int rIdent)const
{
    crackMap_t::const_iterator it = mCrackMap.find(rIdent);
    if (it!=mCrackMap.end())
        return it->second;
    else
    {
    	std::stringstream message;
    	message << "[NuTo::Structure::CrackGetCrackPtr] Crack with identifier " << rIdent << " does not exist." << std::endl;
    	throw MechanicsException(message.str());
    }
}

//! @brief gives the identifier of a crack
//! @param pointer to a crack
//! @return identifier
int NuTo::Structure::CrackGetId(const CrackBase* rCrack)const
{
    for (crackMap_t::const_iterator
            it = mCrackMap.begin(); it!= mCrackMap.end(); it++)
    {
        if (it->second==rCrack)
            return it->first;
    }
    throw MechanicsException("[NuTo::Structure::CrackGetId] Crack does not exist.");
}

//! @brief ... Info routine that prints general information about the cracks
//! @param ... rVerboseLevel describes how detailed the information is
void NuTo::Structure::CrackInfo(int rVerboseLevel)const
{
    std::cout << "number of cracks  : " << mCrackMap.size() << std::endl;
    if (rVerboseLevel>2)
    {
        for (crackMap_t::const_iterator it= mCrackMap.begin(); it!=mCrackMap.end(); it++)
        {
            std::cout << "  Crack " << it->first << std::endl;
            it->second->Info(rVerboseLevel);
            std::cout <<  std::endl;
        }
    }

}

//! @brief ... delete an existing crack
//! @param rIdent ... crack identifier
void NuTo::Structure::CrackDelete(int rIdent)
{
    // find section identifier in map
    crackMap_t::iterator it = this->mCrackMap.find(rIdent);
    if (it == this->mCrackMap.end())
    {
        throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Crack does not exist.");
    }
    else
    {
        this->mCrackMap.erase(it);
    }
}

//! @brief ... create a new crack
//! @param rIdent ... crack identifier
int NuTo::Structure::CrackCreate()
{
	//find unused integer id
    int id(0);
    crackMap_t::iterator it = mCrackMap.find(id);
    while (it!=mCrackMap.end())
    {
        id++;
        it = mCrackMap.find(id);
    }

	// add crack to map

    //! switch dimension of structure
    switch (this->GetDimension())
    {
    case 1:
    	throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Cracks are not implemented for 1D yet!");
    case 2:
    	this->mCrackMap.insert(id, new NuTo::CrackExplicit2D );
    	break;
    case 3:
    	throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Cracks are not implemented for 3D yet!");
    default:
    	throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Wrong dimension of structure found!!!");
    }

	return id;
}

//! @brief ... extends an existing crack
//! @param rIdent ... crack identifier
//! @param rNode ... pointer to the node to be attended to the crack
void NuTo::Structure::CrackPushBack(int rIdent, NuTo::NodeBase* rNode)
{
    switch (this->GetDimension())
    {
    case 2:
		{
			/// find section identifier in map
			boost::ptr_map<int,CrackExplicit2D>::iterator it = this->mCrackMap.find(rIdent);
			if (it == this->mCrackMap.end())
				throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Crack does not exist.");
			else
				(it->second)->PushBack(rNode);
		}
    	break;
    default:
    	throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Not implemented for this dimension of structure!!!");
    }
}

//! @brief ... extends an existing crack
//! @param rIdent ... crack identifier
//! @param rNode ... pointer to the node to be attended to the crack
void NuTo::Structure::CrackPushFront(int rIdent, NuTo::NodeBase* rNode)
{
    switch (this->GetDimension())
    {
    case 2:
		{
			/// find section identifier in map
			boost::ptr_map<int,CrackExplicit2D>::iterator it = this->mCrackMap.find(rIdent);
			if (it == this->mCrackMap.end())
				throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Crack does not exist.");
			else
				(it->second)->PushFront(rNode);
		}
    	break;
    default:
    	throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Not implemented for this dimension of structure!!!");
    }
}


//! @brief ... shortens an existing crack
//! @param rIdent ... crack identifier
void NuTo::Structure::CrackPopBack(int rIdent)
{
    switch (this->GetDimension())
    {
    case 2:
		{
			/// find section identifier in map
			boost::ptr_map<int,CrackExplicit2D>::iterator it = this->mCrackMap.find(rIdent);
			if (it == this->mCrackMap.end())
				throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Crack does not exist.");
			else
				(it->second)->PopBack();
		}
    	break;
    default:
    	throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Not implemented for this dimension of structure!!!");
    }
}

//! @brief ... shortens an existing crack
//! @param rIdent ... crack identifier
void NuTo::Structure::CrackPopFront(int rIdent)
{
    switch (this->GetDimension())
    {
    case 2:
		{
			/// find crack identifier in map
			boost::ptr_map<int,CrackExplicit2D>::iterator it = this->mCrackMap.find(rIdent);
			if (it == this->mCrackMap.end())
				throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Crack does not exist.");
			else
				(it->second)->PopFront();
		}
    	break;
    default:
    	throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Not implemented for this dimension of structure!!!");
    }
}


//! @brief ... store all cracks of a structure in a vector
//! @param rElements ... vector of const crack pointer
void NuTo::Structure::GetCracksTotal(std::vector<const CrackBase*>& rCracks) const
{
	rCracks.reserve(mCrackMap.size());
	crackMap_t::const_iterator CrackIter = this->mCrackMap.begin();
    while (CrackIter != this->mCrackMap.end())
    {
    	rCracks.push_back(CrackIter->second);
    	CrackIter++;
    }
}

//! @brief ... store all cracks of a structure in a vector
//! @param rElements ... vector of crack pointer
void NuTo::Structure::GetCracksTotal(std::vector<CrackBase*>& rCracks)
{
	rCracks.reserve(mCrackMap.size());
	crackMap_t::iterator CrackIter = this->mCrackMap.begin();
    while (CrackIter != this->mCrackMap.end())
    {
    	rCracks.push_back(CrackIter->second);
    	CrackIter++;
    }
}

//! @brief ... merge all cracks to the existing structure
void NuTo::Structure::InitiateCracks()
{
	elementBasePtrSet_t crackedElems;
	NuTo::Structure::InitiateCracks(crackedElems);
}

//! @brief ... merge all cracks to the existing structure
//! @param rCrackedElems (Output) ... vector of cracked elements
void NuTo::Structure::InitiateCracks(elementBasePtrSet_t & rCrackedElems)
{
	// go through all cracks
	BOOST_FOREACH( crackMap_t::value_type thisCrack, mCrackMap )
		InitiateCrack(thisCrack.first, rCrackedElems);
}

//! @brief ... merge specified crack to the existing structure
//! @param rIdent ... crack identifier
//! @param rCrackedElems (Output) ... vector of cracked elements
void NuTo::Structure::InitiateCrack(const int rIdent)
{
	elementBasePtrSet_t crackedElems;
	NuTo::Structure::InitiateCrack(rIdent, crackedElems);
}

//! @brief ... merge specified crack to the existing structure
//! @param rIdent ... crack identifier
//! @param rCrackedElems (Output) ... vector of cracked elements
void NuTo::Structure::InitiateCrack(const int rIdent, elementBasePtrSet_t & rCrackedElems)
{
	//! find crack identifier in map
	NuTo::CrackBase* thisCrack = this->CrackGetCrackPtr(rIdent);

	//! initiate crack into given elements
	elementBasePtrVec_t allElems;
	this->GetElementsTotal ( allElems );

	//! now the crack searching the cracked elements
	elementBasePtrVec_t rThisCrackedElems;
	thisCrack->Initiate(allElems, rThisCrackedElems);

	BOOST_FOREACH( elementBasePtr_t thisElPtr, rThisCrackedElems )
	{
		if(thisElPtr->GetElementDataType()!=NuTo::ElementData::CONSTITUTIVELAWIPCRACK)
	    	throw NuTo::MechanicsException("[NuTo::Structure::InitiateCrack] Not implemented for this type of ElementDataType!!!");

		ElementDataBase* ptrElementData=thisElPtr->GetDataPtr();
		unsigned int thisCrackId=ptrElementData->AddCrack(thisCrack);

		//! @todo check if rCrackedElems contains rThisCrackedElems already, if not: append
		rCrackedElems.insert(thisElPtr);

		std::cout << "Element to be cracked: " << ElementGetId(thisElPtr)  << "; mElementData=" << ptrElementData
				<< "; mIntegrationType=" << ptrElementData->GetIntegrationType()
				<< "; mElementDataType=" << thisElPtr->GetElementDataType()
				<< "; thisCrackId=" << thisCrackId
		<< std::endl;

	}
}
