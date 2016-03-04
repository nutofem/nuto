// $Id$ 
// ElementDataNonlocalBase.cpp
// created Apr 22, 2010 by Joerg F. Unger

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/ElementDataNonlocalBase.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include <assert.h>

NuTo::ElementDataNonlocalBase::ElementDataNonlocalBase() :  NuTo::ElementDataBase::ElementDataBase()
{
}

NuTo::ElementDataNonlocalBase::~ElementDataNonlocalBase()
{
	//std::cout << "NuTo::ElementDataNonlocalBase::~ElementDataNonlocalBase()" << std::endl;
}

const std::vector<const NuTo::ElementBase*>&
  NuTo::ElementDataNonlocalBase::GetNonlocalElements()const
{
    return mNonlocalElements;
}

int NuTo::ElementDataNonlocalBase::GetNumNonlocalElements()const
{
    return mNonlocalElements.size();
}

/*const std::vector<double>&
  NuTo::ElementDataNonlocalBase::GetNonlocalWeights(const ConstitutiveBase* rConstitutive, int rLocalIp, int rNonlocalElement)const
{
    if (rConstitutive==mConstitutive)
    {
	    if (rNonlocalElement>=0 && rNonlocalElement<(int)mNonlocalWeights.size())
    	    return mNonlocalWeights[rNonlocalElement];
	    else
	    	throw MechanicsException("[NuTo::ElementDataNonlocalBase::GetNonlocalElements] Check your nonlocal data.");
    }
    else
    	throw MechanicsException("[NuTo::ElementDataNonlocalBase::GetNonlocalElements] For this constitutive model no nonlocal data is available");
}
*/
//! @brief adds an element to the nonlocal elements
//! @param rConstitutive  constitutive model
//! @return the local element number, the element is either append to the list, or the existing local number is returned
int NuTo::ElementDataNonlocalBase::AddNonlocalElement(const ElementBase* rElement)
{
	for (int theElement=0; theElement<(int)mNonlocalElements.size();theElement++)
	{
		if (mNonlocalElements[theElement]==rElement)
			return theElement;
	}
	mNonlocalElements.push_back(rElement);
	return mNonlocalElements.size()-1;
}


#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ElementDataNonlocalBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ElementDataNonlocalBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ElementDataNonlocalBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ElementDataNonlocalBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ElementDataNonlocalBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ElementDataNonlocalBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);

template<class Archive>
void NuTo::ElementDataNonlocalBase::load(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ElementDataNonlocalBase" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementDataBase);

    int size = 0;
    ar & boost::serialization::make_nvp("mNonlocalElements_size", size);
    std::uintptr_t* mNonlocalElementsAdress = new std::uintptr_t[size];
    ar & boost::serialization::make_nvp("mNonlocalElements", boost::serialization::make_array(mNonlocalElementsAdress, size));
    mNonlocalElements.assign(reinterpret_cast<ElementBase**>(&mNonlocalElementsAdress[0]), reinterpret_cast<ElementBase**>(&mNonlocalElementsAdress[size]));
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ElementDataNonlocalBase" << std::endl;
#endif
}

template<class Archive>
void NuTo::ElementDataNonlocalBase::save(Archive & ar, const unsigned int version) const
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ElementDataNonlocalBase" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementDataBase);

    const std::uintptr_t* mNonlocalElementsAdress = reinterpret_cast<const std::uintptr_t*>(mNonlocalElements.data());
    int size = mNonlocalElements.size();
    ar & boost::serialization::make_nvp("mNonlocalElements_size", size);
    ar & boost::serialization::make_nvp("mNonlocalElements", boost::serialization::make_array(mNonlocalElementsAdress, size));
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ElementDataNonlocalBase" << std::endl;
#endif
}

void NuTo::ElementDataNonlocalBase::SetElementPtrAfterSerialization(const std::map<std::uintptr_t, std::uintptr_t>& mElementMapCast)
{
    for(std::vector<const ElementBase*>::const_iterator it = mNonlocalElements.begin(); it != mNonlocalElements.end(); it++)
    {
        std::uintptr_t temp = reinterpret_cast<std::uintptr_t>(*it);
        std::map<std::uintptr_t, std::uintptr_t>::const_iterator itCast = mElementMapCast.find(temp);
        if(itCast!=mElementMapCast.end())
        {
            ElementBase** tempPtr = const_cast<ElementBase**>(&(*it));
            *tempPtr = reinterpret_cast<ElementBase*>(itCast->second);
        }
        else
            throw MechanicsException("[NuTo::ElementDataNonlocalBase] The ElementBase-Pointer could not be updated.");
    }
}

BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ElementDataNonlocalBase)
#endif // ENABLE_SERIALIZATION
