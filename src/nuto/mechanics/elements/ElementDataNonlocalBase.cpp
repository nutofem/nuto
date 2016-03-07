// $Id$ 
// ElementDataNonlocalBase.cpp
// created Apr 22, 2010 by Joerg F. Unger

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
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ElementDataNonlocalBase)
#endif // ENABLE_SERIALIZATION
