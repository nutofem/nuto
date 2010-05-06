// $ld: $ 
// ElementDataNonlocalBase.cpp
// created Apr 22, 2010 by Joerg F. Unger

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/ElementDataNonlocalBase.h"
#include "nuto/mechanics/elements/ElementWithDataBase.h"
#include <assert.h>

NuTo::ElementDataNonlocalBase::ElementDataNonlocalBase() :  NuTo::ElementDataBase::ElementDataBase()
{
	mConstitutive = 0;
}

NuTo::ElementDataNonlocalBase::~ElementDataNonlocalBase()
{
	//std::cout << "NuTo::ElementDataNonlocalBase::~ElementDataNonlocalBase()" << std::endl;
}

const std::vector<const NuTo::ElementWithDataBase*>&
  NuTo::ElementDataNonlocalBase::GetNonlocalElements(const ConstitutiveBase* rConstitutive)const
{
    if (rConstitutive==mConstitutive)
	    return mNonlocalElements;
    else
    	throw MechanicsException("[NuTo::ElementDataNonlocalBase::GetNonlocalElements] For this constitutive model no nonlocal data is available");
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
int NuTo::ElementDataNonlocalBase::AddNonlocalElement(const ElementWithDataBase* rElement, const ConstitutiveBase* rConstitutive)
{
	printf("\n[NuTo::ElementDataNonlocalBase::AddNonlocalElement %p %p\n]",rConstitutive,mConstitutive);
	if (mConstitutive==0)
		mConstitutive=rConstitutive;
	assert(rConstitutive==mConstitutive);
	for (int theElement=0; theElement<(int)mNonlocalElements.size();theElement++)
	{
		if (mNonlocalElements[theElement]==rElement)
			return theElement;
	}
	mNonlocalElements.push_back(rElement);
	return mNonlocalElements.size()-1;
}
