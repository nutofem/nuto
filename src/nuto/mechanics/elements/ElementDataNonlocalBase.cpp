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

const std::vector<const NuTo::ElementBase*>&
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

//! @brief adds the nonlocal weight to an integration point
//! @param rNonlocalElement element of the nonlocal ip
//! @param rWeight weights for each integration point of the nonlocal element
/*void NuTo::ElementDataNonlocalBase::AddNonlocalElement(int rNonlocalElement)
{
	throw MechanicsException("[NuTo::ElementDataNonlocalBase::AddNonlocalElement] To be implemented.");
}
*/
