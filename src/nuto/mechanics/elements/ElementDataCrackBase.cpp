// $Id$

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/ElementDataCrackBase.h"
#include "nuto/mechanics/elements/ElementBase.h"

NuTo::ElementDataCrackBase::ElementDataCrackBase() :  NuTo::ElementDataBase::ElementDataBase()
{
	isCracked=false;
}

NuTo::ElementDataCrackBase::~ElementDataCrackBase()
{
}

std::vector<NuTo::CrackBase*>& NuTo::ElementDataCrackBase::GetCracks()
{
    return mCracks;
}

int NuTo::ElementDataCrackBase::GetNumCracks()const
{
    return mCracks.size();
}

//! @brief Set the information that the element is already cracked or not
//! @param bool (Input) cracked or not
void NuTo::ElementDataCrackBase::IsCracked(const bool rIsCracked)
{
	isCracked=rIsCracked;
}

//! @brief Give the information if the element is already cracked or not
//! @return bool cracked or not
const bool NuTo::ElementDataCrackBase::IsCracked() const
{
    return isCracked;
}

//! @brief adds a crack to the element
//! @param rCrack  crack
//! @return the local crack number, the crack is either append to the list, or the existing local number is returned
unsigned int NuTo::ElementDataCrackBase::AddCrack(NuTo::CrackBase* rCrack)
{
	for (unsigned int theCrack=0; theCrack<mCracks.size();theCrack++)
	{
		if (mCracks[theCrack]==rCrack)
			return theCrack;
	}
	mCracks.push_back(rCrack);
	return mCracks.size()-1;
}
