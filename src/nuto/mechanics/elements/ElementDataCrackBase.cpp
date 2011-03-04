// $Id$

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/ElementDataCrackBase.h"
#include "nuto/mechanics/elements/ElementBase.h"

NuTo::ElementDataCrackBase::ElementDataCrackBase() :  NuTo::ElementDataBase::ElementDataBase()
{
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
