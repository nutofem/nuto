// $ld: $ 
// IpDataStaticDataNonlocal.cpp
// created Apr 29, 2010 by Joerg F. Unger

#include "nuto/mechanics/elements/IpDataStaticDataNonlocal.h"
#include "nuto/mechanics/elements/ElementWithDataBase.h"
NuTo::IpDataStaticDataNonlocal::IpDataStaticDataNonlocal() :NuTo::IpDataBase::IpDataBase() ,
    NuTo::IpDataStaticDataBase::IpDataStaticDataBase() , NuTo::IpDataNonlocalBase::IpDataNonlocalBase()
{

}

NuTo::IpDataStaticDataNonlocal::~IpDataStaticDataNonlocal()
{
}

void NuTo::IpDataStaticDataNonlocal::Initialize(const ElementWithDataBase* rElement, const ConstitutiveBase* rConstitutive)
{
	if (mStaticData!=0)
		delete mStaticData;
	if (rConstitutive!=0)
	    mStaticData = rElement->AllocateStaticData(rConstitutive);
	else
		mStaticData = 0;
}

