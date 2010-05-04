// $ld: $ 
// IpDataStaticDataNonlocal.cpp
// created Apr 29, 2010 by Joerg F. Unger

#include "nuto/mechanics/elements/IpDataStaticDataNonlocal.h"
#include "nuto/mechanics/elements/ElementWithDataBase.h"

NuTo::IpDataStaticDataNonlocal::~IpDataStaticDataNonlocal()
{
}

void NuTo::IpDataStaticDataNonlocal::Initialize(const ElementWithDataBase* rElement, int rIp)
{
	mStaticData = rElement->AllocateStaticData(rElement->GetConstitutiveLaw(rIp));
}

