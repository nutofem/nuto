// $ld: $ 
// IpDataStaticData.cpp
// created Apr 29, 2010 by Joerg F. Unger

#include "nuto/mechanics/elements/IpDataStaticData.h"
#include "nuto/mechanics/elements/ElementBase.h"

NuTo::IpDataStaticData::IpDataStaticData() : IpDataBase()
{
	mStaticData = 0;
}

NuTo::IpDataStaticData::~IpDataStaticData()
{
}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void NuTo::IpDataStaticData::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IpDataStaticDataBase);
    }
#endif  // ENABLE_SERIALIZATION

void NuTo::IpDataStaticData::Initialize(const ElementBase* rElement, const ConstitutiveBase* rConstitutive)
{
	if (mStaticData!=0)
		delete mStaticData;
	if (rConstitutive!=0)
	    mStaticData = rElement->AllocateStaticData(rConstitutive);
	else
		mStaticData = 0;
}
