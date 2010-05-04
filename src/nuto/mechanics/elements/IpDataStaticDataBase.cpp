// $ld: $ 
// IpDataStaticDataBase.cpp
// created Apr 29, 2010 by Joerg F. Unger
#include "nuto/mechanics/elements/IpDataStaticDataBase.h"

NuTo::IpDataStaticDataBase::IpDataStaticDataBase() : IpDataBase()
{
}

NuTo::IpDataStaticDataBase::~IpDataStaticDataBase()
{
}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void NuTo::IpDataStaticDataBase::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IpDataBase)
           & BOOST_SERIALIZATION_NVP(mStaticData);
    }
#endif  // ENABLE_SERIALIZATION
