#pragma once

#include "nuto/mechanics/elements/IpDataBase.h"
#include "nuto/mechanics/constitutive/staticData/Component.h"

namespace NuTo
{

class IpDataStaticDataBase : public virtual IpDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
	IpDataStaticDataBase();

    virtual ~IpDataStaticDataBase();

    Constitutive::StaticData::Component* GetConstitutiveStaticData() override
    {
        return mStaticData;
    }

    const Constitutive::StaticData::Component* GetConstitutiveStaticData() const override
    {
        return mStaticData;
    }

    void SetConstitutiveStaticData(Constitutive::StaticData::Component* newStaticData) override
    {
        mStaticData = newStaticData;
    }

    IpDataStaticDataBase& GetIpData() override
    {
        return *this;
    }

    const IpDataStaticDataBase& GetIpData() const override
    {
        return *this;
    }

    void Initialize(const ElementBase* rElement, const ConstitutiveBase* rConstitutive) override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

protected:
    Constitutive::StaticData::Component* mStaticData;
};
}
