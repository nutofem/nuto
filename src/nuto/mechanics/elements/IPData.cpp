//
// Created by Thomas Titscher on 10/22/16.
//

#include "nuto/mechanics/elements/IPData.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"

NuTo::IPData::IPData(const IntegrationTypeBase& rIntegrationType)
: mIntegrationType(&rIntegrationType)
{
    mLaws.resize(mIntegrationType->GetNumIntegrationPoints());
}


NuTo::IPData::IPData(const NuTo::IPData& rOther)
{
    *this = rOther;
}

NuTo::IPData::IPData(NuTo::IPData&& rOther)
{
    *this = std::move(rOther);
}

NuTo::IPData& NuTo::IPData::operator=(const NuTo::IPData& rOther)
{
    mIntegrationType = rOther.mIntegrationType;
    mLaws.clear();
    mLaws.reserve(rOther.mLaws.size());
    for (auto& law : rOther.mLaws)
    {
        mLaws.push_back(law->GetConstitutiveLaw().CreateIPLaw());
    }
    return *this;
}

NuTo::IPData& NuTo::IPData::operator=(NuTo::IPData&& rOther)
{
    mIntegrationType = rOther.mIntegrationType;
    mLaws.clear();
    mLaws.reserve(rOther.mLaws.size());
    for (auto& law : rOther.mLaws)
    {
        mLaws.push_back(std::move(law));
    }
    return *this;
}

void NuTo::IPData::SetConstitutiveLaw(NuTo::ConstitutiveBase& rLaw)
{
    for (auto& ipLaw : mLaws)
        ipLaw = rLaw.CreateIPLaw();
}
void NuTo::IPData::SetIntegrationType(const NuTo::IntegrationTypeBase& rIntegrationType)
{
    mIntegrationType = &rIntegrationType;
    if (static_cast<unsigned int>(mIntegrationType->GetNumIntegrationPoints()) == mLaws.size())
        return; // no change required.

    Constitutive::IPConstitutiveLawBase* law = mLaws[0].get();
    mLaws.resize(mIntegrationType->GetNumIntegrationPoints());
    if (law != nullptr)
        SetConstitutiveLaw(law->GetConstitutiveLaw());
}

NuTo::Constitutive::IPConstitutiveLawBase& NuTo::IPData::GetIPConstitutiveLaw(unsigned int rIP)
{
    if (rIP >= mLaws.size())
        throw MechanicsException(__PRETTY_FUNCTION__, "Out of bounds.");

    auto law = mLaws[rIP].get();
    if (law == nullptr)
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law set yet.");
    return *law;
}

const NuTo::Constitutive::IPConstitutiveLawBase& NuTo::IPData::GetIPConstitutiveLaw(unsigned int rIP) const
{
    if (rIP >= mLaws.size())
        throw MechanicsException(__PRETTY_FUNCTION__, "Out of bounds.");

    auto law = mLaws[rIP].get();
    if (law == nullptr)
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law set yet.");
    return *law;
}
bool NuTo::IPData::HasConstitutiveLawAssigned(unsigned int rIP) const
{
    if (rIP >= mLaws.size())
        return false;

    auto law = mLaws[rIP].get();
    return law != nullptr;
}