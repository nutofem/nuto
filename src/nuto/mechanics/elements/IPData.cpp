//
// Created by Thomas Titscher on 10/22/16.
//

#include "nuto/mechanics/elements/IPData.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/base/serializeStream/SerializeStreamIn.h"
#include "nuto/base/serializeStream/SerializeStreamOut.h"

NuTo::IPData::IPData(const IntegrationTypeBase& rIntegrationType)
: mIntegrationType(&rIntegrationType)
{}

void NuTo::IPData::SetConstitutiveLaw(NuTo::ConstitutiveBase& rLaw)
{
    mLaws.clear();
    for (int i = 0; i < mIntegrationType->GetNumIntegrationPoints(); ++i)
        mLaws.push_back(rLaw.CreateIPLaw().release());
}
void NuTo::IPData::SetIntegrationType(const NuTo::IntegrationTypeBase& rIntegrationType)
{
    mIntegrationType = &rIntegrationType;

    if (HasConstitutiveLawAssigned(0))
        SetConstitutiveLaw(mLaws[0].GetConstitutiveLaw());
}

NuTo::Constitutive::IPConstitutiveLawBase& NuTo::IPData::GetIPConstitutiveLaw(unsigned int rIP)
{
    if (HasConstitutiveLawAssigned(rIP))
        return mLaws[rIP];
    throw MechanicsException(__PRETTY_FUNCTION__, "There is no constitutive law at IP " + std::to_string(rIP) + " assigned.");
}

const NuTo::Constitutive::IPConstitutiveLawBase& NuTo::IPData::GetIPConstitutiveLaw(unsigned int rIP) const
{
    if (HasConstitutiveLawAssigned(rIP))
        return mLaws[rIP];
    throw MechanicsException(__PRETTY_FUNCTION__, "There is no constitutive law at IP " + std::to_string(rIP) + " assigned.");
}
bool NuTo::IPData::HasConstitutiveLawAssigned(unsigned int rIP) const
{
    if (rIP >= mLaws.size())
        return false;

    return &mLaws[rIP] != nullptr;
}

void NuTo::IPData::NuToSerializeSave(NuTo::SerializeStreamOut& rStream)
{
    for (auto& ipLaw : mLaws)
        rStream.Serialize(ipLaw);
}

void NuTo::IPData::NuToSerializeLoad(NuTo::SerializeStreamIn& rStream)
{
    for (auto& ipLaw : mLaws)
        rStream.Serialize(ipLaw);
}
