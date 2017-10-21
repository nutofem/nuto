//
// Created by Thomas Titscher on 10/22/16.
//

#include "mechanics/elements/IPData.h"
#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "base/serializeStream/SerializeStreamIn.h"
#include "base/serializeStream/SerializeStreamOut.h"

NuTo::IPData::IPData(const IntegrationTypeBase& rIntegrationType)
    : mIntegrationType(&rIntegrationType)
{
}

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


NuTo::Constitutive::IPConstitutiveLawBase& NuTo::IPData::GetIPConstitutiveLaw(int rIP)
{
    if (HasConstitutiveLawAssigned(rIP))
        return mLaws[rIP];
    throw Exception(__PRETTY_FUNCTION__, "There is no constitutive law at IP " + std::to_string(rIP) + " assigned.");
}


const NuTo::Constitutive::IPConstitutiveLawBase& NuTo::IPData::GetIPConstitutiveLaw(int rIP) const
{
    if (HasConstitutiveLawAssigned(rIP))
        return mLaws[rIP];
    throw Exception(__PRETTY_FUNCTION__, "There is no constitutive law at IP " + std::to_string(rIP) + " assigned.");
}


bool NuTo::IPData::HasConstitutiveLawAssigned(int rIP) const
{
    if (rIP >= mIntegrationType->GetNumIntegrationPoints())
        throw Exception(__PRETTY_FUNCTION__, "The current integration type does not have that many integrationpoints.");
    return !mLaws.empty();
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
