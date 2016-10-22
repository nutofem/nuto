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

}

NuTo::IPData::IPData(NuTo::IPData&& rOther)
{

}

NuTo::IPData& NuTo::IPData::operator=(const NuTo::IPData& rOther)
{
//    return <#initializer#>;
}

NuTo::IPData& NuTo::IPData::operator=(NuTo::IPData&& rOther)
{
//    return <#initializer#>;
}

void NuTo::IPData::SetConstitutiveLaw(NuTo::ConstitutiveBase& rLaw)
{
    for (auto& ipLaw : mLaws)
        ipLaw = rLaw.CreateIPLaw();
}
void NuTo::IPData::SetIntegrationType(const NuTo::IntegrationTypeBase& rIntegrationType)
{
    mIntegrationType = &rIntegrationType;
    // This could be done more efficiently, sure. But let's keep it simple.
    Constitutive::IPConstitutiveLawBase* law = mLaws[0].get();
    mLaws.resize(mIntegrationType->GetNumIntegrationPoints());
    if (law != nullptr)
        SetConstitutiveLaw(law->GetConstitutiveLaw());
}

