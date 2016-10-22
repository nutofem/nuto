//
// Created by Thomas Titscher on 10/22/16.
//
#pragma once
#include <vector>
#include "nuto/mechanics/constitutive/staticData/IPConstitutiveLawBase.h"

namespace NuTo
{
class IntegrationTypeBase;
class ConstitutiveBase;
class IPData
{
public:
    IPData(const IntegrationTypeBase& rIntegrationType);

    IPData(const IPData&  rOther);
    IPData(      IPData&& rOther);

    IPData& operator = (const IPData&  rOther);
    IPData& operator = (      IPData&& rOther);


    void SetConstitutiveLaw(NuTo::ConstitutiveBase& rLaw);

    void SetIntegrationType(const IntegrationTypeBase& rIntegrationType);

    const IntegrationTypeBase& GetIntegrationType() const
    {
        return *mIntegrationType;
    }

private:
    const IntegrationTypeBase* mIntegrationType;
    std::vector<std::unique_ptr<Constitutive::IPConstitutiveLawBase>> mLaws;
};

} // namespace NuTo