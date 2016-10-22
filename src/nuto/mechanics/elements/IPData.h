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

//! @brief class for handling the IPConstitutiveLaw and the IntegrationType.
class IPData
{
public:

    IPData(const IntegrationTypeBase& rIntegrationType);

    IPData(const IPData&  rOther);
    IPData(      IPData&& rOther);

    IPData& operator = (const IPData&  rOther);
    IPData& operator = (      IPData&& rOther);

    //! @brief allocates a new IPConstitiveLaw according to rLaw for each integration point
    //! @param rLaw new constitutive law
    void SetConstitutiveLaw(NuTo::ConstitutiveBase& rLaw);

    //! @brief returns the IPConstitutiveLaw at the given (if allocated)
    //! @param rIP integration point index
    //! @return ip constitutive law base reference
    Constitutive::IPConstitutiveLawBase& GetIPConstitutiveLaw(unsigned int rIP);

    //! @brief returns the IPConstitutiveLaw at the given (if allocated)
    //! @param rIP integration point index
    //! @return ip constitutive law base reference
    const Constitutive::IPConstitutiveLawBase& GetIPConstitutiveLaw(unsigned int rIP) const;


    //! @brief sets a new integration type and performs an update of mLaws
    //! @param rIntegrationType new integration type
    void SetIntegrationType(const IntegrationTypeBase& rIntegrationType);

    //! @brief getter for the integration type
    //! @return integration type reference
    const IntegrationTypeBase& GetIntegrationType() const
    {
        return *mIntegrationType;
    }

private:

    //! @brief integration type pointer
    //! @remark this is only set via references. So no checks for nullptr, please. It is no reference because of slicing
    const IntegrationTypeBase* mIntegrationType;

    //! @brief owning container of IPConstitutiveLawBase - requires special copy/move semantics implemented here.
    std::vector<std::unique_ptr<Constitutive::IPConstitutiveLawBase>> mLaws;
};

} // namespace NuTo