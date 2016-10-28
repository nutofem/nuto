//
// Created by Thomas Titscher on 10/22/16.
//
#pragma once
#include <boost/ptr_container/ptr_vector.hpp>
#include "nuto/mechanics/constitutive/staticData/IPConstitutiveLawBase.h"

namespace NuTo
{
class IntegrationTypeBase;
class ConstitutiveBase;
class SerializeStreamOut;
class SerializeStreamIn;

//! @brief class for handling the IPConstitutiveLaw and the IntegrationType.
class IPData
{
public:

    IPData(const IntegrationTypeBase& rIntegrationType);

    //! @brief allocates a new IPConstitiveLaw according to rLaw for each integration point
    //! @param rLaw new constitutive law
    void SetConstitutiveLaw(NuTo::ConstitutiveBase& rLaw);

    //! @brief returns true, if the constitutive law has been assigned
    //! @param rIP integration point index
    //! @return false if a) rIP is out of bounds, b) constitutive law at rIP is not set yed
    bool HasConstitutiveLawAssigned(int rIP) const;

    //! @brief returns the IPConstitutiveLaw at the given (if allocated)
    //! @param rIP integration point index
    //! @return ip constitutive law base reference
    Constitutive::IPConstitutiveLawBase& GetIPConstitutiveLaw(int rIP);

    //! @brief returns the IPConstitutiveLaw at the given (if allocated)
    //! @param rIP integration point index
    //! @return ip constitutive law base reference
    const Constitutive::IPConstitutiveLawBase& GetIPConstitutiveLaw(int rIP) const;

    //! @brief sets a new integration type and performs an update of mLaws
    //! @param rIntegrationType new integration type
    void SetIntegrationType(const IntegrationTypeBase& rIntegrationType);

    //! @brief getter for the integration type
    //! @return integration type reference
    const IntegrationTypeBase& GetIntegrationType() const
    {
        return *mIntegrationType;
    }

    //! @brief defines the serialization of this class
    //! @param rStream serialize output stream
    void NuToSerializeSave(SerializeStreamOut& rStream);

    //! @brief defines the serialization of this class
    //! @param rStream serialize input stream
    void NuToSerializeLoad(SerializeStreamIn& rStream);

private:
#ifdef ENABLE_SERIALIZATION
    IPData() {};
#endif // ENABLE_SERIALIZATION

    //! @brief integration type pointer
    //! @remark this is only set via references. So no checks for nullptr, please. It is no reference because of slicing
    const IntegrationTypeBase* mIntegrationType;

    //! @brief owning container of IPConstitutiveLawBase. requires new_clone(IPContitutiveLawBase*)
    boost::ptr_vector<Constitutive::IPConstitutiveLawBase> mLaws;
};

} // namespace NuTo
