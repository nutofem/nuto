//
// Created by Thomas Titscher on 10/20/16.
//
#pragma once
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/base/ErrorEnum.h"

namespace NuTo
{
class SerializeStreamIn;
class SerializeStreamOut;

namespace Constitutive
{

//! @brief base class for a combined ConstitutiveLaw - ConstitutiveStaticData structure
class IPConstitutiveLawBase
{
public:

    //! @brief virtual destructor
    virtual ~IPConstitutiveLawBase() = default;

    virtual NuTo::eError Evaluate1D(
        const ConstitutiveInputMap& rConstitutiveInput,
        const ConstitutiveOutputMap& rConstitutiveOutput) = 0;

    //! @brief defines the serialization of this class
    //! @param rStream serialize output stream
    virtual void NuToSerializeSave(SerializeStreamOut& rStream) {/* no members to serialize */};

    //! @brief defines the serialization of this class
    //! @param rStream serialize input stream
    virtual void NuToSerializeLoad(SerializeStreamIn& rStream) {/* no members to serialize */};

protected:

};
} // namespace Constitutive
} // namespace NuTo
