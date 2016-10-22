//
// Created by Thomas Titscher on 10/20/16.
//
#pragma once

#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"

namespace NuTo
{
class ConstitutiveBase;
class SerializeStreamIn;
class SerializeStreamOut;
enum class eError;

namespace Constitutive
{

//! @brief base class for a combined ConstitutiveLaw - ConstitutiveStaticData structure
class IPConstitutiveLawBase
{
public:

    //! @brief virtual destructor
    //! @remark since we don't want to give up move semantics --> rule of 5
    virtual ~IPConstitutiveLawBase() = default;


    IPConstitutiveLawBase() = default;
    IPConstitutiveLawBase(const IPConstitutiveLawBase& ) = default;
    IPConstitutiveLawBase(      IPConstitutiveLawBase&&) = default;
    IPConstitutiveLawBase& operator=(const IPConstitutiveLawBase&)  = default;
    IPConstitutiveLawBase& operator=(      IPConstitutiveLawBase&&) = default;


    virtual ConstitutiveBase& GetConstitutiveLaw() const = 0;

    template <int TDim>
    NuTo::eError Evaluate(
        const ConstitutiveInputMap& rConstitutiveInput,
        const ConstitutiveOutputMap& rConstitutiveOutput)
    {
        static_assert(TDim == 1 or TDim == 2 or TDim == 3, "TDim == 1 or TDim == 2 or TDim == 3 !");
        if (TDim == 1) return Evaluate1D(rConstitutiveInput, rConstitutiveOutput);
        if (TDim == 2) return Evaluate2D(rConstitutiveInput, rConstitutiveOutput);
        if (TDim == 3) return Evaluate3D(rConstitutiveInput, rConstitutiveOutput);
    }

    //! @brief defines the serialization of this class
    //! @param rStream serialize output stream
    virtual void NuToSerializeSave(SerializeStreamOut& rStream) {/* no members to serialize */};

    //! @brief defines the serialization of this class
    //! @param rStream serialize input stream
    virtual void NuToSerializeLoad(SerializeStreamIn& rStream) {/* no members to serialize */};

protected:
    virtual NuTo::eError Evaluate1D(const ConstitutiveInputMap& rConstitutiveInput,
                                    const ConstitutiveOutputMap& rConstitutiveOutput) = 0;
    virtual NuTo::eError Evaluate2D(const ConstitutiveInputMap& rConstitutiveInput,
                                    const ConstitutiveOutputMap& rConstitutiveOutput) = 0;
    virtual NuTo::eError Evaluate3D(const ConstitutiveInputMap& rConstitutiveInput,
                                    const ConstitutiveOutputMap& rConstitutiveOutput) = 0;



};
} // namespace Constitutive
} // namespace NuTo
