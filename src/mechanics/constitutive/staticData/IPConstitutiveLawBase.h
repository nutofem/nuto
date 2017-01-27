//
// Created by Thomas Titscher on 10/20/16.
//
#pragma once

#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"

namespace NuTo
{
class ConstitutiveBase;
class SerializeStreamIn;
class SerializeStreamOut;
enum class eError;

namespace Constitutive
{
class IPAdditiveOutput;
template <typename TLaw> class IPConstitutiveLaw;
namespace StaticData
{
template <typename TLaw> class DataContainer;
}

//! @brief base class for a combined ConstitutiveLaw - ConstitutiveStaticData structure
class IPConstitutiveLawBase
{
    friend class Constitutive::IPAdditiveOutput;


public:

    //! @brief virtual destructor
    //! @remark since we don't want to give up move semantics --> rule of 5
    virtual ~IPConstitutiveLawBase() = default;

    IPConstitutiveLawBase() = default;
    IPConstitutiveLawBase(const IPConstitutiveLawBase& ) = default;
    IPConstitutiveLawBase(      IPConstitutiveLawBase&&) = default;
    IPConstitutiveLawBase& operator=(const IPConstitutiveLawBase&)  = default;
    IPConstitutiveLawBase& operator=(      IPConstitutiveLawBase&&) = default;

    virtual std::unique_ptr<IPConstitutiveLawBase> Clone() const = 0;

    virtual ConstitutiveBase& GetConstitutiveLaw() const = 0;

    template <int TDim>
    void Evaluate(
        const ConstitutiveInputMap& rConstitutiveInput,
        const ConstitutiveOutputMap& rConstitutiveOutput)
    {
        static_assert(TDim == 1 or TDim == 2 or TDim == 3, "TDim == 1 or TDim == 2 or TDim == 3 !");
        if (TDim == 1) Evaluate1D(rConstitutiveInput, rConstitutiveOutput);
        if (TDim == 2) Evaluate2D(rConstitutiveInput, rConstitutiveOutput);
        if (TDim == 3) Evaluate3D(rConstitutiveInput, rConstitutiveOutput);
    }

    //! @brief allocates rNum additional static data
    //! @param rNum number of addtional static data
    virtual void AllocateAdditional(int rNum) = 0;

    //! @brief Puts current static data to previous static data, previous to pre-previous, etc.
    virtual void ShiftToPast() = 0;

    //! @brief Puts previous static data to current static data, pre-previous to previous, etc.
    virtual void ShiftToFuture() = 0;

    //! @brief casts *this to a IPConstitutiveLaw and returns its data
    //! @return static data container of TLaw
    template <typename TLaw>
    StaticData::DataContainer<typename TLaw::StaticDataType>& GetData()
    {
        try
        {
            return dynamic_cast<IPConstitutiveLaw<TLaw>&>(*this).GetStaticData();
        }
        catch (std::bad_cast& e)
        {
            throw MechanicsException(__PRETTY_FUNCTION__, "Wrong ConstitutiveLawType requested.");
        }
    }

    //! @brief defines the serialization of this class
    //! @param rStream serialize output stream
    virtual void NuToSerializeSave(SerializeStreamOut&) {/* no members to serialize */};

    //! @brief defines the serialization of this class
    //! @param rStream serialize input stream
    virtual void NuToSerializeLoad(SerializeStreamIn&) {/* no members to serialize */};

protected:
    virtual void Evaluate1D(const ConstitutiveInputMap& rConstitutiveInput,
                                    const ConstitutiveOutputMap& rConstitutiveOutput) = 0;
    virtual void Evaluate2D(const ConstitutiveInputMap& rConstitutiveInput,
                                    const ConstitutiveOutputMap& rConstitutiveOutput) = 0;
    virtual void Evaluate3D(const ConstitutiveInputMap& rConstitutiveInput,
                                    const ConstitutiveOutputMap& rConstitutiveOutput) = 0;

    //! @brief Searches for a specific IP constitutive law and returns it (Additive laws only)
    //! @param rCLPtr The constitutive law of the IP constitutive law that is requested
    //! @return Searched IP constitutive law - nullptr if law is not found
    virtual IPConstitutiveLawBase* GetSublawIP(ConstitutiveBase* rCLPtr)
    {
        throw MechanicsException(__PRETTY_FUNCTION__,"This function is only used by additive laws!");
    }



};

//! @brief clone methods that enables a boost::ptr_container<this> to copy itself
//! @param rLaw reference to the IPConstitutiveLawBase
//! @return cloned owning raw pointer of rLaw
inline NuTo::Constitutive::IPConstitutiveLawBase* new_clone( const NuTo::Constitutive::IPConstitutiveLawBase& rLaw)
{
    return rLaw.Clone().release();
}


} // namespace Constitutive
} // namespace NuTo

