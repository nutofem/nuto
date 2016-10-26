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
template <typename TLaw> class IPConstitutiveLaw;
namespace StaticData
{
template <typename TLaw> class DataContainer;
}

template<class TLaw>
class HasStaticData
{
private:
    using Yes = char[2];
    using  No = char[1];

    // declare struct with a type StaticDataType
    struct Fallback { typedef void StaticDataType; };
    // if TLaw also has StaticDataType, than Derived has it twice
    struct Derived : TLaw, Fallback { };

    // this can only be instantiated, if U::StaticDataType can be used unambiguosly,
    // i.e. only one StaticDataType in U (will be instantiated with Derived later)
    template<class U>
    static No& test(typename U::StaticDataType*);
    // in case the substitution fails, this test will be instantiated
    template<class U>
    static Yes& test(U*);
public:
    // if test has returntype Yes, we do have staticData
    static constexpr bool value = sizeof(test<Derived>(nullptr)) == sizeof(Yes);
};

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

    virtual std::unique_ptr<IPConstitutiveLawBase> Clone() = 0;

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

    //! @brief allocates rNum additional static data
    //! @param rNum number of addtional static data
    virtual void AllocateAdditional(unsigned int rNum) = 0;

    //! @brief Puts current static data to previous static data, previous to pre-previous, etc.
    virtual void ShiftToPast() = 0;

    //! @brief Puts previous static data to current static data, pre-previous to previous, etc.
    virtual void ShiftToFuture() = 0;

    //! @brief casts *this to a IPConstitutiveLaw and returns its data
    //! @return static data container of TLaw
    template<class TLaw, typename = std::enable_if_t<HasStaticData<TLaw>::value>>
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
    virtual NuTo::eError Evaluate1D(const ConstitutiveInputMap& rConstitutiveInput,
                                    const ConstitutiveOutputMap& rConstitutiveOutput) = 0;
    virtual NuTo::eError Evaluate2D(const ConstitutiveInputMap& rConstitutiveInput,
                                    const ConstitutiveOutputMap& rConstitutiveOutput) = 0;
    virtual NuTo::eError Evaluate3D(const ConstitutiveInputMap& rConstitutiveInput,
                                    const ConstitutiveOutputMap& rConstitutiveOutput) = 0;
};
} // namespace Constitutive
} // namespace NuTo
