//
// Created by Thomas Titscher on 10/20/16.
//
#pragma once

#include "nuto/mechanics/constitutive/staticData/IPConstitutiveLawBase.h"
#include "nuto/mechanics/constitutive/staticData/DataContainer.h"
#include "nuto/base/serializeStream/SerializeStreamIn.h"
#include "nuto/base/serializeStream/SerializeStreamOut.h"
#include <type_traits>

namespace NuTo
{
namespace Constitutive
{

template<typename TLaw>
class IPConstitutiveLaw: public IPConstitutiveLawBase
{
public:

//    static_assert(std::is_base_of<ConstitutiveBase, TLaw>::value,"TLaw must be derived from NuTo::ConstitutiveBase");

    using Type = typename StaticData::DataContainer<TLaw>::Type;
    typedef typename StaticData::DataContainer<TLaw> Data;

    //! @brief constructor
    //! @param rLaw underlying constitutive law
    IPConstitutiveLaw(TLaw& rLaw, const Type& rData)
    : mLaw(rLaw), mData(rData) {}

    virtual std::unique_ptr<IPConstitutiveLawBase> Clone()
    {
        return std::make_unique<IPConstitutiveLaw<TLaw>>(*this);
    }


    ConstitutiveBase& GetConstitutiveLaw() const
    {
        return static_cast<ConstitutiveBase&>(mLaw);
    }

    const StaticData::DataContainer<TLaw>& GetStaticData() const
    {
        return mData;
    }

    StaticData::DataContainer<TLaw>& GetStaticData()
    {
        return mData;
    }

    //! @brief allocates rNum additional static data
    //! @param rNum number of addtional static data
    void AllocateAdditional(unsigned int rNum) override
    {
        mData.AllocateAdditionalData(rNum);
    }

    //! @brief Puts current static data to previous static data, previous to pre-previous, etc.
    void ShiftToPast() override
    {
        mData.ShiftToPast();
    }

    //! @brief Puts previous static data to current static data, pre-previous to previous, etc.
    void ShiftToFuture() override
    {
        mData.ShiftToPast();
    }


    //! @brief defines the serialization of this class
    //! @param rStream serialize output stream
    virtual void NuToSerializeSave(SerializeStreamOut& rStream) override
    {
        IPConstitutiveLawBase::NuToSerializeSave(rStream);
        rStream.Serialize(mData);
    }

    //! @brief defines the serialization of this class
    //! @param rStream serialize input stream
    virtual void NuToSerializeLoad(SerializeStreamIn& rStream) override
    {
        IPConstitutiveLawBase::NuToSerializeLoad(rStream);
        rStream.Serialize(mData);
    }

protected:
    //! @brief Evaluate the constitutive relation in 1D
    //! @param rConstitutiveInput Input to the constitutive law
    //! @param rConstitutiveOutput Output of the constitutive law
    eError Evaluate1D(const ConstitutiveInputMap& rConstitutiveInput,
                      const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        return mLaw.template Evaluate<1>(rConstitutiveInput, rConstitutiveOutput, mData);
    }

    //! @brief Evaluate the constitutive relation in 1D
    //! @param rConstitutiveInput Input to the constitutive law
    //! @param rConstitutiveOutput Output of the constitutive law
    eError Evaluate2D(const ConstitutiveInputMap& rConstitutiveInput,
                      const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        return mLaw.template Evaluate<2>(rConstitutiveInput, rConstitutiveOutput, mData);
    }

    //! @brief Evaluate the constitutive relation in 1D
    //! @param rConstitutiveInput Input to the constitutive law
    //! @param rConstitutiveOutput Output of the constitutive law
    eError Evaluate3D(const ConstitutiveInputMap& rConstitutiveInput,
                      const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        return mLaw.template Evaluate<3>(rConstitutiveInput, rConstitutiveOutput, mData);
    }

private:

    TLaw& mLaw;
    StaticData::DataContainer<TLaw> mData;
};


} // namespace Constitutive
} // namespace NuTo
