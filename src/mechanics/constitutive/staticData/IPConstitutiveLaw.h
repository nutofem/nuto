//
// Created by Thomas Titscher on 10/20/16.
//
#pragma once

#include "mechanics/constitutive/staticData/IPConstitutiveLawBase.h"
#include "mechanics/constitutive/staticData/DataContainer.h"
#include "base/serializeStream/SerializeStreamIn.h"
#include "base/serializeStream/SerializeStreamOut.h"
#include <type_traits>

namespace NuTo
{
namespace Constitutive
{

template<typename TLaw>
class IPConstitutiveLaw: public IPConstitutiveLawBase
{
public:

    static_assert(std::is_base_of<ConstitutiveBase, TLaw>::value,"TLaw must be derived from NuTo::ConstitutiveBase");

    using Type = typename TLaw::StaticDataType;
    typedef typename StaticData::DataContainer<Type> Data;

    //! @brief constructor
    //! @param rLaw underlying constitutive law
    IPConstitutiveLaw(TLaw& rLaw, const Type& rData)
    : mLaw(rLaw), mData(rData) {}

    std::unique_ptr<IPConstitutiveLawBase> Clone() const override
    {
        return std::make_unique<IPConstitutiveLaw<TLaw>>(*this);
    }


    ConstitutiveBase& GetConstitutiveLaw() const override
    {
        return mLaw;
    }

    const StaticData::DataContainer<Type>& GetStaticData() const
    {
        return mData;
    }

    StaticData::DataContainer<Type>& GetStaticData()
    {
        return mData;
    }

    //! @brief allocates rNum additional static data
    //! @param rNum number of addtional static data
    void AllocateAdditional(int rNum) override
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
    StaticData::DataContainer<Type> mData;
};


} // namespace Constitutive
} // namespace NuTo
