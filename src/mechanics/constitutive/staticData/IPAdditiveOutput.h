#pragma once

#include <boost/ptr_container/ptr_vector.hpp>
#include "mechanics/constitutive/staticData/IPConstitutiveLawBase.h"
#include "mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "base/serializeStream/SerializeStreamIn.h"
#include "base/serializeStream/SerializeStreamOut.h"
#include <type_traits>

namespace NuTo
{
class AdditiveOutput;
namespace Constitutive
{

class IPAdditiveOutput: public IPConstitutiveLawBase
{
public:

    //! @brief constructor
    //! @param rLaw underlying constitutive law
    IPAdditiveOutput(AdditiveOutput& rLaw);

    IPAdditiveOutput(const IPAdditiveOutput& rOther);

    std::unique_ptr<IPConstitutiveLawBase> Clone() const override;

    ConstitutiveBase& GetConstitutiveLaw() const override;

    //! @brief allocates rNum additional static data
    //! @param rNum number of addtional static data
    void AllocateAdditional(int rNum) override;

    //! @brief Puts current static data to previous static data, previous to pre-previous, etc.
    void ShiftToPast() override;

    //! @brief Puts previous static data to current static data, pre-previous to previous, etc.
    void ShiftToFuture() override;

    //! @brief defines the serialization of this class
    //! @param rStream serialize output stream
    virtual void NuToSerializeSave(SerializeStreamOut& rStream) override;

    //! @brief defines the serialization of this class
    //! @param rStream serialize input stream
    virtual void NuToSerializeLoad(SerializeStreamIn& rStream) override;

    //! @brief casts *this to a IPConstitutiveLaw and returns its data
    //! @return static data container of TLaw
    template <typename TLaw>
    StaticData::DataContainer<typename TLaw::StaticDataType>& GetSublawData(ConstitutiveBase* rCLPtr)
    {
        try
        {
            IPConstitutiveLawBase* ConstitutiveLaw = GetSublawIP(rCLPtr);
            if(ConstitutiveLaw == nullptr)
                throw Exception(__PRETTY_FUNCTION__,"The requested static data/constitutive law is not attached to this law!");
            return ConstitutiveLaw->GetData<TLaw>();
        }
        catch (std::bad_cast& e)
        {
            throw Exception(__PRETTY_FUNCTION__, "Wrong ConstitutiveLawType requested.");
        }
    }

protected:
    template<int TDim>
    void AdditiveOutputEvaluate(const ConstitutiveInputMap& rConstitutiveInput,
            const ConstitutiveOutputMap& rConstitutiveOutput);



    //! @brief Evaluate the constitutive relation in 1D
    //! @param rConstitutiveInput Input to the constitutive law
    //! @param rConstitutiveOutput Output of the constitutive law
    void Evaluate1D(const ConstitutiveInputMap& rConstitutiveInput,
                      const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        AdditiveOutputEvaluate<1>(rConstitutiveInput, rConstitutiveOutput);
    }

    //! @brief Evaluate the constitutive relation in 1D
    //! @param rConstitutiveInput Input to the constitutive law
    //! @param rConstitutiveOutput Output of the constitutive law
    void Evaluate2D(const ConstitutiveInputMap& rConstitutiveInput,
                      const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        AdditiveOutputEvaluate<2>(rConstitutiveInput, rConstitutiveOutput);
    }

    //! @brief Evaluate the constitutive relation in 1D
    //! @param rConstitutiveInput Input to the constitutive law
    //! @param rConstitutiveOutput Output of the constitutive law
    void Evaluate3D(const ConstitutiveInputMap& rConstitutiveInput,
                      const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        AdditiveOutputEvaluate<3>(rConstitutiveInput, rConstitutiveOutput);
    }

protected:

    //! @brief Searches for a specific IP constitutive law and returns it
    //! @param rCLPtr The constitutive law of the IP constitutive law that is requested
    //! @return Searched IP constitutive law - nullptr if law is not found
    virtual IPConstitutiveLawBase* GetSublawIP(ConstitutiveBase* rCLPtr) override;

private:


    AdditiveOutput& mLaw;
    boost::ptr_vector<IPConstitutiveLawBase> mSublawIPs;
};


} // namespace Constitutive
} // namespace NuTo
