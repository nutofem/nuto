#pragma once

#include <boost/ptr_container/ptr_vector.hpp>
#include "mechanics/constitutive/staticData/IPConstitutiveLawBase.h"
#include "base/serializeStream/SerializeStreamIn.h"
#include "base/serializeStream/SerializeStreamOut.h"

namespace NuTo
{
class AdditiveInputExplicit;
namespace Constitutive
{

class IPAdditiveInputExplicit : public IPConstitutiveLawBase
{
public:
    //! @brief constructor
    //! @param rLaw underlying constitutive law
    IPAdditiveInputExplicit(AdditiveInputExplicit& rLaw);

    IPAdditiveInputExplicit(const IPAdditiveInputExplicit& rOther);

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

protected:
    //! @brief Searches for a specific IP constitutive law and returns it
    //! @param rCLPtr The constitutive law of the IP constitutive law that is requested
    //! @return Searched IP constitutive law - nullptr if law is not found
    virtual IPConstitutiveLawBase* GetSublawIP(ConstitutiveBase* rCLPtr) override;

    template <int TDim>
    void AdditiveInputExplicitEvaluate(const ConstitutiveInputMap& rConstitutiveInput,
                                       const ConstitutiveOutputMap& rConstitutiveOutput);

    template <int TDim>
    void CalculateDerivatives(const ConstitutiveOutputMap& rConstitutiveOutput,
                              std::vector<ConstitutiveOutputMap>& rSublawOutputVec);

    template <int TDim>
    void ApplySublawOutputs(const ConstitutiveInputMap& rMainLawInput, const ConstitutiveOutputMap& rConstitutiveOutput,
                            const ConstitutiveOutputMap& rSublawOutput);

    //! @brief Evaluate the constitutive relation in 1D
    //! @param rConstitutiveInput Input to the constitutive law
    //! @param rConstitutiveOutput Output of the constitutive law
    void Evaluate1D(const ConstitutiveInputMap& rConstitutiveInput,
                    const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        AdditiveInputExplicitEvaluate<1>(rConstitutiveInput, rConstitutiveOutput);
    }

    //! @brief Evaluate the constitutive relation in 1D
    //! @param rConstitutiveInput Input to the constitutive law
    //! @param rConstitutiveOutput Output of the constitutive law
    void Evaluate2D(const ConstitutiveInputMap& rConstitutiveInput,
                    const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        AdditiveInputExplicitEvaluate<2>(rConstitutiveInput, rConstitutiveOutput);
    }

    //! @brief Evaluate the constitutive relation in 1D
    //! @param rConstitutiveInput Input to the constitutive law
    //! @param rConstitutiveOutput Output of the constitutive law
    void Evaluate3D(const ConstitutiveInputMap& rConstitutiveInput,
                    const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        AdditiveInputExplicitEvaluate<3>(rConstitutiveInput, rConstitutiveOutput);
    }

private:
    AdditiveInputExplicit& mLaw;
    std::unique_ptr<IPConstitutiveLawBase> mMainLawIP;
    boost::ptr_vector<IPConstitutiveLawBase> mSublawIPs;
};


} // namespace Constitutive
} // namespace NuTo
