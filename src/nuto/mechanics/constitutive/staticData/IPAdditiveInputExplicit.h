#pragma once

#include <boost/ptr_container/ptr_vector.hpp>
#include "nuto/mechanics/constitutive/staticData/IPConstitutiveLawBase.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "nuto/base/serializeStream/SerializeStreamIn.h"
#include "nuto/base/serializeStream/SerializeStreamOut.h"
#include <type_traits>

namespace NuTo
{
class AdditiveInputExplicit;
namespace Constitutive
{

class IPAdditiveInputExplicit: public IPConstitutiveLawBase
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
    template<int TDim>
    eError AdditiveInputExplicitEvaluate(const ConstitutiveInputMap& rConstitutiveInput,
            const ConstitutiveOutputMap& rConstitutiveOutput);

    template <int TDim>
    void CalculateDerivatives(const ConstitutiveOutputMap& rConstitutiveOutput,
            std::vector<ConstitutiveOutputMap>& rSublawOutputVec);

    template <int TDim>
    void ApplySublawOutputs(const ConstitutiveInputMap& rMainLawInput,
            const ConstitutiveOutputMap& rConstitutiveOutput, const ConstitutiveOutputMap& rSublawOutput);

    //! @brief Evaluate the constitutive relation in 1D
    //! @param rConstitutiveInput Input to the constitutive law
    //! @param rConstitutiveOutput Output of the constitutive law
    eError Evaluate1D(const ConstitutiveInputMap& rConstitutiveInput,
                      const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        return AdditiveInputExplicitEvaluate<1>(rConstitutiveInput, rConstitutiveOutput);
    }

    //! @brief Evaluate the constitutive relation in 1D
    //! @param rConstitutiveInput Input to the constitutive law
    //! @param rConstitutiveOutput Output of the constitutive law
    eError Evaluate2D(const ConstitutiveInputMap& rConstitutiveInput,
                      const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        return AdditiveInputExplicitEvaluate<2>(rConstitutiveInput, rConstitutiveOutput);
    }

    //! @brief Evaluate the constitutive relation in 1D
    //! @param rConstitutiveInput Input to the constitutive law
    //! @param rConstitutiveOutput Output of the constitutive law
    eError Evaluate3D(const ConstitutiveInputMap& rConstitutiveInput,
                      const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        return AdditiveInputExplicitEvaluate<3>(rConstitutiveInput, rConstitutiveOutput);
    }

private:

    AdditiveInputExplicit& mLaw;
    std::unique_ptr<IPConstitutiveLawBase> mMainLawIP;
    boost::ptr_vector<IPConstitutiveLawBase> mSublawIPs;
};


} // namespace Constitutive
} // namespace NuTo
