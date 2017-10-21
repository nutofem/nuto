#pragma once

#include "mechanics/constitutive/laws/AdditiveBase.h"
#include "mechanics/constitutive/staticData/IPAdditiveOutput.h"

namespace NuTo
{

class AdditiveOutput : public AdditiveBase
{
public:
    //! @brief constructor
    AdditiveOutput(const int& rNumTimeDerivatives)
        : AdditiveBase(rNumTimeDerivatives)
    {
    }

    // has no ip static data itself
    std::unique_ptr<Constitutive::IPConstitutiveLawBase> CreateIPLaw() override
    {
        return std::make_unique<Constitutive::IPAdditiveOutput>(*this);
        mStaticDataAllocated = true;
    }


    template <int TDim>
    NuTo::eError Evaluate(const ConstitutiveInputMap& rConstitutiveInput,
                          const ConstitutiveOutputMap& rConstitutiveOutput)
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__,
                              "Additive Law cannot be evaluated. Its IPAdditiveOutputs should be evaluated instead.");
    }


    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    virtual Constitutive::eConstitutiveType GetType() const override;
};
}
