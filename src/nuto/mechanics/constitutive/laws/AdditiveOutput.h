#pragma once

#include "nuto/mechanics/constitutive/laws/AdditiveBase.h"

namespace NuTo
{

class AdditiveOutput : public AdditiveBase
{
public:

    //! @brief constructor
    AdditiveOutput() : AdditiveBase()
    {
        //VHIRTHAMTODO ---> Get number time derivatives during construction (as parameter)
        mComputableDofCombinations.resize(2);
    }

    // has no ip static data itself
    std::unique_ptr<Constitutive::IPConstitutiveLawBase> CreateIPLaw()
    {
        return std::make_unique<Constitutive::IPConstitutiveLawWithoutData<AdditiveOutput>>(*this);
    }


    template <int TDim>
    NuTo::eError Evaluate(const ConstitutiveInputMap& rConstitutiveInput,
                          const ConstitutiveOutputMap& rConstitutiveOutput);


    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    virtual Constitutive::eConstitutiveType GetType() const override;
};
}
