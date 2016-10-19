#pragma once

#include "nuto/mechanics/constitutive/laws/AdditiveBase.h"
#include "nuto/mechanics/constitutive/staticData/Composite.h"

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

    template <int TDim>
    NuTo::eError Evaluate(const ConstitutiveInputMap& rConstitutiveInput,
            const ConstitutiveOutputMap& rConstitutiveOutput,
            Constitutive::StaticData::Component* staticData);

    //! @brief Evaluate the constitutive relation.
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.).
    //! @param rConstitutiveOutput Output of the constitutive law (stress, stiffness, heat flux etc.).
    virtual NuTo::eError Evaluate1D(const ConstitutiveInputMap& rConstitutiveInput,
                                           const ConstitutiveOutputMap& rConstitutiveOutput,
                                           Constitutive::StaticData::Component* staticData) override
    {
        return Evaluate<1>(rConstitutiveInput, rConstitutiveOutput, staticData);
    }

    //! @brief Evaluate the constitutive relation.
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.).
    //! @param rConstitutiveOutput Output of the constitutive law (stress, stiffness, heat flux etc.).
    virtual NuTo::eError Evaluate2D(const ConstitutiveInputMap& rConstitutiveInput,
                                           const ConstitutiveOutputMap& rConstitutiveOutput,
                                           Constitutive::StaticData::Component* staticData) override
    {
        return Evaluate<2>(rConstitutiveInput, rConstitutiveOutput, staticData);
    }

    //! @brief Evaluate the constitutive relation.
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.).
    //! @param rConstitutiveOutput Output of the constitutive law (stress, stiffness, heat flux etc.).
    virtual NuTo::eError Evaluate3D(const ConstitutiveInputMap& rConstitutiveInput,
                                           const ConstitutiveOutputMap& rConstitutiveOutput,
                                           Constitutive::StaticData::Component* staticData) override
    {
        return Evaluate<3>(rConstitutiveInput, rConstitutiveOutput, staticData);
    }

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    virtual Constitutive::eConstitutiveType GetType() const override;
};
}
