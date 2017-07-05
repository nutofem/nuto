#pragma once

#include "mechanics/constitutive/laws/AdditiveBase.h"
#include "mechanics/constitutive/staticData/IPAdditiveInputExplicit.h"

#include <set>
#include <vector>

namespace NuTo
{

class AdditiveInputExplicit : public AdditiveBase
{
public:
    //! @brief constructor
    AdditiveInputExplicit(const int& rNumTimeDerivatives)
        : AdditiveBase(rNumTimeDerivatives)
    {
    }

    // has no ip static data itself
    std::unique_ptr<Constitutive::IPConstitutiveLawBase> CreateIPLaw() override
    {
        return std::make_unique<Constitutive::IPAdditiveInputExplicit>(*this);
        mStaticDataAllocated = true;
    }


    //! @brief ... adds a constitutive law to a model that combines multiple constitutive laws (additive, parallel)
    //! @param rConstitutiveLaw ... additional constitutive law
    //! @param rModiesInput ... enum which defines wich input is modified by a constitutive law.
    virtual void AddConstitutiveLaw(NuTo::ConstitutiveBase& rConstitutiveLaw,
                                    Constitutive::eInput rModiesInput = Constitutive::eInput::NONE) override;

    //! @brief Evaluate the constitutive relation.
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.).
    //! @param rConstitutiveOutput Output of the constitutive law (stress, stiffness, heat flux etc.).
    template <int TDim>
    NuTo::eError Evaluate(const ConstitutiveInputMap&, const ConstitutiveOutputMap&)
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__,
                "Additive Law cannot be evaluated. Their IPAdditiveInputExplicit should be evaluated instead.");
    }

    //! @brief ... determines the constitutive inputs needed to evaluate the constitutive outputs
    //! @param rConstitutiveOutput ... desired constitutive outputs
    //! @param rInterpolationType ... interpolation type to determine additional inputs
    //! @return constitutive inputs needed for the evaluation
    virtual ConstitutiveInputMap GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput,
                                                       const InterpolationType& rInterpolationType) const override;

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    virtual Constitutive::eConstitutiveType GetType() const override;

    //! @brief ... pointer to constitutive law that delivers the output
    NuTo::ConstitutiveBase* mMainLaw = nullptr;

    //! @brief Vector storing which input they modify.
    //! @note Only for @ref AdditiveInputExplicit and @ref AdditiveInputImplicit
    std::vector<Constitutive::eInput> mInputsToModify;

    //! @brief Gets the enum of the sublaw output that is needed to calculate the specified derivative
    //! @param rParameter: Enum of the parameter whose derivative is needed
    //! @param rMainDerivative: The requested global derivative
    //! @return Sublaw derivative enum
    Constitutive::eOutput GetDerivativeEnumSublaw(Constitutive::eOutput rParameter,
                                                  Constitutive::eOutput rMainDerivative) const;

    template <int TDim>
    ConstitutiveOutputMap GetSublawOutputMap(const NuTo::ConstitutiveInputMap& rMainLawInputMap,
                                             const NuTo::ConstitutiveOutputMap& rMainLawOutputMap,
                                             int rSublawIndex) const;
};
}
