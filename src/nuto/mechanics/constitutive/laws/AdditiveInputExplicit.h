#pragma once

#include "nuto/mechanics/constitutive/laws/AdditiveBase.h"

#include <set>
#include <vector>

namespace NuTo
{

class AdditiveInputExplicit : public AdditiveBase
{
public:
    //! @brief Create a new static data object for an integration point.
    //! @return Pointer to the new object.
    Constitutive::StaticData::Component* AllocateStaticData1D(const ElementBase* rElement) const override
    {
        return AllocateStaticData<1>(rElement);
    }

    //! @brief Create a new static data object for an integration point.
    //! @return Pointer to the new object.
    Constitutive::StaticData::Component* AllocateStaticData2D(const ElementBase* rElement) const override
    {
        return AllocateStaticData<2>(rElement);
    }

    //! @brief Create a new static data object for an integration point.
    //! @return Pointer to the new object.
    Constitutive::StaticData::Component* AllocateStaticData3D(const ElementBase* rElement) const override
    {
        return AllocateStaticData<3>(rElement);
    }

    template <int TDim>
    Constitutive::StaticData::Component* AllocateStaticData(const NuTo::ElementBase *rElement) const;

    //! @brief constructor
    AdditiveInputExplicit() : AdditiveBase()
    {
        //VHIRTHAMTODO ---> Get number time derivatives during construction (as parameter)
        mComputableDofCombinations.resize(2); 
    }

    //! @brief ... adds a constitutive law to a model that combines multiple constitutive laws (additive, parallel)
    //! @param rConstitutiveLaw ... additional constitutive law
    //! @param rModiesInput ... enum which defines wich input is modified by a constitutive law.
    virtual void  AddConstitutiveLaw(NuTo::ConstitutiveBase& rConstitutiveLaw, 
            Constitutive::eInput rModiesInput = Constitutive::eInput::NONE) override;

    //! @brief Evaluate the constitutive relation.
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.).
    //! @param rConstitutiveOutput Output of the constitutive law (stress, stiffness, heat flux etc.).
    virtual NuTo::eError Evaluate1D(const ConstitutiveInputMap& rConstitutiveInput,
                                           const ConstitutiveOutputMap& rConstitutiveOutput,
                                           Constitutive::StaticData::Component* staticData) override
    {
        return EvaluateAdditiveInputExplicit<1>(rConstitutiveInput, rConstitutiveOutput, staticData);
    }


    //! @brief Evaluate the constitutive relation.
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.).
    //! @param rConstitutiveOutput Output of the constitutive law (stress, stiffness, heat flux etc.).
    virtual NuTo::eError Evaluate2D(const ConstitutiveInputMap& rConstitutiveInput,
                                           const ConstitutiveOutputMap& rConstitutiveOutput,
                                           Constitutive::StaticData::Component* staticData) override
    {
        return EvaluateAdditiveInputExplicit<2>(rConstitutiveInput, rConstitutiveOutput, staticData);
    }

    //! @brief Evaluate the constitutive relation.
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.).
    //! @param rConstitutiveOutput Output of the constitutive law (stress, stiffness, heat flux etc.).
    virtual NuTo::eError Evaluate3D(const ConstitutiveInputMap& rConstitutiveInput,
                                           const ConstitutiveOutputMap& rConstitutiveOutput,
                                           Constitutive::StaticData::Component* staticData) override
    {
        return EvaluateAdditiveInputExplicit<3>(rConstitutiveInput, rConstitutiveOutput, staticData);
    }

    //! @brief ... determines the constitutive inputs needed to evaluate the constitutive outputs
    //! @param rConstitutiveOutput ... desired constitutive outputs
    //! @param rInterpolationType ... interpolation type to determine additional inputs
    //! @return constitutive inputs needed for the evaluation
    virtual ConstitutiveInputMap GetConstitutiveInputs( const ConstitutiveOutputMap& rConstitutiveOutput,
                                                        const InterpolationType& rInterpolationType) const override;

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    virtual Constitutive::eConstitutiveType GetType() const override;

private:
    //! @brief Applies sublaw dependent modifications to the main laws inputs and the global outputs
    //! @param rMainLawInput: Input map of the main law
    //! @param rConstitutiveOutput: Global output map of this constitutive law
    //! @param rSublawOutput: Output map of a sublaw
    template <int TDim>
    void ApplySublawOutputs(const ConstitutiveInputMap& rMainLawInput,
                            const ConstitutiveOutputMap& rConstitutiveOutput,
                            const ConstitutiveOutputMap& rSublawOutput);

    //! @brief Calculates the derivatives that depend on main law and sublaw outputs
    //! @param rConstitutiveOutput: Global output map of this constitutive law
    //! @param rSublawOutputVec: Vector with all the sublaws output maps
    template <int TDim>
    void CalculateDerivatives(const ConstitutiveOutputMap &rConstitutiveOutput,
                              std::vector<ConstitutiveOutputMap> &rSublawOutputVec);

    //! @brief Evaluate the constitutive relation.
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.).
    //! @param rConstitutiveOutput Output of the constitutive law (stress, stiffness, heat flux etc.).
    //! @param staticData Pointer to the static data.
    template <int TDim>
    eError EvaluateAdditiveInputExplicit(const ConstitutiveInputMap& rConstitutiveInput,
            const ConstitutiveOutputMap& rConstitutiveOutput, Constitutive::StaticData::Component* staticData);

    //! @brief Gets the enum of the sublaw output that is needed to calculate the specified derivative
    //! @param rParameter: Enum of the parameter whose derivative is needed
    //! @param rMainDerivative: The requested global derivative
    //! @return Sublaw derivative enum
    Constitutive::eOutput GetDerivativeEnumSublaw(Constitutive::eOutput rParameter,
                                                          Constitutive::eOutput rMainDerivative) const;


    //! @brief Gets a modified output map for the sublaws depending on the main laws inputs and the specified modified Input (see Member: mModifiedInputs)
    //! @return Modified output map
    template <int TDim>
    ConstitutiveOutputMap GetSublawOutputMap(const ConstitutiveInputMap& rMainLawInputMap,
                                             const NuTo::ConstitutiveOutputMap& rMainLawOutputMap,
                                             unsigned int rSublawIndex) const;

    //! @brief ... pointer to constitutive law that delivers the output
    NuTo::ConstitutiveBase* mMainLaw = nullptr;

    //! @brief Vector storing which input they modify.
    //! @note Only for @ref AdditiveInputExplicit and @ref AdditiveInputImplicit
    std::vector<Constitutive::eInput> mInputsToModify;
};
}
