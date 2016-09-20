#pragma once

#include "nuto/mechanics/constitutive/ConstitutiveBase.h"

namespace NuTo
{
class AdditiveBase : public ConstitutiveBase
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

    //! @brief Adds a constitutive law to a model that combines multiple constitutive laws (additive, parallel)
    //! @param rConstitutiveLaw Constitutive law to be added.
    //! @param rModiesInput Enum which defines which input is modified by a constitutive law.
    virtual void AddConstitutiveLaw(NuTo::ConstitutiveBase& rConstitutiveLaw,
            Constitutive::Input::eInput rModiesInput = Constitutive::Input::NONE);

    //! @brief Check compatibility between element type and type of constitutive relationship.
    //! @param rElementType Element type.
    //! @return `true` if the element is compatible with the constitutive relationship, `false` otherwise.
    virtual bool CheckElementCompatibility(Element::eElementType rElementType) const override;

    //! @brief Check parameters of the constitutive relationship.
    virtual void CheckParameters() const override;

    //! @brief Checks whether a material model has tmp static data, which has to be updated before stress or
    //! stiffness are calculated.
    virtual bool HaveTmpStaticData() const override;

    //! @brief Determines the constitutive inputs needed to evaluate the constitutive outputs.
    //! @param rConstitutiveOutput Desired constitutive outputs.
    //! @param rInterpolationType Interpolation type to determine additional inputs.
    //! @return Constitutive inputs needed for the evaluation.
    virtual ConstitutiveInputMap GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput,
                                                       const InterpolationType& rInterpolationType) const override;

    virtual bool CheckDofCombinationComputable(Node::eDof rDofRow, Node::eDof rDofCol,
            int rTimeDerivative) const;

protected:
    //! @brief Debug variable to avoid that a constitutive law can be attached after allocation of static data.
    mutable bool mStaticDataAllocated = false;

    //! @brief Vector storing the pointers to the sublaws.
    std::vector<NuTo::ConstitutiveBase*> mSublaws;

    //! @brief Vector of all the computable DOF combinations.
    std::vector<std::set<std::pair<Node::eDof,Node::eDof>>> mComputableDofCombinations;

    //! @brief Adds all calculable DOF combinations of an attached constitutive law to an internal storage.
    //! @param rConstitutiveLaw Constitutive law whose DOF combinations are added to additive law.
    void AddCalculableDofCombinations(NuTo::ConstitutiveBase& rConstitutiveLaw);
};
}
