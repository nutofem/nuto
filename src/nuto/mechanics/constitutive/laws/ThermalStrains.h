#pragma once

#include "nuto/mechanics/constitutive/ConstitutiveBase.h"

namespace NuTo
{
class ThermalStrains : public ConstitutiveBase
{
public:
    ThermalStrains() : ConstitutiveBase() {}

    //! @brief Check which DOF combination can be evaluated by this law
    //! @param rDofRow DOF belonging to the row in the assembled system
    //! @param rDofCol DOF belonging to the columns in the assembled system
    //! @param rTimeDerivative The admissable order of the time derivative
    bool CheckDofCombinationComputable(Node::eDof rDofRow,
            Node::eDof rDofCol, int rTimeDerivative) const override;

    //! @brief Check compatibility between element type and the constitutive law
    //! @param rElementType Element type
    bool CheckElementCompatibility( Element::eElementType rElementType) const override;

    //! @brief Determines the constitutive inputs needed to evaluate the constitutive outputs
    //! @param rConstitutiveOutput Desired constitutive outputs
    //! @param rInterpolationType Interpolation type to determine additional inputs
    //! @return Constitutive inputs needed for the evaluation
    ConstitutiveInputMap GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput,
                                               const InterpolationType& rInterpolationType) const override;

    //! @brief Evaluate the constitutive relation in 1D
    //! @param rElement Element
    //! @param rIp Integration point
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput Output of the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate1D(ElementBase* rElement, int rIntegrationPoint,
                                   const ConstitutiveInputMap& rConstitutiveInput,
                                   const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        return Evaluate<1>(rElement, rIntegrationPoint, rConstitutiveInput, rConstitutiveOutput);
    }

    //! @brief Evaluate the constitutive relation in 2D
    //! @param rElement Element
    //! @param rIp Integration point
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput Output of the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate2D(ElementBase* rElement, int rIntegrationPoint,
                                   const ConstitutiveInputMap& rConstitutiveInput,
                                   const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        return Evaluate<2>(rElement, rIntegrationPoint, rConstitutiveInput, rConstitutiveOutput);
    }

    //! @brief Evaluate the constitutive relation in 3D
    //! @param rElement Element
    //! @param rIp Integration point
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput Output of the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate3D(ElementBase* rElement, int rIntegrationPoint,
                                   const ConstitutiveInputMap& rConstitutiveInput,
                                   const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        return Evaluate<3>(rElement, rIntegrationPoint, rConstitutiveInput, rConstitutiveOutput);
    }

    //! @brief Get type of constitutive relationship
    Constitutive::eConstitutiveType GetType() const override
    {
        return NuTo::Constitutive::THERMAL_STRAINS;
    }

    //! @brief True if material model has tmp static data, which has to be
    //! updated before stress or stiffness are calculated
    bool HaveTmpStaticData() const override { return false; }

    //! @brief Check parameters of the constitutive relationship
    // thermal expansion coefficient can be positive or negative, so do nothing
    void CheckParameters() const override {};

    //! @brief Sets a parameter of the constitutive law
    //! @param rIdentifier Enum to identify the requested parameter
    //! @param rValue New value for requested variable
    void SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue) override;

    //! @brief Create new static data object for an integration point.
    //! @return Pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticData1D(const ElementBase* rElement) const override {return nullptr;}

    //! @brief Create new static data object for an integration point.
    //! @return Pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticData2D(const ElementBase* rElement) const override {return nullptr;}

    //! @brief Create new static data object for an integration point.
    //! @return Pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticData3D(const ElementBase* rElement) const override {return nullptr;}

private:

    //! Thermal expansion coefficient \f$ \alpha \f$.
    double mExpansionCoefficient = 0.0;

    template <int TDim>
    NuTo::Error::eError Evaluate(ElementBase* rElement, int rIntegrationPoint,
                                 const ConstitutiveInputMap& rConstitutiveInput,
                                 const ConstitutiveOutputMap& rConstitutiveOutput);

};
}
