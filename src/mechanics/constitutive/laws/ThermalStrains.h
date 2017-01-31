#pragma once

#include "mechanics/constitutive/ConstitutiveBase.h"

namespace NuTo
{
namespace Constitutive
{
class IPConstitutiveLawBase;
} // namespace Constitutive
class ThermalStrains : public ConstitutiveBase
{
public:
    ThermalStrains() : ConstitutiveBase() {}

    //! @brief creates corresponding IPConstitutiveLaw
    std::unique_ptr<Constitutive::IPConstitutiveLawBase> CreateIPLaw() override;

    bool CheckDofCombinationComputable(Node::eDof rDofRow, Node::eDof rDofCol, int rTimeDerivative) const override;

    bool CheckElementCompatibility( Element::eElementType rElementType) const override;

    ConstitutiveInputMap GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput,
                                               const InterpolationType& rInterpolationType) const override;

    //! @brief Evaluate the constitutive relation.
    //! @param rConstitutiveInput Input to the constitutive law (strain, temp gradient etc.).
    //! @param rConstitutiveOutput Output to the constitutive law (stress, stiffness, heat flux etc.).
    template<int TDim>
    NuTo::eError Evaluate(const ConstitutiveInputMap& rConstitutiveInput,
                          const ConstitutiveOutputMap& rConstitutiveOutput);

    //! @brief Get type of constitutive relationship
    Constitutive::eConstitutiveType GetType() const override;

    bool HaveTmpStaticData() const override { return false; }

    // thermal expansion coefficient can be positive or negative, so do nothing
    void CheckParameters() const override {};

    void SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue) override;

    void SetParameterFunction(std::function<std::array<double, 2>(double)> ExpansionFunction) override;
private:
    //! Thermal expansion coefficient \f$ \alpha \f$.
    double mExpansionCoefficient = 0.0;

    std::function<std::array<double, 2>(double)> mNonlinearExpansionFunction;
};
}
