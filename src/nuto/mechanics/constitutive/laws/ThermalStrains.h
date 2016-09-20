#pragma once

#include "nuto/mechanics/constitutive/ConstitutiveBase.h"

namespace NuTo
{
class ThermalStrains : public ConstitutiveBase
{
public:
    ThermalStrains() : ConstitutiveBase() {}

    bool CheckDofCombinationComputable(Node::eDof rDofRow, Node::eDof rDofCol, int rTimeDerivative) const override;

    bool CheckElementCompatibility( Element::eElementType rElementType) const override;

    ConstitutiveInputMap GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput,
                                               const InterpolationType& rInterpolationType) const override;
    template<int TDim>
    NuTo::Error::eError Evaluate(const ConstitutiveInputMap& rConstitutiveInput,
                                 const ConstitutiveOutputMap& rConstitutiveOutput,
                                 Constitutive::StaticData::Component* staticData);

    NuTo::Error::eError Evaluate1D(const ConstitutiveInputMap& rConstitutiveInput,
                                   const ConstitutiveOutputMap& rConstitutiveOutput,
                                   Constitutive::StaticData::Component* staticData) override
    {
        return Evaluate<1>(rConstitutiveInput, rConstitutiveOutput, staticData);
    }

    NuTo::Error::eError Evaluate2D(const ConstitutiveInputMap& rConstitutiveInput,
                                   const ConstitutiveOutputMap& rConstitutiveOutput,
                                   Constitutive::StaticData::Component* staticData) override
    {
        return Evaluate<2>(rConstitutiveInput, rConstitutiveOutput, staticData);
    }

    NuTo::Error::eError Evaluate3D(const ConstitutiveInputMap& rConstitutiveInput,
                                   const ConstitutiveOutputMap& rConstitutiveOutput,
                                   Constitutive::StaticData::Component* staticData) override
    {
        return Evaluate<3>(rConstitutiveInput, rConstitutiveOutput, staticData);
    }

    Constitutive::eConstitutiveType GetType() const override
    {
        return NuTo::Constitutive::THERMAL_STRAINS;
    }

    bool HaveTmpStaticData() const override { return false; }

    // thermal expansion coefficient can be positive or negative, so do nothing
    void CheckParameters() const override {};

    void SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue) override;

    void SetParameterFunction(std::function<std::array<double, 2>(double)> ExpansionFunction) override;
private:
    //! Thermal expansion coefficient \f$ \alpha \f$.
    double mExpansionCoefficient = 0.0;

    std::function<std::array<double, 2>(double)> NonlinearExpansionCoeff;
};
}
