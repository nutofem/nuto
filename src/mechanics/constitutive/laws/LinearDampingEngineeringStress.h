#pragma once

#include "mechanics/constitutive/ConstitutiveBase.h"

namespace NuTo
{
namespace Constitutive
{
class IPConstitutiveLawBase;
} // namespace Constitutive
class LinearDampingEngineeringStress : public ConstitutiveBase
{
public:
    LinearDampingEngineeringStress();

    std::unique_ptr<Constitutive::IPConstitutiveLawBase> CreateIPLaw() override;

    //! @brief ... determines the constitutive inputs needed to evaluate the constitutive outputs
    //! @param rConstitutiveOutput ... desired constitutive outputs
    //! @param rInterpolationType ... interpolation type to determine additional inputs
    //! @return constitutive inputs needed for the evaluation
    virtual ConstitutiveInputMap GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput,
                                                       const InterpolationType& rInterpolationType) const override;

    //! @brief ... determines which submatrices of a multi-doftype problem can be solved by the constitutive law
    //! @param rDofRow ... row dof
    //! @param rDofCol ... column dof
    //! @param rTimeDerivative ... time derivative
    virtual bool CheckDofCombinationComputable(Node::eDof rDofRow, Node::eDof rDofCol,
                                               int rTimeDerivative) const override;

    //! @brief ... check parameters of the constitutive relationship
    //! if one check fails, an exception is thrown
    virtual void CheckParameters() const override;

    //! @brief ... evaluate the constitutive relation in 1D
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    template <int TDim>
    void Evaluate(const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput);

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    virtual Constitutive::eConstitutiveType GetType() const override;

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or
    //! stiffness are calculated)
    //! @return ... see brief explanation
    virtual bool HaveTmpStaticData() const override;


    // Getter / Setter
    // ---------------

    //! @brief ... gets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @return ... value of the requested variable
    virtual double GetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier) const override;

    //! @brief ... sets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @param rValue ... new value for requested variable
    virtual void SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue) override;


private:
    double mDampingCoefficient;
};

} // namespace NuTo
