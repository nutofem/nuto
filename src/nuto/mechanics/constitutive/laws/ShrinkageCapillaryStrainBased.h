#pragma once

#include "nuto/mechanics/constitutive/ConstitutiveBase.h"


//VHIRTHAMTODO Implement missing check routines --- at the moment only temperature
namespace NuTo
{

class ShrinkageCapillaryStrainBased : public ConstitutiveBase
{
public:

    //! @brief constructor
    ShrinkageCapillaryStrainBased()
        : ConstitutiveBase()
    {}

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticData1D(const ElementBase* rElement) const override {return nullptr;}

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticData2D(const ElementBase* rElement) const override  {return nullptr;}

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticData3D(const ElementBase* rElement) const override  {return nullptr;}

    //! @brief ... determines which submatrices of a multi-doftype problem can be solved by the constitutive law
    //! @param rDofRow ... row dof
    //! @param rDofCol ... column dof
    //! @param rTimeDerivative ... time derivative
    virtual bool CheckDofCombinationComputable(Node::eDof rDofRow,
                                                Node::eDof rDofCol,
                                                int rTimeDerivative) const override;

    //! @brief ... check compatibility between element type and type of constitutive relationship
    //! @param rElementType ... element type
    //! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
    virtual bool CheckElementCompatibility( Element::eElementType rElementType) const override
    {
        switch (rElementType)
        {
        case NuTo::Element::CONTINUUMELEMENT:
        case NuTo::Element::CONTINUUMBOUNDARYELEMENT:
        case NuTo::Element::CONTINUUMBOUNDARYELEMENTCONSTRAINEDCONTROLNODE:
            return true;
        default:
            return false;
        }
    }

    //! @brief ... check parameters of the constitutive relationship
    //! if one check fails, an exception is thrwon
    virtual void CheckParameters() const override;

    //! @brief ... evaluate the constitutive relation
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    template <int TDim>
    NuTo::Error::eError EvaluateShrinkageCapillary( ElementBase* rElement,
                                                    int rIp,
                                                    const ConstitutiveInputMap& rConstitutiveInput,
                                                    const ConstitutiveOutputMap& rConstitutiveOutput);


    //! @brief ... evaluate the constitutive relation of every attached constitutive law in 1D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    virtual NuTo::Error::eError Evaluate1D(ElementBase* rElement,
                                           int rIp,
                                           const ConstitutiveInputMap& rConstitutiveInput,
                                           const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        return EvaluateShrinkageCapillary<1>(rElement,
                                             rIp,
                                             rConstitutiveInput,
                                             rConstitutiveOutput);
    }

    //! @brief ... evaluate the constitutive relation of every attached constitutive law in 2D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    virtual NuTo::Error::eError Evaluate2D(ElementBase* rElement,
                                           int rIp,
                                           const ConstitutiveInputMap& rConstitutiveInput,
                                           const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        return EvaluateShrinkageCapillary<2>(rElement,
                                             rIp,
                                             rConstitutiveInput,
                                             rConstitutiveOutput);
    }

    //! @brief ... evaluate the constitutive relation of every attached constitutive law relation in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    virtual NuTo::Error::eError Evaluate3D(ElementBase* rElement,
                                           int rIp,
                                           const ConstitutiveInputMap& rConstitutiveInput,
                                           const ConstitutiveOutputMap& rConstitutiveOutput) override
    {
        return EvaluateShrinkageCapillary<3>(rElement,
                                             rIp,
                                             rConstitutiveInput,
                                             rConstitutiveOutput);
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
    virtual Constitutive::eConstitutiveType GetType() const override
    {
        return NuTo::Constitutive::SHRINKAGE_CAPILLARY_STRAIN_BASED;
    }

    //! @brief ... gets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @return ... value of the requested variable
    virtual double GetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier) const override;

    //! @brief ... sets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @param rValue ... new value for requested variable
    virtual void SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue) override;

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
    //! @return ... see brief explanation
    virtual bool HaveTmpStaticData() const override
    {
        return false;
    }

private:

    //! @brief Temperature in K
    double mTemperature = 293.15;

    //! @brief Macroscopic bulk modulus (whole structure)
    double mMacroscopicBulkModulus;

    //! @brief Solid phase bulk modulus
    double mSolidPhaseModulus;
};

}
