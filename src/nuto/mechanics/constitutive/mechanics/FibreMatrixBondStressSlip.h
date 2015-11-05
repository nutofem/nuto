//============================================================================
// Name        : FibreMatrixBondStressSlip.cpp
// Author      : Philip Huschke
// Version     : 26 Aug 2015
// Copyright   :
// Description : Constitutive law for the interface between fibre and matrix
//============================================================================

#pragma once

#include "nuto/mechanics/constitutive/ConstitutiveBase.h"

namespace NuTo
{

class FibreMatrixBondStressSlip: public ConstitutiveBase
{

public:
    FibreMatrixBondStressSlip();

    //! @brief ... evaluate the constitutive relation in 1D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rUpdateHistory ... update history variables after leaving the routine
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate1D(ElementBase* rElement, int rIp, const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput, std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput) override;

    //! @brief ... evaluate the constitutive relation in 2D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rUpdateHistory ... update history variables after leaving the routine
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate2D(ElementBase* rElement, int rIp, const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput, std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput) override;

    //! @brief ... evaluate the constitutive relation in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rUpdateHistory ... update history variables after leaving the routine
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate3D(ElementBase* rElement, int rIp, const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput, std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput) override;

    //! @brief ... gets a variable of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested variable
    //! @return ... value of the requested variable
    virtual double GetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier) const override;

    //! @brief ... sets a variable of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested variable
    //! @param rValue ... new value for requested variable
    virtual void SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue) override;

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    Constitutive::eConstitutiveType GetType() const override;

    //! @brief ... check parameters of the constitutive relationship
    void CheckParameters() const;

    //! @brief ... check compatibility between element type and type of constitutive relationship
    //! @param rElementType ... element type
    //! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
    bool CheckElementCompatibility(Element::eElementType rElementType) const override;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain2D(const ElementBase* rElement) const override;

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    //! @param rLogger stream for the output
    void Info(unsigned short rVerboseLevel, Logger& rLogger) const override;

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
    //! @return ... see brief explanation
    bool HaveTmpStaticData() const override
    {
        return false;
    }

private:

    //! @brief ... maximum bond stress
    double mMaxBondStress;

    //! @brief ... residual bond stress
    double mResidualBondStress;

    //! @brief ... slip at maximum bond stress
    double mSlipAtMaxBondStress;

    //! @brief ... slip at residual bond stress (bond stress stays constant for greater slips)
    double mSlipAtResidualBondStress;

    //! @brief ... parameter that controls the pre peak shape of the bond stress-slip curve
    double mAlpha;

    //! @brief ... normal penalty stiffness to avoid penetration
    double mNormalStiffness;
};

}

