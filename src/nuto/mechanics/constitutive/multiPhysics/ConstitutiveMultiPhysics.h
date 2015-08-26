#ifndef CONSTITUTIVEMULTIPHYSICS_H
#define CONSTITUTIVEMULTIPHYSICS_H

#include <vector>

#include <nuto/math/FullVector.h>

#include <nuto/mechanics/constitutive/ConstitutiveBase.h>
#include <nuto/mechanics/MechanicsException.h>

namespace NuTo
{


//! @brief ... base class for multi physics
//! @author Volker Hirthammer, BAM
//! @date May 2015
class ConstitutiveMultiPhysics : public ConstitutiveBase
{
public:

    //! @brief ... default constructor
    ConstitutiveMultiPhysics();

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase*                 AllocateStaticDataEngineeringStress_EngineeringStrain1D     (const ElementBase* rElement) const override;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase*                 AllocateStaticDataEngineeringStress_EngineeringStrain2D     (const ElementBase* rElement) const override;

    //! @brief ... check compatibility between element type and type of constitutive relationship
    //! @param rElementType ... element type
    //! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
    virtual bool                                CheckElementCompatibility                                   (NuTo::Element::eElementType rElementType) const override;

    //! @brief ... check parameters of the constitutive relationship
    //! if one check fails, an exception is thrwon
    virtual void                                CheckParameters                                             () const override;

    //! @brief ... evaluate the constitutive relation of every attached constitutive law in 1D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    virtual NuTo::Error::eError                 Evaluate1D                                                  (ElementBase* rElement,
                                                                                                             int rIp,
                                                                                                             const std::map<NuTo::Constitutive::Input::eInput, const NuTo::ConstitutiveInputBase*>& rConstitutiveInput,
                                                                                                             std::map<NuTo::Constitutive::Output::eOutput, NuTo::ConstitutiveOutputBase*>& rConstitutiveOutput) override;

    //! @brief ... evaluate the constitutive relation of every attached constitutive law in 2D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    virtual NuTo::Error::eError                 Evaluate2D                                                  (ElementBase* rElement,
                                                                                                             int rIp,
                                                                                                             const std::map<NuTo::Constitutive::Input::eInput, const NuTo::ConstitutiveInputBase*>& rConstitutiveInput,
                                                                                                             std::map<NuTo::Constitutive::Output::eOutput, NuTo::ConstitutiveOutputBase*>& rConstitutiveOutput) override;

    //! @brief ... evaluate the constitutive relation of every attached constitutive law relation in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    virtual NuTo::Error::eError                 Evaluate3D                                                  (ElementBase* rElement,
                                                                                                             int rIp,
                                                                                                             const std::map<NuTo::Constitutive::Input::eInput, const NuTo::ConstitutiveInputBase*>& rConstitutiveInput,
                                                                                                             std::map<NuTo::Constitutive::Output::eOutput, NuTo::ConstitutiveOutputBase*>& rConstitutiveOutput) override;


    //! @brief ... gets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @return ... value of the requested variable
    virtual double GetParameterDouble (Constitutive::eConstitutiveParameter rIdentifier) const override;


    //! @brief ... sets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @param rValue ... new value for requested variable
    virtual void SetParameterBool(Constitutive::eConstitutiveParameter rIdentifier, bool rValue) override;

    //! @brief ... sets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @param rValue ... new value for requested variable
    virtual void SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue) override;

    //! @brief ... sets a parameter of the constitutive law which is selected by an enum
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @param rValue ... new value for requested variable
    virtual void SetParameterFullVectorDouble(Constitutive::eConstitutiveParameter rIdentifier, NuTo::FullVector<double,Eigen::Dynamic> rValue) override;

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    virtual Constitutive::eConstitutiveType     GetType                                                     () const override;

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
    //! @return ... see brief explanation
    virtual bool                                HaveTmpStaticData                                           () const override;

    //! @brief ... adds a constitutive law to a multi physics model
    //! @param ... additional constitutive law
    virtual void                                MultiPhysicsAddConstitutiveLaw                              (NuTo::ConstitutiveBase* rConstitutiveLaw) override;


private:
    //! @brief ... list of pointers to every included constitutive law
    std::vector<NuTo::ConstitutiveBase*> mConstitutiveLaws;
};


} //namespace NuTo

#endif // CONSTITUTIVEMULTIPHYSICS_H
