#ifndef DRYINGSHRINKAGE_H
#define DRYINGSHRINKAGE_H




#include "nuto/mechanics/constitutive/ConstitutiveBase.h"


namespace NuTo
{



//! @brief Drying shrinkage model
//! @author Volker Hirthammer, BAM
//! @date August 2015
class DryingShrinkage : public ConstitutiveBase
{
public:

    // Constructors / Destructors
    // --------------------------

    //! @brief Constructor
    DryingShrinkage();

    //! @brief Copy constructor
    //! @param rOther: Object that should be copied
    DryingShrinkage(const DryingShrinkage& rOther) = delete;

    //! @brief Move constructor
    //! @param rOther: Object that should be moved
    DryingShrinkage(DryingShrinkage&& rOther) = delete;

    //! @brief Destructor
    ~DryingShrinkage() = default;



    // Operator overloads
    // ------------------

    //! @brief Copy assignment operator
    //! @param rOther: Object that should be copied
    DryingShrinkage&                                operator = (const DryingShrinkage& rOther) = delete;

    //! @brief Move assignment operator
    //! @param rOther: Object that should be moved
    DryingShrinkage&                                operator = (DryingShrinkage&& rOther) = delete;



    // Member functions
    // ----------------

public:

    //! @brief ... check compatibility between element type and type of constitutive relationship
    //! @param rElementType ... element type
    //! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
    virtual bool                                    CheckElementCompatibility                                  (Element::eElementType rElementType) const override;

    //! @brief ... checks if the constitutive law has a specific parameter
    //! @param rIdentifier ... Enum to identify the requested parameter
    //! @return ... true/false
    virtual bool                                    CheckHaveParameter                                          (Constitutive::eConstitutiveParameter rIdentifier) const override;

    //! @brief ... checks if a constitutive law has an specific output
    //! @return ... true/false
    virtual bool                                    CheckOutputTypeCompatibility                                (Constitutive::Output::eOutput rOutputEnum) const override;

    //! @brief ... evaluate the constitutive relation in 2D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError                             Evaluate2D                                                  (ElementBase* rElement,
                                                                                                                 int rIp,
                                                                                                                 const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
                                                                                                                 std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput) override;

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    virtual Constitutive::eConstitutiveType         GetType                                                     () const override;

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
    //! @return ... see brief explanation
    virtual bool                                    HaveTmpStaticData                                           () const override;



protected:

    //! @brief ... check parameters of the constitutive relationship
    //! if one check fails, an exception is thrown
    virtual void                                    CheckParameters                                             () const override;

private:

    //! @brief Atmospheric pressure
    double                      mAtmosphericPressure    = 100000.0;

    //! @brief Density of water
    double                      mDensityWater           = 999.97;

    //! @brief Ideal gas constant
    static constexpr double     mIdealGasConstant       = 8.314459848;

    //! @brief Molar mass of water
    static constexpr double     mMolarMassWater         = 18.01528 / 1000.0;

    //! @brief Porosity
    double                      mPorosity               = 0.25;

    //! @brief Temperature in K
    double                      mTemperature            = 293.15;

};



}// namespace NuTo

#endif // DRYINGSHRINKAGE_H
