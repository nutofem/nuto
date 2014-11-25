#ifndef MOISTURETRANSPORT_H
#define MOISTURETRANSPORT_H

#include "nuto/mechanics/constitutive/ConstitutiveBase.h"


namespace NuTo
{

class MoistureTransport : public ConstitutiveBase
{
public:
    MoistureTransport();

    //! @brief ... check compatibility between element type and type of constitutive relationship
    //! @param rElementType ... element type
    //! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
    virtual bool                                CheckElementCompatibility       (Element::eElementType rElementType) const override;

    //! @brief ... check parameters of the constitutive relationship
    //! if one check fails, an exception is thrwon
    virtual void                                CheckParameters                 () const override;


    //! @brief ... check if water phase density is positive
    //! @param rWaterPhaseDensity ... water phase density
    void                                        CheckWaterPhaseDensity          (double rWaterPhaseDensity) const;


    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    virtual Constitutive::eConstitutiveType     GetType                         () const override;

    //! @brief ... get water phase density
    //! @return ... water phase density
    virtual double                              GetWaterPhaseDensity            () const override;

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
    //! @return ... see brief explanation
    virtual bool                                HaveTmpStaticData               () const override;

    //! @brief ... set water phase density
    //! @param ... water phase density
    void                                        SetWaterPhaseDensity            (double rWaterPhaseDensity);


protected:

    //! @brief ... Water phase density \f$ \rho_w \f$
    double mRhoW = 1.0;

    //! @brief ... Vapor phase density \f$ \rho_v \f$
    double mRhoV = 1.0;

};


}



#endif // MOISTURETRANSPORT_H
