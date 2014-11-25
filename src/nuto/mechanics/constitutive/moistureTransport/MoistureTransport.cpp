#include "nuto/mechanics/MechanicsException.h"

#include "nuto/mechanics/constitutive/moistureTransport/MoistureTransport.h"

NuTo::MoistureTransport::MoistureTransport()
    : ConstitutiveBase()
{
    SetParametersValid();
}


//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool                                        NuTo::MoistureTransport::CheckElementCompatibility          (Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::BRICK8N:
        return false;
    case NuTo::Element::PLANE2D10N:
        return false;
    case NuTo::Element::PLANE2D15N:
        return false;
    case NuTo::Element::PLANE2D3N:
        return false;
    case NuTo::Element::PLANE2D4N:
        return false;
    case NuTo::Element::PLANE2D4NSPECTRALORDER2:
        return false;
    case NuTo::Element::PLANE2D4NSPECTRALORDER3:
        return false;
    case NuTo::Element::PLANE2D4NSPECTRALORDER4:
        return false;
    case NuTo::Element::PLANE2D6N:
        return false;
    case NuTo::Element::TETRAHEDRON4N:
        return false;
    case NuTo::Element::TETRAHEDRON10N:
        return false;
    case NuTo::Element::TRUSS1D2N:
        return true;
    case NuTo::Element::TRUSS1D3N:
        return false;
    default:
        return false;
    }
}


//! @brief ... check parameters of the constitutive relationship
//! if one check fails, an exception is thrwon
void                                        NuTo::MoistureTransport::CheckParameters                    () const
{
    CheckWaterPhaseDensity(mRhoW);
}

//! @brief ... check if water phase density is positive
//! @param rWaterPhaseDensity ... water phase density
void                                        NuTo::MoistureTransport::CheckWaterPhaseDensity             (double rWaterPhaseDensity) const
{
    if (rWaterPhaseDensity <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::::MoistureTransport::CheckWaterPhaseDensity] The water phase density must be a non-negative, non-zero value.");
    }
}

//! @brief ... get type of constitutive relationship
//! @return ... type of constitutive relationship
//! @sa eConstitutiveType
NuTo::Constitutive::eConstitutiveType       NuTo::MoistureTransport::GetType                            () const
{
    return NuTo::Constitutive::MOISTURE_TRANSPORT;
}

//! @brief ... get water phase density
//! @return ... water phase density
double                                      NuTo::MoistureTransport::GetWaterPhaseDensity               () const
{
    return mRhoW;
}

//! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
//! @return ... see brief explanation
bool                                        NuTo::MoistureTransport::HaveTmpStaticData                  () const
{
    return false;
}

//! @brief ... set water phase density
//! @param ... water phase density
void                                        NuTo::MoistureTransport::SetWaterPhaseDensity               (double rWaterPhaseDensity)
{
    CheckWaterPhaseDensity(rWaterPhaseDensity);
    mRhoW = rWaterPhaseDensity;
    SetParametersValid();
}
