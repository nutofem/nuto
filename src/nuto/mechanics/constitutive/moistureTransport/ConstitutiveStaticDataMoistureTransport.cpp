#include "ConstitutiveStaticDataMoistureTransport.h"
#include "nuto/mechanics/MechanicsException.h"

//! @brief constructor
                                                        NuTo::ConstitutiveStaticDataMoistureTransport::ConstitutiveStaticDataMoistureTransport  ()
   : NuTo::ConstitutiveStaticDataBase::ConstitutiveStaticDataBase()
{}

//! @brief copy constructor
                                                        NuTo::ConstitutiveStaticDataMoistureTransport::ConstitutiveStaticDataMoistureTransport  (ConstitutiveStaticDataMoistureTransport const& rOther)
{
    (*this) = rOther;
}

//!@ brief reinterpret as moisture transport
NuTo::ConstitutiveStaticDataMoistureTransport*          NuTo::ConstitutiveStaticDataMoistureTransport::AsMoistureTransport                      ()
{
    return this;
}

//!@ brief reinterpret as moisture transport
const NuTo::ConstitutiveStaticDataMoistureTransport*    NuTo::ConstitutiveStaticDataMoistureTransport::AsMoistureTransport                      ()const
{
    return this;
}

//! @brief check, if the static data is compatible with a given element and a given constitutive model
bool                                                    NuTo::ConstitutiveStaticDataMoistureTransport::CheckConstitutiveCompatibility           (NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType) const
{
    if (rConstitutiveType==NuTo::Constitutive::MOISTURE_TRANSPORT)
    {
        if (rElementType==NuTo::Element::TRUSS1D2N)
            return true;
        else
            return false;
    }
    else
        return false;
}

//! @brief clones (copies) the data
NuTo::ConstitutiveStaticDataMoistureTransport*          NuTo::ConstitutiveStaticDataMoistureTransport::Clone                                    ()const
{
    return new ConstitutiveStaticDataMoistureTransport(*this);
}

//! @brief Gets the actual sorption coefficients
//! @return vector with the actual sorption coefficients
NuTo::FullVector<double,Eigen::Dynamic>                 NuTo::ConstitutiveStaticDataMoistureTransport::GetActualSorptionCoeff                   () const
{
    return mActualSorptionCoeff;
}


//! @brief Gets the relative humidity value of the last timestep
//! @return relative humidity of the last timestep
double                                                  NuTo::ConstitutiveStaticDataMoistureTransport::GetLastRelHumValue                       () const
{
    return mLastRelHumValue;
}

//! @brief Gets the sorption coefficients of the  last timestep
//! @return vector with the sorption coefficients of the last timestep
NuTo::FullVector<double,Eigen::Dynamic>                 NuTo::ConstitutiveStaticDataMoistureTransport::GetLastSorptionCoeff                     () const
{
    return mLastSorptionCoeff;
}

//! @brief Gets the value of the Sorption History
//! @return the sorption behavior of the last timestep (true = desorption, false = adsorption)
bool                                                    NuTo::ConstitutiveStaticDataMoistureTransport::GetSorptionHistoryDesorption             () const
{
    return mSorptionHistoryDesorption;
}

//! @brief Sets the actual sorption coefficients
//! @param rLastSorptionCoeff: vector with the actual sorption coefficients
void                                                    NuTo::ConstitutiveStaticDataMoistureTransport::SetActualSorptionCoeff                   (FullVector<double,Eigen::Dynamic> rActualSorptionCoeff)
{
    switch (rActualSorptionCoeff.GetNumRows())
    {
        case 3:
        {
            mActualSorptionCoeff = rActualSorptionCoeff;
            break;
        }
        case 4:
        {
            mActualSorptionCoeff.Resize(3);
            for(int i=0; i<3; i++)
            {
                mActualSorptionCoeff(i) = rActualSorptionCoeff(i+1);
            }
            break;
        }
        default:
        {
            throw NuTo::MechanicsException("[NuTo::ConstitutiveStaticDataMoistureTransport::SetActualSorptionCoeff] The vector for the sorption coefficients must have 3 or 4 rows. --- Polynom of 3th degree --- in case of 4 coefficients the constant term will be deleted");
            break;
        }
    }
}

//! @brief Sets the relative humidity value of the last timestep
//! @param rLastRelHumValue: relative humidity of the last timestep
void                                                    NuTo::ConstitutiveStaticDataMoistureTransport::SetLastRelHumValue                       (double rLastRelHumValue)
{
    mLastRelHumValue = rLastRelHumValue;
}

//! @brief Sets the sorption coefficients of the  last timestep
//! @param rLastSorptionCoeff: vector with the sorption coefficients of the last timestep
void                                                    NuTo::ConstitutiveStaticDataMoistureTransport::SetLastSorptionCoeff                     (FullVector<double,Eigen::Dynamic> rLastSorptionCoeff)
{
    switch (rLastSorptionCoeff.GetNumRows())
    {
        case 3:
        {
            mLastSorptionCoeff = rLastSorptionCoeff;
            break;
        }
        case 4:
        {
            mLastSorptionCoeff.Resize(3);
            for(int i=0; i<3; i++)
            {
                mLastSorptionCoeff(i) = rLastSorptionCoeff(i+1);
            }
            break;
        }
        default:
        {
            throw NuTo::MechanicsException("[NuTo::ConstitutiveStaticDataMoistureTransport::SetLastSorptionCoeff] The vector for the sorption coefficients must have 3 or 4 rows. --- Polynom of 3th degree --- in case of 4 coefficients the constant term will be deleted");
            break;
        }
    }
}

//! @brief Sets the value of the Sorption History
//! @param rSorptionHistoryDesorption: sets the sorption behavior of the last timestep (true = desorption, false = adsorption)
void                                                    NuTo::ConstitutiveStaticDataMoistureTransport::SetSorptionHistoryDesorption             (bool rSorptionHistoryDesorption)
{
    mSorptionHistoryDesorption = rSorptionHistoryDesorption;
}
