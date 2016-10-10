#pragma once

//VHIRTHAMTODO replace with std::array
#include "nuto/math/FullVector.h"
#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataBase.h"

//! @brief ... base class, storing the static data (history variables) of a constitutive relationship
//! @author Volker Hirthammer, BAM
//! @date February 2015
namespace NuTo
{

class ConstitutiveStaticDataMoistureTransport : public ConstitutiveStaticDataBase
{
    friend class MoistureTransport;

public:
    //! @brief constructor
    ConstitutiveStaticDataMoistureTransport ();

    //! @brief copy constructor
    ConstitutiveStaticDataMoistureTransport (ConstitutiveStaticDataMoistureTransport const& rOther);

    //!@ brief reinterpret as moisture transport
    virtual ConstitutiveStaticDataMoistureTransport* AsMoistureTransport() override;

    //!@ brief reinterpret as moisture transport
    virtual const ConstitutiveStaticDataMoistureTransport* AsMoistureTransport()const override;

    //! @brief check, if the static data is compatible with a given element and a given constitutive model
    virtual bool CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType) const override;

    //! @brief clones (copies) the data
    ConstitutiveStaticDataMoistureTransport* Clone()const override;

    //! @brief Gets the actual sorption coefficients
    //! @return vector with the actual sorption coefficients
    FullVector<double,Eigen::Dynamic> GetCurrentSorptionCoeff() const;

    //! @brief Gets the relative humidity value of the last timestep
    //! @return relative humidity of the last timestep
    double GetLastRelHumValue() const;

    //! @brief Sets the sorption coefficients of the  last timestep
    //! @param rLastSorptionCoeff: vector with the sorption coefficients of the last timestep
    FullVector<double,Eigen::Dynamic> GetLastSorptionCoeff() const;

    //! @brief Gets the value of the Sorption History
    //! @return the sorption behavior of the last timestep (true = desorption, false = adsorption)  false = adsorption)
    bool GetSorptionHistoryDesorption() const;

    //! @brief Sets the actual sorption coefficients
    //! @param rLastSorptionCoeff: vector with the actual sorption coefficients
    void SetCurrentSorptionCoeff(FullVector<double,Eigen::Dynamic> rCurrentSorptionCoeff);

    //! @brief Sets the relative humidity value of the last timestep
    //! @param rLastRelHumValue: relative humidity of the last timestep
    void SetLastRelHumValue(double rLastRelHumValue);

    //! @brief Sets the sorption coefficients of the  last timestep
    //! @param rLastSorptionCoeff: vector with the sorption coefficients of the last timestep
    void SetLastSorptionCoeff(FullVector<double,Eigen::Dynamic> rLastSorptionCoeff);

    //! @brief Sets the value of the Sorption History
    //! @param rSorptionHistoryDesorption: sets the sorption behavior of the last timestep (true = desorption, false = adsorption)
    void SetSorptionHistoryDesorption(bool rSorptionHistoryDesorption);

protected:

    bool    mSorptionHistoryDesorption      = true;
    double  mLastRelHumValue                = 1.0;
    double  mLastJunctionPoint              = 0.0;
    double  mCurrentJunctionPoint           = 0.0;
    FullVector<double,Eigen::Dynamic> mCurrentSorptionCoeff     {{0.0, 0.0, 0.0}};
    FullVector<double,Eigen::Dynamic> mLastSorptionCoeff        {{0.0, 0.0, 0.0}};
};

}

