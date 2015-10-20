#pragma once

#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"

//! @brief Class stores the maximum slip in the loading history of the element
//! @author Philip Huschke, BAM
//! @date October 2015
namespace NuTo
{
class IpDataStaticDataBase;

class ConstitutiveStaticDataBondStressSlip: public ConstitutiveStaticDataBase
{
public:
    //! @brief constructor
    ConstitutiveStaticDataBondStressSlip();

    //! @brief assignment operator
    ConstitutiveStaticDataBondStressSlip& operator=(ConstitutiveStaticDataBondStressSlip const& rOther);

    ConstitutiveStaticDataBondStressSlip* Clone() const override;

    //! @brief check, if the static data is compatible with a given element and a given constitutive model
    virtual bool CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType) const override;

    //!@ brief reinterpret as bond stress slip static data
    virtual ConstitutiveStaticDataBondStressSlip* AsBondStressSlip() override;

    //!@ brief reinterpret as bond stress slip static data
    virtual const ConstitutiveStaticDataBondStressSlip* AsBondStressSlip() const override;

    //! @brief sets the slip
    void SetSlip(double rSlip);

    //! @brief gets the slip
    const double GetSlip() const;

private:
    //! @brief maximum slip in loading history
    double mSlip;

};

}

