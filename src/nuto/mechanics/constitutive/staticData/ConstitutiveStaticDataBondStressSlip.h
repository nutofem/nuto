#pragma once

#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataBase.h"

//! @brief Class stores the maximum slip in the loading history of the element
//! @author Philip Huschke, BAM
//! @date October 2015
namespace NuTo
{
class IpDataStaticDataBase;

class ConstitutiveStaticDataBondStressSlip: public ConstitutiveStaticDataBase
{

#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

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

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION


private:
    //! @brief maximum slip in loading history
    double mSlip;

};

}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveStaticDataBondStressSlip)
#endif // ENABLE_SERIALIZATION
