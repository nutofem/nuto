#pragma once

#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataBase.h"

//! @brief ... base class, storing the static data (history variables) of a constitutive relationship
//! @author Philip Huschke, BAM
//! @date June 2016
namespace NuTo
{

class ConstitutiveStaticDataElasticEnergyDensity : public ConstitutiveStaticDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
	//! @brief constructor
    ConstitutiveStaticDataElasticEnergyDensity();

    //! @brief copy constructor
    ConstitutiveStaticDataElasticEnergyDensity(ConstitutiveStaticDataElasticEnergyDensity const& rOther) = default;

    //! @brief clones (copies) the data
    ConstitutiveStaticDataElasticEnergyDensity* Clone() const override
    {
        return new ConstitutiveStaticDataElasticEnergyDensity(*this);
    }

    //! @brief assignment operator
    ConstitutiveStaticDataElasticEnergyDensity& operator= (ConstitutiveStaticDataElasticEnergyDensity const& rOther) = default;

    //! @brief check, if the static data is compatible with a given element and a given constitutive model
    virtual bool CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const;

    //!@ brief reinterpret as nonlocal damage2d static data
    ConstitutiveStaticDataElasticEnergyDensity* AsElasticEnergyDensity() override
    {
        return this;
    }


    //!@ brief reinterpret as nonlocal damage2d static data
    const ConstitutiveStaticDataElasticEnergyDensity* AsElasticEnergyDensity()const override
    {
        return this;
    }

    void SetMaxElasticEnergyDensity(const double rMaxElasticEnergyDensity)
    {
        mMaxElasticEnergyDensity = rMaxElasticEnergyDensity;
    }

    double GetMaxElasticEnergyDensity() const
    {
        return mMaxElasticEnergyDensity;
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION



protected:
    //! @brief maximal value of the elastic energy density
    double mMaxElasticEnergyDensity;

};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveStaticDataElasticEnergyDensity)
#endif // ENABLE_SERIALIZATION
