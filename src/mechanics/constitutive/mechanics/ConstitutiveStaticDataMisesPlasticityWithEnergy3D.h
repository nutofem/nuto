// $Id: ConstitutiveStaticDataMisesPlasticityWithEnergy3D.h 87 2009-11-06 10:35:39Z unger3 $

#ifndef CONSTITUTIVESTATICDATAMISESPLASTICITYWITHENERGY3D_H
#define CONSTITUTIVESTATICDATAMISESPLASTICITYWITHENERGY3D_H

#include "mechanics/constitutive/mechanics/ConstitutiveStaticDataMisesPlasticity3D.h"
#include "mechanics/constitutive/mechanics/ConstitutiveStaticDataPrevEngineeringStressStrain3D.h"

//! @brief ... class, storing the static data for the Mises plasticity3D including energy updates
//! this is only required, if the energy should be calculated, otherwise the base class ConstitutiveStaticDataMisesPlasticity3D
//! is sufficient
//! @author Jörg F. Unger, ISM
//! @date December 2009
namespace NuTo
{
class MisesPlasticityEngineeringStress;

class ConstitutiveStaticDataMisesPlasticityWithEnergy3D :
	public ConstitutiveStaticDataMisesPlasticity ,
    public ConstitutiveStaticDataPrevEngineeringStressStrain3D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class MisesPlasticityEngineeringStress;
public:
	//! @brief constructor
    ConstitutiveStaticDataMisesPlasticityWithEnergy3D();

    //! @brief copy constructor
    ConstitutiveStaticDataMisesPlasticityWithEnergy3D(ConstitutiveStaticDataMisesPlasticityWithEnergy3D const& rOther)
    {
        (*this) = rOther;
    }

    //! @brief clones (copies) the data
    virtual ConstitutiveStaticDataMisesPlasticityWithEnergy3D* Clone()const
    {
    	return new ConstitutiveStaticDataMisesPlasticityWithEnergy3D(*this);
    }

    //! @brief assignment operator
    ConstitutiveStaticDataMisesPlasticityWithEnergy3D& operator= (ConstitutiveStaticDataMisesPlasticityWithEnergy3D const& rOther);

    //! @brief check, if the static data is compatible with a given element and a given constitutive model
    virtual bool CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveStaticDataMisesPlasticityWithEnergy3D)
namespace boost{
template<>
struct is_virtual_base_of<NuTo::ConstitutiveStaticDataMisesPlasticity, NuTo::ConstitutiveStaticDataMisesPlasticityWithEnergy3D>: public mpl::true_ {};
template<>
struct is_virtual_base_of<NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D, NuTo::ConstitutiveStaticDataMisesPlasticityWithEnergy3D>: public mpl::true_ {};
}
#endif // ENABLE_SERIALIZATION

#endif // CONSTITUTIVESTATICDATAMISESPLASTICITYWITHENERGY3D_H
