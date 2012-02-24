// $Id: $

#ifndef ConstitutiveStaticDataLatticeConcrete2D_H
#define ConstitutiveStaticDataLatticeConcrete2D_H

#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"

//! @brief ... base class, storing the static data (history variables) of a constitutive relationship
//! @author Jörg F. Unger, ISM
//! @date December 2009
namespace NuTo
{
class ConstitutiveMisesPlasticity;
class IpDataStaticDataBase;
class Lattice2D;

class ConstitutiveStaticDataLatticeConcrete2D : virtual public ConstitutiveStaticDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
    friend class ConstitutiveLatticeConcrete;
    friend class Lattice2D;
#endif // ENABLE_SERIALIZATION
public:
	//! @brief constructor
	ConstitutiveStaticDataLatticeConcrete2D();

    //! @brief copy constructor
	ConstitutiveStaticDataLatticeConcrete2D(ConstitutiveStaticDataLatticeConcrete2D const& rOther)
    {
        (*this) = rOther;
    }

    //! @brief clones (copies) the data
    virtual ConstitutiveStaticDataLatticeConcrete2D* Clone()const
    {
    	return new ConstitutiveStaticDataLatticeConcrete2D(*this);
    }

    //! @brief check, if the static data is compatible with a given element and a given constitutive model
    virtual bool CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const;

    //! @brief assignment operator
    ConstitutiveStaticDataLatticeConcrete2D& operator= (ConstitutiveStaticDataLatticeConcrete2D const& rOther);

    //!@ brief reinterpret as lattice concrete 2D static data
    ConstitutiveStaticDataLatticeConcrete2D* AsConstitutiveStaticDataLatticeConcrete2D();

    //!@ brief reinterpret as lattice concrete 2D static data
    const ConstitutiveStaticDataLatticeConcrete2D* AsConstitutiveStaticDataLatticeConcrete2D()const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief accumulated plastic strain (is not always equivalent to epsilon_p)
    double mEpsilonMax;
};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveStaticDataLatticeConcrete2D)
#endif // ENABLE_SERIALIZATION

#endif // ConstitutiveStaticDataLatticeConcrete2D_H
