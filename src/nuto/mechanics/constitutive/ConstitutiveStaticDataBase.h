// $Id$

#ifndef CONSTITUTIVESTATICDATABASE_H_
#define CONSTITUTIVESTATICDATABASE_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif // ENABLE_SERIALIZATION

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION


//! @brief ... base class, storing the static data (history variables) of a constitutive relationship
//! @author Stefan Eckardt, ISM
//! @date November 2009
namespace NuTo
{
class ConstitutiveStaticDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //!@ brief constructor
    ConstitutiveStaticDataBase()
    {}

    //!@ brief destructor (virtual, in order to make the class a polymorphic type)
    virtual ~ConstitutiveStaticDataBase()
    {};

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
};

}

#endif // CONSTITUTIVESTATICDATABASE_H_ 
