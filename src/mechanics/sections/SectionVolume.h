// $Id$

#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "mechanics/sections/SectionBase.h"

namespace NuTo
{
//! @author Stefan Eckardt, ISM
//! @date October 2009
//! @brief ... general three-dimensional section
class SectionVolume : public SectionBase
{
#ifdef ENABLE_SERIALIZATION
   friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief ... constructor
   SectionVolume();

    //! @brief ... get the section type
    //! @return ... section type
    virtual eSectionType GetType() const;

    //! @brief ... print information about the section
    //! @param rVerboseLevel ... verbosity of the information
    virtual void Info(unsigned short rVerboseLevel) const;



#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize SectionVolume" << std::endl;
#endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SectionBase);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize SectionVolume" << std::endl;
#endif
    }
#endif // ENABLE_SERIALIZATION

private:
};

}// namespace NuTo


#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::SectionVolume)
#endif

