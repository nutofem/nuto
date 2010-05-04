// $Id$

#ifndef SECTIONBASE_H
#define SECTIONBASE_H

#include <iostream>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION
#include <boost/noncopyable.hpp>

#include "nuto/mechanics/sections/SectionEnum.h"

namespace NuTo
{
//! @author Stefan Eckardt, ISM
//! @date October 2009
//! @brief ... standard abstract base class for sections
class SectionBase : boost::noncopyable
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief ... get the section type
    //! @return ... section type
    virtual Section::eSectionType GetType() const = 0;

    //! @brief ... get the cross-section area of the section
    //! @return ... section cross-section area
    virtual double GetArea() const;

    //! @brief ... set the cross-section area of the section
    //! @param rArea ... cross-section area
    virtual void SetArea(double rArea);

    //! @brief ... get the thickness of the section
    //! @return ... section thickness
    virtual double GetThickness() const;

    //! @brief ... get the thickness of the section
    //! @return ... section thickness
    virtual void SetThickness(double rThickness);

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
    }
#endif // ENABLE_SERIALIZATION
private:
};

}

#endif // SECTIONBASE_H
