// $Id$

#ifndef SECTIONTRUSS_H
#define SECTIONTRUSS_H

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/sections/SectionBase.h"

namespace NuTo
{
//! @author Stefan Eckardt, ISM
//! @date October 2009
//! @brief ... general one-dimensional section
class SectionTruss : public SectionBase
{
#ifdef ENABLE_SERIALIZATION
   friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief ... constructor
   SectionTruss();

    //! @brief ... get the cross-section area of the section
    //! @return ... section cross-section area
    virtual double GetArea() const;

    //! @brief ... calculates a x-dependent area factor based on the mAreaParameters data
    //! @param rXCoordinate ... x coordinate
    //! @return ... x-dependent area factor
    double GetAreaFactor(double rXCoordinate) const;

    //! @brief ... set the cross-section area of the section
    //! @param rArea ... cross-section area
    virtual void SetArea(double rArea);

    //! @brief ... set the area parameters
    //! @param rAreaParameter ... area parameters [0]-xWeakSpot [1]-lWeakSpot [2]-alpha [3]-exponent
    void SetAreaParameters(double* rAreaParameters);

    //! @brief ... get the section type
    //! @return ... section type
    virtual Section::eSectionType GetType() const;

    //! @brief ... print information about the section
    //! @param rVerboseLevel ... verbosity of the information
    virtual void Info(unsigned short rVerboseLevel) const;

    //! @brief ... cast the base pointer to an SectionTruss, otherwise throws an exception
    SectionTruss* AsSectionTruss();

    //! @brief ... cast the base pointer to an SectionTruss, otherwise throws an exception
    const SectionTruss* AsSectionTruss() const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SectionBase)
           & BOOST_SERIALIZATION_NVP(this->mArea)
           & BOOST_SERIALIZATION_NVP(mAreaParameters);
    }
#endif // ENABLE_SERIALIZATION
private:
    //! @brief ... cross-section area
    double mArea;
    double* mAreaParameters;
};

}

#endif // SECTIONTRUSS_H
