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

#include "nuto/mechanics/sections/SectionBase.h"
#include <vector>

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
    virtual eSectionType GetType() const;

    //! @brief ... print information about the section
    //! @param rVerboseLevel ... verbosity of the information
    virtual void Info(unsigned short rVerboseLevel) const;

    //! @brief ... cast the base pointer to an SectionTruss, otherwise throws an exception
    SectionTruss* AsSectionTruss() override;

    //! @brief ... cast the base pointer to an SectionTruss, otherwise throws an exception
    const SectionTruss* AsSectionTruss() const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

#endif // ENABLE_SERIALIZATION
private:
    //! @brief ... cross-section area
    double mArea;
    std::vector<double> mAreaParameters;
};

} // namespace

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::SectionTruss)
#endif // ENABLE_SERIALIZATION



