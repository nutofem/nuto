// $Id$

#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/sections/SectionBase.h"

namespace NuTo
{
//! @author Stefan Eckardt, ISM
//! @date November 2009
//! @brief ... section for two-dimensional elements
class SectionPlane: public NuTo::SectionBase
{
#ifdef ENABLE_SERIALIZATION
   friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief ... constructor
    SectionPlane(eSectionType rSectionType);

    //! @brief ... get the section thickness
    //! @return ... section thickness
    virtual double GetThickness() const;

    //! @brief ... set the thickness of the section
    //! @param rThickness ... section thickness
    virtual void SetThickness(double rThickness);

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
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

private:
    //! @brief just for serialization
    SectionPlane(){}
    //! @brief ... section thickness
    double        mThickness;
    eSectionType   mSectionType;
};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::SectionPlane)
#endif // ENABLE_SERIALIZATION

