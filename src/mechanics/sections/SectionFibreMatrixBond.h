// $Id$

#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "mechanics/sections/SectionBase.h"

namespace NuTo
{
//! @author Stefan Eckardt, ISM
//! @date November 2015
//! @brief ... section for the fibre-matrix bond element
class SectionFibreMatrixBond: public NuTo::SectionBase
{
#ifdef ENABLE_SERIALIZATION
   friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief ... constructor
   SectionFibreMatrixBond();

    //! @brief ... get the circumference of the fibre
    //! @return ... circumference of the fibre
    virtual double GetCircumference() const override;

    //! @brief ... set the tcircumference of the fibre
    //! @param rCircumference ... circumference of the fibre
    virtual void SetCircumference(double rCircumference) override;

    //! @brief ... get the section type
    //! @return ... section type
    virtual eSectionType GetType() const override;

    //! @brief ... print information about the section
    //! @param rVerboseLevel ... verbosity of the information
    virtual void Info(unsigned short rVerboseLevel) const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

private:

    //! @brief ... circumference of the fibre
    double        mCircumference;
};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::SectionFibreMatrixBond)
#endif // ENABLE_SERIALIZATION

