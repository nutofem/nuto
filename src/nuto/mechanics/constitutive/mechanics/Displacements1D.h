// $Id$

#ifndef TEMPERATURE_H_
#define TEMPERATURE_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveInputBase.h"

namespace NuTo
{
class LinearElastic;
//! @brief ... temperature
//! @author Jörg F. Unger, BAM
//! @date July 2012
class Displacements1D: public ConstitutiveInputBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class LinearElastic;
public:
    //! @brief ... constructor
    Displacements1D();

    //! @brief ... constructor
    virtual ~Displacements1D(){};

    //! @brief ... get the values
    const double* GetData() const;

    //! @brief ... get the values
    double* GetData();

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    void Info(unsigned short rVerboseLevel) const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

private:
    double mDisplacements;
};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::Displacements1D)
#endif // ENABLE_SERIALIZATION

#endif // TEMPERATURE_H_
