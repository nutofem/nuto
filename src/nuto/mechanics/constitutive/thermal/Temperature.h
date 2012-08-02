// $Id$

#ifndef TEMPERATURE_H_
#define TEMPERATURE_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

namespace NuTo
{
class LinearElastic;
//! @brief ... temperature
//! @author Jörg F. Unger, BAM
//! @date July 2012
class Temperature
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class LinearElastic;
public:
    //! @brief ... constructor
    Temperature();

    //! @brief ... constructor
    virtual ~Temperature(){};

    //! @brief ... get the temperature
    //! @sa mTemperature
    double GetTemperature() const;

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
    double mTemperature;
};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::Temperature)
#endif // ENABLE_SERIALIZATION

#endif // TEMPERATURE_H_
