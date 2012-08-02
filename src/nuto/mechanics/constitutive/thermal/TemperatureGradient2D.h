// $Id$

#ifndef TemperatureGradient2D_H
#define TemperatureGradient2D_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveInputBase.h"

namespace NuTo
{
//! @brief ... temperature gradient
//! @author Jörg F. Unger, BAM
//! @date June 2012
class TemperatureGradient2D : public ConstitutiveInputBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class LinearElastic;
    friend class Solid;
public:
    //! @brief ... constructor
    //! @param pStructure ... structure
    //! @param pElement ... element
    //! @param pIntegrationPoint ... integration point
    TemperatureGradient2D();

    //! @brief ... get number of strain components
    //! @return ... number of strain components
    unsigned int GetNumberOfComponents() const;

    //! @brief ... get Engineering Strain
    //! @return ... Engineering Strain (exx)
    //! @sa mDeformationGradient
    const double* GetData() const;


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
private:
    //! @brief ... temperature gradient
    double mTemperatureGradient[2];
};

}

#endif // TemperatureGradient2D_H
