// $Id$

#ifndef HeatFlux3D_H
#define HeatFlux3D_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveOutputBase.h"

namespace NuTo
{

//! @brief ... heat flux
//! @author Jörg F. Unger, BAM
//! @date June 2012
class HeatFlux3D: public ConstitutiveOutputBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief ... constructor
    //! @param pStructure ... structure
    //! @param pElement ... element
    //! @param pIntegrationPoint ... integration point
    HeatFlux3D();

    //! @brief ... return heat flux
    HeatFlux3D& GetHeatFlux3D()
    {
    	return *this;
    }

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
    //! @brief ... heat flux
    double mHeatFlux[3];
};

}

#endif // HeatFlux3D_H
