// $Id: EngineeringStrain1D.h $

#ifndef GREENSTRAIN1D_H
#define GREENSTRAIN1D_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

namespace NuTo
{
class DeformationGradient1D;
class DeformationGradient2D;
class DeformationGradient3D;
class LinearElastic;

//! @brief ... three-dimensional deformation gradient
//! @author Jörg F. Unger, ISM
//! @date November 2009
class GreenLagrangeStrain1D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class LinearElastic;
    friend class DeformationGradient1D;
    friend class DeformationGradient2D;
    friend class DeformationGradient3D;
public:
    //! @brief ... constructor
    //! @param pStructure ... structure
    //! @param pElement ... element
    //! @param pIntegrationPoint ... integration point
    GreenLagrangeStrain1D();

    //! @brief copy constructor
    GreenLagrangeStrain1D(const DeformationGradient1D& rDeformationGradient);

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
    //! @brief ... Green Lagrange strain
    //! The green strain is stored as :
    //! (exx)
    double mGreenLagrangeStrain;
};

}

#endif // GREENSTRAIN1D_H
