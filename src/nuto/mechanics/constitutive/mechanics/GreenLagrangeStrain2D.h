// $Id: EngineeringStrain2D.h $

#ifndef GREENSTRAIN2D_H
#define GREENSTRAIN2D_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

namespace NuTo
{
class DeformationGradient1D;
class DeformationGradient2D;
class LinearElastic;
//! @brief ... three-dimensional deformation gradient
//! @author Jörg F. Unger, ISM
//! @date November 2009
class GreenLagrangeStrain2D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class LinearElastic;
public:
    //! @brief ... constructor
    //! @param pStructure ... structure
    //! @param pElement ... element
    //! @param pIntegrationPoint ... integration point
    GreenLagrangeStrain2D();

    //! @brief ... copy constructor
    GreenLagrangeStrain2D(const DeformationGradient2D& rDeformationGradient);

    //! @brief ... get number of strain components
    //! @return ... number of strain components
    unsigned int GetNumberOfComponents() const;

    //! @brief ... get green Lagrange Strain
    //! @return ... (exx,eyy,gxy)
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
    //! @brief ... green strain
    //! The green strain is stored as :
    //! (exx,eyy,gxy)
    double mGreenLagrangeStrain[3];
};

}

#endif // GREENSTRAIN2D_H
