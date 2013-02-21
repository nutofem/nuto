// $Id: EngineeringStrain1D.h $

#ifndef ENGINEERINGSTRAIN1D_H
#define ENGINEERINGSTRAIN1D_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveOutputBase.h"

namespace NuTo
{
class DeformationGradient1D;
class DeformationGradient2D;
class DeformationGradient3D;
class LinearElastic;
class ConstitutiveMisesPlasticity;

//! @brief ... three-dimensional deformation gradient
//! @author Jörg F. Unger, ISM
//! @date November 2009
class EngineeringStrain1D: public ConstitutiveOutputBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class LinearElastic;
    friend class LinearElasticEngineeringStress;
    friend class ConstitutiveMisesPlasticity;
    friend class ConstitutiveEngineeringStressStrain;
    friend class GradientDamagePlasticityEngineeringStress;
    friend class NonlocalDamagePlasticity;
    friend class Multiscale;
    friend class DeformationGradient1D;
    friend class DeformationGradient2D;
    friend class DeformationGradient3D;
    public:
    //! @brief ... constructor
    //! @param pStructure ... structure
    //! @param pElement ... element
    //! @param pIntegrationPoint ... integration point
    EngineeringStrain1D();

    //! @brief copy constructor
    EngineeringStrain1D(const DeformationGradient1D& rDeformationGradient);

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
    //! @brief ... engineering strain
    //! The engineering strain is stored as :
    //! (exx)
    double mEngineeringStrain;
};

}

#endif // ENGINEERINGSTRAIN1D_H
