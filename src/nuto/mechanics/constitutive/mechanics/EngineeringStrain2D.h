// $Id: EngineeringStrain2D.h $

#ifndef ENGINEERINGSTRAIN2D_H
#define ENGINEERINGSTRAIN2D_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/math/FullVector.h"
#include "nuto/mechanics/constitutive/ConstitutiveInputBase.h"
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
class EngineeringStrain2D: public ConstitutiveOutputBase, public ConstitutiveInputBase, public FullVector<double,3>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class ConstitutiveEngineeringStressStrain;
    friend class ConstraintLinearGlobalCrackAngle;
    friend class ConstraintLinearDisplacementsPeriodic2D;
    friend class ConstraintLinearFineScaleDisplacementsPeriodic2D;
    friend class ConstraintLinearGlobalTotalStrain;
    friend class ConstitutiveMisesPlasticity;
    friend class ConstraintNonlinearGlobalCrackAngle2D;
    friend class ConstitutiveStaticDataMultiscale2DPlaneStrain;
    friend class ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain;
    friend class GradientDamagePlasticityEngineeringStress;
    friend class LinearElastic;
    friend class LinearElasticEngineeringStress;
    friend class Multiscale;
    friend class NonlocalDamagePlasticityEngineeringStress;
    friend class StructureMultiscale;
    friend class DeformationGradient1D;
    friend class DeformationGradient2D;
    friend class DeformationGradient3D;
    public:
    //! @brief ... constructor
    //! @param pStructure ... structure
    //! @param pElement ... element
    //! @param pIntegrationPoint ... integration point
    EngineeringStrain2D();

    //! @brief ... copy constructor
    EngineeringStrain2D(const DeformationGradient2D& rDeformationGradient);

    //! @brief ... copy constructor
    EngineeringStrain2D(const EngineeringStrain2D& rEngineeringStrain2D);

    //! @brief ... get number of strain components
    //! @return ... number of strain components
    unsigned int GetNumberOfComponents() const;

    //! @brief ... get Engineering Strain
    //! @return ... Engineering Strain (exx,eyy,gxy)
    //! @sa mDeformationGradient
    const double* GetData() const;

    //! @brief ... assignment constructor
	//! @param  rOther ... copied element
    template<typename OtherDerived>
    EngineeringStrain2D& operator=( const Eigen::MatrixBase <OtherDerived>& other)
	{
    	this->FullVector<double,3>::operator=(other);
    	return *this;
	}

    //! @brief ... set Engineering Strain
    //! @return ... Engineering Strain (exx,eyy,gxy)
    //! @sa mDeformationGradient
    void SetData(double rEngineeringStrain[3]);

    NuTo::EngineeringStrain2D& GetEngineeringStrain2D()
    {
    	return *this;
    }

    const NuTo::EngineeringStrain2D& GetEngineeringStrain2D()const
    {
    	return *this;
    }


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
    //! (exx,eyy,gxy)
    //double mEngineeringStrain[3];
};

}

#endif // ENGINEERINGSTRAIN2D_H
