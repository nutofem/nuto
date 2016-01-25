// $Id: EngineeringStrain1D.h $

#ifndef ENGINEERINGSTRAIN1D_H
#define ENGINEERINGSTRAIN1D_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/math/FullVector.h"
#include "nuto/mechanics/constitutive/ConstitutiveInputBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveOutputBase.h"
#include "nuto/mechanics/constitutive/mechanics/LocalEqStrain.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"

#include <math.h>

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
class EngineeringStrain1D: public ConstitutiveOutputBase, public ConstitutiveInputBase, public FullVector<double,1>
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
    friend class StrainGradientDamagePlasticityEngineeringStress;
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

    //! @brief ... set Engineering Strain
    //! @return ... Engineering Strain (exx)
    //! @sa mDeformationGradient
    void SetData(const double rData[1]);

    //! @brief ... assignment constructor
	//! @param  rOther ... copied element
    template<typename OtherDerived>
    EngineeringStrain1D& operator=( const Eigen::MatrixBase <OtherDerived>& other)
	{
    	this->FullVector<double,1>::operator=(other);
    	return *this;
	}

    NuTo::EngineeringStrain1D& GetEngineeringStrain1D()
    {
    	return *this;
    }

    const NuTo::EngineeringStrain1D& GetEngineeringStrain1D()const
    {
    	return *this;
    }

    //! @brief ... calculates the norm of the strain tensor in 1D case = equivalent strain
    NuTo::LocalEqStrain Norm() const
    {
    	NuTo::LocalEqStrain localEqStrain;
    	localEqStrain[0] = std::abs((*this)[0]);
    	return localEqStrain;
    }

    //! @brief ... calculates the norm derivative of the strain tensor in 1D case = derivative of the equivalent strain by strain
    NuTo::ConstitutiveTangentLocal<1, 1> DNorm() const
    {
    	NuTo::ConstitutiveTangentLocal<1, 1> tangent;
    	NuTo::LocalEqStrain localEqStrain((*this).Norm());
    	double localEqStrainDerivative(1.);
        if (localEqStrain[0] < 0)
        {
            localEqStrainDerivative = -1.;
        }
        tangent(0, 0) = localEqStrainDerivative;
        return tangent;
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
    //! (exx)
    //double mEngineeringStrain;
};

}

#endif // ENGINEERINGSTRAIN1D_H
