// $Id$

#ifndef ENGINEERINGSTRESS1D_H_
#define ENGINEERINGSTRESS1D_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/math/FullVector.h"
#include "nuto/mechanics/constitutive/ConstitutiveOutputBase.h"

namespace NuTo
{
class LinearElastic;
class ConstitutiveMisesPlasticity;
class ConstitutiveEngineeringStressStrain;
//! @brief ... one dimensional Engineering stress
/*!
 *  In the one-dimensional case, the Engineering stress tensor \f$\boldsymbol{\sigma}\f$ reads
 *  \f[
 *    \boldsymbol{\sigma} = \begin{bmatrix}
 *      \sigma_{xx}
 *    \end{bmatrix}.
 *  \f]
 */
//! @author Stefan Eckardt, ISM
//! @date November 2009
class EngineeringStress1D: public ConstitutiveOutputBase, public FullVector<double,1>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class ConstitutiveMisesPlasticity;
    friend class ConstitutiveEngineeringStressStrain;
    friend class GradientDamagePlasticityEngineeringStress;
    friend class LinearElastic;
    friend class LinearElasticEngineeringStress;
    friend class NonlocalDamagePlasticity;
public:
    //! @brief ... constructor
    EngineeringStress1D();

    //! @brief ... constructor
    virtual ~EngineeringStress1D(){};

    //! @brief ... get number of components stored in this object
    //! @return ... number of components stored in this object
    unsigned int GetNumberOfComponents() const;

    //! @brief ... get the components of the Engineering stress tensor
    //! @return ... components of the Engineering stress tensor (stored in an array)
    //! @sa mEngineeringStress
    const double* GetData() const;

    EngineeringStress1D& GetEngineeringStress1D()override
    {
    	return *this;
    }

    //! @brief ... assignment constructor
	//! @param  rOther ... copied element
    template<typename OtherDerived>
    EngineeringStress1D& operator=( const Eigen::MatrixBase <OtherDerived>& other)
	{
    	this->FullVector<double,1>::operator=(other);
    	return *this;
	}

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
};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::EngineeringStress1D)
#endif // ENABLE_SERIALIZATION

#endif // ENGINEERINGSTRESS1D_H_
