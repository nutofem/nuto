// $Id$

#ifndef ENGINEERINGSTRESS3D_H_
#define ENGINEERINGSTRESS3D_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/math/FullVector.h"
#include "nuto/mechanics/constitutive/ConstitutiveOutputBase.h"
namespace NuTo
{
class LinearElastic;
class LinearElasticEngineeringStress;
class MisesPlasticityEngineeringStress;
class ConstitutiveEngineeringStressStrain;

//! @brief ... three-dimensional Engineering stress
/*!
 *  In the three-dimensional case, the Engineering stress tensor \f$\boldsymbol{\sigma}\f$ reads
 *  \f[
 *    \boldsymbol{\sigma} = \begin{bmatrix}
 *      \sigma_{xx} & \sigma_{xy} & \sigma_{xz}\\
 *      \sigma_{yx} & \sigma_{yy} & \sigma_{yz}\\
 *      \sigma_{zx} & \sigma_{zy} & \sigma_{zz}
 *    \end{bmatrix}.
 *  \f]
 *  Due to the symmetry (\f$ \sigma_{xy} = \sigma_{yx}, \sigma_{xz} = \sigma_{zx}, \sigma_{zy} = \sigma_{yz}\f$) of the Engineering stress tensor, the components can be stored in a vector
 *  \f[
 *     \boldsymbol{\sigma} = \begin{bmatrix}
 *        \sigma_{xx}\\
 *        \sigma_{yy}\\
 *        \sigma_{zz}\\
 *        \sigma_{xy}\\
 *        \sigma_{yz}\\
 *        \sigma_{zx}
 *    \end{bmatrix}.
 *  \f]
 */
//! @author Stefan Eckardt, ISM
//! @date November 2009
class EngineeringStress3D: public ConstitutiveOutputBase, public FullVector<double,6>
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class LinearElastic;
    friend class LinearElasticEngineeringStress;
    friend class MisesPlasticityEngineeringStress;
    friend class ConstitutiveEngineeringStressStrain;
    friend class GradientDamagePlasticityEngineeringStress;
    friend class NonlocalDamagePlasticityEngineeringStress;
    friend class Multiscale;
public:
    //! @brief ... constructor
    EngineeringStress3D();

    //! @brief ... get number of components stored in this object
    //! @return ... number of components stored in this object
    unsigned int GetNumberOfComponents() const;

    //! @brief ... get the components of the Engineering stress tensor
    //! @return ... components of the Engineering stress tensor (stored in an array)
    //! @sa mEngineeringStress
    const double* GetData() const;

    //! @brief ... sets the components of the Engineering stress tensor
    //! @param ... components of the Engineering stress tensor (stored in an array)
    //! @sa mEngineeringStress
    void SetData(double rData[6]);

    EngineeringStress3D& GetEngineeringStress3D()override
    {
    	return *this;
    }

    //! @brief ... assignment constructor
	//! @param  rOther ... copied element
    template<typename OtherDerived>
    EngineeringStress3D& operator=( const Eigen::MatrixBase <OtherDerived>& other)
	{
    	this->FullVector<double,6>::operator=(other);
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

#endif // ENGINEERINGSTRESS3D_H_
