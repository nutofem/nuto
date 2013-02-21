// $Id$

#ifndef ENGINEERINGSTRESS2D_H_
#define ENGINEERINGSTRESS2D_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveOutputBase.h"

namespace NuTo
{
class LinearElastic;
class ConstitutiveMisesPlasticity;
class ConstitutiveEngineeringStressStrain;

//! @brief ... two-dimensional Engineering stress
/*!
 *  In the two-dimensional case, the Engineering stress tensor \f$\boldsymbol{\sigma}\f$ reads
 *  \f[
 *    \boldsymbol{\sigma} = \begin{bmatrix}
 *      \sigma_{xx} & \sigma_{xy}\\
 *      \sigma_{yx} & \sigma_{yy}
 *    \end{bmatrix}.
 *  \f]
 *  Due to the symmetry (\f$ \sigma_{xy} = \sigma_{yx}\f$) of the Engineering stress tensor, the components can be stored in a vector
 *  \f[
 *     \boldsymbol{\sigma} = \begin{bmatrix}
 *        \sigma_{xx}\\
 *        \sigma_{yy}\\
 *        \sigma_{xy}
 *    \end{bmatrix}.
 *  \f]
 */
//! @author Stefan Eckardt, ISM
//! @date November 2009
class EngineeringStress2D: public ConstitutiveOutputBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class ConstitutiveMisesPlasticity;
    friend class ConstitutiveEngineeringStressStrain;
    friend class ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain;
    friend class GradientDamagePlasticityEngineeringStress;
    friend class LinearElastic;
    friend class LinearElasticEngineeringStress;
    friend class NonlocalDamagePlasticityEngineeringStress;
    friend class StructureMultiscale;
public:
    //! @brief ... constructor
    EngineeringStress2D();

    //! @brief ... copy constructor
    EngineeringStress2D(const EngineeringStress2D& rOther);

    //! @brief ... get number of components stored in this object
    //! @return ... number of components stored in this object
    unsigned int GetNumberOfComponents() const;

    EngineeringStress2D& GetEngineeringStress2D() override
    {
    	return *this;
    }

    //! @brief ... get the components of the Engineering stress tensor
    //! @return ... components of the Engineering stress tensor (stored in an array)
    //! @sa mEngineeringStress
    const double* GetData() const;

    EngineeringStress2D operator+ ( const EngineeringStress2D &other ) const
    {
        EngineeringStress2D result;
        result.mEngineeringStress[0] = mEngineeringStress[0] + other.mEngineeringStress[0];
        result.mEngineeringStress[1] = mEngineeringStress[1] + other.mEngineeringStress[1];
        result.mEngineeringStress[2] = mEngineeringStress[2] + other.mEngineeringStress[2];
        return result;
    }

    EngineeringStress2D operator- ( const EngineeringStress2D &other ) const
     {
         EngineeringStress2D result;
         result.mEngineeringStress[0] = mEngineeringStress[0] - other.mEngineeringStress[0];
         result.mEngineeringStress[1] = mEngineeringStress[1] - other.mEngineeringStress[1];
         result.mEngineeringStress[2] = mEngineeringStress[2] - other.mEngineeringStress[2];
         return result;
     }

    EngineeringStress2D operator* ( double scalar ) const
     {
         EngineeringStress2D result;
         result.mEngineeringStress[0] = scalar*mEngineeringStress[0];
         result.mEngineeringStress[1] = scalar*mEngineeringStress[1];
         result.mEngineeringStress[2] = scalar*mEngineeringStress[2];
         return result;
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
    //! @brief ... components of the Engineering stress tensor
    /*!
     *  The components of the Engineering stress tensor are stored in vector notation: \f$ \left[\sigma_{xx}, \sigma_{yy}, \sigma_{xy}\right]. \f$
     */
    double mEngineeringStress[3];
};

}

#endif // ENGINEERINGSTRESS2D_H_
