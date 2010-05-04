// $Id$

#ifndef ENGINEERINGSTRESS3D_H_
#define ENGINEERINGSTRESS3D_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

namespace NuTo
{
class LinearElastic;
class ConstitutiveMisesPlasticity;
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
class EngineeringStress3D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class LinearElastic;
    friend class ConstitutiveMisesPlasticity;
    friend class ConstitutiveEngineeringStressStrain;
    friend class NonlocalDamagePlasticity;
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

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    void Info(unsigned short rVerboseLevel) const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(this->mEngineeringStress);
    }
#endif // ENABLE_SERIALIZATION

private:
    //! @brief ... components of the Engineering stress tensor
    /*!
     *  The components of the Engineering stress tensor are stored in vector notation: \f$ \left[\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sigma_{xy}, \sigma_{yz}, \sigma_{zx} \right]. \f$
     */
    double mEngineeringStress[6];
};

}

#endif // ENGINEERINGSTRESS3D_H_
