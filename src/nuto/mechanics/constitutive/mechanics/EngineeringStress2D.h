// $Id$

#ifndef ENGINEERINGSTRESS2D_H_
#define ENGINEERINGSTRESS2D_H_

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
class EngineeringStress2D
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
    EngineeringStress2D();

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
     *  The components of the Engineering stress tensor are stored in vector notation: \f$ \left[\sigma_{xx}, \sigma_{yy}, \sigma_{xy}\right]. \f$
     */
    double mEngineeringStress[3];
};

}

#endif // ENGINEERINGSTRESS2D_H_
