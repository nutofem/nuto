// $Id$

#ifndef ENGINEERINGSTRESS1D_H_
#define ENGINEERINGSTRESS1D_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

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
class EngineeringStress1D
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
     *  The components of the Engineering stress tensor are stored in vector notation: \f$ \left[\sigma_{xx}\right]. \f$
     */
    double mEngineeringStress;
};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::EngineeringStress1D)
#endif // ENABLE_SERIALIZATION

#endif // ENGINEERINGSTRESS1D_H_
