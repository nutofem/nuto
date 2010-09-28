// $Id$

#ifndef SECONDPIOLAKIRCHHOFFSTRESS1D_H_
#define SECONDPIOLAKIRCHHOFFSTRESS1D_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

namespace NuTo
{
class LinearElastic;

//! @brief ... one-dimensional second Piola-Kirchhoff stress
/*!
 *  In the one-dimensional case, the second Piola-Kirchhoff stress tensor \f$\boldsymbol{S}\f$ reads
 *  \f[
 *    \boldsymbol{S} = \begin{bmatrix}
 *      S_{xx}
 *    \end{bmatrix}.
 *  \f]
 */
//! @author Stefan Eckardt, ISM
//! @date November 2009
class SecondPiolaKirchhoffStress1D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class NuTo::LinearElastic;
public:
    //! @brief ... constructor
    SecondPiolaKirchhoffStress1D();

    //! @brief ... get number of components stored in this object
    //! @return ... number of components stored in this object
    unsigned int GetNumberOfComponents() const;

    //! @brief ... get the components of the second Piola-Kirchhoff stress tensor
    //! @return ... components of the second Piola-Kirchhoff stress tensor (stored in an array)
    //! @sa mCauchyStress
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
    //! @brief ... components of the second Piola-Kirchhoff stress tensor
    /*!
     *  The components of the second Piola-Kirchhoff stress tensor are stored in vector notation: \f$ \left[S_{xx}\right]. \f$
     */
    double mSecondPiolaKirchhoffStress;
};

}

#endif // SECONDPIOLAKIRCHHOFFSTRESS1D_H_
