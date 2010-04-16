// $Id$

#ifndef SECONDPIOLAKIRCHHOFFSTRESS3D_H_
#define SECONDPIOLAKIRCHHOFFSTRESS3D_H_

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

//! @brief ... three-dimensional second Piola-Kirchhoff stress
/*!
 *  In the three-dimensional case, the second Piola-Kirchhoff stress tensor \f$\boldsymbol{S}\f$ reads
 *  \f[
 *    \boldsymbol{S} = \begin{bmatrix}
 *      S_{xx} & S_{xy} & S_{xz}\\
 *      S_{yx} & S_{yy} & S_{yz}\\
 *      S_{zx} & S_{zy} & S_{zz}
 *    \end{bmatrix}.
 *  \f]
 *  Due to the symmetry (\f$ S_{xy} = S_{yx}, S_{yz} = S_{zy}, S_{xz} = S_{zx} \f$) of the second Piola-Kirchhoff stress tensor, the components can be stored in a vector
 *  \f[
 *     \boldsymbol{S} = \begin{bmatrix}
 *        S_{xx}\\
 *        S_{yy}\\
 *        S_{zz}\\
 *        S_{xy}\\
 *        S_{yz}\\
 *        S_{zx}
 *    \end{bmatrix}.
 *  \f]
 */
//! @author Stefan Eckardt, ISM
//! @date November 2009
class SecondPiolaKirchhoffStress3D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class NuTo::LinearElastic;
public:
    //! @brief ... constructor
    SecondPiolaKirchhoffStress3D();

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
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(this->mSecondPiolaKirchhoffStress);
    }
#endif // ENABLE_SERIALIZATION

private:
    //! @brief ... components of the second Piola-Kirchhoff stress tensor
    /*!
     *  The components of the second Piola-Kirchhoff stress tensor are stored in vector notation: \f$ \left[S_{xx}, S_{yy}, S_{zz}, S_{xy}, S_{yz}, S_{zx}\right]. \f$
     */
    double mSecondPiolaKirchhoffStress[6];
};

}

#endif // SECONDPIOLAKIRCHHOFFSTRESS3D_H_
