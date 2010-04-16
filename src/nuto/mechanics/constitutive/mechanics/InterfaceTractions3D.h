// $Id$

#ifndef INTERFACETRACTIONS3D_H_
#define INTERFACETRACTIONS3D_H_

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
//! @brief ... three-dimensional interface tractions
/*!
 *  In the three-dimensional case the vector of interface tractions reads
 *  \f[
 *      \boldsymbol{T} = \begin{bmatrix} T_N \\ T_{T1}\\ T_{T2} \end{bmatrix},
 *  \f]
 *  where \f$  T_N \f$ is the interface traction in normal direction,
 *  and \f$  T_{T1} \f$ and  \f$  T_{T2} \f$ are the interface tractions in tangential direction.
 */
//! @author Stefan Eckardt, ISM
//! @date November 2009
class InterfaceTractions3D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief ... constructor
    InterfaceTractions3D();

    //! @brief ... get number of components
    //! @return ... number of components
    unsigned int GetNumberOfComponents() const;

    //! @brief ... get interface tractions
    //! @return ... interface tractions
    //! @sa mInterfaceTractions
    const double* GetData() const;

    //! @brief ... print information about the strain object
    //! @param rVerboseLevel ... verbosity of the information
    void Info(unsigned short rVerboseLevel) const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(this->mInterfaceTractions);
    }
#endif // ENABLE_SERIALIZATION
private:
    //! @brief ... array of interface tractions
    /*!
     *  The components of the interface traction vector are stored as follows: \f$ \left[ T_N, T_{T1}, T_{T2} \right] \f$
     */
    double mInterfaceTractions[3];
};
}
#endif // INTERFACETRACTIONS3D_H_
