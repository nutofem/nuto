// $Id: $

#ifndef LATTICESTRAINS3D_H_
#define LATTICESTRAINS3D_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

namespace NuTo
{
//! @brief ... three-dimensional lattice strain
/*!
 * which includes the normal and two shear part
 */
//! @author Joerg F. Unger, NU
//! @date Jan 2012
class Lattice2D;
class ConstitutiveLatticeConcrete;
class ConstitutiveLatticeStressStrain;
class LatticeStrain3D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class NuTo::Lattice2D;
    friend class NuTo::ConstitutiveLatticeStressStrain;
    friend class NuTo::ConstitutiveLatticeConcrete;
public:
    //! @brief ... constructor
    LatticeStrain3D();

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
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
private:
    //! @brief ... array of interface tractions
    /*!
     *  The components of the interface traction vector are stored as follows: \f$ \left[ T_N, T_{T1}, T_{T2} \right] \f$
     */
    double mLatticeStrain[3];
};
}

#endif // LATTICESTRAINS3D_H_
