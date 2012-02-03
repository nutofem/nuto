// $Id: InterfaceTractions2D.h 316 2010-09-28 19:40:50Z unger3 $

#ifndef LATTICESTRAINS2D_H_
#define LATTICESTRAINS2D_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

namespace NuTo
{
//! @brief ... two-dimensional lattice strain
/*!
 * which includes the volumetric, deviatoric and shear part
 */
//! @author Joerg F. Unger, NU
//! @date Jan 2012
class Lattice2D;
class ConstitutiveLatticeStressStrain;
class ConstitutiveLatticeConcrete;
class LatticeStrain2D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class NuTo::Lattice2D;
    friend class NuTo::ConstitutiveLatticeStressStrain;
    friend class NuTo::ConstitutiveLatticeConcrete;
public:
    //! @brief ... constructor
    LatticeStrain2D();

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
     *  The components of the interface traction vector are stored as follows: \f$ \left[ T_N, T_T \right] \f$
     */
    double mLatticeStrain[2];
};
}

#endif // LATTICESTRAINS2D_H_
