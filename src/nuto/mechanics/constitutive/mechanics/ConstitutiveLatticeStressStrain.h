// $Id:$

#ifndef ConstitutiveLatticeStressStrain_H_
#define ConstitutiveLatticeStressStrain_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"

namespace NuTo
{
// forward declarations
class ElementBase;
class ConstitutiveTangentBase;
class LatticeStrain2D;
class LatticeStrain3D;
class LatticeStress2D;
class LatticeStress3D;

//! @brief ... base class for the mechanical constitutive relationship using Lattice stress and strains
//! @author Jörg F. Unger, ISM
//! @date December 2009
class ConstitutiveLatticeStressStrain : public virtual ConstitutiveBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:

    //! @brief ... constructor
    ConstitutiveLatticeStressStrain();

    //  Lattice strain /////////////////////////////////////
    //! @brief ... convert an Lattice strain from 2D into 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain2D ... Lattice strain 2D (input)
    //! @param rLatticeStrain3D ... Lattice strain 3D (output)
    virtual Error::eError GetLatticeStrain(const ElementBase* rElement, int rIp,
                                      const LatticeStrain2D& rLatticeStrain2D, LatticeStrain3D& rLatticeStrain3D) const;


    //  Lattice strain /////////////////////////////////////
    //! @brief ... calculate Lattice plastic strain from lattice strain in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... lattice strain
    //! @param rLatticePlasticStrain ... Lattice strain
    virtual Error::eError GetLatticePlasticStrain(const ElementBase* rElement, int rIp,
                                      const LatticeStrain2D& rLatticeStrain, LatticeStrain3D& rLatticePlasticStrain) const=0;

    //  Lattice strain /////////////////////////////////////
    //! @brief ... calculate Lattice plastic strain from lattice strain in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... lattice strain
    //! @param rLatticePlasticStrain ... plastic Lattice strain
    virtual Error::eError GetLatticePlasticStrain(const ElementBase* rElement, int rIp,
                                      const LatticeStrain3D& rLatticeStrain, LatticeStrain3D& rLatticePlasticStrain) const=0;

    // Lattice stress - Lattice strain /////////////////////////////////////
    //! @brief ... calculate Lattice stress from Lattice strain (which is calculated from the lattice strain)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... lattice strain
    //! @param rLatticeStress ... Lattice stress
    virtual Error::eError GetLatticeStressFromLatticeStrain(const ElementBase* rElement, int rIp,
            const LatticeStrain2D& rLatticeStrain, LatticeStress2D& rLatticeStress) const=0;

    // Lattice stress - Lattice strain /////////////////////////////////////
    //! @brief ... calculate Lattice stress from Lattice strain (which is calculated from the lattice strain)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... lattice strain
    //! @param rLatticeStress ... Lattice stress
    virtual Error::eError GetLatticeStressFromLatticeStrain(const ElementBase* rElement, int rIp,
            const LatticeStrain2D& rLatticeStrain, LatticeStress3D& rLatticeStress) const=0;

    // Lattice stress - Lattice strain /////////////////////////////////////
    //! @brief ... calculate Lattice stress from Lattice strain (which is calculated from the lattice strain)
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... lattice strain
    //! @param rCauchyStress ... Cauchy stress
    virtual Error::eError GetLatticeStressFromLatticeStrain(const ElementBase* rElement, int rIp,
            const LatticeStrain3D& rLatticeStrain, LatticeStress3D& rLatticeStress) const=0;

    //  Damage /////////////////////////////////////
    //! @brief ... calculate isotropic damage from lattice strain in 2D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... lattice strain
    //! @param rDamage ... damage variable
    virtual Error::eError GetDamage(const ElementBase* rElement, int rIp,
                                      const LatticeStrain2D& rLatticeStrain, double& rDamage) const=0;

    //  Damage /////////////////////////////////////
    //! @brief ... calculate isotropic damage from lattice strain in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... lattice strain
    //! @param rDamage ... damage variable
    virtual Error::eError GetDamage(const ElementBase* rElement, int rIp,
                                      const LatticeStrain3D& rLatticeStrain, double& rDamage) const=0;

    //! @brief ... calculate the tangent (derivative of the Lattice stresses with respect to the Lattice strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... lattice strain
    //! @param rTangent ... tangent
    virtual Error::eError GetTangent_LatticeStress_LatticeStrain(const ElementBase* rElement, int rIp,
            const LatticeStrain2D& rLatticeStrain,
            ConstitutiveTangentBase* rTangent) const=0;

    //! @brief ... calculate the tangent (derivative of the Lattice stresses with respect to the Lattice strains) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... lattice strain
    //! @param rTangent ... tangent
    virtual Error::eError GetTangent_LatticeStress_LatticeStrain(const ElementBase* rElement, int rIp,
            const LatticeStrain3D& rLatticeStrain,
            ConstitutiveTangentBase* rTangent) const=0;

    //! @brief ... update static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... lattice strain
    virtual Error::eError UpdateStaticData_LatticeStress_LatticeStrain(ElementBase* rElement, int rIp,
            const LatticeStrain2D& rLatticeStrain) const=0;

    //! @brief ... update static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... lattice strain
    virtual Error::eError UpdateStaticData_LatticeStress_LatticeStrain(ElementBase* rElement, int rIp,
            const LatticeStrain3D& rLatticeStrain) const=0;

    //! @brief ... update tmp static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... lattice strain
    virtual Error::eError UpdateTmpStaticData_LatticeStress_LatticeStrain(ElementBase* rElement, int rIp,
            const LatticeStrain2D& rLatticeStrain) const=0;

    //! @brief ... update tmp static data (history variables) of the constitutive relationship
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... lattice strain
    virtual Error::eError UpdateTmpStaticData_LatticeStress_LatticeStrain(ElementBase* rElement, int rIp,
            const LatticeStrain3D& rLatticeStrain) const=0;

    //! @brief ... checks, if a model has to be switched from linear to nonlinear, and then performs the adaption
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... lattice strain
    virtual Error::eError MultiscaleSwitchToNonlinear(ElementBase* rElement, int rIp,
            const LatticeStrain2D& rLatticeStrain)const;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    virtual ConstitutiveStaticDataBase* AllocateStaticDataLatticeStress_LatticeStrain2D(const ElementBase* rElement) const=0;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    virtual ConstitutiveStaticDataBase* AllocateStaticDataLatticeStress_LatticeStrain3D(const ElementBase* rElement) const=0;

    //! @brief ... calculate the total energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... lattice strain
    virtual Error::eError GetTotalEnergy_LatticeStress_LatticeStrain(const ElementBase* rElement, int rIp,
            const LatticeStrain2D& rLatticeStrain, double& energy) const;

    //! @brief ... calculate the total energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... lattice strain
    virtual Error::eError GetTotalEnergy_LatticeStress_LatticeStrain(const ElementBase* rElement, int rIp,
            const LatticeStrain3D& rLatticeStrain, double& energy) const;

    //! @brief ... calculate the elastic energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... lattice strain
    virtual Error::eError GetElasticEnergy_LatticeStress_LatticeStrain(const ElementBase* rElement, int rIp,
            const LatticeStrain2D& rLatticeStrain, double& energy) const;

    //! @brief ... calculate the elastic energy density
    //! @param rStructure ... structure
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rLatticeStrain ... lattice strain
    virtual Error::eError GetElasticEnergy_LatticeStress_LatticeStrain(const ElementBase* rElement, int rIp,
            const LatticeStrain3D& rLatticeStrain, double& energy) const;

    //! @brief ... avoid dynamic cast
    //! @return ... see brief explanation
    virtual const ConstitutiveLatticeStressStrain* AsConstitutiveLatticeStressStrain()const;

    //! @brief ... avoid dynamic cast
    //! @return ... see brief explanation
    virtual ConstitutiveLatticeStressStrain* AsConstitutiveLatticeStressStrain();

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
protected:
    //! @brief ... if true, the static data should be derived from StaticDataPrevStressStrain
    bool mEnergyFlag;

};

}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveLatticeStressStrain)
#endif // ENABLE_SERIALIZATION


#endif // ConstitutiveLatticeStressStrain_H_
