// $Id: ConstitutiveStaticDataNonlocalDamagePlasticity3D.h 87 2009-11-06 10:35:39Z unger3 $

#ifndef CONSTITUTIVESTATICDATANONLOCALDAMAGEPLASTICITY2DPLANESTRAIN_H
#define CONSTITUTIVESTATICDATANONLOCALDAMAGEPLASTICITY2DPLANESTRAIN_H

#include "mechanics/constitutive/mechanics/ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain.h"

//! @brief ... base class, storing the static data (history variables) of a constitutive relationship
//! @author Jörg F. Unger, ISM
//! @date December 2009
namespace NuTo
{
class IpDataStaticDataBase;

class ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain : public ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class NonlocalDamagePlasticityEngineeringStress;
public:
	//! @brief constructor
    ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain();

    //! @brief copy constructor
    ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain(ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain const& rOther)
    {
        (*this) = rOther;
    }

    //! @brief clones (copies) the data
    ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain* Clone()const;

    //! @brief assignment operator
    ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain& operator= (ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain const& rOther);

    //! @brief check, if the static data is compatible with a given element and a given constitutive model
    virtual bool CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const;

    //!@ brief reinterpret as nonlocal damage2d static data
    ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain* AsNonlocalDamagePlasticity2DPlaneStrain();

    //!@ brief reinterpret as nonlocal damage2d static data
    const ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain* AsNonlocalDamagePlasticity2DPlaneStrain()const;


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief local accumulated plastic strain scaled by nonlocal radius
    double mKappa;

    //! @brief plastic strain
    double mEpsilonP[4];

    //! @brief tmp static data derivative of local plastic strain with respect to local total strain (row wise storage 4x4)
    double mTmpdEpsilonPdEpsilon[16];

    //! @brief tmp static data for accumulated plastic strain scaled by equivalente length
    double mTmpKappa;

    //! @brief tmp static data plastic strain
    double mTmpEpsilonP[4];

    //! @brief tmp static data equivalent length
    double mTmpLeq;

    //! @brief tmp static data derivative of eq length with respect to local strain
    double mTmpdLeqdEpsilon[4];
};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain)
#endif // ENABLE_SERIALIZATION

#endif // CONSTITUTIVESTATICDATANONLOCALDAMAGEPLASTICITY3D_H
