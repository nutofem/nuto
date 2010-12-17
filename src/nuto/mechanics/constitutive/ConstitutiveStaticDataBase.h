// $Id$

#ifndef CONSTITUTIVESTATICDATABASE_H_
#define CONSTITUTIVESTATICDATABASE_H_

#include <string>

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

//! @brief ... base class, storing the static data (history variables) of a constitutive relationship
//! @author Stefan Eckardt, ISM
//! @date November 2009
namespace NuTo
{
class ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain;
class ConstitutiveStaticDataMultiscale2DPlaneStrain;

class ConstitutiveStaticDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //!@ brief constructor
    ConstitutiveStaticDataBase()
    {}

    //!@ brief destructor (virtual, in order to make the class a polymorphic type)
    virtual ~ConstitutiveStaticDataBase()
    {};

    //! @brief sets the fine scale model (deserialization from a binary file)
    virtual void SetFineScaleModel(std::string rFileName);

    //!@ brief reinterpret as nonlocal damage2d static data
    virtual ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain* AsNonlocalDamagePlasticity2DPlaneStrain();

    //!@ brief reinterpret as nonlocal damage2d static data
    virtual const ConstitutiveStaticDataNonlocalDamagePlasticity2DPlaneStrain* AsNonlocalDamagePlasticity2DPlaneStrain()const;

    //!@ brief reinterpret as nonlocal damage2d static data
    virtual ConstitutiveStaticDataMultiscale2DPlaneStrain* AsMultiscale2DPlaneStrain();

    //!@ brief reinterpret as nonlocal damage2d static data
    virtual const ConstitutiveStaticDataMultiscale2DPlaneStrain* AsMultiscale2DPlaneStrain()const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
};

}

#endif // CONSTITUTIVESTATICDATABASE_H_ 
