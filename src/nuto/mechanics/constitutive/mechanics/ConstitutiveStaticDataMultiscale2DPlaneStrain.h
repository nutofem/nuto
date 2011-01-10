// $Id: ConstitutiveStaticDataMultiscale3D.h 87 2009-11-06 10:35:39Z unger3 $

#ifndef CONSTITUTIVESTATICDATAMULTISCALE2DPLANESTRAIN_H
#define CONSTITUTIVESTATICDATAMULTISCALE2DPLANESTRAIN_H

#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain.h"

//! @brief ... base class, storing the static data (history variables) of a constitutive relationship
//! @author Jörg F. Unger, ISM
//! @date December 2009
namespace NuTo
{
class StructureIp;
class ConstitutiveStaticDataMultiscale2DPlaneStrain : public ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class Multiscale;
public:
	//! @brief constructor
    ConstitutiveStaticDataMultiscale2DPlaneStrain();

    //! @brief reinterpret as nonlocal damage2d static data
    ConstitutiveStaticDataMultiscale2DPlaneStrain* AsMultiscale2DPlaneStrain();

    //! @brief reinterpret as nonlocal damage2d static data
    const ConstitutiveStaticDataMultiscale2DPlaneStrain* AsMultiscale2DPlaneStrain()const;

    //! @brief return structure
    StructureIp* GetFineScaleStructure();

    //! @brief return structure
    const StructureIp* GetFineScaleStructure()const;

    //! @brief sets the fine scale model (deserialization from a binary file)
    void SetFineScaleModel(std::string rFileName);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief fine scale structure repesentative for a macroscopic integration point
    StructureIp* mStructure;
};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveStaticDataMultiscale2DPlaneStrain)
#endif // ENABLE_SERIALIZATION

#endif // CONSTITUTIVESTATICDATAMULTISCALE2DPLANESTRAIN_H
