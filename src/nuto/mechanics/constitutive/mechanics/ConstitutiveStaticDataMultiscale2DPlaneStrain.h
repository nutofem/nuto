// $Id: ConstitutiveStaticDataMultiscale3D.h 87 2009-11-06 10:35:39Z unger3 $

#ifndef CONSTITUTIVESTATICDATAMULTISCALE2DPLANESTRAIN_H
#define CONSTITUTIVESTATICDATAMULTISCALE2DPLANESTRAIN_H

#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain.h"

//! @brief ... base class, storing the static data (history variables) of a constitutive relationship
//! @author Jörg F. Unger, ISM
//! @date December 2009
namespace NuTo
{
class StructureMultiscale;
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
    StructureMultiscale* GetFineScaleStructure();

    //! @brief return structure
    const StructureMultiscale* GetFineScaleStructure()const;

    //! @brief return the previous hom strain
    const EngineeringStrain2D& GetPrevHomStrain()const;

    //! @brief set the previous hom strain
    void SetPrevHomStrain(EngineeringStrain2D rHomStrain);

    //! @brief sets the fine scale model (deserialization from a binary file)
    void SetFineScaleModel(std::string rFileName, double rMacroLength);

    //! @brief sets the fine scale parameters
    void SetFineScaleParameter(const std::string& rName, double rParameter);

    //! @brief sets the fine scale parameters
    void SetFineScaleParameter(const std::string& rName, std::string rParameter);

    //! @brief in case the fine scale model has not been initialized,
    //! an initial linear elastic model is used
    //! with this routine, the transition to the actual fine scale model is used
    //! with the initialization of the crack angle based on the previous elastic solution
    void UseNonlinearSolution();

    //return if the solution is either linear elastic or from the fine scale model
    bool NonlinearSolutionOn()const
    {
        return mNonlinearSolutionOn;
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief fine scale structure representative for a macroscopic integration point
    StructureMultiscale* mStructure;

    //! @brief previous homogeneous strain
    EngineeringStrain2D mPrevHomStrain;

    bool mNonlinearSolutionOn;
};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveStaticDataMultiscale2DPlaneStrain)
#endif // ENABLE_SERIALIZATION

#endif // CONSTITUTIVESTATICDATAMULTISCALE2DPLANESTRAIN_H
