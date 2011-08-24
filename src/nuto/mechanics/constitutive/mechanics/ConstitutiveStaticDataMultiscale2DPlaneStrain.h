// $Id: ConstitutiveStaticDataMultiscale3D.h 87 2009-11-06 10:35:39Z unger3 $

#ifndef CONSTITUTIVESTATICDATAMULTISCALE2DPLANESTRAIN_H
#define CONSTITUTIVESTATICDATAMULTISCALE2DPLANESTRAIN_H

#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain.h"

//! @brief ... base class, storing the static data (history variables) of a constitutive relationship
//! @author Jörg F. Unger, ISM
//! @date December 2009
namespace NuTo
{
class StructureMultiscale;
class IpDataStaticDataBase;
template<class T> class  FullMatrix;
class ConstitutiveStaticDataMultiscale2DPlaneStrain : public ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class Multiscale;
public:
	//! @brief constructor
    ConstitutiveStaticDataMultiscale2DPlaneStrain();

    //! @brief copy constructor
    ConstitutiveStaticDataMultiscale2DPlaneStrain(ConstitutiveStaticDataMultiscale2DPlaneStrain const& rOther)
    {
        (*this) = rOther;
    }

    //! @brief clones (copies) the data
    virtual ConstitutiveStaticDataMultiscale2DPlaneStrain* Clone()const
    {
    	return new ConstitutiveStaticDataMultiscale2DPlaneStrain(*this);
    }

    //! @brief assignment operator
    ConstitutiveStaticDataMultiscale2DPlaneStrain& operator= (ConstitutiveStaticDataMultiscale2DPlaneStrain const& rOther);

    //! @brief check, if the static data is compatible with a given element and a given constitutive model
    virtual bool CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType)const;

    //! @brief reinterpret as nonlocal damage2d static data
    ConstitutiveStaticDataMultiscale2DPlaneStrain* AsMultiscale2DPlaneStrain();

    //! @brief reinterpret as nonlocal damage2d static data
    const ConstitutiveStaticDataMultiscale2DPlaneStrain* AsMultiscale2DPlaneStrain()const;


    //! @brief return structure
    StructureMultiscale* GetFineScaleStructure();

    //! @brief return structure
    const StructureMultiscale* GetFineScaleStructure()const;

/*    //! @brief return the previous hom strain
    const EngineeringStrain2D& GetPrevHomStrain()const;

    //! @brief set the previous hom strain
    void SetPrevHomStrain(EngineeringStrain2D rHomStrain);
*/
    //! @brief sets the fine scale model (deserialization from a binary file)
    void SetFineScaleModel(std::string rFileName, double rMacroLength, double rCenter[2], std::string rIpName);

    //! @brief sets the fine scale parameters
    void SetFineScaleParameter(const std::string& rName, double rParameter);

    //! @brief sets the fine scale parameters
    void SetFineScaleParameter(const std::string& rName, std::string rParameter);

#ifdef ENABLE_VISUALIZE
    //! @brief Visualize for all integration points the fine scale structure
    //! either visualize the damage zone (rVisualizeDamage=true) or the homogogeneous zone (rVisualizeDamage=false)
    void VisualizeIpMultiscale(VisualizeUnstructuredGrid& rVisualize,
    		const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat, bool rVisualizeDamage)const;
#endif

    //! @brief in case the fine scale model has not been initialized,
    //! an initial linear elastic model is used
    //! with this routine, the transition to the actual fine scale model is used
    //! with the initialization of the crack angle based on the previous elastic solution
    //void UseNonlinearSolution();

    //return if the solution is linear elastic with the homogenized stiffness
    bool LinearSolution()const
    {
        return (mSolutionPhase==Constitutive::HOMOGENIZED_LINEAR_ELASTIC);
    }

    //return if the solution is nonlinear without a crack enrichment
    bool NonLinearSolutionNoCrack()const
    {
    	return (mSolutionPhase==Constitutive::NONLINEAR_NO_CRACK);
    }

    //return if the solution is nonlinear with a crack enrichment
    bool NonLinearSolutionCracked()const
    {
    	return (mSolutionPhase==Constitutive::NONLINEAR_CRACKED);
    }

    //set the type of the solution phase
    void SetSolutionPhase(Constitutive::eSolutionPhaseType rSolutionPhase)
    {
    	mSolutionPhase = rSolutionPhase;
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

    Constitutive::eSolutionPhaseType mSolutionPhase;
};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveStaticDataMultiscale2DPlaneStrain)
#endif // ENABLE_SERIALIZATION

#endif // CONSTITUTIVESTATICDATAMULTISCALE2DPLANESTRAIN_H
