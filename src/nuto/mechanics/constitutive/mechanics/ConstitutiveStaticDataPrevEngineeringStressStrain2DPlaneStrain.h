// $Id: ConstitutiveStaticPrevEngineeringStressStrain2DPlaneStrain.h 87 2009-11-06 10:35:39Z unger3 $

#ifndef CONSTITUTIVESTATICDATAPrevEngineeringStressStrain2D_H
#define CONSTITUTIVESTATICDATAPrevEngineeringStressStrain2D_H

#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress2D.h"

//! @brief ... class, storing the static data for the Mises plasticity3D including energy updates
//! this is only required, if the energy should be calculated, otherwise the base class ConstitutiveStaticDataMisesPlasticity3D
//! is sufficient
//! @author Jörg F. Unger, ISM
//! @date December 2009
namespace NuTo
{
class ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain : virtual public ConstitutiveStaticDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
	//! @brief constructor
    ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain();

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

    //! brief set the previous stress
    void SetPrevStress(const EngineeringStress2D& rPrevSigma);

    //! brief returns the previous stress
    inline const EngineeringStress2D& GetPrevStress()const
    {
        return mPrevSigma;
    }

    //! brief set the previous stress
    void SetPrevStrain(const EngineeringStrain2D& rPrevStrain);

    //! brief returns the previous stress
    inline const EngineeringStrain2D& GetPrevStrain()const
    {
        return  mPrevStrain;
    }

    //! brief returns the previous total energy
    inline double GetPrevTotalEnergy()const
    {
        return  mPrevTotalEnergy;
    }

    //! brief set the previous total energy
    inline void SetPrevTotalEnergy(double rPrevTotalEnergy)
    {
    	mPrevTotalEnergy = rPrevTotalEnergy;
    }

    //! brief returns the previous elastic energy
    inline double GetPrevElasticEnergy()const
    {
        return  mPrevElasticEnergy;
    }

    //! brief set the previous elastic energy
    inline void SetPrevElasticEnergy(double rPrevElasticEnergy)
    {
    	mPrevElasticEnergy = rPrevElasticEnergy;
    }

protected:
    //! @brief previous stress
    EngineeringStress2D mPrevSigma;

    //! @brief previous strain
    EngineeringStrain2D mPrevStrain;

    //! @brief previous total energy
    double mPrevTotalEnergy;

    //! @brief previous elastic energy
    double mPrevElasticEnergy;
};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain)
#endif // ENABLE_SERIALIZATION

#endif // CONSTITUTIVESTATICDATAPrevEngineeringStressStrain3D_H
