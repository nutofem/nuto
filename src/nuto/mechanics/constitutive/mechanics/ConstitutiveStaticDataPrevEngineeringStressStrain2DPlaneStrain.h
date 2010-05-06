// $Id: ConstitutiveStaticPrevEngineeringStressStrain2DPlaneStrain.h 87 2009-11-06 10:35:39Z unger3 $

#ifndef CONSTITUTIVESTATICDATAPrevEngineeringStressStrain2D_H
#define CONSTITUTIVESTATICDATAPrevEngineeringStressStrain2D_H

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"

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
    void SetPrevStress(const double rPrevSigma[4]);

    //! brief returns the previous stress
    inline const double* GetPrevStress()const
    {
        return  mPrevSigma;
    }

    //! brief set the previous stress
    void SetPrevStrain(const double rPrevStrain[3]);

    //! brief returns the previous stress
    inline const double* GetPrevStrain()const
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
    double mPrevSigma[4];

    //! @brief previous strain
    double mPrevStrain[3];

    //! @brief previous total energy
    double mPrevTotalEnergy;

    //! @brief previous elastic energy
    double mPrevElasticEnergy;
};

}

#endif // CONSTITUTIVESTATICDATAPrevEngineeringStressStrain3D_H
