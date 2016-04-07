#pragma once

#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"

//! @brief ... class, storing the static data for the Mises plasticity3D including energy updates
//! this is only required, if the energy should be calculated, otherwise the base class ConstitutiveStaticDataMisesPlasticity3D
//! is sufficient
//! @author Jörg F. Unger, ISM
//! @date December 2009
namespace NuTo
{
template <int TDim>
class ConstitutiveStaticDataPrevEngineeringStressStrain : virtual public ConstitutiveStaticDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief constructor
    ConstitutiveStaticDataPrevEngineeringStressStrain() :
        mPrevTotalEnergy(0), mPrevElasticEnergy(0)
    {
        mPrevStress.setZero();
        mPrevStrain.setZero();
    }

    NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain& operator= (const NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain&  rOther) = default;
    NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain& operator= (      NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain&& rOther) = default;

    //! brief set the previous stress
    void SetPrevStress(const EngineeringStress<TDim>& rPrevStress)
    {
        mPrevStress = rPrevStress;
    }

    //! brief returns the previous stress
    const EngineeringStress<TDim>& GetPrevStress() const
    {
        return mPrevStress;
    }

    //! brief set the previous stress
    void SetPrevStrain(const EngineeringStrain<TDim>& rPrevStrain)
    {
        mPrevStrain = rPrevStrain;
    }

    //! brief returns the previous stress
    const EngineeringStrain<TDim>& GetPrevStrain() const
    {
        return mPrevStrain;
    }

    //! brief returns the previous total energy
    double GetPrevTotalEnergy() const
    {
        return mPrevTotalEnergy;
    }

    //! brief set the previous total energy
    void SetPrevTotalEnergy(double rPrevTotalEnergy)
    {
        mPrevTotalEnergy = rPrevTotalEnergy;
    }

    //! brief returns the previous elastic energy
    double GetPrevElasticEnergy() const
    {
        return mPrevElasticEnergy;
    }

    //! brief set the previous elastic energy
    void SetPrevElasticEnergy(double rPrevElasticEnergy)
    {
        mPrevElasticEnergy = rPrevElasticEnergy;
    }


protected:
    //! @brief previous stress
    EngineeringStress<TDim> mPrevStress;

    //! @brief previous strain
    EngineeringStrain<TDim> mPrevStrain;

    //! @brief previous total energy
    double mPrevTotalEnergy;

    //! @brief previous elastic energy
    double mPrevElasticEnergy;

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain3D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveStaticDataPrevEngineeringStressStrain" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataBase)
       & BOOST_SERIALIZATION_NVP(mPrevSigma)
       & BOOST_SERIALIZATION_NVP(mPrevStrain)
       & BOOST_SERIALIZATION_NVP(mPrevTotalEnergy)
       & BOOST_SERIALIZATION_NVP(mPrevElasticEnergy);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveStaticDataPrevEngineeringStressStrain" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain)
#endif // ENABLE_SERIALIZATION

};

}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstitutiveStaticDataPrevEngineeringStressStrain)
#endif // ENABLE_SERIALIZATION

