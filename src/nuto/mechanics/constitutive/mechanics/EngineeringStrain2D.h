// $Id: EngineeringStrain2D.h $

#ifndef ENGINEERINGSTRAIN2D_H
#define ENGINEERINGSTRAIN2D_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

namespace NuTo
{
class DeformationGradient1D;
class DeformationGradient2D;
class LinearElastic;
class ConstitutiveMisesPlasticity;

//! @brief ... three-dimensional deformation gradient
//! @author Jörg F. Unger, ISM
//! @date November 2009
class EngineeringStrain2D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
    friend class ConstitutiveEngineeringStressStrain;
    friend class ConstraintLinearGlobalCrackAngle;
    friend class ConstraintLinearDisplacementsPeriodic2D;
    friend class ConstraintLinearGlobalTotalStrain;
    friend class ConstitutiveMisesPlasticity;
    friend class ConstraintNonlinearGlobalCrackAngle2D;
    friend class ConstitutiveStaticDataMultiscale2DPlaneStrain;
    friend class ConstitutiveStaticDataPrevEngineeringStressStrain2DPlaneStrain;
    friend class LinearElastic;
    friend class Multiscale;
    friend class NonlocalDamagePlasticity;
    friend class StructureMultiscale;
    public:
    //! @brief ... constructor
    //! @param pStructure ... structure
    //! @param pElement ... element
    //! @param pIntegrationPoint ... integration point
    EngineeringStrain2D();

    //! @brief ... copy constructor
    EngineeringStrain2D(const DeformationGradient2D& rDeformationGradient);

    //! @brief ... copy constructor
    EngineeringStrain2D(const EngineeringStrain2D& rEngineeringStrain2D);

    //! @brief ... get number of strain components
    //! @return ... number of strain components
    unsigned int GetNumberOfComponents() const;

    //! @brief ... get Engineering Strain
    //! @return ... Engineering Strain (exx,eyy,gxy)
    //! @sa mDeformationGradient
    const double* GetData() const;

    EngineeringStrain2D operator+ ( const EngineeringStrain2D &other ) const
    {
        EngineeringStrain2D result;
        result.mEngineeringStrain[0] = mEngineeringStrain[0] + other.mEngineeringStrain[0];
        result.mEngineeringStrain[1] = mEngineeringStrain[1] + other.mEngineeringStrain[1];
        result.mEngineeringStrain[2] = mEngineeringStrain[2] + other.mEngineeringStrain[2];
        return result;
    }

    EngineeringStrain2D operator- ( const EngineeringStrain2D &other ) const
    {
        EngineeringStrain2D result;
        result.mEngineeringStrain[0] = mEngineeringStrain[0] - other.mEngineeringStrain[0];
        result.mEngineeringStrain[1] = mEngineeringStrain[1] - other.mEngineeringStrain[1];
        result.mEngineeringStrain[2] = mEngineeringStrain[2] - other.mEngineeringStrain[2];
        return result;
    }

    EngineeringStrain2D operator* ( double scalar ) const
    {
        EngineeringStrain2D result;
        result.mEngineeringStrain[0] = scalar*mEngineeringStrain[0];
        result.mEngineeringStrain[1] = scalar*mEngineeringStrain[1];
        result.mEngineeringStrain[2] = scalar*mEngineeringStrain[2];
        return result;
    }

    //! @brief ... set Engineering Strain
    //! @return ... Engineering Strain (exx,eyy,gxy)
    //! @sa mDeformationGradient
    void SetData(double rEngineeringStrain[3]);


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
private:
    //! @brief ... engineering strain
    //! The engineering strain is stored as :
    //! (exx,eyy,gxy)
    double mEngineeringStrain[3];
};

}

#endif // ENGINEERINGSTRAIN2D_H
