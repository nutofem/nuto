// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/ptr_container/serialize_ptr_map.hpp>
#endif // ENABLE_SERIALIZATION

#include <boost/spirit/include/classic_core.hpp>

#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/structures/unstructured/StructureIp.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/mechanics/elements/ElementBase.h"

#include <ANN/ANN.h>
#include <set>

#define PRINTRESULT true
#define MAXNUMNEWTONITERATIONS 20
#define MIN_DELTA_STRAIN_FACTOR 1e-7
//! @brief constructor
//! @param mDimension number of nodes
NuTo::StructureIp::StructureIp ( int rDimension )  : Structure ( rDimension )
{
    if (rDimension!=2)
        throw MechanicsException("[NuTo::StructureIp::StructureIp] The concurrent multiscale model is only implemented for 2D.");
    mCrackAngle[0] = 0.;
    mCrackAngle[1] = 0.;
    mCrackAngle[2] = 0.;
    mDOFCrackAngle[0] = -1;
    mDOFCrackAngle[1] = -1;
    mDOFCrackAngle[2] = -1;
    mCrackOpening[0] = 0.;
    mCrackOpening[1] = 0.;
    mCrackOpening[2] = 0.;
    mDOFCrackOpening[0] = -1;
    mDOFCrackOpening[1] = -1;
    mDOFCrackOpening[2] = -1;
}

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::StructureIp::Info()const
{
    Structure::Info();
    //add crack angle and crack orientation

}

#ifdef ENABLE_SERIALIZATION
template void NuTo::StructureIp::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::StructureIp::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::StructureIp::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::StructureIp::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::StructureIp::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::StructureIp::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::StructureIp::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialization of structureIp" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Structure)
       & BOOST_SERIALIZATION_NVP(mCrackAngle)
       & BOOST_SERIALIZATION_NVP(mDOFCrackAngle)
       & BOOST_SERIALIZATION_NVP(mCrackOpening)
       & BOOST_SERIALIZATION_NVP(mDOFCrackOpening);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialization of structureIp" << std::endl;
#endif
}

//! @brief ... save the object to a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::StructureIp::Save (const std::string &filename, std::string rType )const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    try
    {
        //transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), (int(*)(int))toupper);

        // open file
        std::ofstream ofs ( filename.c_str(), std::ios_base::binary );
        if(! ofs.is_open())
        {
            throw MechanicsException("[NuTo::StructureIp::Save] Error opening file.");
        }

        // write data to file
        std::string typeIdString(this->GetTypeId());
        if (rType=="BINARY")
        {
            boost::archive::binary_oarchive oba ( ofs, std::ios::binary );
            oba & boost::serialization::make_nvp ("Object_type", typeIdString );
            oba & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_oarchive oxa ( ofs, std::ios::binary );
            std::string tmpString(this->GetTypeId());
            oxa & boost::serialization::make_nvp ("Object_type", typeIdString );
            oxa & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_oarchive ota ( ofs, std::ios::binary );
            ota & boost::serialization::make_nvp("Object_type", typeIdString );
            ota & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else
        {
            throw MechanicsException ( "[NuTo::StructureIp::Save] File type not implemented." );
        }

        // close file
        ofs.close();
    }
    catch ( boost::archive::archive_exception e )
    {
        std::string s ( std::string ( "[NuTo::StructureIp::Save]File save exception in boost - " ) + std::string ( e.what() ) );
        throw MechanicsException ( s );
    }
    catch ( MechanicsException &e )
    {
        throw e;
    }
    catch ( std::exception &e )
    {
        throw MechanicsException ( e.what() );
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::Structure::Save] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}


//! @brief ... restore the object from a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::StructureIp::Restore (const std::string &filename, std::string rType )
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    try
    {
        //transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), (int(*)(int))toupper);

        // open file
        std::ifstream ifs ( filename.c_str(), std::ios_base::binary );
        if(! ifs.is_open())
        {
            throw MechanicsException("[NuTo::StructureIp::Restore] Error opening file.");
        }

        std::string typeIdString;
        if (rType=="BINARY")
        {
            boost::archive::binary_iarchive oba ( ifs, std::ios::binary );
            oba & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if ( typeIdString != this->GetTypeId() )
            {
                throw MechanicsException ( "[NuTo::StructureIp::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
            }
            oba & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if ( typeIdString != this->GetTypeId() )
            {
                throw MechanicsException ( "[NuTo::StructureIp::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
            }
            oxa & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_iarchive ota ( ifs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if ( typeIdString != this->GetTypeId() )
            {
                throw MechanicsException ( "[NuTo::StructureIp::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
            }
            ota & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else
        {
            throw MechanicsException ( "[NuTo::StructureIp::Restore]File type not implemented" );
        }
        // close file
        ifs.close();
    }
    catch ( boost::archive::archive_exception e )
    {
        std::string s ( std::string ( "[NuTo::StructureIp::Restore] File save exception in boost - " ) + std::string ( e.what() ) );
        throw MechanicsException ( s );
    }
    catch ( MechanicsException &e )
    {
        throw e;
    }
    catch ( std::exception &e )
    {
        throw MechanicsException ( e.what() );
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::StructureIp::Restore] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}
#endif// ENABLE_SERIALIZATION


//************ constitutive routines    ***********
//**  defined in structures/StructureIpConstitutive.cpp *********
//*************************************************
//  Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering plastic strain from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
void NuTo::StructureIp::GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
                                  const DeformationGradient1D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const
{
    throw MechanicsException("[NuTo::StructureIp::GetEngineeringPlasticStrain] not implemented.");
}
//  Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering plastic strain from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
void NuTo::StructureIp::GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
                                  const DeformationGradient2D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const
{
    throw MechanicsException("[NuTo::StructureIp::GetEngineeringPlasticStrain] not implemented.");
}

//  Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering plastic strain from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
void NuTo::StructureIp::GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
                                  const DeformationGradient3D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const
{
    throw MechanicsException("[NuTo::StructureIp::GetEngineeringPlasticStrain] not implemented.");
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStress ... Engineering stress
void NuTo::StructureIp::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient, EngineeringStress1D& rEngineeringStress) const
{
    throw MechanicsException("[NuTo::StructureIp::GetEngineeringStressFromEngineeringStrain] not implemented for 1D.");
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStress ... Engineering stress
void NuTo::StructureIp::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const
{
    throw MechanicsException("[NuTo::StructureIp::GetEngineeringStressFromEngineeringStrain] not implemented for 1D.");
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStress ... Engineering stress
void NuTo::StructureIp::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient, EngineeringStress2D& rEngineeringStress) const
{
     // calculate engineering strain
    EngineeringStrain2D engineeringStrain;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);

    //apply boundary conditions and solve
    //myStructureFineScale.ConstraintPeriodicSetCrackOpening(constraintPeriodic,curCrackOpening);
    if (mSolveForLocalization)
        throw MechanicsException("[NuTo::StructureIp::GetEngineeringStressFromEngineeringStrain] not implemented.");

    //solve
    double tolerance(1e-6);
    // this is a somehow weird situation, since for a normal material law nothing should be changed
    // since the material law is a full structure whose bc change, this can either be implemented with a cast (as I did)
    // or by removing the const flag from all material routines (which I do not consider is good)
    (const_cast<NuTo::StructureIp*>(&*this))->Solve(engineeringStrain,tolerance);

    //calculate average stress
    NuTo::FullMatrix<double> averageStress;
    ElementTotalGetAverageStress(mlX*mlY,averageStress);
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStress ... Engineering stress
void NuTo::StructureIp::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const
{
    throw MechanicsException("[NuTo::StructureIp::GetEngineeringStressFromEngineeringStrain] not implemented.");
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rCauchyStress ... Cauchy stress
void NuTo::StructureIp::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const
{
    throw MechanicsException("[NuTo::StructureIp::GetEngineeringStressFromEngineeringStrain] not implemented.");
}

//  Damage /////////////////////////////////////
//! @brief ... calculate isotropic damage from deformation gradient in 1D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDamage ... damage variable
void NuTo::StructureIp::GetDamage(const ElementBase* rElement, int rIp,
                                  const DeformationGradient1D& rDeformationGradient, double& rDamage) const
{
    throw MechanicsException("[NuTo::StructureIp::GetDamage] not implemented.");
}

//  Damage /////////////////////////////////////
//! @brief ... calculate isotropic damage from deformation gradient in 2D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDamage ... damage variable
void NuTo::StructureIp::GetDamage(const ElementBase* rElement, int rIp,
                                  const DeformationGradient2D& rDeformationGradient, double& rDamage) const
{
    throw MechanicsException("[NuTo::StructureIp::GetDamage] not implemented.");
}

//  Damage /////////////////////////////////////
//! @brief ... calculate isotropic damage from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDamage ... damage variable
void NuTo::StructureIp::GetDamage(const ElementBase* rElement, int rIp,
                                  const DeformationGradient3D& rDeformationGradient, double& rDamage) const
{
    throw MechanicsException("[NuTo::StructureIp::GetDamage] not implemented.");
}

//! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rTangent ... tangent
void NuTo::StructureIp::GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient,
        ConstitutiveTangentBase* rTangent) const
{
    throw MechanicsException("[NuTo::StructureIp::GetTangent_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rTangent ... tangent
void NuTo::StructureIp::GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient,
        ConstitutiveTangentBase* rTangent) const
{
    throw MechanicsException("[NuTo::StructureIp::GetTangent_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rTangent ... tangent
void NuTo::StructureIp::GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient,
        ConstitutiveTangentBase* rTangent) const
{
    throw MechanicsException("[NuTo::StructureIp::GetTangent_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::StructureIp::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::StructureIp::UpdateStaticData_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::StructureIp::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::StructureIp::UpdateStaticData_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::StructureIp::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::StructureIp::UpdateStaticData_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... update tmp static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::StructureIp::UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::StructureIp::UpdateTmpStaticData_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... update tmp static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::StructureIp::UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::StructureIp::UpdateTmpStaticData_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... update tmp static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::StructureIp::UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::StructureIp::UpdateTmpStaticData_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::StructureIp::AllocateStaticDataEngineeringStress_EngineeringStrain1D(const ElementBase* rElement) const
{
    throw MechanicsException("[NuTo::StructureIp::AllocateStaticDataEngineeringStress_EngineeringStrain1D] not implemented.");
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::StructureIp::AllocateStaticDataEngineeringStress_EngineeringStrain2D(const ElementBase* rElement) const
{
    throw MechanicsException("[NuTo::StructureIp::AllocateStaticDataEngineeringStress_EngineeringStrain1D] not implemented.");
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::StructureIp::AllocateStaticDataEngineeringStress_EngineeringStrain3D(const ElementBase* rElement) const
{
    throw MechanicsException("[NuTo::StructureIp::AllocateStaticDataEngineeringStress_EngineeringStrain1D] not implemented.");
}

//! @brief ... calculates the difference of the elastic strain between the current state and the previous update
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
void NuTo::StructureIp::GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient, EngineeringStrain1D& rDeltaElasticEngineeringStrain) const
{
    throw MechanicsException("[NuTo::StructureIp::GetDeltaElasticEngineeringStrain] not implemented.");
}

//! @brief ... calculates the difference of the elastic strain between the current state and the previous update
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
void NuTo::StructureIp::GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient, EngineeringStrain2D& rDeltaElasticEngineeringStrain) const
{
    throw MechanicsException("[NuTo::StructureIp::GetDeltaElasticEngineeringStrain] not implemented.");
}

//! @brief ... calculates the difference of the elastic strain between the current state and the previous update
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
void NuTo::StructureIp::GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient, EngineeringStrain3D& rDeltaElasticEngineeringStrain) const
{
    throw MechanicsException("[NuTo::StructureIp::GetDeltaElasticEngineeringStrain] not implemented.");
}

//! @brief ... check parameters of the constitutive relationship
void NuTo::StructureIp::CheckParameters()const
{
    throw MechanicsException("[NuTo::StructureIp::CheckParameters] not implemented.");
}

//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool NuTo::StructureIp::CheckElementCompatibility(Element::eElementType rElementType) const
{
    throw MechanicsException("[NuTo::StructureIp::CheckElementCompatibility] not implemented.");
}

//! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
//! @return ... see brief explanation
bool NuTo::StructureIp::HaveTmpStaticData() const
{
    throw MechanicsException("[NuTo::StructureIp::HaveTmpStaticData] not implemented.");
}

//! @brief ... get type of constitutive relationship
//! @return ... type of constitutive relationship
//! @sa eConstitutiveType
NuTo::Constitutive::eConstitutiveType NuTo::StructureIp::GetType() const
{
    throw MechanicsException("[NuTo::StructureIp::GetType] not implemented.");
}


void NuTo::StructureIp::Solve(NuTo::EngineeringStrain2D rEngineeringStrain, double rTolerance)
{
    // start analysis
    double deltaStrainFactor(1.0);
    double curStrainFactor(1.0);

    NuTo::FullMatrix<double> totalEngineeringStrain(3,1);
    totalEngineeringStrain(0,0) = rEngineeringStrain.mEngineeringStrain[0];
    totalEngineeringStrain(1,0) = rEngineeringStrain.mEngineeringStrain[1];
    totalEngineeringStrain(2,0) = rEngineeringStrain.mEngineeringStrain[2];

    NuTo::FullMatrix<double> deltaEngineeringStrain(3,1);
    deltaEngineeringStrain(0,0) = rEngineeringStrain.mEngineeringStrain[0]-mPrevStrain[0];
    deltaEngineeringStrain(1,0) = rEngineeringStrain.mEngineeringStrain[1]-mPrevStrain[1];
    deltaEngineeringStrain(2,0) = rEngineeringStrain.mEngineeringStrain[2]-mPrevStrain[2];

    NuTo::FullMatrix<double> curEngineeringStrain(3,1);
    curEngineeringStrain(0,0) = mPrevStrain[0] + deltaEngineeringStrain(0,0) * curStrainFactor;
    curEngineeringStrain(1,0) = mPrevStrain[1] + deltaEngineeringStrain(1,0) * curStrainFactor;
    curEngineeringStrain(2,0) = mPrevStrain[2] + deltaEngineeringStrain(2,0) * curStrainFactor;

    NuTo::FullMatrix<double> curCrackOpening(3,1);
    //update displacement of boundary (disp controlled)
    this->ConstraintPeriodicSetStrain(mConstraintPeriodic,curEngineeringStrain);

    //update conre mat
    this->NodeBuildGlobalDofs();

    //update tmpstatic data with zero displacements
    this->ElementTotalUpdateTmpStaticData();

    //init some auxiliary variables
    NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixCSRVector2;
    NuTo::FullMatrix<double> dispForceVector;
    NuTo::FullMatrix<double> intForceVector;
    NuTo::FullMatrix<double> extForceVector;
    NuTo::FullMatrix<double> rhsVector;

    //allocate solver
    NuTo::SparseDirectSolverMUMPS mySolver;
    mySolver.SetShowTime(false);

    //calculate stiffness
    this->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);

    // build global external load vector and RHS vector
    this->BuildGlobalExternalLoadVector(extForceVector);
    rhsVector = extForceVector + dispForceVector;

    //calculate absolute tolerance for matrix entries to be not considered as zero
    double maxValue, minValue, ToleranceZeroStiffness;
    stiffnessMatrixCSRVector2.Max(maxValue);
    stiffnessMatrixCSRVector2.Min(minValue);
    //std::cout << "min and max " << minValue << " , " << maxValue << std::endl;

    ToleranceZeroStiffness = (1e-14) * (fabs(maxValue)>fabs(minValue) ?  fabs(maxValue) : fabs(minValue));
    this->SetToleranceStiffnessEntries(ToleranceZeroStiffness);
    int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
    int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
    //std::cout << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << std::endl;

    //update displacements of all nodes according to the new conre mat
    {
        NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
        NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
        this->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
        this->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
        this->ElementTotalUpdateTmpStaticData();
    }

    //repeat until max displacement is reached
    bool convergenceStatusLoadSteps(false);
    while (!convergenceStatusLoadSteps)
    {

        double normResidual(1);
        double maxResidual(1);
        int numNewtonIterations(0);
        double normRHS(1.);
        double alpha(1.);
        int convergenceStatus(0);
        //0 - not converged, continue Newton iteration
        //1 - converged
        //2 - stop iteration, decrease load step
        while(convergenceStatus==0)
        {
            numNewtonIterations++;

            if (numNewtonIterations>MAXNUMNEWTONITERATIONS)
            {
                if (PRINTRESULT)
                {
                    std::cout << "numNewtonIterations (" << numNewtonIterations << ") > MAXNUMNEWTONITERATIONS (" << MAXNUMNEWTONITERATIONS << ")" << std::endl;
                }
                convergenceStatus = 2; //decrease load step
                break;
            }

            normRHS = rhsVector.Norm();

            // solve
            NuTo::FullMatrix<double> deltaDisplacementsActiveDOFs;
            NuTo::FullMatrix<double> oldDisplacementsActiveDOFs;
            NuTo::FullMatrix<double> displacementsActiveDOFs;
            NuTo::FullMatrix<double> displacementsDependentDOFs;
            NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrixCSR(stiffnessMatrixCSRVector2);
            stiffnessMatrixCSR.SetOneBasedIndexing();
            mySolver.Solve(stiffnessMatrixCSR, rhsVector, deltaDisplacementsActiveDOFs);

            // write displacements to node
            this->NodeExtractDofValues(oldDisplacementsActiveDOFs, displacementsDependentDOFs);

            //perform a linesearch
            alpha = 1.;
            do
            {
                //add new displacement state
                displacementsActiveDOFs = oldDisplacementsActiveDOFs + deltaDisplacementsActiveDOFs*alpha;
                this->NodeMergeActiveDofValues(displacementsActiveDOFs);
                this->ElementTotalUpdateTmpStaticData();

                // calculate residual
                this->BuildGlobalGradientInternalPotentialVector(intForceVector);
                rhsVector = extForceVector - intForceVector;
                normResidual = rhsVector.Norm();
                //std::cout << "alpha " << alpha << ", normResidual " << normResidual << ", normResidualInit "<< normRHS << ", normRHS*(1-0.5*alpha) " << normRHS*(1-0.5*alpha) << std::endl;
                alpha*=0.5;
            }
            while(alpha>1e-3 && normResidual>normRHS*(1-0.5*alpha) && normResidual>1e-5);
            if (normResidual>normRHS*(1-0.5*alpha) && normResidual>1e-5)
            {
                convergenceStatus=2;
                break;
            }

            maxResidual = rhsVector.Max();

            //std::cout << std::endl << "Newton iteration " << numNewtonIterations << ", final alpha " << 2*alpha << ", normResidual " << normResidual<< ", maxResidual " << maxResidual<<std::endl;

            //check convergence
            if (normResidual<rTolerance || maxResidual<rTolerance)
            {
                if (PRINTRESULT)
                {
                    std::cout << "Convergence after " << numNewtonIterations << " Newton iterations, curStrainFactor " << curStrainFactor << ", deltaStrainFactor "<< deltaStrainFactor << std::endl<< std::endl;
                }
                convergenceStatus=1;
                break;
            }

            //convergence status == 0 (continue Newton iteration)
            //build new stiffness matrix
            this->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
            int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
            int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
            std::cout << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << std::endl;
        }

        if (deltaStrainFactor==0)
            throw NuTo::MechanicsException("[NuTo::StructureIp::Solve] No convergence, delta strain factor < 1e-7");

        if (convergenceStatus==1)
        {
            throw MechanicsException("[NuTo::StructureIp::Solve] Do not update the static data without restoring it afterwards");
            this->ElementTotalUpdateStaticData();
            // visualize results

        #ifdef ENABLE_VISUALIZE
            std::cout << " store element id and ip in output file" << std::endl;
            this->ExportVtkDataFile(mIPName+std::string(".vtk"));
        #endif
        #ifdef ENABLE_SERIALIZATION
            this->Save(mIPName+std::string("bin"),"BINARY");
        #endif // ENABLE_SERIALIZATION

            if (curStrainFactor==1)
                convergenceStatusLoadSteps=true;
            else
            {
                //eventually increase load step
                if (numNewtonIterations<MAXNUMNEWTONITERATIONS/3)
                {
                    deltaStrainFactor*=1.5;
                }

                //increase displacement
                curStrainFactor+=deltaStrainFactor;
                if (curStrainFactor>1)
                {
                    deltaStrainFactor -= curStrainFactor -1.;
                    curStrainFactor=1;
                }

                curEngineeringStrain(0,0) = mPrevStrain[0] + deltaEngineeringStrain(0,0) * curStrainFactor;
                curEngineeringStrain(1,0) = mPrevStrain[1] + deltaEngineeringStrain(1,0) * curStrainFactor;
                curEngineeringStrain(2,0) = mPrevStrain[2] + deltaEngineeringStrain(2,0) * curStrainFactor;

                //old stiffness matrix is used in first step of next load increment in order to prevent spurious problems at the boundary
                //std::cout << "press enter to next load increment, delta strain factor " << deltaStrainFactor << " max delta strain factor " <<  maxDeltaStrainFactor << std::endl << std::endl;
                //char cDummy[100]="";
                //std::cin.getline(cDummy, 100);
            }
        }
        else
        {
            assert(convergenceStatus==2);
            //calculate stiffness of previous loadstep (used as initial stiffness in the next load step)
            //this is done within the loop in order to ensure, that for the first step the stiffness matrix of the previous step is used
            //otherwise, the additional boundary displacements will result in an artifical localization in elements at the boundary
            curStrainFactor-=deltaStrainFactor;
            curEngineeringStrain(0,0) = mPrevStrain[0] + deltaEngineeringStrain(0,0) * curStrainFactor;
            curEngineeringStrain(1,0) = mPrevStrain[1] + deltaEngineeringStrain(1,0) * curStrainFactor;
            curEngineeringStrain(2,0) = mPrevStrain[2] + deltaEngineeringStrain(2,0) * curStrainFactor;

            this->ConstraintPeriodicSetStrain(mConstraintPeriodic,curEngineeringStrain);

            // build global dof numbering
            this->NodeBuildGlobalDofs();

            //set previous converged displacements
            NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
            NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
            this->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
            this->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
            this->ElementTotalUpdateTmpStaticData();

            // calculate previous residual (should be almost zero)
            this->BuildGlobalGradientInternalPotentialVector(intForceVector);

            //decrease load step
            deltaStrainFactor*=0.5;
            curStrainFactor+=deltaStrainFactor;

            //check for minimum delta (this mostly indicates an error in the software
            if (deltaStrainFactor<MIN_DELTA_STRAIN_FACTOR)
            {
                deltaStrainFactor = 0;
                //throw NuTo::MechanicsException("Example ConcurrentMultiscale : No convergence, delta strain factor < 1e-7");
            }

            //std::cout << "press enter to reduce load increment" << std::endl;
            //char cDummy[100]="";
            //std::cin.getline(cDummy, 100);;
        }

        if (!convergenceStatusLoadSteps)
        {
            //update new displacement of RHS
            this->ConstraintPeriodicSetStrain(mConstraintPeriodic,curEngineeringStrain);

            // build global dof numbering
            this->NodeBuildGlobalDofs();

            //update stiffness in order to calculate new dispForceVector
            this->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
            int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
            int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
            //std::cout << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << std::endl;

            //update rhs vector for next Newton iteration
            rhsVector = dispForceVector + extForceVector - intForceVector;

            //update displacements of all nodes according to the new conre mat
            NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
            NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
            this->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
            this->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
            this->ElementTotalUpdateTmpStaticData();
        }
    }
}

#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::StructureIp)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
