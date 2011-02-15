// $Id: ConstraintGlobalCrackAngle.cpp 314 2010-09-27 16:31:43Z unger3 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constraints/ConstraintLinearGlobalCrackAngle.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/structures/unstructured/StructureIp.h"
#include "nuto/math/FullMatrix.h"

// constructor
NuTo::ConstraintLinearGlobalCrackAngle::ConstraintLinearGlobalCrackAngle(const StructureIp* rStructure):
        ConstraintLinear()
{
    mStructure = rStructure;
}

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintLinearGlobalCrackAngle::GetNumLinearConstraints()const
{
    return 1;
}

//! @brief adds the constraint equations to the matrix
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
//! @param rRHS right hand side of the constraint equation
void NuTo::ConstraintLinearGlobalCrackAngle::AddToConstraintMatrix(int& curConstraintEquation,
        NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix,
        NuTo::FullMatrix<double>& rRHS)const
{
    // set right hand side value
    // calcute principal direction of total strain
    EngineeringStrain2D totalStrain(mStructure->GetTotalStrain());
    if (totalStrain.mEngineeringStrain[0]-totalStrain.mEngineeringStrain[1]!=0.)
        rRHS(curConstraintEquation,0) = atan(totalStrain.mEngineeringStrain[2]/(totalStrain.mEngineeringStrain[0]-totalStrain.mEngineeringStrain[1]));
    else
        rRHS(curConstraintEquation,0) = 0.;
    rRHS(curConstraintEquation,0) = M_PI*0.49;
    rConstraintMatrix.AddEntry(curConstraintEquation,mStructure->GetDofCrackAngle(), 1);

    // increase constraint equation number
    curConstraintEquation++;
}

//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void NuTo::ConstraintLinearGlobalCrackAngle::Info(unsigned short rVerboseLevel) const
{
    double anglePrescribed(0);
    // calcute principal direction of total strain
    EngineeringStrain2D totalStrain(mStructure->GetTotalStrain());
    if (totalStrain.mEngineeringStrain[0]-totalStrain.mEngineeringStrain[1]!=0.)
        anglePrescribed = atan(totalStrain.mEngineeringStrain[2]/(totalStrain.mEngineeringStrain[0]-totalStrain.mEngineeringStrain[1]));
    else
        anglePrescribed = 0.;
    std::cout << "NuTo::ConstraintLinearGlobalCrackAngle : prescribed angle " <<  anglePrescribed << std::endl;
}


#ifdef ENABLE_SERIALIZATION
// serialize
template void NuTo::ConstraintLinearGlobalCrackAngle::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearGlobalCrackAngle::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearGlobalCrackAngle::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearGlobalCrackAngle::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearGlobalCrackAngle::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearGlobalCrackAngle::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintLinearGlobalCrackAngle::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLinearGlobalCrackAngle" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintBase)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintLinear)
       & BOOST_SERIALIZATION_NVP(mStructure);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLinearGlobalCrackAngle" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintLinearGlobalCrackAngle)
#endif // ENABLE_SERIALIZATION
