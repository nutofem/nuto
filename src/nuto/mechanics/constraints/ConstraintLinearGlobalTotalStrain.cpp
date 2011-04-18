// $Id: ConstraintLinearGlobalTotalStrain.cpp 314 2010-09-27 16:31:43Z unger3 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include <iostream>

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constraints/ConstraintLinearGlobalTotalStrain.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/structures/unstructured/StructureMultiscale.h"
#include "nuto/math/FullMatrix.h"

// constructor
NuTo::ConstraintLinearGlobalTotalStrain::ConstraintLinearGlobalTotalStrain(const StructureMultiscale* rStructure, const EngineeringStrain2D& rStrain):
        ConstraintLinear()
{
    mStructure = rStructure;
    mStrain = rStrain;
}

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintLinearGlobalTotalStrain::GetNumLinearConstraints()const
{
    return 3;
}

//! @brief sets the rhs of the constraint equation
void NuTo::ConstraintLinearGlobalTotalStrain::SetRHS(const EngineeringStrain2D& rStrain)
{
    mStrain = rStrain;
}

//! @brief adds the constraint equations to the matrix
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
//! @param rRHS right hand side of the constraint equation
void NuTo::ConstraintLinearGlobalTotalStrain::AddToConstraintMatrix(int& curConstraintEquation,
        NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix,
        NuTo::FullMatrix<double>& rRHS)const
{
    rRHS(curConstraintEquation,0) = mStrain.mEngineeringStrain[0];
    rConstraintMatrix.AddEntry(curConstraintEquation,mStructure->GetDofGlobalTotalStrain2D()[0],1);
    curConstraintEquation++;

    rRHS(curConstraintEquation,0) = mStrain.mEngineeringStrain[1];
    rConstraintMatrix.AddEntry(curConstraintEquation,mStructure->GetDofGlobalTotalStrain2D()[1],1);
    curConstraintEquation++;

    rRHS(curConstraintEquation,0) = mStrain.mEngineeringStrain[2];
    rConstraintMatrix.AddEntry(curConstraintEquation,mStructure->GetDofGlobalTotalStrain2D()[2],1);
    curConstraintEquation++;
}

#ifdef ENABLE_SERIALIZATION
// serialize
template void NuTo::ConstraintLinearGlobalTotalStrain::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearGlobalTotalStrain::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearGlobalTotalStrain::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearGlobalTotalStrain::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearGlobalTotalStrain::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearGlobalTotalStrain::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintLinearGlobalTotalStrain::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLinearGlobalTotalStrain" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintLinear)
       & BOOST_SERIALIZATION_NVP(mStructure)
       & BOOST_SERIALIZATION_NVP(mStrain);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLinearGlobalTotalStrain" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintLinearGlobalTotalStrain)
#endif // ENABLE_SERIALIZATION
