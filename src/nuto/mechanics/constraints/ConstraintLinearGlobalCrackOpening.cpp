// $Id: ConstraintGlobalCrackOpening.cpp 314 2010-09-27 16:31:43Z unger3 $

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
#include "nuto/mechanics/constraints/ConstraintLinearGlobalCrackOpening.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/structures/unstructured/StructureMultiscale.h"
#include "nuto/math/FullMatrix.h"

// constructor
NuTo::ConstraintLinearGlobalCrackOpening::ConstraintLinearGlobalCrackOpening(const StructureMultiscale* rStructure, const NuTo::FullMatrix<double>& rDirection, double rValue):
        ConstraintLinear()
{
    mStructure = rStructure;
    if (rDirection.GetNumColumns()!=1 || rDirection.GetNumRows()!=2)
        throw MechanicsException("[NuTo::ConstraintLinearGlobalCrackOpening::ConstraintLinearGlobalCrackOpening] Dimension of the direction matrix must be equal to the dimension of the structure.");

    memcpy(mDirection,rDirection.mEigenMatrix.data(),2*sizeof(double));
    //normalize the direction
    double norm = sqrt(mDirection[0]*mDirection[0]+mDirection[1]*mDirection[1]);
    if (norm < 1e-14)
    {
        throw MechanicsException("[NuTo::ConstraintLinearGlobalCrackOpening::ConstraintLinearGlobalCrackOpening] direction vector has zero length.");
    }
    double invNorm = 1./norm;
    mDirection[0]*=invNorm;
    mDirection[1]*=invNorm;
    mRHS = rValue;
}

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintLinearGlobalCrackOpening::GetNumLinearConstraints()const
{
    return 1;
}

//!@brief sets/modifies the right hand side of the constraint equation
//!@param rRHS new right hand side
void NuTo::ConstraintLinearGlobalCrackOpening::SetRHS(double rRHS)
{
	mRHS=rRHS;
}

//! @brief adds the constraint equations to the matrix
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
//! @param rRHS right hand side of the constraint equation
void NuTo::ConstraintLinearGlobalCrackOpening::AddToConstraintMatrix(int& curConstraintEquation,
        NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix,
        NuTo::FullMatrix<double>& rRHS)const
{
    rRHS(curConstraintEquation,0) = mRHS;

    if (fabs(mDirection[0])>1e-18)
        rConstraintMatrix.AddEntry(curConstraintEquation,mStructure->GetDofGlobalCrackOpening2D()[0],mDirection[0]);
    if (fabs(mDirection[1])>1e-18)
        rConstraintMatrix.AddEntry(curConstraintEquation,mStructure->GetDofGlobalCrackOpening2D()[1],mDirection[1]);

    curConstraintEquation++;
}

#ifdef ENABLE_SERIALIZATION
// serialize
template void NuTo::ConstraintLinearGlobalCrackOpening::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearGlobalCrackOpening::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearGlobalCrackOpening::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearGlobalCrackOpening::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearGlobalCrackOpening::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearGlobalCrackOpening::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintLinearGlobalCrackOpening::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLinearGlobalCrackOpening" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintLinear)
       & BOOST_SERIALIZATION_NVP(mStructure)
       & BOOST_SERIALIZATION_NVP(mRHS)
       & BOOST_SERIALIZATION_NVP(mDirection);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLinearGlobalCrackOpening" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintLinearGlobalCrackOpening)
#endif // ENABLE_SERIALIZATION
