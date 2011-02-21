// $Id: ConstraintGlobalCrackAngle.cpp 314 2010-09-27 16:31:43Z unger3 $

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
#include "nuto/mechanics/constraints/ConstraintNonlinearGlobalCrackOpeningTangential2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/structures/unstructured/StructureIp.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRVector2Symmetric.h"

// constructor
NuTo::ConstraintNonlinearGlobalCrackOpeningTangential2D::ConstraintNonlinearGlobalCrackOpeningTangential2D(const StructureIp* rStructure,
        double rScalingFactor, double rPenaltyStiffness) : ConstraintNonlinear()
{
    mStructure = rStructure;
    mScalingFactor = rScalingFactor;
    mPenaltyStiffness = rPenaltyStiffness;
}

NuTo::ConstraintNonlinearGlobalCrackOpeningTangential2D* NuTo::ConstraintNonlinearGlobalCrackOpeningTangential2D::AsConstraintNonlinearGlobalCrackOpeningTangential2D()
{
    return this;
}

const NuTo::ConstraintNonlinearGlobalCrackOpeningTangential2D* NuTo::ConstraintNonlinearGlobalCrackOpeningTangential2D::AsConstraintNonlinearGlobalCrackOpeningTangential2D()const
{
    return this;
}

void NuTo::ConstraintNonlinearGlobalCrackOpeningTangential2D::SetPenaltyStiffness(double rPenaltyStiffness)
{
    mPenaltyStiffness = rPenaltyStiffness;
}

void NuTo::ConstraintNonlinearGlobalCrackOpeningTangential2D::SetScalingFactor(double rScalingFactor)
{
    mScalingFactor = rScalingFactor;
}

//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void NuTo::ConstraintNonlinearGlobalCrackOpeningTangential2D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "NuTo::ConstraintNonlinearGlobalCrackOpeningTangential2D " << std::endl;
}

//! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the stiffness matrix (Hessian) of the total potential
//! @param rResult ... coefficient matrix
//! @param rGlobalDofs ... row and column numbers in global system
void NuTo::ConstraintNonlinearGlobalCrackOpeningTangential2D::CalculateCoefficientMatrix_0(NuTo::SparseMatrixCSRVector2Symmetric<double>& rResult,
        std::vector<int>& rGlobalDofs)const
{
    //crack angle
    int dof(1);
    rResult.Resize(dof,dof);
    rGlobalDofs.resize(dof,1);

    rGlobalDofs[0] = mStructure->GetDofGlobalCrackOpening2D()[0];
    double ut(mStructure->GetGlobalCrackOpening2D()[0]);
    double ratioSquare(ut/mScalingFactor);
    ratioSquare*=ratioSquare;
    rResult.AddEntry(0,0, mPenaltyStiffness*exp(-ratioSquare));
}

//! @brief calculates the gradient of the internal potential
//! for a mechanical problem, this corresponds to the internal force vector
void NuTo::ConstraintNonlinearGlobalCrackOpeningTangential2D::CalculateGradientInternalPotential(NuTo::FullMatrix<double>& rResult,
        std::vector<int>& rGlobalDofs)const
{
    //crack angle
    int dof(1);
    rResult.Resize(dof,1);
    rGlobalDofs.resize(dof,1);

    rGlobalDofs[0] = mStructure->GetDofGlobalCrackOpening2D()[0];
    double ut(mStructure->GetGlobalCrackOpening2D()[0]);
    double ratioSquare(ut/mScalingFactor);
    ratioSquare*=ratioSquare;

    //calculate angle orthogonal to second principal stress
    rResult(0,0) =  mPenaltyStiffness*mScalingFactor*sqrt(M_PI)*0.5*erf(ut/mScalingFactor);
}

//! @brief calculates the internal potential
double NuTo::ConstraintNonlinearGlobalCrackOpeningTangential2D::CalculateTotalPotential()const
{
    //calculate angle orthogonal to second principal stress
    double ut(mStructure->GetGlobalCrackOpening2D()[0]);
    double ratio(ut/mScalingFactor);
    double ratioSquare=ratio*ratio;
    return  mPenaltyStiffness*mScalingFactor*mScalingFactor*sqrt(M_PI)*0.5*
            (erf(ratio)*ratio + exp(-ratioSquare)/sqrt(M_PI));

}


#ifdef ENABLE_SERIALIZATION
// serialize
template void NuTo::ConstraintNonlinearGlobalCrackOpeningTangential2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintNonlinearGlobalCrackOpeningTangential2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintNonlinearGlobalCrackOpeningTangential2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintNonlinearGlobalCrackOpeningTangential2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintNonlinearGlobalCrackOpeningTangential2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintNonlinearGlobalCrackOpeningTangential2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintNonlinearGlobalCrackOpeningTangential2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintNonlinearGlobalCrackOpeningTangential2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintBase)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintNonlinear)
       & BOOST_SERIALIZATION_NVP(mStructure)
       & BOOST_SERIALIZATION_NVP(mScalingFactor)
       & BOOST_SERIALIZATION_NVP(mPenaltyStiffness);

#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintNonlinearGlobalCrackOpeningTangential2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintNonlinearGlobalCrackOpeningTangential2D)
#endif // ENABLE_SERIALIZATION
