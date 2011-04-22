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
#include "nuto/mechanics/constraints/ConstraintNonlinearGlobalCrackAngle2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/structures/unstructured/StructureMultiscale.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRVector2Symmetric.h"

// constructor
NuTo::ConstraintNonlinearGlobalCrackAngle2D::ConstraintNonlinearGlobalCrackAngle2D(const StructureMultiscale* rStructure,
        double rPenaltyStiffness, double rScalingFactor):
        ConstraintNonlinear()
{
    mStructure = rStructure;
    mPenaltyStiffness = rPenaltyStiffness;
    mScalingFactor = rScalingFactor;
}

NuTo::ConstraintNonlinearGlobalCrackAngle2D* NuTo::ConstraintNonlinearGlobalCrackAngle2D::AsConstraintNonlinearGlobalCrackAngle2D()
{
    return this;
}

const NuTo::ConstraintNonlinearGlobalCrackAngle2D* NuTo::ConstraintNonlinearGlobalCrackAngle2D::AsConstraintNonlinearGlobalCrackAngle2D()const
{
    return this;
}

void NuTo::ConstraintNonlinearGlobalCrackAngle2D::SetPenaltyStiffness(double rPenaltyStiffness)
{
    mPenaltyStiffness = rPenaltyStiffness;
}

void NuTo::ConstraintNonlinearGlobalCrackAngle2D::SetScalingFactor(double rScalingFactor)
{
    mScalingFactor = rScalingFactor;
}

//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void NuTo::ConstraintNonlinearGlobalCrackAngle2D::Info(unsigned short rVerboseLevel) const
{
    double anglePrescribed(0);
    // calcute principal direction of total strain
    EngineeringStrain2D totalStrain(mStructure->GetTotalEngineeringStrain());
    if (totalStrain.mEngineeringStrain[0]-totalStrain.mEngineeringStrain[1]!=0.)
        anglePrescribed = atan(totalStrain.mEngineeringStrain[2]/(totalStrain.mEngineeringStrain[0]-totalStrain.mEngineeringStrain[1]));
    else
    {
        if (totalStrain.mEngineeringStrain[2]>0)
            anglePrescribed = 1.5*M_PI;
        else
            anglePrescribed = 0.5*M_PI;
    }

    //std::cout << "NuTo::ConstraintLinearGlobalCrackAngle : prescribed angle " <<  anglePrescribed << std::endl;
}

//! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the stiffness matrix (Hessian) of the total potential
//! @param rResult ... coefficient matrix
//! @param rGlobalDofs ... row and column numbers in global system
void NuTo::ConstraintNonlinearGlobalCrackAngle2D::CalculateCoefficientMatrix_0(NuTo::SparseMatrixCSRVector2Symmetric<double>& rResult,
        std::vector<int>& rGlobalDofs)const
{
    //crack angle
    int dof(4);
    rResult.Resize(dof,dof);
    rGlobalDofs.resize(dof,1);

    rGlobalDofs[0] = mStructure->GetDofCrackAngle();

    //calculate angle orthogonal to second principal stress
    double delta_alpha = mStructure->CalculateDeltaCrackAngleElastic();
    double ratioSquare(delta_alpha/mScalingFactor);
    ratioSquare*=ratioSquare;
    rResult.AddEntry(0,0, mPenaltyStiffness*exp(-ratioSquare));
}

//! @brief calculates the gradient of the internal potential
//! for a mechanical problem, this corresponds to the internal force vector
void NuTo::ConstraintNonlinearGlobalCrackAngle2D::CalculateGradientInternalPotential(NuTo::FullMatrix<double>& rResult,
        std::vector<int>& rGlobalDofs)const
{
    //crack angle
    int dof(1);
    rResult.Resize(dof,1);
    rGlobalDofs.resize(dof,1);

    rGlobalDofs[0] = mStructure->GetDofCrackAngle();

    //calculate angle orthogonal to second principal stress
    double delta_alpha = mStructure->CalculateDeltaCrackAngleElastic();
    rResult(0,0) =  mPenaltyStiffness*mScalingFactor*sqrt(M_PI)*0.5*erf(delta_alpha/mScalingFactor);
}

//! @brief calculates the internal potential
double NuTo::ConstraintNonlinearGlobalCrackAngle2D::CalculateTotalPotential()const
{
    //calculate angle orthogonal to second principal stress
    double delta_alpha = mStructure->CalculateDeltaCrackAngleElastic();
    double ratio(delta_alpha/mScalingFactor);
    double ratioSquare=ratio*ratio;
    return  mPenaltyStiffness*mScalingFactor*mScalingFactor*sqrt(M_PI)*0.5*
            (erf(ratio)*ratio + exp(-ratioSquare)/sqrt(M_PI));
}


#ifdef ENABLE_SERIALIZATION
// serialize
template void NuTo::ConstraintNonlinearGlobalCrackAngle2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintNonlinearGlobalCrackAngle2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintNonlinearGlobalCrackAngle2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintNonlinearGlobalCrackAngle2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintNonlinearGlobalCrackAngle2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintNonlinearGlobalCrackAngle2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintNonlinearGlobalCrackAngle2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintNonLinearGlobalCrackAngle2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintNonlinear)
       & BOOST_SERIALIZATION_NVP(mPenaltyStiffness)
       & BOOST_SERIALIZATION_NVP(mScalingFactor)
       & BOOST_SERIALIZATION_NVP(mStructure);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintNonLinearGlobalCrackAngle2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintNonlinearGlobalCrackAngle2D)
#endif // ENABLE_SERIALIZATION
