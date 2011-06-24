// $Id: ConstraintLagrangeGlobalCrackOpening2D.cpp -1   $

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
#include "nuto/mechanics/constraints/ConstraintLagrangeGlobalCrackOpening2D.h"
#include "nuto/mechanics/structures/unstructured/StructureMultiscale.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRVector2Symmetric.h"

NuTo::ConstraintLagrangeGlobalCrackOpening2D::ConstraintLagrangeGlobalCrackOpening2D(const NuTo::StructureMultiscale* rStructure, double rPenaltyStiffness) :
        ConstraintLagrange(NuTo::Constraint::EQUAL, rPenaltyStiffness)
{
    mStructure = rStructure;
    mLagrangeValue = 0;
    mLagrangeDOF = -1;
}

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintLagrangeGlobalCrackOpening2D::GetNumLagrangeMultipliers()const
{
    return 1;
}

//! @brief returns the Lagrange Multiplier
//! first col Lagrange, second column slack variables
void NuTo::ConstraintLagrangeGlobalCrackOpening2D::GetLagrangeMultiplier(FullMatrix<double>& rLagrangeMultiplier)const
{
    rLagrangeMultiplier.Resize(1,1);
    rLagrangeMultiplier(0,0) = mLagrangeValue;
}

//! @brief returns the Lagrange Multiplier dofs
//! first col Lagrangedofs
void NuTo::ConstraintLagrangeGlobalCrackOpening2D::GetDofsLagrangeMultiplier(FullMatrix<int>& rLagrangeMultiplier)const
{
    rLagrangeMultiplier.Resize(1,1);
    rLagrangeMultiplier(0,0) = mLagrangeDOF;
}

//! @brief sets the global dofs
//! @param rDOF current maximum DOF, this variable is increased within the routine
void NuTo::ConstraintLagrangeGlobalCrackOpening2D::SetGlobalDofs(int& rDOF)
{
    mLagrangeDOF = rDOF++;
}

//! @brief write dof values to constraints (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::ConstraintLagrangeGlobalCrackOpening2D::SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);

    double value;
    if (mLagrangeDOF >= rActiveDofValues.GetNumRows())
    {
        int dof = mLagrangeDOF - rActiveDofValues.GetNumRows();
        assert(dof < rDependentDofValues.GetNumRows());
        value = rDependentDofValues(dof,0);
    }
    else
    {
        value = rActiveDofValues(mLagrangeDOF,0);
    }
    this->mLagrangeValue = value;
}

//! @brief extract dof values from the constraints (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::ConstraintLagrangeGlobalCrackOpening2D::GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);

    if (mLagrangeDOF >= rActiveDofValues.GetNumRows())
    {
        int dof = mLagrangeDOF - rActiveDofValues.GetNumRows();
        assert(dof < rDependentDofValues.GetNumRows());
        rDependentDofValues(dof,0) = mLagrangeValue;
    }
    else
    {
        rActiveDofValues(mLagrangeDOF,0) = mLagrangeValue;
    }
}

//! @brief renumber the global dofs according to predefined ordering
//! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
void NuTo::ConstraintLagrangeGlobalCrackOpening2D::RenumberGlobalDofs(const std::vector<int>& rMappingInitialToNewOrdering)
{
    mLagrangeDOF = rMappingInitialToNewOrdering[mLagrangeDOF];
}

//! @brief cast to Lagrange constraint - the corresponding dofs are eliminated in the global system
NuTo::ConstraintLagrange* NuTo::ConstraintLagrangeGlobalCrackOpening2D::AsConstraintLagrange()
{
    return this;
}

//! @brief cast to Lagrange constraint - the corresponding dofs are eliminated in the global system
const NuTo::ConstraintLagrange* NuTo::ConstraintLagrangeGlobalCrackOpening2D::AsConstraintLagrange()const
{
    return this;
}

//! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the stiffness matrix (Hessian) of the total potential
//! @param rResult ... coefficient matrix
//! @param rGlobalDofs ... row and column numbers in global system
void NuTo::ConstraintLagrangeGlobalCrackOpening2D::CalculateCoefficientMatrix_0(NuTo::SparseMatrixCSRVector2Symmetric<double>& rResult,
        std::vector<int>& rGlobalDofs)const
{
    //Lagrange mult + normal dof
    int dof(2);
    rResult.Resize(dof,dof);
    rGlobalDofs.resize(dof,1);

    boost::array<int,2> dofCrackOpening(mStructure->GetDofGlobalCrackOpening2D());
    boost::array<double,2> crackOpening(mStructure->GetGlobalCrackOpening2D());

    rGlobalDofs[0] = mLagrangeDOF;
    rGlobalDofs[1] = dofCrackOpening[1];

    if (-crackOpening[1]>-mLagrangeValue/mPenalty)
    {
        //derivative with respect to un and lambda
        rResult.AddEntry(0,1,-mStructure->GetScalingFactorCrackOpening());
        //derivative with respect to un^2
        rResult.AddEntry(1,1,mPenalty*mStructure->GetScalingFactorCrackOpening()*mStructure->GetScalingFactorCrackOpening());
    }
    else
    {
        //derivative with respect to lambda^2
        rResult.AddEntry(0,0,-1./mPenalty);
    }
/*
    //check result matrix by calling the gradient routine
    NuTo::FullMatrix<double> rResultFull(rResult);
    NuTo::FullMatrix<double> rResultCDF(rResult);
    NuTo::FullMatrix<double> gradient1,gradient2;
    double delta(-1e-8);
    CalculateGradientInternalPotential(gradient1,rGlobalDofs);

    const_cast<ConstraintLagrangeGlobalCrackOpening2D*>(this)->mLagrangeValue+=delta;
    CalculateGradientInternalPotential(gradient2,rGlobalDofs);
    rResultCDF.SetColumn(0,(gradient2-gradient1)*(1./delta));
    const_cast<ConstraintLagrangeGlobalCrackOpening2D*>(this)->mLagrangeValue-=delta;

    crackOpening[1]+=delta;
    mStructure->SetGlobalCrackOpening(crackOpening);
    CalculateGradientInternalPotential(gradient2,rGlobalDofs);
    rResultCDF.SetColumn(1,(gradient2-gradient1)*(mStructure->GetScalingFactorCrackOpening()/delta));
    crackOpening[1]-=delta;
    mStructure->SetGlobalCrackOpening(crackOpening);


    if ((rResultFull-rResultCDF).Abs().Max()<1e-2)
    {
        std::cout << "ConstraintLagrangeGlobalCrackOpening2D constraint crack opening "  << std::endl;
        std::cout<< "stiffness algo " << std::endl;
        rResultFull.Info(12,10);
        std::cout<< "stiffness cdf " << std::endl;
        rResultCDF.Info(12,10);
    }
    else
    {
        std::cout << "hessian constraint crack opening is fine " << std::endl;
    }
*/
}

//! @brief calculates the gradient of the internal potential
//! for a mechanical problem, this corresponds to the internal force vector
void NuTo::ConstraintLagrangeGlobalCrackOpening2D::CalculateGradientInternalPotential(NuTo::FullMatrix<double>& rResult,
        std::vector<int>& rGlobalDofs)const
{
    int dof(2);
    rResult.Resize(dof,1);
    rGlobalDofs.resize(dof);

    boost::array<int,2> dofCrackOpening(mStructure->GetDofGlobalCrackOpening2D());
    boost::array<double,2> crackOpening(mStructure->GetGlobalCrackOpening2D());

    rGlobalDofs[0] = mLagrangeDOF;
    rGlobalDofs[1] = dofCrackOpening[1];

    if (-crackOpening[1]>-mLagrangeValue/mPenalty)
    {
        //derivative with respect to lambda
        rResult(0,0)=-crackOpening[1];
        //derivative with respect to un
        rResult(1,0)=(-mLagrangeValue-mPenalty*(-crackOpening[1]))*mStructure->GetScalingFactorCrackOpening();
    }
    else
    {
        //derivative with respect to lambda
        rResult(0,0)=-mLagrangeValue/mPenalty;
    }

/*
    //check result matrix by calling the potential routine
    NuTo::FullMatrix<double> rResultFull(rResult);
    NuTo::FullMatrix<double> rResultCDF(rResult);
    double energy1,energy2;
    double delta(-1e-8);
    energy1 = CalculateTotalPotential();

    const_cast<ConstraintLagrangeGlobalCrackOpening2D*>(this)->mLagrangeValue+=delta;
    energy2 = CalculateTotalPotential();
    rResultCDF(0,0) = (energy2-energy1)*(1./delta);
    const_cast<ConstraintLagrangeGlobalCrackOpening2D*>(this)->mLagrangeValue-=delta;

    crackOpening[1]+=delta;
    mStructure->SetGlobalCrackOpening(crackOpening);
    energy2 = CalculateTotalPotential();
    rResultCDF(1,0) = (energy2-energy1)*(mStructure->GetScalingFactorCrackOpening()/delta);
    crackOpening[1]-=delta;
    mStructure->SetGlobalCrackOpening(crackOpening);

    //if ((rResultFull-rResultCDF).Abs().Max()>1e-2)
    if (true)
    {
        std::cout<< "constraint crack opening gradient algo " << std::endl;
        rResultFull.Info(12,10);
        std::cout<< "constraint crack opening gradient cdf " << std::endl;
        rResultCDF.Info(12,10);
    }
    else
    {
        std::cout << "gradient constraint crack opening is fine " << std::endl;
    }
*/
}

//! @brief calculates the internal potential
double NuTo::ConstraintLagrangeGlobalCrackOpening2D::CalculateTotalPotential()const
{
    boost::array<double,2> crackOpening(mStructure->GetGlobalCrackOpening2D());

    double g(-crackOpening[1]);

    if (g>-mLagrangeValue/mPenalty)
    {
        return g*(mLagrangeValue+0.5*mPenalty*g);
    }
    else
    {
        g=-mLagrangeValue/mPenalty;
        return g*(mLagrangeValue+0.5*mPenalty*g);
    }
}

//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void NuTo::ConstraintLagrangeGlobalCrackOpening2D::Info(unsigned short rVerboseLevel) const
{
    boost::array<int,2> dofCrackOpening(mStructure->GetDofGlobalCrackOpening2D());
    boost::array<double,2> crackOpening(mStructure->GetGlobalCrackOpening2D());

    double g(-crackOpening[1]);
    mStructure->GetLogger() << "ConstraintLagrangeGlobalCrackOpening2D : lambda " << mLagrangeValue << "(" <<  mLagrangeDOF << ")" << " g " << g <<
            " crackOpening " << crackOpening[0] << " " << crackOpening[1] << "\n";
    if (g>-mLagrangeValue/mPenalty)
    {
        std::cout << "\t\t constraint crack opening active " << std::endl;
    }
    else
    {
        std::cout << "\t\t constraint crack opening inactive " << std::endl;
    }
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ConstraintLagrangeGlobalCrackOpening2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrangeGlobalCrackOpening2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrangeGlobalCrackOpening2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrangeGlobalCrackOpening2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrangeGlobalCrackOpening2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrangeGlobalCrackOpening2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintLagrangeGlobalCrackOpening2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLagrangeGlobalCrackOpening2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintBase)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintLagrange)
       & BOOST_SERIALIZATION_NVP(mStructure)
       & BOOST_SERIALIZATION_NVP(mLagrangeValue)
       & BOOST_SERIALIZATION_NVP(mLagrangeDOF);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLagrangeGlobalCrackOpening2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintLagrangeGlobalCrackOpening2D)
#endif // ENABLE_SERIALIZATION
