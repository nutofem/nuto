// $Id: ConstraintLagrangeGlobalCrackAngle2D.cpp -1   $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constraints/ConstraintLagrangeGlobalCrackAngle2D.h"
#include "nuto/mechanics/structures/unstructured/StructureIp.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRVector2Symmetric.h"

#define UMAX 10.1 //0.03mm

NuTo::ConstraintLagrangeGlobalCrackAngle2D::ConstraintLagrangeGlobalCrackAngle2D(const NuTo::StructureIp* rStructure) :
        ConstraintLagrange(NuTo::Constraint::EQUAL)
{
    mStructure = rStructure;
    mLagrangeValue = 0;
    mLagrangeDOF = -1;
    std::cout<< "scaling factor for ConstraintLagrangeGlobalCrackAngle2D has to be adapted" << std::endl;
    mScaling = 2.*M_PI/(UMAX);
    mScaling = 0.;
}

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintLagrangeGlobalCrackAngle2D::GetNumLagrangeMultipliers()const
{
    return 1;
}

//! @brief returns the Lagrange Multiplier
//! first col Lagrange, second column slack variables
void NuTo::ConstraintLagrangeGlobalCrackAngle2D::GetLagrangeMultiplier(FullMatrix<double>& rLagrangeMultiplier)const
{
    rLagrangeMultiplier.Resize(1,1);
    rLagrangeMultiplier(0,0) = mLagrangeValue;
}

//! @brief sets the global dofs
//! @param rDOF current maximum DOF, this variable is increased within the routine
void NuTo::ConstraintLagrangeGlobalCrackAngle2D::SetGlobalDofs(int& rDOF)
{
    mLagrangeDOF = rDOF++;
}

//! @brief write dof values to constraints (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::ConstraintLagrangeGlobalCrackAngle2D::SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
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
void NuTo::ConstraintLagrangeGlobalCrackAngle2D::GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
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
void NuTo::ConstraintLagrangeGlobalCrackAngle2D::RenumberGlobalDofs(const std::vector<int>& rMappingInitialToNewOrdering)
{
    mLagrangeDOF = rMappingInitialToNewOrdering[mLagrangeDOF];
}

//! @brief cast to Lagrange constraint - the corresponding dofs are eliminated in the global system
NuTo::ConstraintLagrange* NuTo::ConstraintLagrangeGlobalCrackAngle2D::AsConstraintLagrange()
{
    return this;
}

//! @brief cast to Lagrange constraint - the corresponding dofs are eliminated in the global system
const NuTo::ConstraintLagrange* NuTo::ConstraintLagrangeGlobalCrackAngle2D::AsConstraintLagrange()const
{
    return this;
}

//! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the stiffness matrix (Hessian) of the total potential
//! @param rResult ... coefficient matrix
//! @param rGlobalDofs ... row and column numbers in global system
void NuTo::ConstraintLagrangeGlobalCrackAngle2D::CalculateCoefficientMatrix_0(NuTo::SparseMatrixCSRVector2Symmetric<double>& rResult,
        std::vector<int>& rGlobalDofs)const
{
    //Lagrange mult + alpha+ 2 crack opening
    int dof(4);
    rResult.Resize(dof,dof);
    rGlobalDofs.resize(dof,1);

    boost::array<int,2> dofCrackOpening(mStructure->GetDofGlobalCrackOpening2D());
    boost::array<double,2> crackOpening(mStructure->GetGlobalCrackOpening2D());

    rGlobalDofs[0] = mLagrangeDOF;
    rGlobalDofs[1] = mStructure->GetDofCrackAngle();
    rGlobalDofs[2] = dofCrackOpening[0];
    rGlobalDofs[3] = dofCrackOpening[1];

    double alpha(mStructure->GetCrackAngle());
    //calculate angle orthogonal to second principal stress
    double alpha_2 = mStructure->CalculateCrackAngleElastic();

    double delta_alpha = alpha-alpha_2;
    while (fabs(delta_alpha>M_PI))
    {
        if (delta_alpha>0)
        {
            delta_alpha-=M_PI;
            alpha_2+=M_PI;
        }
        else
        {
            delta_alpha+=M_PI;
            alpha_2-=M_PI;
        }
    }
    //now delta_alpha is in [0..PI]

    if (delta_alpha>M_PI-delta_alpha)
    {
        alpha_2 += M_PI;
    }

    double abs_delta_alpha = alpha-alpha_2;
    std::cout << abs_delta_alpha << " " << 0.5*M_PI << std::endl;
    assert(fabs(abs_delta_alpha)<=0.50001*M_PI);
    double sign_delta_alpha(1.);
    if (delta_alpha<0)
    {
        abs_delta_alpha=fabs(abs_delta_alpha);
        sign_delta_alpha = -1;
    }
    double signU0 = crackOpening[0] < 0 ? -1 : 1;
    double signU1 = crackOpening[1] < 0 ? -1 : 1;
    double g = abs_delta_alpha-mScaling*(signU0*crackOpening[0]+signU1*crackOpening[1]);
    //double g = abs_delta_alpha-mScaling*(signU1*crackOpening[1]);

    //std::cout << "alpha_2 " << alpha_2 << " alpha "<< alpha << " abs_delta_alpha " <<  abs_delta_alpha << " sign " << sign_delta_alpha << std::endl;

    if (g>=-mLagrangeValue/mPenalty)
    {
        //std::cout << "constraint crack angle active(stiffness) " << g << std::endl;
        //derivative with respect to lambda and alpha
        rResult.AddEntry(0,1,sign_delta_alpha);
        //derivative with respect to lambda and u0
        rResult.AddEntry(0,2,-mScaling*signU0);
        //derivative with respect to lambda and u1
        rResult.AddEntry(0,3,-mScaling*signU1);
        //derivative with respect to alpha^2
        rResult.AddEntry(1,1,mPenalty);
        //derivative with respect to alpha and u0
        rResult.AddEntry(1,2,-mPenalty*mScaling*sign_delta_alpha*signU0);
        //derivative with respect to alpha and u1
        rResult.AddEntry(1,3,-mPenalty*mScaling*sign_delta_alpha*signU1);
        //derivative with respect to u0^2
        rResult.AddEntry(2,2,mPenalty*mScaling*mScaling);
        //derivative with respect to u0 and u1
        rResult.AddEntry(2,3,mPenalty*mScaling*mScaling*signU0*signU1);
        //derivative with respect to uy^2
        rResult.AddEntry(3,3,mPenalty*mScaling*mScaling);
    }
    else
    {
        //std::cout << "constraint crack angle inactive(stiffness) " << g << std::endl;
        //derivative with respect to lambda^2
        rResult.AddEntry(0,0,-1./mPenalty);
    }

    FullMatrix<double> rResultFull(rResult);

    //check result matrix by calling the gradient routine
    NuTo::FullMatrix<double> rResultCDF(rResult);
    NuTo::FullMatrix<double> gradient1,gradient2;
    double delta(1e-8);
    CalculateGradientInternalPotential(gradient1,rGlobalDofs);

    const_cast<ConstraintLagrangeGlobalCrackAngle2D*>(this)->mLagrangeValue+=delta;
    CalculateGradientInternalPotential(gradient2,rGlobalDofs);
    rResultCDF.SetColumn(0,(gradient2-gradient1)*(1./delta));
    const_cast<ConstraintLagrangeGlobalCrackAngle2D*>(this)->mLagrangeValue-=delta;

    alpha+=delta;
    mStructure->SetCrackAngle(alpha);
    CalculateGradientInternalPotential(gradient2,rGlobalDofs);
    rResultCDF.SetColumn(1,(gradient2-gradient1)*(1./delta));
    alpha-=delta;
    mStructure->SetCrackAngle(alpha);

    crackOpening[0]+=delta;
    mStructure->SetGlobalCrackOpening(crackOpening);
    CalculateGradientInternalPotential(gradient2,rGlobalDofs);
    rResultCDF.SetColumn(2,(gradient2-gradient1)*(1./delta));
    crackOpening[0]-=delta;
    mStructure->SetGlobalCrackOpening(crackOpening);

    crackOpening[1]+=delta;
    mStructure->SetGlobalCrackOpening(crackOpening);
    CalculateGradientInternalPotential(gradient2,rGlobalDofs);
    rResultCDF.SetColumn(3,(gradient2-gradient1)*(1./delta));
    crackOpening[1]-=delta;
    mStructure->SetGlobalCrackOpening(crackOpening);

    if ((rResultFull-rResultCDF).Abs().Max()>1e-2)
    {
        std::cout << "constraint matrix alpha" << std::endl;
        rResultFull.Info(12,5);
        std::cout<< "stiffness cdf " << std::endl;
        rResultCDF.Info(10,3);
    }
    else
    {
        std::cout << "constraint crack angle is fine " << std::endl;
    }
}

//! @brief calculates the gradient of the internal potential
//! for a mechanical problem, this corresponds to the internal force vector
void NuTo::ConstraintLagrangeGlobalCrackAngle2D::CalculateGradientInternalPotential(NuTo::FullMatrix<double>& rResult,
        std::vector<int>& rGlobalDofs)const
{
    int dof(4);
    rResult.Resize(dof,1);
    rGlobalDofs.resize(dof);

    boost::array<int,2> dofCrackOpening(mStructure->GetDofGlobalCrackOpening2D());
    boost::array<double,2> crackOpening(mStructure->GetGlobalCrackOpening2D());

    rGlobalDofs[0] = mLagrangeDOF;
    rGlobalDofs[1] = mStructure->GetDofCrackAngle();
    rGlobalDofs[2] = dofCrackOpening[0];
    rGlobalDofs[3] = dofCrackOpening[1];

    double alpha(mStructure->GetCrackAngle());
    //calculate angle orthogonal to second principal stress
    double alpha_2 = mStructure->CalculateCrackAngleElastic();

    double delta_alpha = alpha-alpha_2;
    while (fabs(delta_alpha>M_PI))
    {
        if (delta_alpha>0)
        {
            delta_alpha-=M_PI;
            alpha_2+=M_PI;
        }
        else
        {
            delta_alpha+=M_PI;
            alpha_2-=M_PI;
        }
    }
    //now delta_alpha is in [0..PI]

    if (delta_alpha>0.5*M_PI)
    {
        alpha_2 += M_PI;
    }

    double abs_delta_alpha = alpha-alpha_2;
    assert(fabs(abs_delta_alpha)<0.50001*M_PI);
    double sign_delta_alpha(1.);
    if (delta_alpha<0)
    {
        abs_delta_alpha=fabs(abs_delta_alpha);
        sign_delta_alpha = -1;
    }
    double signU0 = crackOpening[0] < 0 ? -1 : 1;
    double signU1 = crackOpening[1] < 0 ? -1 : 1;
    double g = abs_delta_alpha-mScaling*(signU0*crackOpening[0]+signU1*crackOpening[1]);
    //double g = abs_delta_alpha-mScaling*(signU1*crackOpening[1]);

    if (g>=-mLagrangeValue/mPenalty)
    {
        //std::cout << "constraint crack angle active(gradient) " << g << "  " << mLagrangeValue << std::endl;
        //derivative with respect to lambda
        rResult(0,0)= g;
        //derivative with respect to alpha
        rResult(1,0)= sign_delta_alpha*(mLagrangeValue + mPenalty*g);
        //derivative with respect to u0 (tangential)
        rResult(2,0)= -mScaling*signU0*(mLagrangeValue+mPenalty*g);
        //derivative with respect to u1 (normal)
        rResult(3,0)= -mScaling*signU1*(mLagrangeValue+mPenalty*g);
    }
    else
    {
        //derivative with respect to lambda
        rResult(0,0)=-mLagrangeValue/mPenalty;
    }

    //check result matrix by calling the potential routine
    NuTo::FullMatrix<double> rResultFull(rResult);
    NuTo::FullMatrix<double> rResultCDF(rResult);
    double energy1,energy2;
    double delta(1e-8);
    energy1 = CalculateTotalPotential();

    const_cast<ConstraintLagrangeGlobalCrackAngle2D*>(this)->mLagrangeValue+=delta;
    energy2 = CalculateTotalPotential();
    rResultCDF(0,0) = (energy2-energy1)*(1./delta);
    const_cast<ConstraintLagrangeGlobalCrackAngle2D*>(this)->mLagrangeValue-=delta;

    alpha+=delta;
    mStructure->SetCrackAngle(alpha);
    energy2 = CalculateTotalPotential();
    rResultCDF(1,0) = (energy2-energy1)*(1./delta);
    alpha-=delta;
    mStructure->SetCrackAngle(alpha);

    crackOpening[0]+=delta;
    mStructure->SetGlobalCrackOpening(crackOpening);
    energy2 = CalculateTotalPotential();
    rResultCDF(2,0) = (energy2-energy1)*(1./delta);
    crackOpening[0]-=delta;
    mStructure->SetGlobalCrackOpening(crackOpening);

    crackOpening[1]+=delta;
    mStructure->SetGlobalCrackOpening(crackOpening);
    energy2 = CalculateTotalPotential();
    rResultCDF(3,0) = (energy2-energy1)*(1./delta);
    crackOpening[1]-=delta;
    mStructure->SetGlobalCrackOpening(crackOpening);

    if ((rResultCDF-rResultFull).Abs().Max()>1e-2)
    {
        std::cout<< "constraint crack angle gradient algo " << std::endl;
        rResultFull.Info(10,5);

        std::cout<< "constraint crack angle gradient cdf " << std::endl;
        rResultCDF.Info(10,5);
    }
}

//! @brief calculates the internal potential
double NuTo::ConstraintLagrangeGlobalCrackAngle2D::CalculateTotalPotential()const
{
    boost::array<int,2> dofCrackOpening(mStructure->GetDofGlobalCrackOpening2D());
    boost::array<double,2> crackOpening(mStructure->GetGlobalCrackOpening2D());

    double alpha(mStructure->GetCrackAngle());
    //calculate angle orthogonal to second principal stress
    double alpha_2 = mStructure->CalculateCrackAngleElastic();

    double delta_alpha = alpha-alpha_2;
    while (fabs(delta_alpha>M_PI))
    {
        if (delta_alpha>0)
        {
            delta_alpha-=M_PI;
            alpha_2+=M_PI;
        }
        else
        {
            delta_alpha+=M_PI;
            alpha_2-=M_PI;
        }
    }
    //now delta_alpha is in [0..PI]

    if (delta_alpha>0.5*M_PI)
    {
        alpha_2 += M_PI;
    }

    double abs_delta_alpha = alpha-alpha_2;
    assert(fabs(abs_delta_alpha)<0.50001*M_PI);
    double sign_delta_alpha(1.);
    if (delta_alpha<0)
    {
        abs_delta_alpha=fabs(abs_delta_alpha);
        sign_delta_alpha = -1;
    }
    double signU0 = crackOpening[0] < 0 ? -1 : 1;
    double signU1 = crackOpening[1] < 0 ? -1 : 1;
    double g = abs_delta_alpha-mScaling*(signU0*crackOpening[0]+signU1*crackOpening[1]);
//    double g = abs_delta_alpha-mScaling*(crackOpening[1]);

    if (g>=-mLagrangeValue/mPenalty)
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
void NuTo::ConstraintLagrangeGlobalCrackAngle2D::Info(unsigned short rVerboseLevel) const
{
    boost::array<int,2> dofCrackOpening(mStructure->GetDofGlobalCrackOpening2D());
    boost::array<double,2> crackOpening(mStructure->GetGlobalCrackOpening2D());

    double alpha(mStructure->GetCrackAngle());

    double alpha_2 = mStructure->CalculateCrackAngleElastic();

    double delta_alpha = alpha-alpha_2;
    while (fabs(delta_alpha>M_PI))
    {
        if (delta_alpha>0)
        {
            delta_alpha-=M_PI;
            alpha_2+=M_PI;
        }
        else
        {
            delta_alpha+=M_PI;
            alpha_2-=M_PI;
        }
    }
    //now delta_alpha is in [0..PI]
    double signU0 = crackOpening[0] < 0 ? -1 : 1;
    double signU1 = crackOpening[1] < 0 ? -1 : 1;

    if (delta_alpha>0.5*M_PI)
    {
        alpha_2 += M_PI;
    }

    double abs_delta_alpha = alpha-alpha_2;
    assert(fabs(abs_delta_alpha)<0.50001*M_PI);
    double sign_delta_alpha(1.);
    if (delta_alpha<0)
    {
        abs_delta_alpha=fabs(abs_delta_alpha);
        sign_delta_alpha = -1;
    }

    double g = abs_delta_alpha-mScaling*(signU0*crackOpening[0]+signU1*crackOpening[1]);
    //    double g = abs_delta_alpha-mScaling*(crackOpening[1]);

    std::cout << "CrackAngle : lambda " << mLagrangeValue << "(" <<  mLagrangeDOF << ")" <<" g " << g << std::endl;
    if (g>=-mLagrangeValue/mPenalty)
    {
        std::cout << "constraint crack angle active ";
    }
    else
    {
        std::cout << "constraint crack angle inactive ";
    }
    std::cout << "alpha_2 " << alpha_2 << " alpha "<< alpha << " abs_delta_alpha " <<  abs_delta_alpha << " sign " << sign_delta_alpha << std::endl;
}


#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ConstraintLagrangeGlobalCrackAngle2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrangeGlobalCrackAngle2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrangeGlobalCrackAngle2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrangeGlobalCrackAngle2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrangeGlobalCrackAngle2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrangeGlobalCrackAngle2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintLagrangeGlobalCrackAngle2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLagrangeGlobalCrackAngle2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintLagrange)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintBase)
       & BOOST_SERIALIZATION_NVP(mStructure)
       & BOOST_SERIALIZATION_NVP(mLagrangeValue)
       & BOOST_SERIALIZATION_NVP(mLagrangeDOF)
       & BOOST_SERIALIZATION_NVP(mScaling);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLagrangeGlobalCrackAngle2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintLagrangeGlobalCrackAngle2D)
#endif // ENABLE_SERIALIZATION
