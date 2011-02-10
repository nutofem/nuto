// $Id: ConstraintLagrangeGlobalCrackOpening2D.cpp -1   $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constraints/ConstraintLagrangeGlobalCrackOpening2D.h"
#include "nuto/mechanics/structures/unstructured/StructureIp.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRVector2Symmetric.h"

NuTo::ConstraintLagrangeGlobalCrackOpening2D::ConstraintLagrangeGlobalCrackOpening2D(const NuTo::StructureIp* rStructure) :
        ConstraintLagrange(NuTo::Constraint::EQUAL) , ConstraintBase()
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
    //Lagrange mult + alpha + 2 disp Dofs
    int dof(4);
    rResult.Resize(dof,dof);
    rGlobalDofs.resize(dof,1);

    boost::array<int,2> dofCrackOpening(mStructure->GetDofGlobalCrackOpening2D());
    boost::array<double,2> crackOpening(mStructure->GetGlobalCrackOpening2D());

    rGlobalDofs[0] = mLagrangeDOF;
    rGlobalDofs[1] = mStructure->GetDofCrackAngle();
    rGlobalDofs[2] = dofCrackOpening[0];
    rGlobalDofs[3] = dofCrackOpening[1];

    //calculate normal crack opening
    double crackAngle(mStructure->GetCrackAngle());
    double cosAlpha(cos(crackAngle));
    double sinAlpha(sin(crackAngle));

    //std::cout << "normal crack opening" << -sinAlpha*crackOpening[0]+cosAlpha*crackOpening[1] << std::endl;

    double g(sinAlpha*crackOpening[0]-cosAlpha*crackOpening[1]);
    double dgda(cosAlpha*crackOpening[0]+sinAlpha*crackOpening[1]);

    if (g>-mLagrangeValue/mPenalty)
    {
        //derivative with respect to alpha and lambda
        rResult.AddEntry(0,1,dgda);
        //derivative with respect to ux and lambda
        rResult.AddEntry(0,2,sinAlpha);
        //derivative with respect to uy and lambda
        rResult.AddEntry(0,3,-cosAlpha);
        //derivative with respect to alpha and alpha
        rResult.AddEntry(1,1,-mLagrangeValue*g+mPenalty*(dgda*dgda-g*g));
        //derivative with respect to ux and alpha
        rResult.AddEntry(1,2,mLagrangeValue*cosAlpha + mPenalty*(sinAlpha*dgda+g*cosAlpha));
        //derivative with respect to uy and alpha
        rResult.AddEntry(1,3,mLagrangeValue*sinAlpha + mPenalty*(-cosAlpha*dgda+g*sinAlpha));
        //derivative with respect to ux^2
        rResult.AddEntry(2,2,mPenalty*sinAlpha*sinAlpha);
        //derivative with respect to ux and uy
        rResult.AddEntry(2,3,-mPenalty*sinAlpha*cosAlpha);
        //derivative with respect to uy^2
        rResult.AddEntry(3,3,mPenalty*cosAlpha*cosAlpha);
    }
    else
    {
        //derivative with respect to lambda^2
        rResult.AddEntry(0,0,-1./mPenalty);
    }

    //check result matrix by calling the gradient routine
    NuTo::FullMatrix<double> rResultFull(rResult);
    NuTo::FullMatrix<double> rResultCDF(rResult);
    NuTo::FullMatrix<double> gradient1,gradient2;
    double delta(1e-8);
    CalculateGradientInternalPotential(gradient1,rGlobalDofs);

    const_cast<ConstraintLagrangeGlobalCrackOpening2D*>(this)->mLagrangeValue+=delta;
    CalculateGradientInternalPotential(gradient2,rGlobalDofs);
    rResultCDF.SetColumn(0,(gradient2-gradient1)*(1./delta));
    const_cast<ConstraintLagrangeGlobalCrackOpening2D*>(this)->mLagrangeValue-=delta;

    crackAngle+=delta;
    mStructure->SetCrackAngle(crackAngle);
    CalculateGradientInternalPotential(gradient2,rGlobalDofs);
    rResultCDF.SetColumn(1,(gradient2-gradient1)*(1./delta));
    crackAngle-=delta;
    mStructure->SetCrackAngle(crackAngle);

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

    if ((rResultFull-rResultCDF).Abs().Max()<1e-2)
    {
        std::cout << "constraint crack opening "  << std::endl;
        std::cout<< "stiffness algo " << std::endl;
        rResultFull.Info(10,3);
        std::cout<< "stiffness cdf " << std::endl;
        rResultCDF.Info(10,3);
    }
}

//! @brief calculates the gradient of the internal potential
//! for a mechanical problem, this corresponds to the internal force vector
void NuTo::ConstraintLagrangeGlobalCrackOpening2D::CalculateGradientInternalPotential(NuTo::FullMatrix<double>& rResult,
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

    //calculate normal crack opening
    double crackAngle(mStructure->GetCrackAngle());
    double cosAlpha(cos(crackAngle));
    double sinAlpha(sin(crackAngle));

    //normal crack opening is -g, but since nco>=0 is the constraint, <=> -nco<=0
    double g(sinAlpha*crackOpening[0]-cosAlpha*crackOpening[1]);
    double dgda(cosAlpha*crackOpening[0]+sinAlpha*crackOpening[1]);

    if (g>-mLagrangeValue/mPenalty)
    {
        //std::cout << "constraint crack opening active(gradient) " << g << "  " << mLagrangeValue << std::endl;
        //derivative with respect to lambda
        rResult(0,0)= g;
        //derivative with respect to alpha
        rResult(1,0)= mLagrangeValue*dgda+mPenalty*g*dgda;
        //derivative with respect to ux
        rResult(2,0)= mLagrangeValue*sinAlpha+mPenalty*sinAlpha*g;
        //derivative with respect to uy
        rResult(3,0)= -mLagrangeValue*cosAlpha-mPenalty*cosAlpha*g;
    }
    else
    {
        //std::cout << "constraint crack opening inactive(gradient) " << g << "  " << mLagrangeValue << std::endl;
        //derivative with respect to lambda
        rResult(0,0)=-mLagrangeValue/mPenalty;
    }

    //check result matrix by calling the potential routine
    NuTo::FullMatrix<double> rResultFull(rResult);
    NuTo::FullMatrix<double> rResultCDF(rResult);
    double energy1,energy2;
    double delta(1e-8);
    energy1 = CalculateTotalPotential();

    const_cast<ConstraintLagrangeGlobalCrackOpening2D*>(this)->mLagrangeValue+=delta;
    energy2 = CalculateTotalPotential();
    rResultCDF(0,0) = (energy2-energy1)*(1./delta);
    const_cast<ConstraintLagrangeGlobalCrackOpening2D*>(this)->mLagrangeValue-=delta;

    crackAngle+=delta;
    mStructure->SetCrackAngle(crackAngle);
    energy2 = CalculateTotalPotential();
    rResultCDF(1,0) = (energy2-energy1)*(1./delta);
    crackAngle-=delta;
    mStructure->SetCrackAngle(crackAngle);

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

    if ((rResultFull-rResultCDF).Abs().Max()>1e-2)
    {
        std::cout<< "constraint crack opening gradient algo " << std::endl;
        rResultFull.Info(10,5);
        std::cout<< "constraint crack opening gradient cdf " << std::endl;
        rResultCDF.Info(10,5);
    }
}

//! @brief calculates the internal potential
double NuTo::ConstraintLagrangeGlobalCrackOpening2D::CalculateTotalPotential()const
{
    boost::array<int,2> dofCrackOpening(mStructure->GetDofGlobalCrackOpening2D());
    boost::array<double,2> crackOpening(mStructure->GetGlobalCrackOpening2D());

    //calculate normal crack opening
    double crackAngle(mStructure->GetCrackAngle());
    double cosAlpha(cos(crackAngle));
    double sinAlpha(sin(crackAngle));

    double g(sinAlpha*crackOpening[0]-cosAlpha*crackOpening[1]);

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

    //calculate normal crack opening
    double crackAngle(mStructure->GetCrackAngle());
    double cosAlpha(cos(crackAngle));
    double sinAlpha(sin(crackAngle));

    double g(sinAlpha*crackOpening[0]-cosAlpha*crackOpening[1]);
    std::cout << "CrackOpening : lambda " << mLagrangeValue << " g " << g <<
            " crackOpening " << crackOpening[0] << " " << crackOpening[1] <<
            "crackAngle  " << crackAngle <<
            std::endl;
    if (g>-mLagrangeValue/mPenalty)
    {
        std::cout << "constraint crack opening active " << std::endl;
    }
    else
    {
        std::cout << "constraint crack opening inactive " << std::endl;
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
