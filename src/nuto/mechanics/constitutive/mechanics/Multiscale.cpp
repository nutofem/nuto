// $Id: Multiscale.cpp 342 2010-10-18 12:39:08Z arnold2 $
// Multiscale.cpp
// created Apr 26, 2010 by Joerg F. Unger
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include <eigen2/Eigen/LU>
#include <eigen2/Eigen/Array>

#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentNonlocal3x3.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataMultiscale2DPlaneStrain.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"
#include "nuto/mechanics/constitutive/mechanics/Multiscale.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/SecondPiolaKirchhoffStress1D.h"
#include "nuto/mechanics/constitutive/mechanics/SecondPiolaKirchhoffStress2D.h"
#include "nuto/mechanics/constitutive/mechanics/SecondPiolaKirchhoffStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/GreenLagrangeStrain3D.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/IpDataBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionEnum.h"
#include "nuto/mechanics/structures/unstructured/StructureIp.h"

#define sqrt3 1.732050808
#define MAX_OMEGA 0.999
//#define ENABLE_DEBUG

NuTo::Multiscale::Multiscale() : ConstitutiveEngineeringStressStrain()
{
    SetParametersValid();
    mTolerance=1e-6;
}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void NuTo::Multiscale::serialize(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
       std::cout << "start serialize Multiscale" << std::endl;
#endif
       ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveEngineeringStressStrain)
           & BOOST_SERIALIZATION_NVP(mTolerance);
#ifdef DEBUG_SERIALIZATION
       std::cout << "finish serialize Multiscale" << std::endl;
#endif
    }
    BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Multiscale)
#endif // ENABLE_SERIALIZATION

//************ constitutive routines    ***********
//**  defined in structures/StructureIpConstitutive.cpp *********
//*************************************************
//  Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering plastic strain from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
void NuTo::Multiscale::GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
                                  const DeformationGradient1D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const
{
    throw MechanicsException("[NuTo::Multiscale::GetEngineeringPlasticStrain] not implemented.");
}
//  Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering plastic strain from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
void NuTo::Multiscale::GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
                                  const DeformationGradient2D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const
{
    throw MechanicsException("[NuTo::Multiscale::GetEngineeringPlasticStrain] not implemented.");
}

//  Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering plastic strain from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStrain ... engineering strain
void NuTo::Multiscale::GetEngineeringPlasticStrain(const ElementBase* rElement, int rIp,
                                  const DeformationGradient3D& rDeformationGradient, EngineeringStrain3D& rEngineeringPlasticStrain) const
{
    throw MechanicsException("[NuTo::Multiscale::GetEngineeringPlasticStrain] not implemented.");
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStress ... Engineering stress
void NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient, EngineeringStress1D& rEngineeringStress) const
{
    throw MechanicsException("[NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain] not implemented for 1D.");
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStress ... Engineering stress
void NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const
{
    throw MechanicsException("[NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain] not implemented for 1D.");
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStress ... Engineering stress
void NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient, EngineeringStress2D& rEngineeringStress) const
{
    // this is a somehow weird situation, since for a normal material law nothing should be changed
    // since the material law is a full structure whose bc change, this can either be implemented with a cast (as I did)
    // or by removing the const flag from all material routines (which I do not consider is good)
    std::stringstream previousState;
    bool updatePerformed;
    const ConstitutiveStaticDataMultiscale2DPlaneStrain *staticData = (rElement->GetStaticData(rIp))->AsMultiscale2DPlaneStrain();

    Solve(rElement, rIp, rDeformationGradient,mTolerance,previousState,updatePerformed);
    const StructureIp *fineScaleStructure = staticData->GetFineScaleStructure();


    //calculate average stress
    NuTo::FullMatrix<double> averageStress;

    double mlX, mlY;
    mlX = fineScaleStructure->GetDimensionX();
    mlY = fineScaleStructure->GetDimensionY();
    fineScaleStructure->ElementTotalGetAverageStress(mlX*mlY,averageStress);
    rEngineeringStress.mEngineeringStress[0] = averageStress(0,0);
    rEngineeringStress.mEngineeringStress[1] = averageStress(1,0);
    rEngineeringStress.mEngineeringStress[2] = averageStress(3,0);

    if (updatePerformed)
    {
        //restore previous state
        //don't be astonished, a const pointer does only mean, you are not allowed to change the pointer
        // but the object it points to can be changed - weird, but that's the way it is
//TODO test, if I need the const_cast
        const_cast<StructureIp*>(fineScaleStructure)->Restore(previousState);
   }
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rEngineeringStress ... Engineering stress
void NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const
{
    throw MechanicsException("[NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain] not implemented.");
}

// Engineering stress - Engineering strain /////////////////////////////////////
//! @brief ... calculate engineering stress from engineering strain (which is calculated from the deformation gradient)
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rCauchyStress ... Cauchy stress
void NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient, EngineeringStress3D& rEngineeringStress) const
{
    throw MechanicsException("[NuTo::Multiscale::GetEngineeringStressFromEngineeringStrain] not implemented.");
}

//  Damage /////////////////////////////////////
//! @brief ... calculate isotropic damage from deformation gradient in 1D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDamage ... damage variable
void NuTo::Multiscale::GetDamage(const ElementBase* rElement, int rIp,
                                  const DeformationGradient1D& rDeformationGradient, double& rDamage) const
{
    throw MechanicsException("[NuTo::Multiscale::GetDamage] not implemented.");
}

//  Damage /////////////////////////////////////
//! @brief ... calculate isotropic damage from deformation gradient in 2D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDamage ... damage variable
void NuTo::Multiscale::GetDamage(const ElementBase* rElement, int rIp,
                                  const DeformationGradient2D& rDeformationGradient, double& rDamage) const
{
    throw MechanicsException("[NuTo::Multiscale::GetDamage] not implemented.");
}

//  Damage /////////////////////////////////////
//! @brief ... calculate isotropic damage from deformation gradient in 3D
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDamage ... damage variable
void NuTo::Multiscale::GetDamage(const ElementBase* rElement, int rIp,
                                  const DeformationGradient3D& rDeformationGradient, double& rDamage) const
{
    throw MechanicsException("[NuTo::Multiscale::GetDamage] not implemented.");
}

//! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rTangent ... tangent
void NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient,
        ConstitutiveTangentBase* rTangent) const
{
    throw MechanicsException("[NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rTangent ... tangent
void NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient,
        ConstitutiveTangentBase* rTangent) const
{
#ifdef HAVE_MUMPS
    //solve
    std::stringstream previousState;
    bool updatePerformed;
    const ConstitutiveStaticDataMultiscale2DPlaneStrain *staticData = (rElement->GetStaticData(rIp))->AsMultiscale2DPlaneStrain();

    const StructureIp *fineScaleStructure = staticData->GetFineScaleStructure();
    std::cout << "deformation gradient artificially modified" << std::endl;
    const_cast<DeformationGradient2D &> (rDeformationGradient).mDeformationGradient[0]=1.0;
    const_cast<DeformationGradient2D &> (rDeformationGradient).mDeformationGradient[1]=0.01;
    const_cast<DeformationGradient2D &> (rDeformationGradient).mDeformationGradient[2]=0.01;
    const_cast<DeformationGradient2D &> (rDeformationGradient).mDeformationGradient[3]=1.0;

    Solve(rElement, rIp, rDeformationGradient,mTolerance,previousState,updatePerformed);

    // use Schur complement to calculate the stiffness
    SparseMatrixCSRVector2General<double>
        matrixJJ(fineScaleStructure->GetNumActiveDofs(), fineScaleStructure->GetNumActiveDofs()),
        matrixJK(fineScaleStructure->GetNumActiveDofs(), fineScaleStructure->GetNumDofs() - fineScaleStructure->GetNumActiveDofs()),
        matrixJM(fineScaleStructure->GetNumActiveDofs(), 3),
        matrixKJ(fineScaleStructure->GetNumDofs() - fineScaleStructure->GetNumActiveDofs(), fineScaleStructure->GetNumActiveDofs()),
        matrixKK(fineScaleStructure->GetNumDofs() - fineScaleStructure->GetNumActiveDofs(), fineScaleStructure->GetNumDofs() - fineScaleStructure->GetNumActiveDofs()),
        matrixKM(fineScaleStructure->GetNumDofs() - fineScaleStructure->GetNumActiveDofs(), 3),
        matrixMJ(3, fineScaleStructure->GetNumActiveDofs()),
        matrixMK(3, fineScaleStructure->GetNumDofs() - fineScaleStructure->GetNumActiveDofs()),
        matrixMM(3, 3);

    fineScaleStructure->BuildGlobalCoefficientSubMatrices0General(
            matrixJJ, matrixJK, matrixJM,
            matrixKJ, matrixKK, matrixKM,
            matrixMJ, matrixMK, matrixMM);

    const_cast<StructureIp*> (fineScaleStructure)->ElementTotalUpdateTmpStaticData();
    double totEnergy = fineScaleStructure->ElementTotalGetTotalEnergy();
    const_cast<StructureIp*> (fineScaleStructure)->CalculateHomogeneousEngineeringStrain();
    //calculate average strain
    NuTo::FullMatrix<double> averageStrain;
    double mlX, mlY;
    mlX = fineScaleStructure->GetDimensionX();
    mlY = fineScaleStructure->GetDimensionY();
    fineScaleStructure->ElementTotalGetAverageStrain(mlX*mlY,averageStrain);
    std::cout<<"average strain " << std::endl;
    averageStrain.Info(12,4);
    fineScaleStructure->NodeInfo(10);
    fineScaleStructure->ExportVtkDataFile("/home/unger3/develop/nuto_build/examples/c++/FineScaleConcurrentMultiscale.vtk");

/*    {
    std::cout << "matrixJJ algo "<< std::endl;
    NuTo::FullMatrix<double> matrixJJFull(matrixJJ);
    matrixJJFull.Info(12,3);

    std::cout << "matrixKJ algo "<< std::endl;
    NuTo::FullMatrix<double> matrixKJFull(matrixKJ);
    matrixKJFull.Info(12,3);

    std::cout << "matrixMJ algo "<< std::endl;
    NuTo::FullMatrix<double> matrixMJFull(matrixMJ);
    matrixMJFull.Info(12,3);

    std::cout << "matrixJM algo "<< std::endl;
    NuTo::FullMatrix<double> matrixJMFull(matrixJM);
    matrixJMFull.Info(12,3);

    std::cout << "matrixKM algo "<< std::endl;
    NuTo::FullMatrix<double> matrixKMFull(matrixKM);
    matrixKMFull.Info(12,3);

    std::cout << "matrixMM algo "<< std::endl;
    NuTo::FullMatrix<double> matrixMMFull(matrixMM);
    matrixMMFull.Info(12,3);
}
*/

    // build global matrix
    SparseMatrixCSRVector2General<double> constraintMatrixVector2 (fineScaleStructure->GetConstraintMatrix());
    SparseMatrixCSRVector2General<double> transConstraintMatrixVector2 (constraintMatrixVector2.Transpose());

    matrixJJ -= transConstraintMatrixVector2 * matrixKJ + matrixJK * constraintMatrixVector2;
    matrixJM -= transConstraintMatrixVector2 * matrixKM;
    matrixMJ -= matrixMK * constraintMatrixVector2;

    matrixJJ += transConstraintMatrixVector2 * matrixKK * constraintMatrixVector2;

    //join the matrix
    matrixJJ.ConcatenateColumns(matrixJM);
    matrixMJ.ConcatenateColumns(matrixMM);
    matrixJJ.ConcatenateRows(matrixMJ);

    /*
    {
        //check the Global Matrix
        NuTo::FullMatrix<double> matrixJJCDF(fineScaleStructure->GetNumActiveDofs(),fineScaleStructure->GetNumActiveDofs());
        NuTo::FullMatrix<double> matrixKJCDF(fineScaleStructure->GetNumDofs() - fineScaleStructure->GetNumActiveDofs(),fineScaleStructure->GetNumActiveDofs());
        NuTo::FullMatrix<double> matrixMJCDF(3,fineScaleStructure->GetNumActiveDofs());
        NuTo::FullMatrix<double> matrixJMCDF(fineScaleStructure->GetNumActiveDofs(),3);
        NuTo::FullMatrix<double> matrixKMCDF(fineScaleStructure->GetNumDofs() - fineScaleStructure->GetNumActiveDofs(),3);
        NuTo::FullMatrix<double> matrixMMCDF(3,3);
        //std::cout<<"stiffness matrix" << std::endl;
        //stiffnessMatrixCSRVector2Full.Info(10,3);
        double interval(1e-7);
        NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
        NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
        StructureIp *fineScaleStructureNonConst = const_cast<StructureIp*> (fineScaleStructure);

        fineScaleStructureNonConst->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
        fineScaleStructureNonConst->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
        fineScaleStructureNonConst->ElementTotalUpdateTmpStaticData();

        NuTo::FullMatrix<double> activeGrad1(displacementsActiveDOFsCheck),
                                 dependentGrad1(displacementsDependentDOFsCheck),
                                 epsilonMGrad1(3,1),
                                 activeGrad2(displacementsActiveDOFsCheck),
                                 dependentGrad2(displacementsDependentDOFsCheck),
                                 epsilonMGrad2(3,1);
        fineScaleStructureNonConst->BuildGlobalGradientInternalPotentialSubVectors(activeGrad1,dependentGrad1,&epsilonMGrad1);
        for (int count=0; count<displacementsActiveDOFsCheck.GetNumRows(); count++)
        {
            displacementsActiveDOFsCheck(count,0)+=interval;
            fineScaleStructureNonConst->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
            fineScaleStructureNonConst->ElementTotalUpdateTmpStaticData();
            fineScaleStructureNonConst->BuildGlobalGradientInternalPotentialSubVectors(activeGrad2,dependentGrad2,&epsilonMGrad2);
            matrixJJCDF.SetColumn(count,(activeGrad2-activeGrad1)*(1./interval));
            matrixKJCDF.SetColumn(count,(dependentGrad2-dependentGrad1)*(1./interval));
            matrixMJCDF.SetColumn(count,(epsilonMGrad2-epsilonMGrad1)*(1./interval));

            displacementsActiveDOFsCheck(count,0)-=interval;
        }
        fineScaleStructureNonConst->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);

        std::cout << "matrixJJ cdf "<< std::endl;
        matrixJJCDF.Info(12,3);

        std::cout << "matrixKJ cdf "<< std::endl;
        matrixKJCDF.Info(12,3);

        std::cout << "matrixMJ cdf "<< std::endl;
        matrixMJCDF.Info(12,3);

        fineScaleStructureNonConst->BuildGlobalGradientInternalPotentialSubVectors(activeGrad1,dependentGrad1,&epsilonMGrad1);
        NuTo::EngineeringStrain2D strain(fineScaleStructure->GetTotalEngineeringStrain());

        interval = 0.1;
        std::cout << "active grad 1" << std::endl;
        activeGrad1.Trans().Info(12,3);
        std::cout << "dependent grad 1" << std::endl;
        dependentGrad1.Trans().Info(12,3);
        std::cout << "epsilonMGrad 1" << std::endl;
        epsilonMGrad1.Trans().Info(12,3);
        for (int count=0; count<3; count++)
        {
            strain.mEngineeringStrain[count]+=interval;
            fineScaleStructureNonConst->SetTotalEngineeringStrain(strain);
            fineScaleStructureNonConst->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
            fineScaleStructureNonConst->ElementTotalUpdateTmpStaticData();
            fineScaleStructureNonConst->BuildGlobalGradientInternalPotentialSubVectors(activeGrad2,dependentGrad2,&epsilonMGrad2);
            matrixJMCDF.SetColumn(count,(activeGrad2-activeGrad1)*(1./interval));
            matrixKMCDF.SetColumn(count,(dependentGrad2-dependentGrad1)*(1./interval));
            matrixMMCDF.SetColumn(count,(epsilonMGrad2-epsilonMGrad1)*(1./interval));

            strain.mEngineeringStrain[count]-=interval;
            fineScaleStructureNonConst->SetTotalEngineeringStrain(strain);

            std::cout << "active grad 2" << std::endl;
            activeGrad2.Trans().Info(12,3);
            std::cout << "dependent grad 2" << std::endl;
            dependentGrad2.Trans().Info(12,3);
            std::cout << "epsilonMGrad 2" << std::endl;
            epsilonMGrad2.Trans().Info(12,3);
        }
        fineScaleStructureNonConst->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
        std::cout << "matrixJM cdf "<< std::endl;
        matrixJMCDF.Info(12,3);

        std::cout << "matrixKM cdf "<< std::endl;
        matrixKMCDF.Info(12,3);

        std::cout << "matrixMM cdf "<< std::endl;
        matrixMMCDF.Info(12,3);

        exit(0);

    }
*/
    {
    std::cout << "stiffness before solution" << std::endl;
    NuTo::FullMatrix<double> matrixJJFull(matrixJJ);
    matrixJJFull.Info(12,3);
    NuTo::FullMatrix<double> eigenvalues, eigenvectors;
    matrixJJFull.EigenValuesSymmetric(eigenvalues);
    matrixJJFull.EigenVectorsSymmetric(eigenvectors);
    std::cout << "eigenvalues " << std::endl;
    eigenvalues.Trans().Info(12,3);
    std::cout << "eigenvectors " << std::endl;
    eigenvectors.Info(12,3);
    }

    //calculate schur complement
    NuTo::FullMatrix<int> schurIndicesMatrix(3,1);
    NuTo::FullMatrix<double> stiffness(3,3);
    //attention, the index is in the zero based indexing system
    schurIndicesMatrix(0,0) = matrixJJ.GetNumRows()-3;
    schurIndicesMatrix(1,0) = matrixJJ.GetNumRows()-2;
    schurIndicesMatrix(2,0) = matrixJJ.GetNumRows()-1;

    NuTo::SparseDirectSolverMUMPS mumps;
    NuTo::SparseMatrixCSRGeneral<double> stiffnessFineScale(matrixJJ);
    stiffnessFineScale.SetOneBasedIndexing();
    mumps.SchurComplement(stiffnessFineScale,schurIndicesMatrix,stiffness);

    //scale with the dimension of the structure (area)
    stiffness*=1./(fineScaleStructure->GetDimensionX()*fineScaleStructure->GetDimensionY());
    *(rTangent->AsConstitutiveTangentLocal3x3()) = stiffness;

    if (updatePerformed)
    {
        //restore previous state
        //don't be astonished, a const pointer does only mean, you are not allowed to change the pointer
        // but the object it points to can be changed - weird, but that's the way it is
        const_cast<StructureIp*>(fineScaleStructure)->Restore(previousState);
    }

    //just for test purpose, calculate stiffness via resforces
    EngineeringStress2D stress1, stress2;
    GetEngineeringStressFromEngineeringStrain(rElement, rIp, rDeformationGradient, stress1);
    double delta(1e-8);
    NuTo::FullMatrix<double> stiffnessCDF(3,3);
    for (int count=0; count<3; count++)
    {
       DeformationGradient2D deformationGradient(rDeformationGradient);
       switch(count)
        {
        case 0:
            deformationGradient.mDeformationGradient[0]+=delta;
            break;
        case 1:
            deformationGradient.mDeformationGradient[3]+=delta;
            break;
        case 2:
            deformationGradient.mDeformationGradient[1]+=0.5*delta;
            deformationGradient.mDeformationGradient[2]+=0.5*delta;
            break;
        default:
            throw MechanicsException("");
        }
        GetEngineeringStressFromEngineeringStrain(rElement, rIp, deformationGradient, stress2);

        for (int count2=0; count2<3; count2++)
        {
            stiffnessCDF(count2,count) = (stress2.GetData()[count2]- stress1.GetData()[count2])/delta;
        }
    }
    std::cout << "Schur stiffness algo" << std::endl;
    stiffness.Info(12,3);
    std::cout << "Schur stiffness cdf" << std::endl;
    stiffnessCDF.Info(12,3);

    exit(0);
#else
    throw MechanicsException("[NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain] MUMPS solver required to calculate Schur complement.");
#endif //HAVE_MUMPS
}

//! @brief ... calculate the tangent (derivative of the Engineering stresses with respect to the engineering strains) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rTangent ... tangent
void NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient,
        ConstitutiveTangentBase* rTangent) const
{
    throw MechanicsException("[NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient) const
{

    double tolerance(1e-6);
    // this is a somehow weird situation, since for a normal material law nothing should be changed
    // since the material law is a full structure whose bc change, this can either be implemented with a cast (as I did)
    // or by removing the const flag from all material routines (which I do not consider is good)
    std::stringstream previousState;
    bool updatePerformed;
    const ConstitutiveStaticDataMultiscale2DPlaneStrain *staticData = (rElement->GetStaticData(rIp))->AsMultiscale2DPlaneStrain();

    Solve(rElement, rIp, rDeformationGradient,tolerance,previousState,updatePerformed);

    StructureIp* fineScaleStructure =  const_cast<StructureIp*>(staticData->GetFineScaleStructure());
    fineScaleStructure->ElementTotalUpdateStaticData();
    fineScaleStructure->SetPrevCrackAngle(fineScaleStructure->GetCrackAngle());
}

//! @brief ... update static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::Multiscale::UpdateStaticData_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... update tmp static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::Multiscale::UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::Multiscale::UpdateTmpStaticData_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... update tmp static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::Multiscale::UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient) const
{
    return;
}

//! @brief ... update tmp static data (history variables) of the constitutive relationship
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
void NuTo::Multiscale::UpdateTmpStaticData_EngineeringStress_EngineeringStrain(ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::Multiscale::UpdateTmpStaticData_EngineeringStress_EngineeringStrain] not implemented.");
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::Multiscale::AllocateStaticDataEngineeringStress_EngineeringStrain1D(const ElementBase* rElement) const
{
    throw MechanicsException("[NuTo::Multiscale::AllocateStaticDataEngineeringStress_EngineeringStrain1D] not implemented.");
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::Multiscale::AllocateStaticDataEngineeringStress_EngineeringStrain2D(const ElementBase* rElement) const
{
    if (rElement->GetSection()==0)
        throw MechanicsException("[NuTo::Multiscale::AllocateStaticDataEngineeringStress_EngineeringStrain1D] Section required to distinguish between plane stress and plane strain and thickness information.");
    if (rElement->GetSection()->GetType()==NuTo::Section::PLANE_STRESS)
        throw MechanicsException("[NuTo::Multiscale::AllocateStaticDataEngineeringStress_EngineeringStrain1D] Nonlocal damage plasticity model not implemented for plane stress.");
    else
        return new ConstitutiveStaticDataMultiscale2DPlaneStrain();
}

//! @brief ... create new static data object for an integration point
//! @return ... pointer to static data object
NuTo::ConstitutiveStaticDataBase* NuTo::Multiscale::AllocateStaticDataEngineeringStress_EngineeringStrain3D(const ElementBase* rElement) const
{
    throw MechanicsException("[NuTo::Multiscale::AllocateStaticDataEngineeringStress_EngineeringStrain3D] not implemented.");
}

//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain] not yet implemented.");
    /*
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
           //throw an exception giving information related to the wrong parameter
        CheckParameters();
        //if there is no exception thrown there is a problem with the source code
        //since every time a material parameter is changed, the parametes should be checked
        throw MechanicsException("[NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }

    // calculate engineering strain
    EngineeringStrain1D engineeringStrain;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);
    return 0.5 * engineeringStrain.mEngineeringStrain * this->mE * engineeringStrain.mEngineeringStrain;
*/
}


//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain] not yet implemented.");
/*
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
           //throw an exception giving information related to the wrong parameter
        CheckParameters();
        //if there is no exception thrown there is a problem with the source code
        //since every time a material parameter is changed, the parametes should be checked
        throw MechanicsException("[NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }
    // calculate engineering strain
    EngineeringStrain2D engineeringStrain;
    EngineeringStress2D engineeringStress;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);

    const SectionBase* theSection(rElement->GetSection());
    if (theSection==0)
        throw MechanicsException("[NuTo::Multiscale::GetTangent_EngineeringStress_EngineeringStrain] No section defined for element.");
    if (theSection->GetType()==Section::PLANE_STRAIN)
    {
        // calculate coefficients of the material matrix
        double C11, C12, C33;
        this->CalculateCoefficients3D(C11, C12, C33);

        // calculate Engineering stress
        engineeringStress.mEngineeringStress[0] = C11 * engineeringStrain.mEngineeringStrain[0] + C12 * engineeringStrain.mEngineeringStrain[1];
        engineeringStress.mEngineeringStress[1] = C11 * engineeringStrain.mEngineeringStrain[1] + C12 * engineeringStrain.mEngineeringStrain[0];
        engineeringStress.mEngineeringStress[2] = C33 * engineeringStrain.mEngineeringStrain[2] ;
    }
    else
    {
        throw MechanicsException("[NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain] Plane stress is to be implemented.");
    }
    return 0.5*(
            engineeringStrain.mEngineeringStrain[0]*engineeringStress.mEngineeringStress[0]
           +engineeringStrain.mEngineeringStrain[1]*engineeringStress.mEngineeringStress[1]
           +engineeringStrain.mEngineeringStrain[2]*engineeringStress.mEngineeringStress[2]);
*/
}


//! @brief ... calculate the total energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain] not yet implemented.");
/*
    // check if parameters are valid
    if (this->mParametersValid == false)
    {
           //throw an exception giving information related to the wrong parameter
        CheckParameters();
        //if there is no exception thrown there is a problem with the source code
        //since every time a material parameter is changed, the parametes should be checked
        throw MechanicsException("[NuTo::Multiscale::GetTotalEnergy_EngineeringStress_EngineeringStrain] Check the material parameters.");
    }
    // calculate engineering strain
    EngineeringStrain3D engineeringStrain;
    EngineeringStress3D engineeringStress;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);

    // calculate coefficients of the material matrix
    double C11, C12, C44;
    this->CalculateCoefficients3D(C11, C12, C44);

    // calculate Engineering stress
    engineeringStress.mEngineeringStress[0] = C11 * engineeringStrain.mEngineeringStrain[0] + C12 * (engineeringStrain.mEngineeringStrain[1]+engineeringStrain.mEngineeringStrain[2]);
    engineeringStress.mEngineeringStress[1] = C11 * engineeringStrain.mEngineeringStrain[1] + C12 * (engineeringStrain.mEngineeringStrain[0]+engineeringStrain.mEngineeringStrain[2]);
    engineeringStress.mEngineeringStress[2] = C11 * engineeringStrain.mEngineeringStrain[2] + C12 * (engineeringStrain.mEngineeringStrain[0]+engineeringStrain.mEngineeringStrain[1]);
    engineeringStress.mEngineeringStress[3] = C44 * engineeringStrain.mEngineeringStrain[3] ;
    engineeringStress.mEngineeringStress[4] = C44 * engineeringStrain.mEngineeringStrain[4] ;
    engineeringStress.mEngineeringStress[5] = C44 * engineeringStrain.mEngineeringStrain[5] ;

    return 0.5*(
            engineeringStrain.mEngineeringStrain[0]*engineeringStress.mEngineeringStress[0]
           +engineeringStrain.mEngineeringStrain[1]*engineeringStress.mEngineeringStress[1]
           +engineeringStrain.mEngineeringStrain[2]*engineeringStress.mEngineeringStress[2]
           +engineeringStrain.mEngineeringStrain[3]*engineeringStress.mEngineeringStress[3]
           +engineeringStrain.mEngineeringStrain[4]*engineeringStress.mEngineeringStress[4]
           +engineeringStrain.mEngineeringStrain[5]*engineeringStress.mEngineeringStress[5]);
*/
}


//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::Multiscale::GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::Multiscale::GetElasticEnergy_EngineeringStress_EngineeringStrain] not yet implemented.");
}


//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::Multiscale::GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient) const
{
    throw MechanicsException("[NuTo::Multiscale::GetElasticEnergy_EngineeringStress_EngineeringStrain] not yet implemented.");
}


//! @brief ... calculate the elastic energy density
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
double NuTo::Multiscale::GetElasticEnergy_EngineeringStress_EngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient) const
{
    return GetTotalEnergy_EngineeringStress_EngineeringStrain(rElement, rIp, rDeformationGradient);
}


//! @brief ... calculates the difference of the elastic strain between the current state and the previous update
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
void NuTo::Multiscale::GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient1D& rDeformationGradient, EngineeringStrain1D& rDeltaElasticEngineeringStrain) const
{
    throw MechanicsException("[NuTo::Multiscale::GetElasticEnergy_EngineeringStress_EngineeringStrain] not yet implemented.");
}

//! @brief ... calculates the difference of the elastic strain between the current state and the previous update
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
void NuTo::Multiscale::GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient2D& rDeformationGradient, EngineeringStrain2D& rDeltaElasticEngineeringStrain) const
{
    throw MechanicsException("[NuTo::Multiscale::GetDeltaElasticEngineeringStrain] not implemented.");
}

//! @brief ... calculates the difference of the elastic strain between the current state and the previous update
//! @param rStructure ... structure
//! @param rElement ... element
//! @param rIp ... integration point
//! @param rDeformationGradient ... deformation gradient
//! @param rDeltaElasticEngineeringStrain ... delta elastic engineering strain (return value)
void NuTo::Multiscale::GetDeltaElasticEngineeringStrain(const ElementBase* rElement, int rIp,
        const DeformationGradient3D& rDeformationGradient, EngineeringStrain3D& rDeltaElasticEngineeringStrain) const
{
    throw MechanicsException("[NuTo::Multiscale::GetDeltaElasticEngineeringStrain] not implemented.");
}

#define MAXNUMNEWTONITERATIONS 100
#define PRINTRESULT true
#define MIN_DELTA_STRAIN_FACTOR 1e-7
void NuTo::Multiscale::Solve(const ElementBase* rElement, int rIp, const NuTo::DeformationGradient2D &rDeformationGradient, double rTolerance,
        std::stringstream& rStringStreamBeforeSolve, bool& rStringStreamBeforeSolveWritten)const
{
    const ConstitutiveStaticDataMultiscale2DPlaneStrain *staticData = (rElement->GetStaticData(rIp))->AsMultiscale2DPlaneStrain();
    StructureIp *fineScaleStructure = const_cast<StructureIp*>(staticData->GetFineScaleStructure());

    // calculate engineering strain
    EngineeringStrain2D engineeringStrain;
    rDeformationGradient.GetEngineeringStrain(engineeringStrain);

    // write data to stringstream
    rStringStreamBeforeSolveWritten = false;
    //if there is a problem with memory this could be replaced by writing it to a file
    //std::cout << "sizeof string " << stringStreamBeforeSolve.str().length() << std::endl;

    // start analysis
    double deltaStrainFactor(1.0);
    double curStrainFactor(1.0);

    const ConstitutiveStaticDataMultiscale2DPlaneStrain *oldStaticData = rElement->GetStaticData(rIp)->AsMultiscale2DPlaneStrain();
    EngineeringStrain2D prevStrain(oldStaticData->GetPrevStrain());

    EngineeringStrain2D totalEngineeringStrain(engineeringStrain);

    EngineeringStrain2D deltaEngineeringStrain;
    deltaEngineeringStrain = engineeringStrain-prevStrain;

    //update conre mat and set displacement according to linear constraints
    fineScaleStructure->NodeBuildGlobalDofs();
    {
        NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
        NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
        fineScaleStructure->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
        fineScaleStructure->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
        fineScaleStructure->ElementTotalUpdateTmpStaticData();
    }
    //fineScaleStructure->ExportVtkDataFile("/home/unger3/develop/nuto_build/examples/c++/FineScaleConcurrentMultiscaleCrack.vtk");
    //exit(0);

    //set the total strain and calculate from the existing crack opening the homogeneous strain
    fineScaleStructure->SetTotalEngineeringStrain(prevStrain);

    //update tmpstatic data with zero displacements
    fineScaleStructure->ElementTotalUpdateTmpStaticData();

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
    fineScaleStructure->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);

//Check the stiffness matrix
CheckStiffness(fineScaleStructure);

    //set the total strain and calculate from the existing crack opening the homogeneous strain
    EngineeringStrain2D curEngineeringStrain(prevStrain+deltaEngineeringStrain*curStrainFactor);
    fineScaleStructure->SetTotalEngineeringStrain(curEngineeringStrain);

    //update conre mat
    fineScaleStructure->NodeBuildGlobalDofs();

    //update displacements of all nodes according to the new conre mat
    {
        NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
        NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
        fineScaleStructure->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
        fineScaleStructure->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
        fineScaleStructure->ElementTotalUpdateTmpStaticData();
    }

    // build global external load vector and RHS vector
    fineScaleStructure->BuildGlobalExternalLoadVector(extForceVector);
    fineScaleStructure->BuildGlobalGradientInternalPotentialVector(intForceVector);
    // calculate residual
    std::cout << "initial intForceVector "  << std::endl;
    intForceVector.Trans().Info(10,3);

    rhsVector = extForceVector + dispForceVector - intForceVector;

    //calculate absolute tolerance for matrix entries to be not considered as zero
    double maxValue, minValue, ToleranceZeroStiffness;
    stiffnessMatrixCSRVector2.Max(maxValue);
    stiffnessMatrixCSRVector2.Min(minValue);
    //std::cout << "min and max " << minValue << " , " << maxValue << std::endl;

    ToleranceZeroStiffness = (1e-14) * (fabs(maxValue)>fabs(minValue) ?  fabs(maxValue) : fabs(minValue));
    fineScaleStructure->SetToleranceStiffnessEntries(ToleranceZeroStiffness);
    int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
    int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
    //std::cout << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << std::endl;


    //repeat until max displacement is reached
    bool convergenceStatusLoadSteps(false);
    int loadstep(1);
    NuTo::FullMatrix<double> displacementsActiveDOFsLastConverged,displacementsDependentDOFsLastConverged;
    fineScaleStructure->NodeExtractDofValues(displacementsActiveDOFsLastConverged,displacementsDependentDOFsLastConverged);
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
            NuTo::FullMatrix<double> stiffnessMatrixFull(stiffnessMatrixCSR);
            std::cout << " stiffness" << std::endl;
            stiffnessMatrixFull.Info(10,3);
            NuTo::FullMatrix<double> eigenValues;
            stiffnessMatrixFull.EigenValuesSymmetric(eigenValues);
            NuTo::FullMatrix<double> eigenVectors;
            stiffnessMatrixFull.EigenVectorsSymmetric(eigenVectors);
            std::cout << "DOFs : alpha "<< fineScaleStructure->GetDofCrackAngle() << " crack opening "
                      << fineScaleStructure->GetDofGlobalCrackOpening2D()[0] << " "
                      << fineScaleStructure->GetDofGlobalCrackOpening2D()[1] << std::endl;
            std::cout << " eigenvalues" << std::endl;
            eigenValues.Trans().Info(12,7);
            std::cout << " eigenvectors" << std::endl;
            eigenVectors.Info(10,3);
            //std::cout << " rhsVector" << std::endl;
            //rhsVector.Trans().Info(10,3);
            std::cout << " delta_disp" << std::endl;
            deltaDisplacementsActiveDOFs.Trans().Info(10,3);
            double deltaAlpha=deltaDisplacementsActiveDOFs(fineScaleStructure->GetDofCrackAngle(),0);
            std::cout << " delta_alpha " << deltaAlpha <<std::endl;

            // write displacements to node
            fineScaleStructure->NodeExtractDofValues(oldDisplacementsActiveDOFs, displacementsDependentDOFs);

            double alphaCrackAngle;
            if (fabs(deltaAlpha)>M_PI/180.*30.) //5 degree
            {
                alphaCrackAngle=M_PI/180.*30./fabs(deltaAlpha);
                //scaling of delta disp so that delta_alpha==5degree
            }
            else
                alphaCrackAngle = 1;
            //perform a linesearch
            alpha = alphaCrackAngle;
            do
            {
                //add new displacement state
                displacementsActiveDOFs = oldDisplacementsActiveDOFs + deltaDisplacementsActiveDOFs*alpha;

                std::cout << " displacementsActiveDOFs" << std::endl;
                displacementsActiveDOFs.Trans().Info(10,3);
                fineScaleStructure->NodeMergeActiveDofValues(displacementsActiveDOFs);
                fineScaleStructure->ElementTotalUpdateTmpStaticData();

                // calculate residual
                fineScaleStructure->BuildGlobalGradientInternalPotentialVector(intForceVector);
                std::cout << "intForceVector "  << std::endl;
                intForceVector.Trans().Info(10,3);

                rhsVector = extForceVector - intForceVector;
                normResidual = rhsVector.Norm();

                alpha*=0.5;
            }
            while(alpha>1e-2*alphaCrackAngle && normResidual>normRHS*(1-0.5*alpha) && normResidual>rTolerance);

            double energyElement(fineScaleStructure->ElementTotalGetTotalEnergy());
            double energyConstraint(fineScaleStructure->ConstraintTotalGetTotalEnergy());
            EngineeringStrain2D strainHom(fineScaleStructure->GetHomogeneousEngineeringStrain());
            std::cout << "cutback factor " << alpha*2 << ", normResidual " << normResidual << ", normResidualInit "<< normRHS << ", normRHS*(1-0.5*alpha) " << normRHS*(1-0.5*alpha)
                      << " energy " << energyConstraint+energyElement <<"(" << energyElement << ","<< energyConstraint <<")"
                      << " hom strain " << strainHom.GetData()[0] << " "
                      << strainHom.GetData()[1] << " " << strainHom.GetData()[2] << std::endl;
            double crackAngle(displacementsActiveDOFs(fineScaleStructure->GetDofCrackAngle(),0));
            while (crackAngle>2.*M_PI)
                crackAngle-=2.*M_PI;
            while (crackAngle<0)
                crackAngle+=2.*M_PI;
            std::cout << "crack angle "<< crackAngle*180/M_PI << " crack opening " << std::endl;
/*                      << displacementsActiveDOFs(fineScaleStructure->GetDofGlobalCrackOpening2D()[0],0) << " "
                      << displacementsActiveDOFs(fineScaleStructure->GetDofGlobalCrackOpening2D()[1],0) << std::endl;
*/
            std::stringstream ss;
            ss << numNewtonIterations;
            fineScaleStructure->ExportVtkDataFile(std::string("/home/unger3/develop/nuto_build/examples/c++/FineScaleConcurrentMultiscale") + ss.str() + std::string(".vtk"));

            if (normResidual>normRHS*(1-0.5*alpha) && normResidual>rTolerance)
            {
                convergenceStatus=2;
                break;
            }

            maxResidual = rhsVector.Abs().Max();

            //std::cout << std::endl << "Newton iteration " << numNewtonIterations << ", final alpha " << 2*alpha << ", normResidual " << normResidual<< ", maxResidual " << maxResidual<<std::endl;

            //check convergence
            if (normResidual<rTolerance || maxResidual<rTolerance)
            {
                if (PRINTRESULT)
                {
                    std::cout <<fineScaleStructure->GetIPName() << " Convergence after " << numNewtonIterations << " Newton iterations, curStrainFactor " << curStrainFactor << ", deltaStrainFactor "<< deltaStrainFactor << std::endl<< std::endl;
                    fineScaleStructure->ConstraintInfo(10);
                }
                convergenceStatus=1;
                break;
            }

            //convergence status == 0 (continue Newton iteration)
            //build new stiffness matrix
            fineScaleStructure->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);

//check stiffness
CheckStiffness(fineScaleStructure);
            int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
            int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
            std::cout << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << std::endl;
        }

        if (deltaStrainFactor<1e-7)
            throw NuTo::MechanicsException("[NuTo::Multiscale::Solve] No convergence, delta strain factor < 1e-7");

        if (convergenceStatus==1)
        {
            // visualize results

//#ifdef ENABLE_VISUALIZE
//            std::cout << " store element id and ip in output file" << std::endl;
//            this->ExportVtkDataFile(mIPName+std::string(".vtk"));
//#endif
//#ifdef ENABLE_SERIALIZATION
//            this->Save(mIPName+std::string("bin"),"BINARY");
//#endif // ENABLE_SERIALIZATION

            fineScaleStructure->NodeExtractDofValues(displacementsActiveDOFsLastConverged,displacementsDependentDOFsLastConverged);
            if (curStrainFactor==1)
            {
                convergenceStatusLoadSteps=true;
            }
            else
            {
                if (rStringStreamBeforeSolveWritten==false)
                {
                    fineScaleStructure->Save(rStringStreamBeforeSolve);
                    rStringStreamBeforeSolveWritten=true;
                }

                // the update is only required to allow for a stepwise solution procedure in the fine scale model
                // a final update is only required for an update on the macroscale, otherwise,the original state has
                // to be reconstructed.
                fineScaleStructure->ElementTotalUpdateStaticData();

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

                curEngineeringStrain = prevStrain + deltaEngineeringStrain * curStrainFactor;

                //old stiffness matrix is used in first step of next load increment in order to prevent spurious problems at the boundary
                //std::cout << "press enter to next load increment, delta strain factor " << deltaStrainFactor << " max delta strain factor " <<  maxDeltaStrainFactor << std::endl << std::endl;
                //char cDummy[100]="";
                //std::cin.getline(cDummy, 100);
            }
            loadstep++;
        }
        else
        {
            assert(convergenceStatus==2);
            //calculate stiffness of previous loadstep (used as initial stiffness in the next load step)
            //this is done within the loop in order to ensure, that for the first step the stiffness matrix of the previous step is used
            //otherwise, the additional boundary displacements will result in an artifical localization in elements at the boundary
            curStrainFactor-=deltaStrainFactor;
            curEngineeringStrain = prevStrain + deltaEngineeringStrain * curStrainFactor;

            //set the total strain and calculate from the existing crack opening and the total strain the homogeneous strain
            fineScaleStructure->SetTotalEngineeringStrain(curEngineeringStrain);

            // build global dof numbering
            fineScaleStructure->NodeBuildGlobalDofs();

            //set previous converged displacements
            fineScaleStructure->NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);
            fineScaleStructure->ElementTotalUpdateTmpStaticData();

            //decrease load step
            deltaStrainFactor*=0.5;
            curStrainFactor+=deltaStrainFactor;
            curEngineeringStrain = prevStrain + deltaEngineeringStrain * curStrainFactor;

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
            //set the total strain and calculate from the existing crack opening and the total strain the homogeneous strain
            fineScaleStructure->SetTotalEngineeringStrain(curEngineeringStrain);

            // build global dof numbering
            fineScaleStructure->NodeBuildGlobalDofs();

            //update stiffness in order to calculate new dispForceVector
            fineScaleStructure->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
            int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
            int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
            //std::cout << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << std::endl;

            //update rhs vector for next Newton iteration
            rhsVector = dispForceVector + extForceVector - intForceVector;

            //update displacements of all nodes according to the new conre mat
            NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
            NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
            fineScaleStructure->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
            fineScaleStructure->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
            fineScaleStructure->ElementTotalUpdateTmpStaticData();

            // calculate initial residual for next load step
            fineScaleStructure->BuildGlobalGradientInternalPotentialVector(intForceVector);

            //update rhs vector for next Newton iteration
            rhsVector = dispForceVector + extForceVector - intForceVector;
        }
    }
}
///////////////////////////////////////////////////////////////////////////

//! @brief ... get type of constitutive relationship
//! @return ... type of constitutive relationship
//! @sa eConstitutiveType
NuTo::Constitutive::eConstitutiveType NuTo::Multiscale::GetType() const
{
    return NuTo::Constitutive::MULTISCALE;
}

//! @brief ... check compatibility between element type and type of constitutive relationship
//! @param rElementType ... element type
//! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
bool NuTo::Multiscale::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::PLANE2D3N:
        return true;
    case NuTo::Element::PLANE2D4N:
        return true;
    case NuTo::Element::PLANE2D6N:
        return true;
    default:
        return false;
    }
}

//! @brief ... print information about the object
//! @param rVerboseLevel ... verbosity of the information
void NuTo::Multiscale::Info(unsigned short rVerboseLevel) const
{
}

// check parameters
void NuTo::Multiscale::CheckParameters()const
{
}

//! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
//! @return ... see brief explanation
bool NuTo::Multiscale::HaveTmpStaticData() const
{
    return false;
}

//! @brief ... returns true, if a material model has is nonlocal (stiffness is of dynamic size, nonlocal averaging)
//! @return ... see brief explanation
bool NuTo::Multiscale::IsNonlocalModel()const
{
    return false;
}

bool NuTo::Multiscale::CheckStiffness(NuTo::StructureIp* rFineScaleStructure)const
{
    std::cout << "test of stiffness still included " << std::endl;
    NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixCSRVector2;
    NuTo::FullMatrix<double> dispForceVector;

    //recalculate stiffness
    rFineScaleStructure->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
    rFineScaleStructure->ConstraintInfo(10);

    NuTo::FullMatrix<double> stiffnessMatrixCSRVector2Full(stiffnessMatrixCSRVector2);
    //std::cout<<"stiffness matrix" << std::endl;
    //stiffnessMatrixCSRVector2Full.Info(10,3);
    double interval(1e-9);
    NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
    NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
    NuTo::FullMatrix<double> stiffnessMatrixCSRVector2_CDF(stiffnessMatrixCSRVector2.GetNumRows(), stiffnessMatrixCSRVector2.GetNumColumns());
    NuTo::FullMatrix<double> intForceVector1, intForceVector2, intForceVectorCDF(stiffnessMatrixCSRVector2.GetNumRows(),1);
    double energy1,energy2;
    rFineScaleStructure->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
    rFineScaleStructure->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
    rFineScaleStructure->ElementTotalUpdateTmpStaticData();
    rFineScaleStructure->BuildGlobalGradientInternalPotentialVector(intForceVector1);
    std::cout << "check stiffness:: intForceVector1"<< std::endl;
    intForceVector1.Trans().Info(12,3);
    energy1 = rFineScaleStructure->ElementTotalGetTotalEnergy();
    energy1 += rFineScaleStructure->ConstraintTotalGetTotalEnergy();
    for (int count=0; count<displacementsActiveDOFsCheck.GetNumRows(); count++)
    {
        displacementsActiveDOFsCheck(count,0)+=interval;
        rFineScaleStructure->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
        rFineScaleStructure->ElementTotalUpdateTmpStaticData();
        rFineScaleStructure->BuildGlobalGradientInternalPotentialVector(intForceVector2);
        std::cout << "check stiffness:: intForceVector2"<< std::endl;
        intForceVector2.Trans().Info(12,3);
        rFineScaleStructure->ConstraintInfo(10);
        energy2 = rFineScaleStructure->ElementTotalGetTotalEnergy();
        energy2 += rFineScaleStructure->ConstraintTotalGetTotalEnergy();
        stiffnessMatrixCSRVector2_CDF.SetColumn(count,(intForceVector2-intForceVector1)*(1./interval));
        intForceVectorCDF(count,0) = (energy2-energy1)/interval;
        displacementsActiveDOFsCheck(count,0)-=interval;
    }
    rFineScaleStructure->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
    if ((stiffnessMatrixCSRVector2_CDF-stiffnessMatrixCSRVector2Full).Abs().Max()>1e-3)
    {
        std::cout << "globalStiffnessMatrix algo" << std::endl;
        stiffnessMatrixCSRVector2Full.Info(10,3);
        std::cout<< std::endl << "globalStiffnessMatrix cdf" << std::endl;
        stiffnessMatrixCSRVector2_CDF.Info(10,3);
        std::cout<< std::endl << "error" << std::endl;
        (stiffnessMatrixCSRVector2_CDF-stiffnessMatrixCSRVector2Full).Info(10);
        std::cout << "maximum error is " << (stiffnessMatrixCSRVector2_CDF-stiffnessMatrixCSRVector2Full).Abs().Max() << std::endl;
        std::cout<< std::endl << "intForceVector algo" << std::endl;
        intForceVector1.Trans().Info(10);
        std::cout<< std::endl << "intForceVector cdf" << std::endl;
        intForceVectorCDF.Trans().Info(10);
        //throw MechanicsException("[NuTo::Multiscale::Solve] Stiffness matrix is not correct.");
        std::cout << "stiffness is wrong!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<< std::endl;
        return false;
    }
    else
    {
        std::cout << "stiffness is OK "<< std::endl;
        return true;
    }

}
