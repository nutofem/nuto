#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal1x1.h"
#include "nuto/mechanics/elements/Truss.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveEngineeringStressStrain.h"
#include "nuto/mechanics/sections/SectionBase.h"

//! @brief constructor
NuTo::Truss::Truss(const StructureBase* rStructure, ElementDataBase::eElementDataType rElementDataType, IntegrationTypeBase::eIntegrationType rIntegrationType) :
        NuTo::ElementWithDataBase::ElementWithDataBase(rStructure, rElementDataType, rIntegrationType)
{
    mSection = 0;
}

//! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the stiffness matrix
void NuTo::Truss::CalculateCoefficientMatrix_0(NuTo::FullMatrix<double>& rCoefficientMatrix,
        std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn, bool& rSymmetry)const
{
    //calculate local coordinates
    std::vector<double> localNodeCoord(GetNumLocalDofs());
    CalculateLocalCoordinates(localNodeCoord);

    //calculate local displacements
    std::vector<double> localNodeDisp(GetNumLocalDofs());
    CalculateLocalDisplacements(localNodeDisp);

    //allocate space for local ip coordinates
    double localIPCoord;

    //allocate space for local shape functions
    std::vector<double> derivativeShapeFunctions(GetLocalDimension()*GetNumShapeFunctions());

    //allocate deformation gradient
    DeformationGradient1D deformationGradient;

    //allocate deformation gradient
    ConstitutiveTangentLocal1x1 tangent;

    //material pointer
    const ConstitutiveEngineeringStressStrain *constitutivePtr;

    //allocate and initialize result matrix
    rCoefficientMatrix.Resize(GetNumLocalDofs(),GetNumLocalDofs());
    bool areAllIpsSymmetric=(true);
    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, localIPCoord);

        CalculateDerivativeShapeFunctions(localIPCoord, derivativeShapeFunctions);

        // determine deformation gradient from the local Displacements and the derivative of the shape functions
        // this is not included in the AddIpStiffness to avoid reallocation of the deformation gradient for each IP
        CalculateDeformationGradient(derivativeShapeFunctions, localNodeCoord, localNodeDisp, deformationGradient);

        //call material law to calculate tangent
        constitutivePtr = dynamic_cast<const ConstitutiveEngineeringStressStrain*>(GetConstitutiveLaw(theIP));
        if (constitutivePtr==0)
            throw MechanicsException("[NuTo::Solid::GetEngineeringStress] Constitutive law can not deal with engineering stresses and strains");
        constitutivePtr->GetTangent_EngineeringStress_EngineeringStrain(this, theIP,
                deformationGradient, tangent);
        areAllIpsSymmetric &= tangent.GetSymmetry();

        // calculate local stiffness matrix
        // don't forget to include determinant of the Jacobian and area
        // theoretically, the factor  is
        // detJ * area*BtCB, but B includes 1/detJ, which finally gives:
        double factor (mSection->GetArea()
                       /DetJacobian(derivativeShapeFunctions,localNodeCoord)
                       *(mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP)));

        AddDetJBtCB(derivativeShapeFunctions,tangent, factor, rCoefficientMatrix);
    }

    // eventually blow local matrix to full matrix - only relevant for
    // truss in 2D and 3D
    BlowLocalMatrixToGlobal(rCoefficientMatrix);

    // symmetry flag
    rSymmetry = areAllIpsSymmetric;

    // calculate list of global dofs related to the entries in the element stiffness matrix
    this->CalculateGlobalRowDofs(rGlobalDofsRow);
    this->CalculateGlobalColumnDofs(rGlobalDofsColumn);
}

//! @brief calculates the gradient of the internal potential
//! for a mechanical problem, this corresponds to the internal force vector
void NuTo::Truss::CalculateGradientInternalPotential(NuTo::FullMatrix<double>& rResult,
        std::vector<int>& rGlobalDofs)const
{
    //calculate local coordinates
    std::vector<double> localNodeCoord(GetNumLocalDofs());
    CalculateLocalCoordinates(localNodeCoord);

    //calculate local displacements
    std::vector<double> localNodeDisp(GetNumLocalDofs());
    CalculateLocalDisplacements(localNodeDisp);

    //allocate space for local ip coordinates
    double localIPCoord;

    //allocate space for local shape functions
    std::vector<double> derivativeShapeFunctions(GetLocalDimension()*GetNumShapeFunctions());

    //allocate deformation gradient
    DeformationGradient1D deformationGradient;

    //allocate global engineering stress
    EngineeringStress1D engineeringStress;

    //material pointer
    const ConstitutiveEngineeringStressStrain *constitutivePtr;

    //allocate and initialize result matrix
    rResult.Resize(GetNumLocalDofs(),1);
    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, localIPCoord);

        CalculateDerivativeShapeFunctions(localIPCoord, derivativeShapeFunctions);

        // determine deformation gradient from the local Displacements and the derivative of the shape functions
        // this is not included in the AddIpStiffness to avoid reallocation of the deformation gradient for each IP
        CalculateDeformationGradient(derivativeShapeFunctions, localNodeCoord, localNodeDisp, deformationGradient);

        //call material law to calculate stress
        constitutivePtr = dynamic_cast<const ConstitutiveEngineeringStressStrain*>(GetConstitutiveLaw(theIP));
        if (constitutivePtr==0)
            throw MechanicsException("[NuTo::Truss::CalculateGradientInternalPotential] Constitutive law can not deal with engineering stresses and strains");
        constitutivePtr->GetEngineeringStressFromEngineeringStrain(this, theIP,
                deformationGradient, engineeringStress);

        // calculate local stiffness matrix
        // theoretically, the factor  is
        // detJ * area*BtSigma, but B includes 1/detJ, which finally gives:
        double factor (mSection->GetArea()
                       *(mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP)));

        AddDetJBtSigma(derivativeShapeFunctions,engineeringStress, factor, rResult);
    }

    // eventually blow local matrix to full matrix - only relevant for
    // truss in 2D and 3D
    BlowLocalVectorToGlobal(rResult);

    // calculate list of global dofs related to the entries in the element stiffness matrix
    this->CalculateGlobalRowDofs(rGlobalDofs);
}
//! @brief Update the static data of an element
void NuTo::Truss::UpdateStaticData()
{
    //calculate local coordinates
    std::vector<double> localNodeCoord(GetNumLocalDofs());
    CalculateLocalCoordinates(localNodeCoord);

    //calculate local displacements
    std::vector<double> localNodeDisp(GetNumLocalDofs());
    CalculateLocalDisplacements(localNodeDisp);

    //allocate space for local ip coordinates
    double localIPCoord;

    //allocate space for local shape functions
    std::vector<double> derivativeShapeFunctions(GetLocalDimension()*GetNumShapeFunctions());

    //allocate deformation gradient
    DeformationGradient1D deformationGradient;

    //material pointer
    const ConstitutiveEngineeringStressStrain *constitutivePtr;

    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, localIPCoord);

        CalculateDerivativeShapeFunctions(localIPCoord, derivativeShapeFunctions);

        // determine deformation gradient from the local Displacements and the derivative of the shape functions
        // this is not included in the AddIpStiffness to avoid reallocation of the deformation gradient for each IP
        CalculateDeformationGradient(derivativeShapeFunctions, localNodeCoord, localNodeDisp, deformationGradient);

        //call material law to update static data
        constitutivePtr = dynamic_cast<const ConstitutiveEngineeringStressStrain*>(GetConstitutiveLaw(theIP));
        if (constitutivePtr==0)
            throw MechanicsException("[NuTo::Truss::CalculateGradientInternalPotential] Constitutive law can not deal with engineering stresses and strains");
        constitutivePtr->UpdateStaticData_EngineeringStress_EngineeringStrain(this, theIP, deformationGradient);
    }
}


//! @brief calculates deformation gradient1D
//! @param rRerivativeShapeFunctions derivatives of the shape functions
//! @param rLocalDisp local displacements
//! @param rConstitutiveInput (return value)
void NuTo::Truss::CalculateDeformationGradient(const std::vector<double>& rDerivativeShapeFunctions,
        const std::vector<double>& rLocalCoord,
        const std::vector<double>& rLocalDisp,
        DeformationGradient1D& rDeformationGradient)const
{
    assert((int)rLocalDisp.size()==GetNumNodes() && (int)rDerivativeShapeFunctions.size()==GetNumNodes());
    rDeformationGradient.mDeformationGradient = 0;

    //normally, the inverse Jacobian should be calculated, but for a truss element, it is sufficient to use the inverse of the Jacobian determinant
    double factor(1./DetJacobian(rDerivativeShapeFunctions, rLocalCoord));
    for (int count=0; count<GetNumNodes(); count++)
    {
        rDeformationGradient.mDeformationGradient+=rLocalDisp[count]*rDerivativeShapeFunctions[count];
    }
    rDeformationGradient.mDeformationGradient*=factor;
    rDeformationGradient.mDeformationGradient+=1.;
}

//! @brief calculates the coefficient matrix for the 1-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the damping matrix
void NuTo::Truss::CalculateCoefficientMatrix_1(NuTo::FullMatrix<double>& rResult,
        std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn)const
{
    throw MechanicsException("[NuTo::Truss::CalculateCoefficientMatrix_1] to be implemented.");
}

//! @brief calculates the coefficient matrix for the 2-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the Mass matrix
void NuTo::Truss::CalculateCoefficientMatrix_2(NuTo::FullMatrix<double>& rResult,
        std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn)const
{
    throw MechanicsException("[NuTo::Truss::CalculateCoefficientMatrix_2] to be implemented.");
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point
//! @param rCoordinates local coordinates (return value)
void NuTo::Truss::GetLocalIntegrationPointCoordinates(int rIpNum, double& rCoordinates)const
{
    this->mElementData->GetIntegrationType()->GetLocalIntegrationPointCoordinates1D(rIpNum, rCoordinates);
    return;
}


//! @brief calculates the engineering strain
//! @param rEngineerungStrain engineering strain (return value, always 6xnumIp matrix)
void NuTo::Truss::GetEngineeringStrain(FullMatrix<double>& rEngineeringStrain)const
{
    //calculate local coordinates
    std::vector<double> localNodeCoord(GetNumLocalDofs());
    CalculateLocalCoordinates(localNodeCoord);

    //calculate local displacements
    std::vector<double> localNodeDisp(GetNumLocalDofs());
    CalculateLocalDisplacements(localNodeDisp);

    //allocate space for local ip coordinates
    double localIPCoord;

    //allocate space for local shape functions
    std::vector<double> derivativeShapeFunctions(GetLocalDimension()*GetNumShapeFunctions());

    //allocate deformation gradient
    DeformationGradient1D deformationGradient;

    //allocate global engineering strain
    EngineeringStrain3D engineeringStrain;

    //material pointer
    const ConstitutiveEngineeringStressStrain *constitutivePtr;

    //allocate and initialize result matrix
    rEngineeringStrain.Resize(6,GetNumIntegrationPoints());
    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, localIPCoord);

        CalculateDerivativeShapeFunctions(localIPCoord, derivativeShapeFunctions);

        // determine deformation gradient from the local Displacements and the derivative of the shape functions
        CalculateDeformationGradient(derivativeShapeFunctions, localNodeCoord, localNodeDisp, deformationGradient);

        //call material law to calculate engineering strain
        constitutivePtr = dynamic_cast<const ConstitutiveEngineeringStressStrain*>(GetConstitutiveLaw(theIP));
        if (constitutivePtr==0)
            throw MechanicsException("[NuTo::Truss::GetEngineeringStrain] Constitutive law can not deal with engineering stresses and strains");
        constitutivePtr->GetEngineeringStrain(this, theIP, deformationGradient, engineeringStrain);

        //copy to FullMatrix
        memcpy(&(rEngineeringStrain.mEigenMatrix.data()[theIP*6]),engineeringStrain.GetData(),6*sizeof(double));
    }
}
//! @brief calculates the engineering plastic strain
//! @param rEngineerungStrain engineering strain (return value, always 6xnumIp matrix)
void NuTo::Truss::GetEngineeringPlasticStrain(FullMatrix<double>& rEngineeringPlasticStrain)const
{
	throw MechanicsException("[NuTo::Truss::GetEngineeringPlasticStrain] To be implemented.");
}

//! @brief calculates the engineering stress
//! @param rEngineeringStress engineering stress (return value, always 6xnumIp matrix)
void NuTo::Truss::GetEngineeringStress(FullMatrix<double>& rEngineeringStress)const
{
    //calculate local coordinates
    std::vector<double> localNodeCoord(GetNumLocalDofs());
    CalculateLocalCoordinates(localNodeCoord);

    //calculate local displacements
    std::vector<double> localNodeDisp(GetNumLocalDofs());
    CalculateLocalDisplacements(localNodeDisp);

    //allocate space for local ip coordinates
    double localIPCoord;

    //allocate space for local shape functions
    std::vector<double> derivativeShapeFunctions(GetLocalDimension()*GetNumShapeFunctions());

    //allocate deformation gradient
    DeformationGradient1D deformationGradient;

    //allocate global engineering strain
    EngineeringStress3D engineeringStress;

    //material pointer
    const ConstitutiveEngineeringStressStrain *constitutivePtr;

    //allocate and initialize result matrix
    rEngineeringStress.Resize(6,GetNumIntegrationPoints());
    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, localIPCoord);

        CalculateDerivativeShapeFunctions(localIPCoord, derivativeShapeFunctions);

        // determine deformation gradient from the local Displacements and the derivative of the shape functions
        CalculateDeformationGradient(derivativeShapeFunctions, localNodeCoord, localNodeDisp, deformationGradient);

        //call material law to calculate engineering strain
        constitutivePtr = dynamic_cast<const ConstitutiveEngineeringStressStrain*>(GetConstitutiveLaw(theIP));
        if (constitutivePtr==0)
            throw MechanicsException("[NuTo::Truss::GetEngineeringStress] Constitutive law can not deal with engineering stresses and strains");
        constitutivePtr->GetEngineeringStressFromEngineeringStrain(this, theIP, deformationGradient, engineeringStress);

        //copy to FullMatrix
        memcpy(&(rEngineeringStress.mEigenMatrix.data()[theIP*6]),engineeringStress.GetData(),6*sizeof(double));
    }
}
//! @brief sets the section of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need a section
//! @param rSection pointer to section
//! @return pointer to constitutive law
void NuTo::Truss::SetSection(const SectionBase* rSection)
{
    mSection = rSection;
}

//! @brief returns a pointer to the section of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need a section
//! @return pointer to section
const NuTo::SectionBase* NuTo::Truss::GetSection()const
{
    return mSection;
}

//! @brief Allocates static data for an integration point of an element
//! @param rConstitutiveLaw constitutive law, which is called to allocate the static data object
//! actually, both - the element type and the constitutive law are required to determine the static data object actually required
NuTo::ConstitutiveStaticDataBase* NuTo::Truss::AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw)const
{
    const ConstitutiveEngineeringStressStrain* constitutivePtr = dynamic_cast<const ConstitutiveEngineeringStressStrain*>(rConstitutiveLaw);
    if (constitutivePtr==0)
        throw MechanicsException("[NuTo::Solid::AllocateStaticData] Constitutive law can not deal with engineering stresses and strains");
    return constitutivePtr->AllocateStaticDataEngineeringStress_EngineeringStrain1D(this);
}

//! @brief returns determinant of the Jacobian
//! @param derivativeShapeFunctions derivatives of the shape functions
//! @param localCoord local coordinates
//! @return determinant of the Jacobian
double NuTo::Truss::DetJacobian(const std::vector<double>& derivativeShapeFunctions,const std::vector<double>& localCoord)const
{
    assert((int)localCoord.size()==GetNumNodes() && (int)derivativeShapeFunctions.size()==GetNumNodes());
    double detJ(0);
    for (int count=0; count<GetNumNodes(); count++)
        detJ+=derivativeShapeFunctions[count]*localCoord[count];
    return detJ;
}

//! @brief adds to a matrix the product B^tCB, where B contains the derivatives of the shape functions and C is the constitutive tangent
//! eventually include also area/width of an element
//! @param rDerivativeShapeFunctions derivatives of the shape functions
//! @param ConstitutiveTangentBase constitutive tangent matrix
//! @param rFactor factor including area, determinant of Jacobian and IP weight
//! @param rCoefficientMatrix to be added to
void NuTo::Truss::AddDetJBtCB(const std::vector<double>& rDerivativeShapeFunctions,
                              const ConstitutiveTangentLocal1x1& rConstitutiveTangent, double rFactor,
                              FullMatrix<double>& rCoefficientMatrix)const
{
    rFactor *=rConstitutiveTangent.GetData()[0];
    for (int node1=0; node1<GetNumNodes(); node1++)
    {
        for (int node2=0; node2<GetNumNodes(); node2++)
        {
            rCoefficientMatrix(node1,node2)+=rFactor*rDerivativeShapeFunctions[node1]*rDerivativeShapeFunctions[node2];
        }
    }
}
//! @brief adds up the internal force vector
//! @param derivativeShapeFunctions derivatives of the shape functions
//! @param rEngineeringStress stress
//! @param factor factor including det Jacobian area and integration point weight
//! @param rResult resforce vector
void NuTo::Truss::AddDetJBtSigma(const std::vector<double>& rDerivativeShapeFunctions,
                                 const EngineeringStress1D& rEngineeringStress, double rFactor, FullMatrix<double>& rResult)const
{
    rFactor*=rEngineeringStress.GetData()[0];
    for (int node1=0; node1<GetNumNodes(); node1++)
    {
        rResult(node1,0)+=rFactor*rDerivativeShapeFunctions[node1];
    }
}

