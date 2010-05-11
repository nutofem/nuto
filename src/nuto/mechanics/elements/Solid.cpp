#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal6x6.h"
#include "nuto/mechanics/elements/Solid.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveEngineeringStressStrain.h"

//! @brief constructor
NuTo::Solid::Solid(const StructureBase* rStructure, ElementData::eElementDataType rElementDataType,
		IntegrationType::eIntegrationType rIntegrationType, IpData::eIpDataType rIpDataType) :
        NuTo::ElementWithDataBase::ElementWithDataBase(rStructure, rElementDataType, rIntegrationType, rIpDataType)
{
}

//! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the stiffness matrix
void NuTo::Solid::CalculateCoefficientMatrix_0(NuTo::FullMatrix<double>& rCoefficientMatrix,
        std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn, bool& rSymmetry)const
{
    //calculate coordinates
    std::vector<double> nodeCoord(GetNumDofs());
    CalculateCoordinates(nodeCoord);

    //calculate local displacements
    std::vector<double> nodeDisp(GetNumDofs());
    CalculateDisplacements(nodeDisp);

    //allocate space for local ip coordinates
    double localIPCoord[3];

    //allocate space for derivatives of shape functions
    std::vector<double> derivativeShapeFunctionsLocal(GetNumDofs());
    std::vector<double> derivativeShapeFunctionsGlobal(GetNumDofs());

    //allocate deformation gradient
    DeformationGradient3D deformationGradient;

    //allocate tangent
    ConstitutiveTangentLocal6x6 tangent;

    //InvJacobian and determinant of Jacobian
    double invJacobian[9], detJac;

    //material pointer
    const ConstitutiveEngineeringStressStrain *constitutivePtr;

    //allocate and initialize result matrix
    rCoefficientMatrix.Resize(GetNumDofs(),GetNumDofs());
    bool areAllIpsSymmetric=(true);
    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, localIPCoord);

        CalculateDerivativeShapeFunctionsLocal(localIPCoord, derivativeShapeFunctionsLocal);

        CalculateJacobian(derivativeShapeFunctionsLocal,nodeCoord, invJacobian, detJac);

        CalculateDerivativeShapeFunctionsGlobal(derivativeShapeFunctionsLocal,invJacobian,
                                                derivativeShapeFunctionsGlobal);

        // determine deformation gradient from the local Displacements and the derivative of the shape functions
        // this is not included in the AddIpStiffness to avoid reallocation of the deformation gradient for each IP
        CalculateDeformationGradient(derivativeShapeFunctionsGlobal, nodeDisp, deformationGradient);

        //call material law to calculate tangent
        constitutivePtr = dynamic_cast<const ConstitutiveEngineeringStressStrain*>(GetConstitutiveLaw(theIP));
        if (constitutivePtr==0)
            throw MechanicsException("[NuTo::Solid::CalculateCoefficientMatrix_0] Constitutive law can not deal with engineering stresses and strains");
        constitutivePtr->GetTangent_EngineeringStress_EngineeringStrain(this, theIP,
                deformationGradient, tangent);
        areAllIpsSymmetric &= tangent.GetSymmetry();

        // calculate element stiffness matrix
        // don't forget to include determinant of the Jacobian and area
        double factor(fabs(detJac*(mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP))));
        AddDetJBtCB(derivativeShapeFunctionsGlobal, tangent, factor, rCoefficientMatrix);
    }
    rSymmetry = areAllIpsSymmetric;

    // calculate list of global dofs related to the entries in the element stiffness matrix
    this->CalculateGlobalRowDofs(rGlobalDofsRow);
    this->CalculateGlobalColumnDofs(rGlobalDofsColumn);
}

//! @brief adds to a matrix the product B^tCB, where B contains the derivatives of the shape functions and C is the constitutive tangent
//! eventually include also area/width of an element
//! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param ConstitutiveTangentBase constitutive tangent matrix
//! @param rFactor factor including determinant of Jacobian and IP weight
//! @param rCoefficientMatrix to be added to
void NuTo::Solid::AddDetJBtCB(const std::vector<double>& rDerivativeShapeFunctionsGlobal,
                              const ConstitutiveTangentLocal6x6& rConstitutiveTangent, double rFactor,
                              FullMatrix<double>& rCoefficientMatrix)const
{
    assert(rCoefficientMatrix.GetNumRows()==3*GetNumNodes() && rCoefficientMatrix.GetNumColumns()==3*GetNumNodes());
    assert((int)rDerivativeShapeFunctionsGlobal.size()==3*GetNumNodes());
    const double *C = rConstitutiveTangent.GetData();
    double x1,x2,y1,y2,z1,z2,x2x1,y2x1,z2x1,x2y1,y2y1,z2y1,x2z1,y2z1,z2z1;
    for (int theNode1=0; theNode1<GetNumNodes(); theNode1++)
    {
        int node1mul3 = 3*theNode1;
        int node1mul3plus1 = node1mul3+1;
        int node1mul3plus2 = node1mul3plus1+1;

        x1 = rFactor * rDerivativeShapeFunctionsGlobal[node1mul3];
        y1 = rFactor * rDerivativeShapeFunctionsGlobal[node1mul3plus1];
        z1 = rFactor * rDerivativeShapeFunctionsGlobal[node1mul3plus2];
        for (int theNode2=0; theNode2<GetNumNodes(); theNode2++)
        {
            int node2mul3 = 3*theNode2;
            int node2mul3plus1 = node2mul3+1;
            int node2mul3plus2 = node2mul3plus1+1;

            x2 = rDerivativeShapeFunctionsGlobal[node2mul3];
            y2 = rDerivativeShapeFunctionsGlobal[node2mul3plus1];
            z2 = rDerivativeShapeFunctionsGlobal[node2mul3plus2];

            x2x1 = x2*x1;
            y2x1 = y2*x1;
            z2x1 = z2*x1;
            x2y1 = x2*y1;
            y2y1 = y2*y1;
            z2y1 = z2*y1;
            x2z1 = x2*z1;
            y2z1 = y2*z1;
            z2z1 = z2*z1;

            rCoefficientMatrix(node1mul3,node2mul3)          +=x2x1*C[0] +x2y1*C[3] +x2z1*C[5] +y2x1*C[18]+y2y1*C[21]+y2z1*C[23]+z2x1*C[30]+z2y1*C[33]+z2z1*C[35];
            rCoefficientMatrix(node1mul3,node2mul3plus1)     +=y2x1*C[6] +y2y1*C[9] +y2z1*C[11]+x2x1*C[18]+x2y1*C[21]+x2z1*C[23]+z2x1*C[24]+z2y1*C[27]+z2z1*C[29];
            rCoefficientMatrix(node1mul3,node2mul3plus2)     +=z2x1*C[12]+z2y1*C[15]+z2z1*C[17]+y2x1*C[24]+y2y1*C[27]+y2z1*C[29]+x2x1*C[30]+x2y1*C[33]+x2z1*C[35];
            rCoefficientMatrix(node1mul3plus1,node2mul3)     +=x2y1*C[1] +x2x1*C[3] +x2z1*C[4] +y2y1*C[19]+y2x1*C[21]+y2z1*C[22]+z2y1*C[31]+z2x1*C[33]+z2z1*C[34];
            rCoefficientMatrix(node1mul3plus1,node2mul3plus1)+=y2y1*C[7] +y2x1*C[9] +y2z1*C[10]+x2y1*C[19]+x2x1*C[21]+x2z1*C[22]+z2y1*C[25]+z2x1*C[27]+z2z1*C[28];
            rCoefficientMatrix(node1mul3plus1,node2mul3plus2)+=z2y1*C[13]+z2x1*C[15]+z2z1*C[16]+y2y1*C[25]+y2x1*C[27]+y2z1*C[28]+x2y1*C[31]+x2x1*C[33]+x2z1*C[34];
            rCoefficientMatrix(node1mul3plus2,node2mul3)     +=x2z1*C[2] +x2y1*C[4] +x2x1*C[5] +y2z1*C[20]+y2y1*C[22]+y2x1*C[23]+z2z1*C[32]+z2y1*C[34]+z2x1*C[35];
            rCoefficientMatrix(node1mul3plus2,node2mul3plus1)+=y2z1*C[8] +y2y1*C[10]+y2x1*C[11]+x2z1*C[20]+x2y1*C[22]+x2x1*C[23]+z2z1*C[26]+z2y1*C[28]+z2x1*C[29];
            rCoefficientMatrix(node1mul3plus2,node2mul3plus2)+=z2z1*C[14]+z2y1*C[16]+z2x1*C[17]+y2z1*C[26]+y2y1*C[28]+y2x1*C[29]+x2z1*C[32]+x2y1*C[34]+x2x1*C[35];
        }
    }
}

//! @brief adds up the internal force vector
//! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param rEngineeringStress stress
//! @param factor factor including det Jacobian area and integration point weight
//! @param rResult resforce vector
void NuTo::Solid::AddDetJBtSigma(const std::vector<double>& rDerivativeShapeFunctionsGlobal,
                                 const EngineeringStress3D& rEngineeringStress,
                                 double rFactor,
                                 FullMatrix<double>& rResult)const
{
    assert(rResult.GetNumRows()==3*GetNumNodes() && rResult.GetNumColumns()==1);
    const double *s = rEngineeringStress.GetData();
    double x1,y1,z1;
    for (int theNode1=0; theNode1<GetNumNodes(); theNode1++)
    {
        int node1mul3 = 3*theNode1;
        int node1mul3plus1 = node1mul3+1;
        int node1mul3plus2 = node1mul3plus1+1;

        x1 = rFactor * rDerivativeShapeFunctionsGlobal[node1mul3];
        y1 = rFactor * rDerivativeShapeFunctionsGlobal[node1mul3plus1];
        z1 = rFactor * rDerivativeShapeFunctionsGlobal[node1mul3plus2];

        rResult(node1mul3,0)     +=x1*s[0]+y1*s[3]+z1*s[5];
        rResult(node1mul3plus1,0)+=y1*s[1]+x1*s[3]+z1*s[4];
        rResult(node1mul3plus2,0)+=z1*s[2]+y1*s[4]+x1*s[5];
    }
}

//! @brief Calculates the the inverse of the Jacobian and its determinant
//! @param rDerivativeShapeFunctions Derivatives of the shape functions (dN1dx, dN1dy, dN1dz, dN2dx, ..
//! @param rNodeCoordinates Node coordinates (X1,Y1,Z1,X2,Y2,Z2,...
//! @param rInvJacobian inverse Jacobian matrix (return value)
//! @param rDetJac determinant of the Jacobian (return value)
void NuTo::Solid::CalculateJacobian(const std::vector<double>& rDerivativeShapeFunctions,
                                    const std::vector<double>& rNodeCoordinates,
                                    double rInvJacobian[9],
                                    double& rDetJac)const
{
    /*	   jacobian
    	   j0, j1, j2
    	   j3, j4, j5
    	   j6, j7, j8*/

    assert((int)rDerivativeShapeFunctions.size()==3*GetNumNodes() && (int)rNodeCoordinates.size()==3*GetNumNodes());
    double  j0(0.),j1(0.),j2(0.),j3(0.),j4(0.),j5(0.),j6(0.),j7(0.),j8(0.),
    j48_57,j27_18,j15_24,x,y,z;

    int theDerivative(0);
    for (int count = 0; count < GetNumNodes(); count++)
    {
        x = rNodeCoordinates[theDerivative];
        y = rNodeCoordinates[theDerivative+1];
        z = rNodeCoordinates[theDerivative+2];

        j0 += rDerivativeShapeFunctions[theDerivative] * x;
        j3 += rDerivativeShapeFunctions[theDerivative] * y;
        j6 += rDerivativeShapeFunctions[theDerivative] * z;
        theDerivative++;

        j1 += rDerivativeShapeFunctions[theDerivative] * x;
        j4 += rDerivativeShapeFunctions[theDerivative] * y;
        j7 += rDerivativeShapeFunctions[theDerivative] * z;
        theDerivative++;

        j2 += rDerivativeShapeFunctions[theDerivative] * x;
        j5 += rDerivativeShapeFunctions[theDerivative] * y;
        j8 += rDerivativeShapeFunctions[theDerivative] * z;
        theDerivative++;
    }

    j48_57 = j4 * j8 - j5 * j7;
    j27_18 = j2 * j7 - j1 * j8;
    j15_24 = j1 * j5 - j2 * j4;
    /**********************************/
    rDetJac = j0 * j48_57 +j3 * j27_18 + j6 * j15_24;

    if (rDetJac==0)
        throw MechanicsException("[NuTo::Solid::CalculateJacobian] Determinant of the Jacobian is zero, no inversion possible.");

    if (rInvJacobian!=0)
    {
        double invDeterminant(1./rDetJac);
        rInvJacobian[0]=j48_57*invDeterminant;
        rInvJacobian[1]=j27_18*invDeterminant;
        rInvJacobian[2]=j15_24*invDeterminant;
        rInvJacobian[3]=(j5*j6-j3*j8)*invDeterminant;
        rInvJacobian[4]=(j0*j8-j2*j6)*invDeterminant;
        rInvJacobian[5]=(j2*j3-j0*j5)*invDeterminant;
        rInvJacobian[6]=(j3*j7-j4*j6)*invDeterminant;
        rInvJacobian[7]=(j1*j6-j0*j7)*invDeterminant;
        rInvJacobian[8]=(j0*j4-j1*j3)*invDeterminant;
    }
}

//! @brief calculates the derivative of the shape functions with respect to global coordinates
//! @param std::vector<double>& rDerivativeShapeFunctions derivatives of the shape functions
//! @param rJacInv inverse of the Jacobian
//! @param rDerivativeShapeFunctionsGlobal derivaties of the shape functions with respect to global coordinates
//! size should already be correct, but can be checked with an assert
//! first all the directions for a single node, and then for the next node
void NuTo::Solid::CalculateDerivativeShapeFunctionsGlobal(const std::vector<double>& rDerivativeShapeFunctionsLocal, const double rJacInv[9], std::vector<double>& rDerivativeShapeFunctionsGlobal)const
{
    assert(rDerivativeShapeFunctionsLocal.size()==rDerivativeShapeFunctionsGlobal.size());
    for (int count=0; count<GetNumNodes(); count++)
    {
        int mul3count = 3*count;
        int mul3countplus1 = mul3count+1;
        int mul3countplus2 = mul3countplus1+1;
        rDerivativeShapeFunctionsGlobal[mul3count] =
            rDerivativeShapeFunctionsLocal[mul3count]     *rJacInv[0]+
            rDerivativeShapeFunctionsLocal[mul3countplus1]*rJacInv[3]+
            rDerivativeShapeFunctionsLocal[mul3countplus2]*rJacInv[6];
        rDerivativeShapeFunctionsGlobal[mul3countplus1] =
            rDerivativeShapeFunctionsLocal[mul3count]     *rJacInv[1]+
            rDerivativeShapeFunctionsLocal[mul3countplus1]*rJacInv[4]+
            rDerivativeShapeFunctionsLocal[mul3countplus2]*rJacInv[7];
        rDerivativeShapeFunctionsGlobal[mul3countplus2] =
            rDerivativeShapeFunctionsLocal[mul3count]     *rJacInv[2]+
            rDerivativeShapeFunctionsLocal[mul3countplus1]*rJacInv[5]+
            rDerivativeShapeFunctionsLocal[mul3countplus2]*rJacInv[8];
    }
}
//! @brief calculates the gradient of the internal potential
//! for a mechanical problem, this corresponds to the internal force vector
void NuTo::Solid::CalculateGradientInternalPotential(NuTo::FullMatrix<double>& rResult,
        std::vector<int>& rGlobalDofs)const
{
    //calculate coordinates
    std::vector<double> nodeCoord(GetNumDofs());
    CalculateCoordinates(nodeCoord);

    //calculate local displacements
    std::vector<double> nodeDisp(GetNumDofs());
    CalculateDisplacements(nodeDisp);

    //allocate space for local ip coordinates
    double localIPCoord[3];

    //allocate space for derivatives of shape functions
    std::vector<double> derivativeShapeFunctionsLocal(GetNumDofs());
    std::vector<double> derivativeShapeFunctionsGlobal(GetNumDofs());

    //allocate deformation gradient
    DeformationGradient3D deformationGradient;

    //allocate global engineering stress
    EngineeringStress3D engineeringStress;

    //InvJacobian and determinant of Jacobian
    double invJacobian[9], detJac;

    //material pointer
    const ConstitutiveEngineeringStressStrain *constitutivePtr;

    //allocate and initialize result matrix
    rResult.Resize(GetNumDofs(),1);
    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, localIPCoord);

        CalculateDerivativeShapeFunctionsLocal(localIPCoord, derivativeShapeFunctionsLocal);

        CalculateJacobian(derivativeShapeFunctionsLocal,nodeCoord, invJacobian, detJac);

        CalculateDerivativeShapeFunctionsGlobal(derivativeShapeFunctionsLocal,invJacobian,
                                                derivativeShapeFunctionsGlobal);

        // determine deformation gradient from the local Displacements and the derivative of the shape functions
        // this is not included in the AddIpStiffness to avoid reallocation of the deformation gradient for each IP
        CalculateDeformationGradient(derivativeShapeFunctionsGlobal, nodeDisp, deformationGradient);

        //call constitutive law and calculate stress
        constitutivePtr = dynamic_cast<const ConstitutiveEngineeringStressStrain*>(GetConstitutiveLaw(theIP));
        if (constitutivePtr==0)
            throw MechanicsException("[NuTo::Solid::CalculateGradientInternalPotential] Constitutive law can not deal with engineering stresses and strains");
        constitutivePtr->GetEngineeringStressFromEngineeringStrain(this, theIP,
                deformationGradient, engineeringStress);

        //add to local resforce vector
        double factor(fabs(detJac*(mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP))));
        AddDetJBtSigma(derivativeShapeFunctionsGlobal,engineeringStress, factor, rResult);
    }
    // calculate list of global dofs related to the entries in the element stiffness matrix
    this->CalculateGlobalRowDofs(rGlobalDofs);
}

//! @brief Update the static data of an element
void NuTo::Solid::UpdateStaticData()
{
    //calculate coordinates
    std::vector<double> nodeCoord(GetNumDofs());
    CalculateCoordinates(nodeCoord);

    //calculate local displacements
    std::vector<double> nodeDisp(GetNumDofs());
    CalculateDisplacements(nodeDisp);

    //allocate space for local ip coordinates
    double localIPCoord[3];

    //allocate space for derivatives of shape functions
    std::vector<double> derivativeShapeFunctionsLocal(GetNumDofs());
    std::vector<double> derivativeShapeFunctionsGlobal(GetNumDofs());

    //allocate deformation gradient
    DeformationGradient3D deformationGradient;

    //InvJacobian and determinant of Jacobian
    double invJacobian[9], detJac;

    //material pointer
    const ConstitutiveEngineeringStressStrain *constitutivePtr;

    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, localIPCoord);

        CalculateDerivativeShapeFunctionsLocal(localIPCoord, derivativeShapeFunctionsLocal);

        CalculateJacobian(derivativeShapeFunctionsLocal,nodeCoord, invJacobian, detJac);

        CalculateDerivativeShapeFunctionsGlobal(derivativeShapeFunctionsLocal,invJacobian,
                                                derivativeShapeFunctionsGlobal);

        // determine deformation gradient from the local Displacements and the derivative of the shape functions
        // this is not included in the AddIpStiffness to avoid reallocation of the deformation gradient for each IP
        CalculateDeformationGradient(derivativeShapeFunctionsGlobal, nodeDisp, deformationGradient);

        //call material law to update static data
        constitutivePtr = dynamic_cast<const ConstitutiveEngineeringStressStrain*>(GetConstitutiveLaw(theIP));
        if (constitutivePtr==0)
            throw MechanicsException("[NuTo::Solid::UpdateStaticData] Constitutive law can not deal with engineering stresses and strains");
        constitutivePtr->UpdateStaticData_EngineeringStress_EngineeringStrain(this, theIP, deformationGradient);
    }
}

//! @brief calculates the deformation gradient in 3D
//! @param rRerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param rLocalDisp local displacements
//! @param rDeformationGradient (return value)
void NuTo::Solid::CalculateDeformationGradient(const std::vector<double>& rDerivativeShapeFunctionsGlobal,
        const std::vector<double>& rLocalDisp,
        DeformationGradient3D& rDeformationGradient)const
{
    assert((int)rLocalDisp.size()==3*GetNumNodes() && (int)rDerivativeShapeFunctionsGlobal.size()==3*GetNumNodes());

    rDeformationGradient.mDeformationGradient[0] = 1.;
    rDeformationGradient.mDeformationGradient[1] = 0.;
    rDeformationGradient.mDeformationGradient[2] = 0.;
    rDeformationGradient.mDeformationGradient[3] = 0.;
    rDeformationGradient.mDeformationGradient[4] = 1.;
    rDeformationGradient.mDeformationGradient[5] = 0.;
    rDeformationGradient.mDeformationGradient[6] = 0.;
    rDeformationGradient.mDeformationGradient[7] = 0.;
    rDeformationGradient.mDeformationGradient[8] = 1.;

    int theDisp(0);
    double dNdX,dNdY,dNdZ;
    for (int count=0; count<GetNumNodes(); count++)
    {
        dNdX = rDerivativeShapeFunctionsGlobal[theDisp];
        dNdY = rDerivativeShapeFunctionsGlobal[theDisp+1];
        dNdZ = rDerivativeShapeFunctionsGlobal[theDisp+2];

        rDeformationGradient.mDeformationGradient[0]+=rLocalDisp[theDisp]* dNdX;
        rDeformationGradient.mDeformationGradient[1]+=rLocalDisp[theDisp]* dNdY;
        rDeformationGradient.mDeformationGradient[2]+=rLocalDisp[theDisp]* dNdZ;
        theDisp++;
        rDeformationGradient.mDeformationGradient[3]+=rLocalDisp[theDisp]* dNdX;
        rDeformationGradient.mDeformationGradient[4]+=rLocalDisp[theDisp]* dNdY;
        rDeformationGradient.mDeformationGradient[5]+=rLocalDisp[theDisp]* dNdZ;
        theDisp++;
        rDeformationGradient.mDeformationGradient[6]+=rLocalDisp[theDisp]* dNdX;
        rDeformationGradient.mDeformationGradient[7]+=rLocalDisp[theDisp]* dNdY;
        rDeformationGradient.mDeformationGradient[8]+=rLocalDisp[theDisp]* dNdZ;
        theDisp++;
    }
}

//! @brief calculates the coefficient matrix for the 1-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the damping matrix
void NuTo::Solid::CalculateCoefficientMatrix_1(NuTo::FullMatrix<double>& rResult,
        std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn)const
{
    throw MechanicsException("[NuTo::Solid::CalculateCoefficientMatrix_1] to be implemented.");
}

//! @brief calculates the coefficient matrix for the 2-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the Mass matrix
void NuTo::Solid::CalculateCoefficientMatrix_2(NuTo::FullMatrix<double>& rResult,
        std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn)const
{
    throw MechanicsException("[NuTo::Solid::CalculateCoefficientMatrix_2] to be implemented.");
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point
//! @param rCoordinates local coordinates (return value)
void NuTo::Solid::GetLocalIntegrationPointCoordinates(int rIpNum, double rCoordinates[3])const
{
    this->mElementData->GetIntegrationType()->GetLocalIntegrationPointCoordinates3D(rIpNum, rCoordinates);
    return;
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point
//! @param rCoordinates local coordinates (return value)
void  NuTo::Solid::GetGlobalIntegrationPointCoordinates(int rIpNum, double rCoordinates[3])const
{
    double localCoordinates[3];
    double nodeCoordinates[3];
    this->mElementData->GetIntegrationType()->GetLocalIntegrationPointCoordinates3D(rIpNum, localCoordinates);
    std::vector<double> shapeFunctions(GetNumNodes());
    CalculateShapeFunctions(localCoordinates, shapeFunctions);
    rCoordinates[0] = 0.;
    rCoordinates[1] = 0.;
    rCoordinates[2] = 0.;

    nodeCoordinates[0] = 0;
    nodeCoordinates[1] = 0;
    nodeCoordinates[2] = 0;
    for (int theNode=0; theNode<GetNumNodes(); theNode++)
    {
        GetNode(theNode)->GetCoordinates3D(nodeCoordinates);
        rCoordinates[0]+=shapeFunctions[theNode]*nodeCoordinates[0];
        rCoordinates[1]+=shapeFunctions[theNode]*nodeCoordinates[1];
        rCoordinates[2]+=shapeFunctions[theNode]*nodeCoordinates[2];
    }
    return;
}


//! @brief calculates the engineering strain
//! @param rEngineerungStrain engineering strain (return value, always 6xnumIp matrix)
void NuTo::Solid::GetEngineeringStrain(FullMatrix<double>& rEngineeringStrain)const
{
    //calculate coordinates
    std::vector<double> nodeCoord(GetNumDofs());
    CalculateCoordinates(nodeCoord);

    //calculate local displacements
    std::vector<double> nodeDisp(GetNumDofs());
    CalculateDisplacements(nodeDisp);

    //allocate space for local ip coordinates
    double localIPCoord[3];

    //allocate space for derivatives of shape functions
    std::vector<double> derivativeShapeFunctionsLocal(GetNumDofs());
    std::vector<double> derivativeShapeFunctionsGlobal(GetNumDofs());

    //allocate deformation gradient
    DeformationGradient3D deformationGradient;

    //allocate global engineering stress
    EngineeringStrain3D engineeringStrain;

    //InvJacobian and determinant of Jacobian
    double invJacobian[9], detJac;

    //material pointer
    const ConstitutiveEngineeringStressStrain *constitutivePtr;

    //allocate and initialize result matrix
    rEngineeringStrain.Resize(6,GetNumIntegrationPoints());
    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, localIPCoord);

        CalculateDerivativeShapeFunctionsLocal(localIPCoord, derivativeShapeFunctionsLocal);

        CalculateJacobian(derivativeShapeFunctionsLocal,nodeCoord, invJacobian, detJac);

        CalculateDerivativeShapeFunctionsGlobal(derivativeShapeFunctionsLocal,invJacobian,
                                                derivativeShapeFunctionsGlobal);

        // determine deformation gradient from the local Displacements and the derivative of the shape functions
        // this is not included in the AddIpStiffness to avoid reallocation of the deformation gradient for each IP
        CalculateDeformationGradient(derivativeShapeFunctionsGlobal, nodeDisp, deformationGradient);

        //call material law to calculate engineering strain
        constitutivePtr = dynamic_cast<const ConstitutiveEngineeringStressStrain*>(GetConstitutiveLaw(theIP));
        if (constitutivePtr==0)
            throw MechanicsException("[NuTo::Solid::GetEngineeringStrain] Constitutive law can not deal with engineering stresses and strains");
        constitutivePtr->GetEngineeringStrain(this, theIP, deformationGradient, engineeringStrain);

        //copy to FullMatrix
        memcpy(&(rEngineeringStrain.mEigenMatrix.data()[theIP*6]),engineeringStrain.GetData(),6*sizeof(double));
    }
}

//! @brief calculates the engineering plastic strain
//! @param rEngineerungStrain engineering strain (return value, always 6xnumIp matrix)
void NuTo::Solid::GetEngineeringPlasticStrain(FullMatrix<double>& rEngineeringPlasticStrain)const
{
    //calculate coordinates
    std::vector<double> nodeCoord(GetNumDofs());
    CalculateCoordinates(nodeCoord);

    //calculate local displacements
    std::vector<double> nodeDisp(GetNumDofs());
    CalculateDisplacements(nodeDisp);

    //allocate space for local ip coordinates
    double localIPCoord[3];

    //allocate space for derivatives of shape functions
    std::vector<double> derivativeShapeFunctionsLocal(GetNumDofs());
    std::vector<double> derivativeShapeFunctionsGlobal(GetNumDofs());

    //allocate deformation gradient
    DeformationGradient3D deformationGradient;

    //allocate global engineering stress
    EngineeringStrain3D engineeringPlasticStrain;

    //InvJacobian and determinant of Jacobian
    double invJacobian[9], detJac;

    //material pointer
    const ConstitutiveEngineeringStressStrain *constitutivePtr;

    //allocate and initialize result matrix
    rEngineeringPlasticStrain.Resize(6,GetNumIntegrationPoints());
    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, localIPCoord);

        CalculateDerivativeShapeFunctionsLocal(localIPCoord, derivativeShapeFunctionsLocal);

        CalculateJacobian(derivativeShapeFunctionsLocal,nodeCoord, invJacobian, detJac);

        CalculateDerivativeShapeFunctionsGlobal(derivativeShapeFunctionsLocal,invJacobian,
                                                derivativeShapeFunctionsGlobal);

        // determine deformation gradient from the local Displacements and the derivative of the shape functions
        // this is not included in the AddIpStiffness to avoid reallocation of the deformation gradient for each IP
        CalculateDeformationGradient(derivativeShapeFunctionsGlobal, nodeDisp, deformationGradient);

        //call material law to calculate engineering strain
        constitutivePtr = dynamic_cast<const ConstitutiveEngineeringStressStrain*>(GetConstitutiveLaw(theIP));
        if (constitutivePtr==0)
            throw MechanicsException("[NuTo::Solid::GetEngineeringStrain] Constitutive law can not deal with engineering stresses and strains");
        constitutivePtr->GetEngineeringPlasticStrain(this, theIP, deformationGradient, engineeringPlasticStrain);

        //copy to FullMatrix
        memcpy(&(rEngineeringPlasticStrain.mEigenMatrix.data()[theIP*6]),engineeringPlasticStrain.GetData(),6*sizeof(double));
    }
}
//! @brief calculates the engineering stress
//! @param rEngineeringStress engineering stress (return value, always 6xnumIp matrix)
void NuTo::Solid::GetEngineeringStress(FullMatrix<double>& rEngineeringStress)const
{
    //calculate coordinates
    std::vector<double> nodeCoord(GetNumDofs());
    CalculateCoordinates(nodeCoord);

    //calculate local displacements
    std::vector<double> nodeDisp(GetNumDofs());
    CalculateDisplacements(nodeDisp);

    //allocate space for local ip coordinates
    double localIPCoord[3];

    //allocate space for derivatives of shape functions
    std::vector<double> derivativeShapeFunctionsLocal(GetNumDofs());
    std::vector<double> derivativeShapeFunctionsGlobal(GetNumDofs());

    //allocate deformation gradient
    DeformationGradient3D deformationGradient;

    //allocate global engineering stress
    EngineeringStress3D engineeringStress;

    //InvJacobian and determinant of Jacobian
    double invJacobian[9], detJac;

    //material pointer
    const ConstitutiveEngineeringStressStrain *constitutivePtr;

    //allocate and initialize result matrix
    rEngineeringStress.Resize(6,GetNumIntegrationPoints());
    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, localIPCoord);

        CalculateDerivativeShapeFunctionsLocal(localIPCoord, derivativeShapeFunctionsLocal);

        CalculateJacobian(derivativeShapeFunctionsLocal,nodeCoord, invJacobian, detJac);

        CalculateDerivativeShapeFunctionsGlobal(derivativeShapeFunctionsLocal,invJacobian,
                                                derivativeShapeFunctionsGlobal);

        // determine deformation gradient from the local Displacements and the derivative of the shape functions
        // this is not included in the AddIpStiffness to avoid reallocation of the deformation gradient for each IP
        CalculateDeformationGradient(derivativeShapeFunctionsGlobal, nodeDisp, deformationGradient);

        //call material law to calculate engineering strain
        constitutivePtr = dynamic_cast<const ConstitutiveEngineeringStressStrain*>(GetConstitutiveLaw(theIP));
        if (constitutivePtr==0)
            throw MechanicsException("[NuTo::Solid::GetEngineeringStress] Constitutive law can not deal with engineering stresses and strains");
        constitutivePtr->GetEngineeringStressFromEngineeringStrain(this, theIP, deformationGradient, engineeringStress);

        //copy to FullMatrix
        memcpy(&(rEngineeringStress.mEigenMatrix.data()[theIP*6]),engineeringStress.GetData(),6*sizeof(double));
    }
}

//! @brief Allocates static data for an integration point of an element
//! @param rConstitutiveLaw constitutive law, which is called to allocate the static data object
//! actually, both - the element type and the constitutive law are required to determine the static data object actually required
NuTo::ConstitutiveStaticDataBase* NuTo::Solid::AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw)const
{
	const ConstitutiveEngineeringStressStrain *constitutivePtr =
        dynamic_cast<const ConstitutiveEngineeringStressStrain*>(rConstitutiveLaw);
    if (constitutivePtr==0)
        throw MechanicsException("[NuTo::Solid::AllocateStaticData] Constitutive law can not deal with engineering stresses and strains");
    return constitutivePtr->AllocateStaticDataEngineeringStress_EngineeringStrain3D(this);
}

//! @brief stores the coordinates of the nodes
//! @param localCoordinates vector with already correct size allocated
//! this can be checked with an assertation
void NuTo::Solid::CalculateCoordinates(std::vector<double>& rCoordinates)const
{
    assert((int)rCoordinates.size()==3*GetNumNodes());
    for (int count=0; count<GetNumNodes(); count++)
    {
        GetNode(count)->GetCoordinates3D(&(rCoordinates[3*count]));
    }

}

//! @brief stores the displacements of the nodes
//! @param localDisplacements vector with already correct size allocated
//! this can be checked with an assertation
void NuTo::Solid::CalculateDisplacements(std::vector<double>& rDisplacements)const
{
    assert((int)rDisplacements.size()==3*GetNumNodes());
    for (int count=0; count<GetNumNodes(); count++)
    {
        GetNode(count)->GetDisplacements3D(&(rDisplacements[3*count]));
    }
}


// interpolate geometry
void NuTo::Solid::InterpolateCoordinatesFrom3D(double rLocalCoordinates[3], double rGlobalCoordinates[3]) const
{
    // calculate shape functions
    std::vector<double> ShapeFunctions(this->GetNumNodes());
    this->CalculateShapeFunctions(rLocalCoordinates, ShapeFunctions);

    // start interpolation
    rGlobalCoordinates[0] = 0.0;
    rGlobalCoordinates[1] = 0.0;
    rGlobalCoordinates[2] = 0.0;
    for (int NodeCount = 0; NodeCount < this->GetNumNodes(); NodeCount++)
    {
        // get node coordinate
        double NodeCoordinate[3];
        GetNode(NodeCount)->GetCoordinates3D(NodeCoordinate);

        // add node contribution
        rGlobalCoordinates[0] += ShapeFunctions[NodeCount] *  NodeCoordinate[0];
        rGlobalCoordinates[1] += ShapeFunctions[NodeCount] *  NodeCoordinate[1];
        rGlobalCoordinates[2] += ShapeFunctions[NodeCount] *  NodeCoordinate[2];
    }
}

// interpolate displacements
void NuTo::Solid::InterpolateDisplacementsFrom3D(double rLocalCoordinates[3], double rGlobalDisplacements[3]) const
{
    // calculate shape functions
    std::vector<double> ShapeFunctions(this->GetNumNodes());
    this->CalculateShapeFunctions(rLocalCoordinates, ShapeFunctions);

    // start interpolation
    rGlobalDisplacements[0] = 0.0;
    rGlobalDisplacements[1] = 0.0;
    rGlobalDisplacements[2] = 0.0;
    for (int NodeCount = 0; NodeCount < this->GetNumNodes(); NodeCount++)
    {
        // get node displacements
        double NodeDisplacement[3];
        GetNode(NodeCount)->GetDisplacements3D(NodeDisplacement);

        // add node contribution
        rGlobalDisplacements[0] += ShapeFunctions[NodeCount] *  NodeDisplacement[0];
        rGlobalDisplacements[1] += ShapeFunctions[NodeCount] *  NodeDisplacement[1];
        rGlobalDisplacements[2] += ShapeFunctions[NodeCount] *  NodeDisplacement[2];
    }
}


// build global row dofs
void NuTo::Solid::CalculateGlobalRowDofs(std::vector<int>& rGlobalRowDofs) const
{
    rGlobalRowDofs.resize(3 * this->GetNumNodes());
    for (int nodeCount = 0; nodeCount < this->GetNumNodes(); nodeCount++)
    {
        const NodeBase * nodePtr(GetNode(nodeCount));
        rGlobalRowDofs[3 * nodeCount    ] = nodePtr->GetDofDisplacement(0);
        rGlobalRowDofs[3 * nodeCount + 1] = nodePtr->GetDofDisplacement(1);
        rGlobalRowDofs[3 * nodeCount + 2] = nodePtr->GetDofDisplacement(2);
    }
}

// build global column dof
void NuTo::Solid::CalculateGlobalColumnDofs(std::vector<int>& rGlobalColumnDofs) const
{
    this->CalculateGlobalRowDofs(rGlobalColumnDofs);
}


// check element definition
void NuTo::Solid::CheckElement()
{
    // check nodes
    for (int nodeCount = 0; nodeCount < this->GetNumNodes(); nodeCount++)
    {
        if (GetNode(nodeCount)->GetNumCoordinates()!=3)
        {
            throw MechanicsException("[NuTo::Solid::CheckElement] invalid node type (check node definition for coordinates).");
        }
    }

    // check node ordering (element length must be positive) and for changing sign in jacobian determinant
    // calculate coordinates
    std::vector<double> nodeCoord(this->GetNumDofs());
    this->CalculateCoordinates(nodeCoord);

    // check number of integration points
    if (this->GetNumIntegrationPoints() < 1)
    {
        throw MechanicsException("[NuTo::Solid::CheckElement] invalid integration type.");
    }

    // check sign of the jacobian determinant of the first integration point
    double localIPCoord[3];
    this->GetLocalIntegrationPointCoordinates(0, localIPCoord);

    std::vector<double> derivativeShapeFunctionsLocal(GetNumDofs());
    this->CalculateDerivativeShapeFunctionsLocal(localIPCoord, derivativeShapeFunctionsLocal);

    double invJacobian[9], detJacobian;
    this->CalculateJacobian(derivativeShapeFunctionsLocal,nodeCoord, invJacobian, detJacobian);
    // reorder nodes if determinant is negative
    if (detJacobian < 0.0)
    {
        this->ReorderNodes();
        // recalculate node coordinates after renumbering
        this->CalculateCoordinates(nodeCoord);
    }

    // check jacobian determinant for all integration points for positive sign and calculate element volume
    double volume = 0;
    for (int ipCount = 0; ipCount < this->GetNumIntegrationPoints(); ipCount++)
    {
        // calculate jacobian determinant
        this->GetLocalIntegrationPointCoordinates(ipCount, localIPCoord);
        this->CalculateDerivativeShapeFunctionsLocal(localIPCoord, derivativeShapeFunctionsLocal);
        this->CalculateJacobian(derivativeShapeFunctionsLocal,nodeCoord, invJacobian, detJacobian);
        if (detJacobian <= 0)
        {
            throw MechanicsException("[NuTo::Solid::CheckElement] element is not properly defined by this nodes (zero or negative jacobian determinant).");
        }
        volume += this->GetIntegrationPointWeight(ipCount) * detJacobian;
    }

    // check element volume
    if (volume < 1e-14)
    {
        throw MechanicsException("[NuTo::Solid::CheckElement] element with zero volume (check nodes).");
    }
}

//! @brief calculates the volume of an integration point (weight * detJac)
//! @param rVolume  vector for storage of the ip volumes (area in 2D)
void NuTo::Solid::GetIntegrationPointVolume(std::vector<double>& rVolume)const
{
    //calculate coordinates
    std::vector<double> nodeCoord(GetNumDofs());
    CalculateCoordinates(nodeCoord);

    //allocate space for local ip coordinates
    double localIPCoord[3];

    //allocate space for derivatives of shape functions
    std::vector<double> derivativeShapeFunctionsLocal(GetNumDofs());

    //determinant of Jacobian
    double detJac;

    rVolume.resize(GetNumIntegrationPoints());

    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, localIPCoord);

        CalculateDerivativeShapeFunctionsLocal(localIPCoord, derivativeShapeFunctionsLocal);

        CalculateJacobian(derivativeShapeFunctionsLocal,nodeCoord, 0, detJac);

		rVolume[theIP] = detJac * mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP);
    }
}

