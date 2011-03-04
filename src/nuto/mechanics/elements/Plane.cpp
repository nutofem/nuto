// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentNonlocal3x3.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal3x3.h"
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/elements/Plane.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveEngineeringStressStrain.h"
#include "nuto/mechanics/sections/SectionBase.h"

#include "nuto/math/FullMatrix.h"

//! @brief constructor
NuTo::Plane::Plane(const StructureBase* rStructure, ElementData::eElementDataType rElementDataType,
        IntegrationType::eIntegrationType rIntegrationType, IpData::eIpDataType rIpDataType) :
        NuTo::ElementBase::ElementBase(rStructure, rElementDataType, rIntegrationType, rIpDataType)
{
    mSection = 0;
}

//! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the stiffness matrix
void NuTo::Plane::CalculateCoefficientMatrix_0(NuTo::FullMatrix<double>& rCoefficientMatrix,
        std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn, bool& rSymmetry)const
{
    //calculate local coordinates
    std::vector<double> localNodeCoord(GetNumLocalDofs());
    CalculateLocalCoordinates(localNodeCoord);
    //calculate local displacements
    std::vector<double> localNodeDisp(GetNumLocalDofs());
    CalculateLocalDisplacements(localNodeDisp);

    //allocate space for ip coordinates in natural coordinate system (-1,1)
    double naturalIPCoord[2];
    double nonlocalNaturalIPCoord[2];

    //allocate space for derivatives of shape functions in natural coordinate system
    std::vector<double> derivativeShapeFunctionsNatural(2*GetNumShapeFunctions());
    //allocate space for derivatives of shape functions in local coordinate system
    std::vector<double> derivativeShapeFunctionsLocal(2*GetNumShapeFunctions());

    std::vector<double> nonlocalDerivativeShapeFunctionsNatural;
    std::vector<double> nonlocalDerivativeShapeFunctionsLocal;

    //allocate deformation gradient
    DeformationGradient2D deformationGradient;

    //allocate vector of tangent matrices
    int NumNonlocalElements(GetNumNonlocalElements());
    // in case of a local formulation, just use a nonlocal matrix with 1 entry

    // calculate number of DOFS involved in all nonlocal elements
    int NumNonlocalIps;
    if (NumNonlocalElements==0)
    {
        NumNonlocalIps = 1;
        //allocate and initialize result matrix
        rCoefficientMatrix.Resize(GetNumLocalDofs(),GetNumLocalDofs());
    }
    else
    {
        //allocate and initialize result matrix
        const std::vector<const ElementBase*>& nonlocalElements(GetNonlocalElements());
        int NumCols(0);
        NumNonlocalIps=0;
        for (int theNonlocalElement=0; theNonlocalElement<NumNonlocalElements; theNonlocalElement++)
        {
            NumCols += nonlocalElements[theNonlocalElement]->AsPlane()->GetNumLocalDofs();
            NumNonlocalIps += nonlocalElements[theNonlocalElement]->GetNumIntegrationPoints();
        }
        rCoefficientMatrix.Resize(GetNumLocalDofs(),NumCols);
    }
    ConstitutiveTangentNonlocal3x3 NonlocalTangent(NumNonlocalIps);

    //InvJacobian and determinant of Jacobian
    double invJacobian[4], detJac;

    //InvJacobian and determinant of Jacobian for nonlocal model
    double nonlocalInvJacobian[4], nonlocalDetJac;

    //material pointer
    const ConstitutiveEngineeringStressStrain *constitutivePtr;


    bool areAllIpsSymmetric(true);
    bool areAllIpsLocal(true);
    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, naturalIPCoord);

        CalculateDerivativeShapeFunctionsNatural(naturalIPCoord, derivativeShapeFunctionsNatural);

        CalculateJacobian(derivativeShapeFunctionsNatural,localNodeCoord, invJacobian, detJac);

        CalculateDerivativeShapeFunctionsLocal(derivativeShapeFunctionsNatural,invJacobian,
                                                derivativeShapeFunctionsLocal);

        // determine deformation gradient from the local Displacements and the derivative of the shape functions
        // this is not included in the AddIpStiffness to avoid reallocation of the deformation gradient for each IP
        CalculateDeformationGradient(derivativeShapeFunctionsLocal, localNodeDisp, deformationGradient);

        // get material law to calculate tangent
        constitutivePtr = GetConstitutiveLaw(theIP)->AsConstitutiveEngineeringStressStrain();

        //factor for the numerical integration
        assert(mSection->GetThickness()>0);
        double factor(mSection->GetThickness()*detJac*(mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP)));

        if (NumNonlocalElements!=0)
        {
            //Nonlocal Model
            constitutivePtr->GetTangent_EngineeringStress_EngineeringStrain(this, theIP,
                deformationGradient, &NonlocalTangent);


            /*
            NuTo::FullMatrix<double> stiffnessAnalytic(Eigen::Matrix<double,3,3>::Map(NonlocalTangent.GetSubMatrix(0)->GetData(),3,3));
std::cout << "Element " << this << std::endl;
std::cout << "stiffness analytic " << std::endl<< stiffnessAnalytic << std::endl << std::endl;

//test the stiffness
double delta(1e-6);
NuTo::EngineeringStress2D  engineeringStress1,engineeringStress2;
constitutivePtr->GetEngineeringStressFromEngineeringStrain(this, theIP,
        deformationGradient, engineeringStress1);
NuTo::FullMatrix<double> stiffnessCDF(3,3);
for (int countStrain=0; countStrain<3; countStrain++)
{
    switch(countStrain)
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
    }
    constitutivePtr->UpdateTmpStaticData_EngineeringStress_EngineeringStrain(const_cast<NuTo::Plane*>(this), theIP, deformationGradient);

    constitutivePtr->GetEngineeringStressFromEngineeringStrain(this, theIP,
                    deformationGradient, engineeringStress2);

    stiffnessCDF(0,countStrain) = (engineeringStress2.mEngineeringStress[0]-engineeringStress1.mEngineeringStress[0])/delta;
    stiffnessCDF(1,countStrain) = (engineeringStress2.mEngineeringStress[1]-engineeringStress1.mEngineeringStress[1])/delta;
    stiffnessCDF(2,countStrain) = (engineeringStress2.mEngineeringStress[2]-engineeringStress1.mEngineeringStress[2])/delta;

    switch(countStrain)
    {
    case 0:
        deformationGradient.mDeformationGradient[0]-=delta;
        break;
    case 1:
        deformationGradient.mDeformationGradient[3]-=delta;
        break;
    case 2:
        deformationGradient.mDeformationGradient[1]-=0.5*delta;
        deformationGradient.mDeformationGradient[2]-=0.5*delta;
    }
    constitutivePtr->UpdateTmpStaticData_EngineeringStress_EngineeringStrain(const_cast<NuTo::Plane*>(this), theIP, deformationGradient);
}
std::cout << "stiffnessCDF " << std::endl << stiffnessCDF << std::endl << std::endl;
if ((stiffnessAnalytic-stiffnessCDF).Abs().Max()/stiffnessCDF.Abs().Max()>0.1)
{
    std::cout << "error " << (stiffnessAnalytic-stiffnessCDF).Abs().Max()/stiffnessCDF.Abs().Max() << std::endl;
    std::cout << "deformation gradient " << std::endl;
    printf("%20.15g %20.15g %20.15g %20.15g\n",deformationGradient.mDeformationGradient[0], deformationGradient.mDeformationGradient[1],
            deformationGradient.mDeformationGradient[2],deformationGradient.mDeformationGradient[3]);
    throw MechanicsException("[NuTo::Plane::CalculateCoefficientMatrix_0] stiffness is wrong (ip stiffness");
}
exit(0);
*/

            if (NonlocalTangent.IsLocal())
            {
                //same as local model, e.g. in the unloading range or
                areAllIpsSymmetric &= NonlocalTangent.GetSubMatrix(0)->GetSymmetry();
                //nonlocal BMatrix of nonlocal other integration point
                int firstCol(0);
                const std::vector<const ElementBase*>& nonlocalElements(GetNonlocalElements());
                for (unsigned int theNonlocalElement = 0; theNonlocalElement<nonlocalElements.size(); theNonlocalElement++)
                {
                    const Plane* nonlocalElement(nonlocalElements[theNonlocalElement]->AsPlane());

                    if (nonlocalElement!=this)
                    {
                        throw MechanicsException("[NuTo::Plane::CalculateCoefficientMatrix_0] The first nonlocal element for an element should always be itself.");
                    }

                    for (int theNonlocalIp = 0; theNonlocalIp < nonlocalElement->GetNumIntegrationPoints(); theNonlocalIp++)
                    {
                        if (theNonlocalIp!= theIP)
                            continue;

                        // calculate element stiffness matrix
                        // don't forget to include determinant of the Jacobian and area
                        AddDetJBtCB(derivativeShapeFunctionsLocal, derivativeShapeFunctionsLocal, NonlocalTangent.GetSubMatrix(0), factor, rCoefficientMatrix, firstCol);
                        break;

                    }
                    firstCol+=nonlocalElement->GetNumLocalDofs();
                    break;
                }
            }
            else
            {
                //nonlocal BMatrix of nonlocal other integration point
                areAllIpsSymmetric = false;
                areAllIpsLocal = false;
                //sum over nonlocal elements and their ips
                const std::vector<const ElementBase*>& nonlocalElements(GetNonlocalElements());
                int firstCol(0);
                int totalNonlocalIp(0);
                for (unsigned int theNonlocalElement = 0; theNonlocalElement<nonlocalElements.size(); theNonlocalElement++)
                {
                    const Plane* nonlocalElement(nonlocalElements[theNonlocalElement]->AsPlane());
                    //calculate local coordinates
                    std::vector<double> nonlocalLocalNodeCoord(nonlocalElement->GetNumLocalDofs());
                    nonlocalElement->CalculateLocalCoordinates(nonlocalLocalNodeCoord);
/*
                    for (int counter=0; counter<(int)nonlocalLocalNodeCoord.size(); counter++)
                    {
                        std::cout << counter << "  " << nonlocalLocalNodeCoord[counter]<< std::endl;
                    }
*/
                    //get weights
                    const std::vector<double>& weights(GetNonlocalWeights(theIP,theNonlocalElement));

                    for (int theNonlocalIp = 0; theNonlocalIp < nonlocalElement->GetNumIntegrationPoints(); theNonlocalIp++, totalNonlocalIp++)
                    {
                        if (weights[theNonlocalIp]==0.)
                            continue;

                        // get IP coordinates
                        nonlocalElement->GetLocalIntegrationPointCoordinates(theNonlocalIp, nonlocalNaturalIPCoord);

                        // calculate derivatives of shape functions in natural coordinate system
                        nonlocalDerivativeShapeFunctionsNatural.resize(2*nonlocalElement->GetNumShapeFunctions());
                        nonlocalElement->CalculateDerivativeShapeFunctionsNatural(nonlocalNaturalIPCoord, nonlocalDerivativeShapeFunctionsNatural);

                        // calculate Jacobian in order to transform natural to local coordinate system
                        nonlocalElement->CalculateJacobian(nonlocalDerivativeShapeFunctionsNatural,nonlocalLocalNodeCoord, nonlocalInvJacobian, nonlocalDetJac);

                        // calculate derivates of shape functions in local coordinate system
                        nonlocalDerivativeShapeFunctionsLocal.resize(nonlocalDerivativeShapeFunctionsNatural.size());
                        nonlocalElement->CalculateDerivativeShapeFunctionsLocal(nonlocalDerivativeShapeFunctionsNatural,nonlocalInvJacobian,
                                                                nonlocalDerivativeShapeFunctionsLocal);

                        // calculate element stiffness matrix
                        // don't forget to include determinant of the Jacobian and area
                        AddDetJBtCB(derivativeShapeFunctionsLocal, nonlocalDerivativeShapeFunctionsLocal, NonlocalTangent.GetSubMatrix(totalNonlocalIp), factor, rCoefficientMatrix, firstCol);
                    }
                    firstCol+=nonlocalElement->GetNumLocalDofs();
                }
                assert(totalNonlocalIp==NumNonlocalIps);
            }
        }
        else
        {
        constitutivePtr->GetTangent_EngineeringStress_EngineeringStrain(this, theIP,
                            deformationGradient, NonlocalTangent.GetSubMatrix(0));
            areAllIpsSymmetric &= NonlocalTangent.GetSubMatrix(0)->GetSymmetry();

            // calculate element stiffness matrix
            // don't forget to include determinant of the Jacobian and area
            AddDetJBtCB(derivativeShapeFunctionsLocal, NonlocalTangent.GetSubMatrix(0), factor, rCoefficientMatrix);
        }
    }
    rSymmetry = areAllIpsSymmetric;

    // calculate list of global dofs related to the entries in the element stiffness matrix
    this->CalculateGlobalRowDofs(rGlobalDofsRow);
    if (areAllIpsLocal)
    {
        // only the first nonlocal matrix has nonzero entries, the rest is zero, so resize the matrix
        // to avoid adding unnecessary elements to the global matrix
        rCoefficientMatrix = rCoefficientMatrix.GetBlock(0,0,rGlobalDofsRow.size(),rGlobalDofsRow.size());
        rGlobalDofsColumn = rGlobalDofsRow;
        assert((int)rGlobalDofsRow.size()==rCoefficientMatrix.GetNumColumns());
    }
    else
    {
        this->CalculateGlobalColumnDofs(rGlobalDofsColumn);
    }

    assert((int)rGlobalDofsRow.size()==rCoefficientMatrix.GetNumRows());
    assert((int)rGlobalDofsColumn.size()==rCoefficientMatrix.GetNumColumns());
}

//! @brief adds to a matrix the product B^tCBnonlocal, where B contains the derivatives of the shape functions and C is the constitutive tangent and Bnonlocal is the nonlocal B matrix
//! eventually include also area/width of an element
//! @param rLocalDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param rNonlocalDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param ConstitutiveTangentBase constitutive tangent matrix
//! @param rFactor factor including determinant of Jacobian and IP weight
//! @param rCoefficientMatrix to be added to
//! &param rFirstCol first column of the coefficient matrix to be modified (corresponding to the current nonlocal element)
void NuTo::Plane::AddDetJBtCB(const std::vector<double>& rLocalDerivativeShapeFunctionsLocal,const std::vector<double>& rNonlocalDerivativeShapeFunctionsLocal,
                              const ConstitutiveTangentLocal3x3* rConstitutiveTangent, double rFactor,
                              FullMatrix<double>& rCoefficientMatrix, int rFirstCol)const
{
    assert(rCoefficientMatrix.GetNumRows()==2*GetNumShapeFunctions() && rFirstCol + (int)rNonlocalDerivativeShapeFunctionsLocal.size()<=rCoefficientMatrix.GetNumColumns());
    assert((int)rLocalDerivativeShapeFunctionsLocal.size()==2*GetNumShapeFunctions());
    const double *C = rConstitutiveTangent->GetData();
    double x1,x2,y1,y2,x1x2,y2x1,x2y1,y2y1;
    for (int theNode1=0; theNode1<GetNumNodes(); theNode1++)
    {
        int node1mul2 = 2*theNode1;
        int node1mul2plus1 = node1mul2+1;

        x1 = rFactor * rLocalDerivativeShapeFunctionsLocal[node1mul2];
        y1 = rFactor * rLocalDerivativeShapeFunctionsLocal[node1mul2plus1];
        // << is division by two, but faster
        for (unsigned int theNode2=0; theNode2<rNonlocalDerivativeShapeFunctionsLocal.size()>>1; theNode2++)
        {
            int node2mul2 = 2*theNode2;
            int node2mul2plus1 = node2mul2+1;

            x2 = rNonlocalDerivativeShapeFunctionsLocal[node2mul2];
            y2 = rNonlocalDerivativeShapeFunctionsLocal[node2mul2plus1];

            x1x2 = x2*x1;
            y2x1 = y2*x1;
            x2y1 = x2*y1;
            y2y1 = y2*y1;

            rCoefficientMatrix(node1mul2,rFirstCol+node2mul2)          +=x1x2*C[0] +x2y1*C[2] + y2x1*C[6] +y2y1*C[8];
            rCoefficientMatrix(node1mul2,rFirstCol+node2mul2plus1)     +=x1x2*C[6] +x2y1*C[8] + y2x1*C[3] +y2y1*C[5];
            rCoefficientMatrix(node1mul2plus1,rFirstCol+node2mul2)     +=x1x2*C[2] +x2y1*C[1] + y2x1*C[8] +y2y1*C[7];
            rCoefficientMatrix(node1mul2plus1,rFirstCol+node2mul2plus1)+=x1x2*C[8] +x2y1*C[7] + y2x1*C[5] +y2y1*C[4];
        }
    }
}

//! @brief adds to a matrix the product B^tCB, where B contains the derivatives of the shape functions and C is the constitutive tangent
//! eventually include also area/width of an element
//! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param ConstitutiveTangentBase constitutive tangent matrix
//! @param rFactor factor including determinant of Jacobian and IP weight
//! @param rCoefficientMatrix to be added to
void NuTo::Plane::AddDetJBtCB(const std::vector<double>& rDerivativeShapeFunctionsLocal,
                              const ConstitutiveTangentLocal3x3* rConstitutiveTangent, double rFactor,
                              FullMatrix<double>& rCoefficientMatrix)const
{
    assert(rCoefficientMatrix.GetNumRows()==2*GetNumShapeFunctions() && rCoefficientMatrix.GetNumColumns()==2*GetNumShapeFunctions());
    assert((int)rDerivativeShapeFunctionsLocal.size()==2*GetNumShapeFunctions());
    const double *C = rConstitutiveTangent->GetData();
    double x1,x2,y1,y2,x1x2,y2x1,x2y1,y2y1;
    for (int theNode1=0; theNode1<GetNumNodes(); theNode1++)
    {
        int node1mul2 = 2*theNode1;
        int node1mul2plus1 = node1mul2+1;

        x1 = rFactor * rDerivativeShapeFunctionsLocal[node1mul2];
        y1 = rFactor * rDerivativeShapeFunctionsLocal[node1mul2plus1];
        for (int theNode2=0; theNode2<GetNumNodes(); theNode2++)
        {
            int node2mul2 = 2*theNode2;
            int node2mul2plus1 = node2mul2+1;

            x2 = rDerivativeShapeFunctionsLocal[node2mul2];
            y2 = rDerivativeShapeFunctionsLocal[node2mul2plus1];

            x1x2 = x2*x1;
            y2x1 = y2*x1;
            x2y1 = x2*y1;
            y2y1 = y2*y1;

            rCoefficientMatrix(node1mul2,node2mul2)          +=x1x2*C[0] +x2y1*C[2] + y2x1*C[6] +y2y1*C[8];
            rCoefficientMatrix(node1mul2,node2mul2plus1)     +=x1x2*C[6] +x2y1*C[8] + y2x1*C[3] +y2y1*C[5];
            rCoefficientMatrix(node1mul2plus1,node2mul2)     +=x1x2*C[2] +x2y1*C[1] + y2x1*C[8] +y2y1*C[7];
            rCoefficientMatrix(node1mul2plus1,node2mul2plus1)+=x1x2*C[8] +x2y1*C[7] + y2x1*C[5] +y2y1*C[4];
        }
    }
}

//! @brief adds up the internal force vector
//! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param rEngineeringStress stress
//! @param factor factor including det Jacobian area and integration point weight
//! @param rResult resforce vector
void NuTo::Plane::AddDetJBtSigma(const std::vector<double>& rDerivativeShapeFunctionsLocal,
                                 const EngineeringStress2D& rEngineeringStress,
                                 double rFactor,
                                 FullMatrix<double>& rResult)const
{
    assert(rResult.GetNumRows()==2*GetNumShapeFunctions() && rResult.GetNumColumns()==1);
    const double *s = rEngineeringStress.GetData();
    double x1,y1;
    for (int theNode1=0; theNode1<GetNumNodes(); theNode1++)
    {
        int node1mul2 = 2*theNode1;
        int node1mul2plus1 = node1mul2+1;

        x1 = rFactor * rDerivativeShapeFunctionsLocal[node1mul2];
        y1 = rFactor * rDerivativeShapeFunctionsLocal[node1mul2plus1];

        rResult(node1mul2,0)     +=x1*s[0]+y1*s[2];
        rResult(node1mul2plus1,0)+=y1*s[1]+x1*s[2];
    }
}

//! @brief Calculates the the inverse of the Jacobian and its determinant
//! @param rDerivativeShapeFunctions Derivatives of the shape functions (dN1dx, dN1dy, dN1dz, dN2dx, ..
//! @param rNodeCoordinates Node coordinates (X1,Y1,Z1,X2,Y2,Z2,...
//! @param rInvJacobian inverse Jacobian matrix (return value)
//! @param rDetJac determinant of the Jacobian (return value)
void NuTo::Plane::CalculateJacobian(const std::vector<double>& rDerivativeShapeFunctions,
                                    const std::vector<double>& rNodeCoordinates,
                                    double rInvJacobian[4],
                                    double& rDetJac)const
{
    /*       jacobian
           j0, j2,
           j1, j3 */

    assert((int)rDerivativeShapeFunctions.size()==2*GetNumNodes() && (int)rNodeCoordinates.size()==2*GetNumNodes());
    double  j0(0.),j1(0.),j2(0.),j3(0.),x,y;

    int theDerivative(0);
    for (int count = 0; count < GetNumNodes(); count++)
    {
        x = rNodeCoordinates[theDerivative];
        y = rNodeCoordinates[theDerivative+1];

        j0 += rDerivativeShapeFunctions[theDerivative] * x;
        j1 += rDerivativeShapeFunctions[theDerivative] * y;
        theDerivative++;

        j2 += rDerivativeShapeFunctions[theDerivative] * x;
        j3 += rDerivativeShapeFunctions[theDerivative] * y;
        theDerivative++;
    }

    rDetJac = (j0*j3-j1*j2);

    if (rDetJac==0)
        throw MechanicsException("[NuTo::Plane::CalculateJacobian] Determinant of the Jacobian is zero, no inversion possible.");

    if (rInvJacobian!=0)
    {
        double invDeterminant(1./rDetJac);
        rInvJacobian[0]=j3*invDeterminant;
        rInvJacobian[1]=-j2*invDeterminant;
        rInvJacobian[2]=-j1*invDeterminant;
        rInvJacobian[3]=j0*invDeterminant;
    }
}

//! @brief calculates the derivative of the shape functions with respect to global coordinates
//! @param std::vector<double>& rDerivativeShapeFunctions derivatives of the shape functions
//! @param rJacInv inverse of the Jacobian
//! @param rDerivativeShapeFunctionsGlobal derivaties of the shape functions with respect to global coordinates
//! size should already be correct, but can be checked with an assert
//! first all the directions for a single node, and then for the next node
void NuTo::Plane::CalculateDerivativeShapeFunctionsLocal(const std::vector<double>& rDerivativeShapeFunctionsNatural, const double rJacInv[4], std::vector<double>& rDerivativeShapeFunctionsLocal)const
{
    assert(rDerivativeShapeFunctionsLocal.size()==rDerivativeShapeFunctionsNatural.size());
    for (int count=0; count<GetNumNodes(); count++)
    {
        int mul2count = 2*count;
        int mul2countplus1 = mul2count+1;
        rDerivativeShapeFunctionsLocal[mul2count] =
            rDerivativeShapeFunctionsNatural[mul2count]     *rJacInv[0]+
            rDerivativeShapeFunctionsNatural[mul2countplus1]*rJacInv[2];

        rDerivativeShapeFunctionsLocal[mul2countplus1] =
            rDerivativeShapeFunctionsNatural[mul2count]     *rJacInv[1]+
            rDerivativeShapeFunctionsNatural[mul2countplus1]*rJacInv[3];
    }
}

//! @brief calculates the gradient of the internal potential
//! for a mechanical problem, this corresponds to the internal force vector
void NuTo::Plane::CalculateGradientInternalPotential(NuTo::FullMatrix<double>& rResult,
        std::vector<int>& rGlobalDofs)const
{
    //calculate local coordinates
    std::vector<double> localNodeCoord(GetNumLocalDofs());
    CalculateLocalCoordinates(localNodeCoord);

    //calculate local displacements
    std::vector<double> localNodeDisp(GetNumLocalDofs());
    CalculateLocalDisplacements(localNodeDisp);

    //allocate space for ip coordinates in natural coordinate system (-1,1)
    double naturalIPCoord[2];

    //allocate space for derivatives of shape functions in natural coordinate system
    std::vector<double> derivativeShapeFunctionsNatural(2*GetNumShapeFunctions());
    //allocate space for derivatives of shape functions in local coordinate system
    std::vector<double> derivativeShapeFunctionsLocal(2*GetNumShapeFunctions());

    //allocate deformation gradient
    DeformationGradient2D deformationGradient;

    //allocate global engineering stress
    EngineeringStress2D engineeringStress;

    //InvJacobian and determinant of Jacobian
    double invJacobian[4], detJac;

    //material pointer
    const ConstitutiveEngineeringStressStrain *constitutivePtr;

    //allocate and initialize result matrix
    rResult.Resize(GetNumLocalDofs(),1);
    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, naturalIPCoord);

        CalculateDerivativeShapeFunctionsNatural(naturalIPCoord, derivativeShapeFunctionsNatural);

        CalculateJacobian(derivativeShapeFunctionsNatural,localNodeCoord, invJacobian, detJac);

        CalculateDerivativeShapeFunctionsLocal(derivativeShapeFunctionsNatural,invJacobian,
                                                derivativeShapeFunctionsLocal);

        // determine deformation gradient from the local Displacements and the derivative of the shape functions
        // this is not included in the AddIpStiffness to avoid reallocation of the deformation gradient for each IP
        CalculateDeformationGradient(derivativeShapeFunctionsLocal, localNodeDisp, deformationGradient);

        //call constitutive law and calculate stress
        constitutivePtr = dynamic_cast<const ConstitutiveEngineeringStressStrain*>(GetConstitutiveLaw(theIP));
        if (constitutivePtr==0)
            throw MechanicsException("[NuTo::Plane::CalculateGradientInternalPotential] Constitutive law can not deal with engineering stresses and strains");
        constitutivePtr->GetEngineeringStressFromEngineeringStrain(this, theIP,
                deformationGradient, engineeringStress);

        //add to local resforce vector
        assert(mSection->GetThickness()>0);
        double factor(mSection->GetThickness()*detJac*(mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP)));
        AddDetJBtSigma(derivativeShapeFunctionsLocal, engineeringStress, factor, rResult);
    }
    // calculate list of global dofs related to the entries in the element stiffness matrix
    this->CalculateGlobalRowDofs(rGlobalDofs);
}

//! @brief Update the static data of an element
void NuTo::Plane::UpdateStaticData(NuTo::Element::eUpdateType rUpdateType)
{
    //calculate local coordinates
    std::vector<double> localNodeCoord(GetNumLocalDofs());
    CalculateLocalCoordinates(localNodeCoord);

    //calculate local displacements
    std::vector<double> localNodeDisp(GetNumLocalDofs());
    CalculateLocalDisplacements(localNodeDisp);

    //allocate space for ip coordinates in natural coordinate system (-1,1)
    double naturalIPCoord[2];

    //allocate space for derivatives of shape functions in natural coordinate system
    std::vector<double> derivativeShapeFunctionsNatural(2*GetNumShapeFunctions());
    //allocate space for derivatives of shape functions in local coordinate system
    std::vector<double> derivativeShapeFunctionsLocal(2*GetNumShapeFunctions());

    //allocate deformation gradient
    DeformationGradient2D deformationGradient;

    //allocate global engineering stress
    EngineeringStress2D engineeringStress;

    //InvJacobian and determinant of Jacobian
    double invJacobian[4], detJac;

    //material pointer
    const ConstitutiveEngineeringStressStrain *constitutivePtr;

    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, naturalIPCoord);

        CalculateDerivativeShapeFunctionsNatural(naturalIPCoord, derivativeShapeFunctionsNatural);

        CalculateJacobian(derivativeShapeFunctionsNatural,localNodeCoord, invJacobian, detJac);

        CalculateDerivativeShapeFunctionsLocal(derivativeShapeFunctionsNatural,invJacobian,
                                                derivativeShapeFunctionsLocal);

        // determine deformation gradient from the local Displacements and the derivative of the shape functions
        // this is not included in the AddIpStiffness to avoid reallocation of the deformation gradient for each IP
        CalculateDeformationGradient(derivativeShapeFunctionsLocal, localNodeDisp, deformationGradient);

        //call constitutive law and calculate stress
        constitutivePtr = dynamic_cast<const ConstitutiveEngineeringStressStrain*>(GetConstitutiveLaw(theIP));
        if (constitutivePtr==0)
            throw MechanicsException("[NuTo::Plane::UpdateStaticData] Constitutive law can not deal with engineering stresses and strains");
        switch(rUpdateType)
        {
        case NuTo::Element::STATICDATA:
            constitutivePtr->UpdateStaticData_EngineeringStress_EngineeringStrain(this, theIP, deformationGradient);
        break;
        case NuTo::Element::TMPSTATICDATA:
            constitutivePtr->UpdateTmpStaticData_EngineeringStress_EngineeringStrain(this, theIP, deformationGradient);
        break;
        default:
            throw MechanicsException("[NuTo::Plane::UpdateStaticData] update static data flag not known (neither static not tmpstatic data");
        }
    }
}

//! @brief calculates the deformation gradient in 2D
//! @param rRerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param rLocalDisp local displacements
//! @param rDeformationGradient (return value)
void NuTo::Plane::CalculateDeformationGradient(const std::vector<double>& rDerivativeShapeFunctionsLocal,
        const std::vector<double>& rLocalDisp,
        DeformationGradient2D& rDeformationGradient)const
{
    assert((int)rLocalDisp.size()==2*GetNumShapeFunctions() && (int)rDerivativeShapeFunctionsLocal.size()==2*GetNumShapeFunctions());

    rDeformationGradient.mDeformationGradient[0] = 1.;
    rDeformationGradient.mDeformationGradient[1] = 0.;
    rDeformationGradient.mDeformationGradient[2] = 0.;
    rDeformationGradient.mDeformationGradient[3] = 1.;

    int theDisp(0);
    double dNdX,dNdY;
    for (int count=0; count<GetNumNodes(); count++)
    {
        dNdX = rDerivativeShapeFunctionsLocal[theDisp];
        dNdY = rDerivativeShapeFunctionsLocal[theDisp+1];

        rDeformationGradient.mDeformationGradient[0]+=rLocalDisp[theDisp]* dNdX;
        rDeformationGradient.mDeformationGradient[1]+=rLocalDisp[theDisp]* dNdY;
        theDisp++;
        rDeformationGradient.mDeformationGradient[2]+=rLocalDisp[theDisp]* dNdX;
        rDeformationGradient.mDeformationGradient[3]+=rLocalDisp[theDisp]* dNdY;
        theDisp++;
    }
}

//! @brief calculates the coefficient matrix for the 1-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the damping matrix
void NuTo::Plane::CalculateCoefficientMatrix_1(NuTo::FullMatrix<double>& rResult,
        std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn)const
{
    throw MechanicsException("[NuTo::Plane::CalculateCoefficientMatrix_1] to be implemented.");
}

//! @brief calculates the coefficient matrix for the 2-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the Mass matrix
void NuTo::Plane::CalculateCoefficientMatrix_2(NuTo::FullMatrix<double>& rResult,
        std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn)const
{
    throw MechanicsException("[NuTo::Plane::CalculateCoefficientMatrix_2] to be implemented.");
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point
//! @param rCoordinates local coordinates (return value)
void NuTo::Plane::GetLocalIntegrationPointCoordinates(int rIpNum, double rCoordinates[2])const
{
    this->mElementData->GetIntegrationType()->GetLocalIntegrationPointCoordinates2D(rIpNum, rCoordinates);
    return;
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point
//! @param rCoordinates local coordinates (return value)
void  NuTo::Plane::GetGlobalIntegrationPointCoordinates(int rIpNum, double rCoordinates[3])const
{
    double naturalCoordinates[2];
    double nodeCoordinates[3];
    this->mElementData->GetIntegrationType()->GetLocalIntegrationPointCoordinates2D(rIpNum, naturalCoordinates);
    std::vector<double> shapeFunctions(GetNumNodes());
    CalculateShapeFunctions(naturalCoordinates, shapeFunctions);
    rCoordinates[0] = 0.;
    rCoordinates[1] = 0.;
    rCoordinates[2] = 0.;

    nodeCoordinates[0] = 0;
    nodeCoordinates[1] = 0;
    nodeCoordinates[2] = 0;
    for (int theNode=0; theNode<GetNumNodes(); theNode++)
    {
        const NodeBase *nodePtr(GetNode(theNode));
        if (nodePtr->GetNumCoordinates()==2)
            nodePtr->GetCoordinates2D(nodeCoordinates);
        else
            nodePtr->GetCoordinates3D(nodeCoordinates);
        for (int theCoordinate=0; theCoordinate<nodePtr->GetNumCoordinates(); theCoordinate++)
        {
            rCoordinates[theCoordinate]+=shapeFunctions[theNode]*nodeCoordinates[theCoordinate];
        }
    }
    return;
}

// these headers are just for test purpose
//#include <eigen2/Eigen/LU>
//#include <eigen2/Eigen/Array>

//! @brief calculates the integration point data with the current displacements applied
//! @param rIpDataType data type to be stored for each integration point
//! @param rIpData return value with dimension (dim of data type) x (numIp)
void NuTo::Plane::GetIpData(NuTo::IpData::eIpStaticDataType rIpDataType, FullMatrix<double>& rIpData)const
{
    //calculate local coordinates
    std::vector<double> localNodeCoord(GetNumLocalDofs());
    CalculateLocalCoordinates(localNodeCoord);

    //std::cout<< "localNodeCoord " << Eigen::Matrix<double,Eigen::Dynamic,1>::Map(&(localNodeCoord[0]),localNodeCoord.size(),1) << std::endl;

    //calculate local displacements
    std::vector<double> localNodeDisp(GetNumLocalDofs());
    CalculateLocalDisplacements(localNodeDisp);

    //allocate space for ip coordinates in natural coordinate system (-1,1)
    double naturalIPCoord[2];

    //allocate space for derivatives of shape functions in natural coordinate system
    std::vector<double> derivativeShapeFunctionsNatural(2*GetNumShapeFunctions());
    //allocate space for derivatives of shape functions in local coordinate system
    std::vector<double> derivativeShapeFunctionsLocal(2*GetNumShapeFunctions());

    //allocate deformation gradient
    DeformationGradient2D deformationGradient;

    //allocate global engineering stress
    EngineeringStrain3D engineeringStrain;

    //allocate global engineering stress
    EngineeringStress3D engineeringStress;

    //InvJacobian and determinant of Jacobian
    double invJacobian[4], detJac;

    //material pointer
    const ConstitutiveEngineeringStressStrain *constitutivePtr;

    //allocate and initialize result matrix
    switch (rIpDataType)
    {
    case NuTo::IpData::ENGINEERING_STRAIN:
    case NuTo::IpData::ENGINEERING_STRESS:
    case NuTo::IpData::ENGINEERING_PLASTIC_STRAIN:
           rIpData.Resize(6,GetNumIntegrationPoints());
    break;
    case NuTo::IpData::DAMAGE:
           rIpData.Resize(1,GetNumIntegrationPoints());
    break;
    case NuTo::IpData::ELASTIC_ENERGY:
    case NuTo::IpData::TOTAL_ENERGY:
        rIpData.Resize(2,GetNumIntegrationPoints());
    break;
    default:
        throw MechanicsException("[NuTo::Plane::GetIpData] Ip data not implemented.");
    }

    //store the data
    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, naturalIPCoord);

        CalculateDerivativeShapeFunctionsNatural(naturalIPCoord, derivativeShapeFunctionsNatural);
        //std::cout<< "derivativeShapeFunctionsNatural " << Eigen::Matrix<double,Eigen::Dynamic,1>::Map(&(derivativeShapeFunctionsNatural[0]),derivativeShapeFunctionsNatural.size(),1) << std::endl;

        CalculateJacobian(derivativeShapeFunctionsNatural,localNodeCoord, invJacobian, detJac);
        //std::cout<< "invJacobian " << Eigen::Matrix<double,2,2>::Map(&(invJacobian[0]),2,2) << std::endl;

        CalculateDerivativeShapeFunctionsLocal(derivativeShapeFunctionsNatural,invJacobian,
                                                derivativeShapeFunctionsLocal);
        //std::cout<< "derivativeShapeFunctionsLocal " << Eigen::Matrix<double,Eigen::Dynamic,1>::Map(&(derivativeShapeFunctionsLocal[0]),derivativeShapeFunctionsLocal.size(),1) << std::endl;

        // determine deformation gradient from the local Displacements and the derivative of the shape functions
        // this is not included in the AddIpStiffness to avoid reallocation of the deformation gradient for each IP
        CalculateDeformationGradient(derivativeShapeFunctionsLocal, localNodeDisp, deformationGradient);

        //call constitutive law and calculate stress
        constitutivePtr = GetConstitutiveLaw(theIP)->AsConstitutiveEngineeringStressStrain();

        switch (rIpDataType)
        {
        case NuTo::IpData::ENGINEERING_STRAIN:
            constitutivePtr->GetEngineeringStrain(this, theIP, deformationGradient, engineeringStrain);
            memcpy(&(rIpData.mEigenMatrix.data()[theIP*6]),engineeringStrain.GetData(),6*sizeof(double));
        break;
        case NuTo::IpData::ENGINEERING_STRESS:
            constitutivePtr->GetEngineeringStressFromEngineeringStrain(this, theIP, deformationGradient, engineeringStress);
            memcpy(&(rIpData.mEigenMatrix.data()[theIP*6]),engineeringStress.GetData(),6*sizeof(double));
        break;
        case NuTo::IpData::ENGINEERING_PLASTIC_STRAIN:
            constitutivePtr->GetEngineeringPlasticStrain(this, theIP, deformationGradient, engineeringStrain);
            memcpy(&(rIpData.mEigenMatrix.data()[theIP*6]),engineeringStrain.GetData(),6*sizeof(double));
        break;
        case NuTo::IpData::DAMAGE:
            constitutivePtr->GetDamage(this, theIP, deformationGradient, rIpData.mEigenMatrix.data()[theIP]);
        break;
        case NuTo::IpData::ELASTIC_ENERGY:
        {
            rIpData.mEigenMatrix(0,theIP) = constitutivePtr->GetElasticEnergy_EngineeringStress_EngineeringStrain(this, theIP, deformationGradient);
            assert(mSection->GetThickness()>0);
            double factor(mSection->GetThickness()*detJac*(mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP)));
            rIpData.mEigenMatrix(1,theIP) = factor;
        }
        break;
        case NuTo::IpData::TOTAL_ENERGY:
        {
            rIpData.mEigenMatrix(0,theIP) = constitutivePtr->GetTotalEnergy_EngineeringStress_EngineeringStrain(this, theIP, deformationGradient);
            assert(mSection->GetThickness()>0);
            double factor(mSection->GetThickness()*detJac*(mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP)));
            rIpData.mEigenMatrix(1,theIP) = factor;
        }
        break;
        default:
            throw MechanicsException("[NuTo::Plane::GetIpData] Ip data not implemented.");
        }
    }
}


//! @brief Allocates static data for an integration point of an element
//! @param rConstitutiveLaw constitutive law, which is called to allocate the static data object
//! actually, both - the element type and the constitutive law are required to determine the static data object actually required
NuTo::ConstitutiveStaticDataBase* NuTo::Plane::AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw)const
{
    const ConstitutiveEngineeringStressStrain *constitutivePtr =
        dynamic_cast<const ConstitutiveEngineeringStressStrain*>(rConstitutiveLaw);
    if (constitutivePtr==0)
        throw MechanicsException("[NuTo::Plane::AllocateStaticData] Constitutive law can not deal with engineering stresses and strains");
    return constitutivePtr->AllocateStaticDataEngineeringStress_EngineeringStrain2D(this);
}

// interpolate geometry
void NuTo::Plane::InterpolateCoordinatesFrom2D(double rNaturalCoordinates[2], double rGlobalCoordinates[3]) const
{
    // calculate shape functions
    std::vector<double> ShapeFunctions(this->GetNumNodes());
    this->CalculateShapeFunctions(rNaturalCoordinates, ShapeFunctions);

    // start interpolation
    rGlobalCoordinates[0] = 0.0;
    rGlobalCoordinates[1] = 0.0;
    rGlobalCoordinates[2] = 0.0;
    for (int NodeCount = 0; NodeCount < this->GetNumNodes(); NodeCount++)
    {
        // get node coordinate
        double NodeCoordinate[3];
        const NodeBase *nodePtr(GetNode(NodeCount));
        if (nodePtr->GetNumCoordinates()==2)
            nodePtr->GetCoordinates2D(NodeCoordinate);
        else
            nodePtr->GetCoordinates3D(NodeCoordinate);

        // add node contribution
        for (int theCoordinate=0; theCoordinate<nodePtr->GetNumCoordinates(); theCoordinate++)
        {
            rGlobalCoordinates[theCoordinate] += ShapeFunctions[NodeCount] *  NodeCoordinate[theCoordinate];
        }
    }
}

// interpolate displacements
void NuTo::Plane::InterpolateDisplacementsFrom2D(double rNaturalCoordinates[2], double rGlobalDisplacements[3]) const
{
    // calculate shape functions
    std::vector<double> ShapeFunctions(this->GetNumNodes());
    this->CalculateShapeFunctions(rNaturalCoordinates, ShapeFunctions);

    // start interpolation
    rGlobalDisplacements[0] = 0.0;
    rGlobalDisplacements[1] = 0.0;
    rGlobalDisplacements[2] = 0.0;
    for (int NodeCount = 0; NodeCount < this->GetNumNodes(); NodeCount++)
    {
        // get node displacements
        double NodeDisplacement[3];
        const NodeBase *nodePtr(GetNode(NodeCount));
        if (nodePtr->GetNumDisplacements()==2 || nodePtr->GetNumFineScaleDisplacements()==2)
            nodePtr->GetDisplacements2D(NodeDisplacement);
        else
            nodePtr->GetDisplacements3D(NodeDisplacement);

        // add node contribution
        for (int theCoordinate=0; theCoordinate<nodePtr->GetNumCoordinates(); theCoordinate++)
        {
            rGlobalDisplacements[theCoordinate] += ShapeFunctions[NodeCount] *  NodeDisplacement[theCoordinate];
        }
    }
}

// check element definition
void NuTo::Plane::CheckElement()
{
    // check nodes
    for (int nodeCount = 0; nodeCount < this->GetNumNodes(); nodeCount++)
    {
        int numCoordinates(GetNode(nodeCount)->GetNumCoordinates());
        if (numCoordinates<2 || numCoordinates>3)
        {
            throw MechanicsException("[NuTo::Plane::CheckElement] invalid node type (check node definition for coordinates).");
        }
    }

    // check node ordering (element length must be positive) and for changing sign in jacobian determinant
    // calculate coordinates
    std::vector<double> nodeCoord(2*this->GetNumNodes());
    this->CalculateLocalCoordinates(nodeCoord);
    /*for (int count=0; count<GetNumNodes(); count++)
    {
        std::cout << "Node " << count+1 << " with coordinates " << nodeCoord[2*count]<<","
                << nodeCoord[2*count+1]<<std::endl;
    }*/

    // check number of integration points
    if (this->GetNumIntegrationPoints() < 1)
    {
        throw MechanicsException("[NuTo::Plane::CheckElement] invalid integration type.");
    }

    // check sign of the jacobian determinant of the first integration point
    double naturalIPCoord[2];
    this->GetLocalIntegrationPointCoordinates(0, naturalIPCoord);

    std::vector<double> derivativeShapeFunctionsNatural(2*GetNumNodes());
    this->CalculateDerivativeShapeFunctionsNatural(naturalIPCoord, derivativeShapeFunctionsNatural);

    double invJacobian[4], detJacobian;
    this->CalculateJacobian(derivativeShapeFunctionsNatural,nodeCoord, invJacobian, detJacobian);
    // reorder nodes if determinant is negative
    if (detJacobian < 0.0)
    {
        this->ReorderNodes();
        // recalculate node coordinates after renumbering
        this->CalculateLocalCoordinates(nodeCoord);
    }

    // check jacobian determinant for all integration points for positive sign and calculate element volume
    double volume = 0;
    for (int ipCount = 0; ipCount < this->GetNumIntegrationPoints(); ipCount++)
    {
        // calculate jacobian determinant
        this->GetLocalIntegrationPointCoordinates(ipCount, naturalIPCoord);
        this->CalculateDerivativeShapeFunctionsNatural(naturalIPCoord, derivativeShapeFunctionsNatural);
        this->CalculateJacobian(derivativeShapeFunctionsNatural,nodeCoord, invJacobian, detJacobian);
        //std::cout << "Jacobian " << detJacobian << std::endl;
        if (detJacobian <= 0)
        {
            throw MechanicsException("[NuTo::Plane::CheckElement] element is not properly defined by this nodes (zero or negative jacobian determinant).");
        }
        volume += this->GetIntegrationPointWeight(ipCount) * detJacobian;
    }

    // check element volume
    if (volume < 1e-14)
    {
        throw MechanicsException("[NuTo::Plane::CheckElement] element with zero volume (check nodes).");
    }
}

//! @brief sets the section of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need a section
//! @param rSection pointer to section
//! @return pointer to constitutive law
void NuTo::Plane::SetSection(const SectionBase* rSection)
{
    mSection = rSection;
}

//! @brief returns a pointer to the section of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need a section
//! @return pointer to section
const NuTo::SectionBase* NuTo::Plane::GetSection()const
{
    return mSection;
}

//! @brief calculates the volume of an integration point (weight * detJac)
//! @param rVolume  vector for storage of the ip volumes (area in 2D)
void NuTo::Plane::GetIntegrationPointVolume(std::vector<double>& rVolume)const
{
   //calculate local coordinates
    std::vector<double> localNodeCoord(GetNumLocalDofs());
    CalculateLocalCoordinates(localNodeCoord);

    //allocate space for ip coordinates in natural coordinate system (-1,1)
    double naturalIPCoord[2];

    //allocate space for derivatives of shape functions in natural coordinate system
    std::vector<double> derivativeShapeFunctionsNatural(2*GetNumShapeFunctions());

    //determinant of Jacobian
    double detJac;

    rVolume.resize(GetNumIntegrationPoints());

    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, naturalIPCoord);

        CalculateDerivativeShapeFunctionsNatural(naturalIPCoord, derivativeShapeFunctionsNatural);

        CalculateJacobian(derivativeShapeFunctionsNatural,localNodeCoord, 0, detJac);

        //attention in 2D, this is just the area, but that is required for the nonlocal model
        rVolume[theIP] = detJac * mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP);
    }
}

//! @brief cast the base pointer to an ElementPlane, otherwise throws an exception
const NuTo::Plane* NuTo::Plane::AsPlane()const
{
    return this;
}

//! @brief cast the base pointer to an Plane, otherwise throws an exception
NuTo::Plane* NuTo::Plane::AsPlane()
{
    return this;
}
//! @brief sets the fine scale model (deserialization from a binary file)
void NuTo::Plane::SetFineScaleModel(int rIp, std::string rFileName)
{
    mElementData->SetFineScaleModel(rIp,rFileName);
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::Plane::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Plane::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Plane::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Plane::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Plane::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Plane::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Plane::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Plane" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementBase)
       & BOOST_SERIALIZATION_NVP(mSection);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Plane" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Plane)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::Plane)
#endif // ENABLE_SERIALIZATION

