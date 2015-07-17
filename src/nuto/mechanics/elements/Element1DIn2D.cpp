/*
 * Element1DIn2D.cpp
 *
 *  Created on: 13 July 2015
 *      Author: phuschke
 */
#include "nuto/mechanics/elements/Element1D.h"
#include "nuto/mechanics/elements/Element1DIn2D.h"
#include "nuto/mechanics/constitutive/mechanics/Damage.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/thermal/HeatFlux1D.h"
#include "nuto/mechanics/constitutive/thermal/HeatFlux3D.h"
#include "nuto/mechanics/constitutive/thermal/Temperature.h"
#include "nuto/mechanics/constitutive/thermal/TemperatureGradient1D.h"
#include "nuto/mechanics/constitutive/thermal/TemperatureGradient3D.h"
#include "nuto/mechanics/constitutive/mechanics/LocalEqStrain.h"
#include "nuto/mechanics/constitutive/mechanics/NonlocalEqStrain.h"
#include "nuto/mechanics/constitutive/moistureTransport/RelativeHumidity.h"
#include "nuto/mechanics/constitutive/moistureTransport/WaterVolumeFraction.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentNonlocal.h"
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixInt.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixDouble.h"
#include "nuto/mechanics/elements/ElementOutputVectorInt.h"
#include "nuto/mechanics/elements/Plane.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionTruss.h"
#include "nuto/mechanics/structures/StructureBase.h"

#include "nuto/math/FullMatrix.h"

NuTo::Element1DIn2D::Element1DIn2D(const NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase*>& rNodes, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType, InterpolationType* rInterpolationType) :
        NuTo::Element1D::Element1D(rStructure, rNodes, rElementDataType, rIpDataType, rInterpolationType)
{
    mElementLength = CalculateElementLength();
    mRotationMatrix = CalculateRotationMatrix();
}

NuTo::Element::eElementType NuTo::Element1DIn2D::GetEnumType() const
{
    return NuTo::Element::ELEMENT1DIn2D;
}

int NuTo::Element1DIn2D::GetGlobalDimension() const
{
    return 2;
}

int NuTo::Element1DIn2D::GetLocalDimension() const
{
    return 1;
}


Eigen::Matrix2d NuTo::Element1DIn2D::CalculateRotationMatrix()
{
    Eigen::MatrixXd globalNodeCoordinates = ExtractGlobalNodeValues(0, Node::COORDINATES);

    unsigned int numberOfNodes = globalNodeCoordinates.cols();
    assert(globalNodeCoordinates.rows() == 2 && "Need 2D coordinates");

    double x1 = globalNodeCoordinates.at(0, 0);
    double y1 = globalNodeCoordinates.at(1, 0);
    double x2 = globalNodeCoordinates.at(0, numberOfNodes - 1);
    double y2 = globalNodeCoordinates.at(1, numberOfNodes - 1);

    double cos = (x2 - x1) / mElementLength;
    double sin = (y2 - y1) / mElementLength;

    Eigen::Matrix2d rotationMatrix;

    rotationMatrix(0, 0) = cos;
    rotationMatrix(1, 0) = -sin;
    rotationMatrix(0, 1) = sin;
    rotationMatrix(1, 1) = cos;

    return rotationMatrix;
}

const Eigen::MatrixXd NuTo::Element1DIn2D::ExtractNodeValues(int rTimeDerivative, Node::eAttributes rDofType) const
{

    Eigen::MatrixXd globalNodeCoordinates = ExtractGlobalNodeValues(rTimeDerivative, rDofType);
    Eigen::MatrixXd nodeCoordinates;

    nodeCoordinates = mRotationMatrix * globalNodeCoordinates;

    return nodeCoordinates.row(0);
}

const Eigen::MatrixXd NuTo::Element1DIn2D::ExtractGlobalNodeValues(int rTimeDerivative, Node::eAttributes rDofType) const
{
    const InterpolationBase& interpolationTypeDof = GetInterpolationType()->Get(rDofType);

    int numNodes = interpolationTypeDof.GetNumNodes();
    int numDofsPerNode = interpolationTypeDof.GetNumDofsPerNode();

    Eigen::MatrixXd globalNodeValues(numDofsPerNode, numNodes);

    for (int iNode = 0; iNode < numNodes; ++iNode)
    {
        const NodeBase* node = mNodes[interpolationTypeDof.GetNodeIndex(iNode)];

        switch (rDofType)
        {
        case Node::COORDINATES:
            globalNodeValues.block(0, iNode, numDofsPerNode, 1) = node->GetCoordinates2D();
            break;

        case Node::DISPLACEMENTS:
            globalNodeValues.block(0, iNode, numDofsPerNode, 1) = node->GetDisplacements2D(rTimeDerivative);
            break;
        default:
            throw MechanicsException("[NuTo::Element1DIn2D::ExtractNodeValues] Not implemented for " + Node::AttributeToString(rDofType) + ".");
        }
    }
    return globalNodeValues;

}

const Eigen::VectorXd NuTo::Element1DIn2D::InterpolateDofGlobal(const Eigen::VectorXd& rNaturalCoordinates, NuTo::Node::eAttributes rDofType) const
{
    return InterpolateDofGlobal(0, rNaturalCoordinates, rDofType);

}
const Eigen::VectorXd NuTo::Element1DIn2D::InterpolateDofGlobal(int rTimeDerivative, const Eigen::VectorXd& rNaturalCoordinates, NuTo::Node::eAttributes rDofType) const
{

    const InterpolationBase& interpolationType = mInterpolationType->Get(rDofType);
    const Eigen::MatrixXd nodalValues = ExtractGlobalNodeValues(rTimeDerivative, rDofType);
    const Eigen::VectorXd shapeFunctions = interpolationType.CalculateShapeFunctions(rNaturalCoordinates);

    return nodalValues * shapeFunctions;

}


double NuTo::Element1DIn2D::CalculateElementLength()
{

    Eigen::MatrixXd globalNodeCoordinates = ExtractGlobalNodeValues(0, NuTo::Node::eAttributes::COORDINATES);

    assert(globalNodeCoordinates.rows() == 2 && "Need 2D coordinates");

    unsigned int numberOfNodes = globalNodeCoordinates.cols();

    double x1 = globalNodeCoordinates.at(0, 0);
    double y1 = globalNodeCoordinates.at(1, 0);
    double x2 = globalNodeCoordinates.at(0, numberOfNodes - 1);
    double y2 = globalNodeCoordinates.at(1, numberOfNodes - 1);

    return std::sqrt(std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2));
}

double NuTo::Element1DIn2D::GetElementLength() const
{
    return mElementLength;
}

void NuTo::Element1DIn2D::SetElementLength(const double rElementLength)
{
    mElementLength = rElementLength;
}

const Eigen::VectorXi NuTo::Element1DIn2D::CalculateGlobalRowDofs() const
{
    int numActiveDos = mInterpolationType->GetNumActiveDofs();
    Eigen::VectorXi globalRowDofs(numActiveDos);

    for (auto dof : mInterpolationType->GetActiveDofs())
    {
        const InterpolationBase& interpolationType = mInterpolationType->Get(dof);
        int index = interpolationType.GetLocalStartIndex(); //

        for (int iNodeDof = 0; iNodeDof < interpolationType.GetNumNodes(); ++iNodeDof)
        {
            const NodeBase* nodePtr = mNodes[interpolationType.GetNodeIndex(iNodeDof)];
            switch (dof)
            {
            case Node::DISPLACEMENTS:
            {
                globalRowDofs[index++] = nodePtr->GetDofDisplacement(0);
                globalRowDofs[index++] = nodePtr->GetDofDisplacement(1);
            }
                break;
            default:
                throw MechanicsException("[NuTo::Element1DIn2D::CalculateGlobalRowDofs] Not implemented for " + Node::AttributeToString(dof) + ".");

            }
        }
    }
    return globalRowDofs;
}


void NuTo::Element1DIn2D::AddDetJBtCB(const Eigen::MatrixXd& rDerivativeShapeFunctions, const ConstitutiveTangentLocal<1, 1>& rConstitutiveTangent, double rFactor, int rRow, int rCol, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCoefficientMatrix) const
{
    int numNodes = mInterpolationType->Get(Node::DISPLACEMENTS).GetNumNodes();
    assert(rDerivativeShapeFunctions.rows() == numNodes);
    assert(rDerivativeShapeFunctions.cols() == 1);

    Eigen::MatrixXd globalDerivativeShapeFunctions(2 * numNodes, 1);
    for (int iNode = 0; iNode < numNodes; ++iNode)
    {
        globalDerivativeShapeFunctions(2 * iNode + 0, 0) = rDerivativeShapeFunctions.at(iNode, 0);
        globalDerivativeShapeFunctions(2 * iNode + 1, 0) = 0.0;
    }
    Eigen::MatrixXd result = rFactor * rConstitutiveTangent.data()[0] * globalDerivativeShapeFunctions * globalDerivativeShapeFunctions.transpose();

    Eigen::MatrixXd transformationMatrix(2 * numNodes, 2 * numNodes);
    transformationMatrix.setZero();

    for (int i = 0; i < numNodes; ++i)
    {
        transformationMatrix.block(2*i, 2*i, 2, 2) = mRotationMatrix;
    }

    rCoefficientMatrix.block(rRow, rCol, 2 * numNodes, 2 * numNodes) += transformationMatrix.transpose() * result * transformationMatrix;


}

void NuTo::Element1DIn2D::AddDetJBtSigma(const Eigen::MatrixXd& rDerivativeShapeFunctions, const EngineeringStress1D& rEngineeringStress, double rFactor, int rRow, FullVector<double, Eigen::Dynamic>& rResult) const
{
//    assert(rResult.GetNumRows()==2*GetNumNodesField() && rResult.GetNumColumns()==1);
    assert(rResult.GetNumColumns() == 1);
    int numNodes = mInterpolationType->Get(Node::DISPLACEMENTS).GetNumNodes();
    assert(rDerivativeShapeFunctions.rows() == numNodes);
    assert(rDerivativeShapeFunctions.cols() == 1);

    Eigen::MatrixXd globalDerivativeShapeFunctions(2 * numNodes, 1);
    for (int iNode = 0; iNode < numNodes; ++iNode)
    {
        globalDerivativeShapeFunctions(2 * iNode + 0, 0) = rDerivativeShapeFunctions.at(iNode, 0);
        globalDerivativeShapeFunctions(2 * iNode + 1, 0) = 0.0;
    }

    Eigen::MatrixXd transformationMatrix(2 * numNodes, 2 * numNodes);
    transformationMatrix.setZero();

    for (int i = 0; i < numNodes; ++i)
    {
        transformationMatrix.block(2*i, 2*i, 2, 2) = mRotationMatrix;
    }

    Eigen::VectorXd result = rFactor * rEngineeringStress.GetData()[0] * globalDerivativeShapeFunctions;

    rResult.block(rRow, 0, 2 * numNodes, 1) += transformationMatrix.transpose() * result;

}

void NuTo::Element1DIn2D::CheckElement()
{

    int numIntegrationPoints = GetNumIntegrationPoints();

    // check number of integration points
    assert(numIntegrationPoints > 0);

    Eigen::MatrixXd nodeCoordinates = ExtractNodeValues(0, Node::COORDINATES);

    double length = 0;
    for (int iIP = 0; iIP < numIntegrationPoints; ++iIP)
    {
        const Eigen::MatrixXd& derivativeShapeFunctions = mInterpolationType->Get(Node::COORDINATES).GetDerivativeShapeFunctionsNatural(iIP);
        double detJacobian = CalculateJacobian(derivativeShapeFunctions, nodeCoordinates);
        if (detJacobian <= 0)
            throw MechanicsException("[NuTo::Element1DIn2D::CheckElement] Determinant of the Jacobian <= zero, no inversion possible.");
        length += this->GetIntegrationPointWeight(iIP) * detJacobian;
    }

    // check element length
    if (length < 1e-14)
    {
        throw MechanicsException("[NuTo::Element1DIn2D::CheckElement] element with zero length (check nodes).");
    }

}

