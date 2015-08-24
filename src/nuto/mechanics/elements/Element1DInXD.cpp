#include "nuto/mechanics/elements/Element1DInXD.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixInt.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixDouble.h"
#include "nuto/mechanics/elements/ElementOutputVectorInt.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/structures/StructureBase.h"

#include "nuto/math/FullMatrix.h"

NuTo::Element1DInXD::Element1DInXD(const NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase*>& rNodes, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType, InterpolationType* rInterpolationType) :
        NuTo::Element1D::Element1D(rStructure, rNodes, rElementDataType, rIpDataType, rInterpolationType)
{
    mRotationMatrix = CalculateRotationMatrix();
}

NuTo::Element::eElementType NuTo::Element1DInXD::GetEnumType() const
{
    return NuTo::Element::ELEMENT1DINXD;
}


Eigen::MatrixXd NuTo::Element1DInXD::CalculateRotationMatrix()
{
    const Eigen::MatrixXd globalNodeCoordinates = ExtractGlobalNodeValues(0, Node::COORDINATES);

    unsigned int globalDimension = globalNodeCoordinates.rows();

    assert((globalDimension == 2 or globalDimension == 3) and "Need 2d or 3d coordinates");

    // the rotation matrix consists of three basis vectors that define the new coordinate system.
    Eigen::Vector3d basisVector00 = Eigen::Vector3d::Zero();
    Eigen::Vector3d basisVector01 = Eigen::Vector3d::Zero();
    Eigen::Vector3d basisVector02 = Eigen::Vector3d::Zero();

    // basisVector00 is aligned with the truss. Note: Trusses have to be straight!
    basisVector00.block(0, 0, globalDimension, 1) = globalNodeCoordinates.col(1) - globalNodeCoordinates.col(0);

    // check if basisVector00 is linear independent to unitVectorZ and calculate the cross product to determine basisVector01
    if (basisVector00.at(0, 0) > 1e-8 or basisVector00.at(1, 0) > 1e-8)
    {
        const Eigen::Vector3d unitVectorZ = Eigen::Vector3d::UnitZ();
        basisVector01 = basisVector00.cross(unitVectorZ);
    } else
    {
        const Eigen::Vector3d unitVectorX = Eigen::Vector3d::UnitX();
        basisVector01 = basisVector00.cross(unitVectorX);
    }

    // basisVector02 is othogonal to basisVector01 and basisVector00
    basisVector02 = basisVector00.cross(basisVector01);

    Eigen::Matrix3d rotationMatrix3d = Eigen::Matrix3d::Zero();

    basisVector00.normalize();
    basisVector01.normalize();
    basisVector02.normalize();

    rotationMatrix3d.col(0) = basisVector00;
    rotationMatrix3d.col(1) = basisVector01;
    rotationMatrix3d.col(2) = basisVector02;

    return rotationMatrix3d.block(0, 0, globalDimension, globalDimension);
}

Eigen::MatrixXd NuTo::Element1DInXD::CalculateTransformationMatrix(unsigned int rGlobalDimension, unsigned int rNumberOfNodes) const
{

    Eigen::MatrixXd transformationMatrix(rGlobalDimension * rNumberOfNodes, rGlobalDimension * rNumberOfNodes);
    transformationMatrix.setZero();

    for (unsigned int i = 0; i < rNumberOfNodes; ++i)
    {
        transformationMatrix.block(rGlobalDimension * i, rGlobalDimension * i, rGlobalDimension, rGlobalDimension) = mRotationMatrix;
    }

    return transformationMatrix;

}

const Eigen::MatrixXd NuTo::Element1DInXD::ExtractNodeValues(int rTimeDerivative, Node::eAttributes rDofType) const
{
    Eigen::MatrixXd globalNodeCoordinates = ExtractGlobalNodeValues(rTimeDerivative, rDofType);
    Eigen::MatrixXd nodeCoordinates;

    return (mRotationMatrix.transpose() * globalNodeCoordinates).row(0);
}

const Eigen::MatrixXd NuTo::Element1DInXD::ExtractGlobalNodeValues(int rTimeDerivative, Node::eAttributes rDofType) const
{
    const InterpolationBase& interpolationTypeDof = GetInterpolationType()->Get(rDofType);
    int numNodes = interpolationTypeDof.GetNumNodes();
    int numDofsPerNode = GetNumDofsPerNode(rDofType);

    Eigen::MatrixXd globalNodeValues(numDofsPerNode, numNodes);

    for (int iNode = 0; iNode < numNodes; ++iNode)
    {
        const NodeBase* node = mNodes[interpolationTypeDof.GetNodeIndex(iNode)];

        switch (rDofType)
        {
        case Node::COORDINATES:
            globalNodeValues.block(0, iNode, numDofsPerNode, 1) = node->GetCoordinates();
            break;

        case Node::DISPLACEMENTS:
            globalNodeValues.block(0, iNode, numDofsPerNode, 1) = node->GetDisplacements(rTimeDerivative);
            break;
        default:
            throw MechanicsException("[NuTo::Element1DInXD::ExtractNodeValues] Not implemented for " + Node::AttributeToString(rDofType) + ".");
        }
    }

    return globalNodeValues;

}

int NuTo::Element1DInXD::GetNumDofsPerNode(Node::eAttributes rDofType) const
{
    switch (rDofType)
    {
    case NuTo::Node::COORDINATES:
        return GetStructure()->GetDimension();
    case NuTo::Node::DISPLACEMENTS:
        return GetStructure()->GetDimension();
    default:
        throw NuTo::MechanicsException("[NuTo::Element1DInXD::GetNumDofsPerNode] dof type not found.");
    }
}

const Eigen::VectorXd NuTo::Element1DInXD::InterpolateDofGlobal(const Eigen::VectorXd& rNaturalCoordinates, NuTo::Node::eAttributes rDofType) const
{
    return InterpolateDofGlobal(0, rNaturalCoordinates, rDofType);

}
const Eigen::VectorXd NuTo::Element1DInXD::InterpolateDofGlobal(int rTimeDerivative, const Eigen::VectorXd& rNaturalCoordinates, NuTo::Node::eAttributes rDofType) const
{

    const InterpolationBase& interpolationType = mInterpolationType->Get(rDofType);
    const Eigen::MatrixXd nodalValues = ExtractGlobalNodeValues(rTimeDerivative, rDofType);
    const Eigen::VectorXd shapeFunctions = interpolationType.CalculateShapeFunctions(rNaturalCoordinates);

    return nodalValues * shapeFunctions;

}

const Eigen::VectorXi NuTo::Element1DInXD::CalculateGlobalRowDofs() const
{

    int numActiveDos = mInterpolationType->GetNumActiveDofs();

    Eigen::VectorXi globalRowDofs(numActiveDos);

    for (auto dof : mInterpolationType->GetActiveDofs())
    {
        const InterpolationBase& interpolationType = mInterpolationType->Get(dof);
        int index = interpolationType.GetLocalStartIndex();

        for (int iNodeDof = 0; iNodeDof < interpolationType.GetNumNodes(); ++iNodeDof)
        {
            const NodeBase* nodePtr = mNodes[interpolationType.GetNodeIndex(iNodeDof)];
            switch (dof)
            {
            case Node::DISPLACEMENTS:
            {
                if (GetStructure()->GetDimension() == 2)
                {
                    globalRowDofs[index++] = nodePtr->GetDofDisplacement(0);
                    globalRowDofs[index++] = nodePtr->GetDofDisplacement(1);
                } else if (GetStructure()->GetDimension() == 3)
                {
                    globalRowDofs[index++] = nodePtr->GetDofDisplacement(0);
                    globalRowDofs[index++] = nodePtr->GetDofDisplacement(1);
                    globalRowDofs[index++] = nodePtr->GetDofDisplacement(2);
                } else
                {
                    throw MechanicsException("[NuTo::Element1DInXD::CalculateGlobalRowDofs] Global dimension must be either 2 or 3");
                }
            }
                break;
            default:
                throw MechanicsException("[NuTo::Element1DInXD::CalculateGlobalRowDofs] Not implemented for " + Node::AttributeToString(dof) + ".");

            }
        }
    }
    return globalRowDofs;
}

void NuTo::Element1DInXD::AddDetJBtCB(const Eigen::MatrixXd& rDerivativeShapeFunctions, const ConstitutiveTangentLocal<1, 1>& rConstitutiveTangent, double rFactor, int rRow, int rCol, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCoefficientMatrix) const
{
    const unsigned int globalDimension = GetStructure()->GetDimension();
    const unsigned int numberOfNodes = mInterpolationType->Get(Node::DISPLACEMENTS).GetNumNodes();
    assert(rDerivativeShapeFunctions.rows() == numberOfNodes);
    assert(rDerivativeShapeFunctions.cols() == 1);

    Eigen::MatrixXd globalDerivativeShapeFunctions(globalDimension * numberOfNodes, 1);
    globalDerivativeShapeFunctions.setZero();

    for (unsigned int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        globalDerivativeShapeFunctions.row(globalDimension * iNode) = rDerivativeShapeFunctions.row(iNode);
    }

    Eigen::MatrixXd transformationMatrix = CalculateTransformationMatrix(globalDimension, numberOfNodes);

    Eigen::MatrixXd result = rFactor * rConstitutiveTangent.data()[0] * globalDerivativeShapeFunctions * globalDerivativeShapeFunctions.transpose();
    rCoefficientMatrix.block(rRow, rCol, globalDimension * numberOfNodes, globalDimension * numberOfNodes) += transformationMatrix * result * transformationMatrix.transpose();

}

void NuTo::Element1DInXD::AddDetJBtSigma(const Eigen::MatrixXd& rDerivativeShapeFunctions, const EngineeringStress1D& rEngineeringStress, double rFactor, int rRow, FullVector<double, Eigen::Dynamic>& rResult) const
{

    assert(rResult.GetNumColumns() == 1);
    assert(rDerivativeShapeFunctions.cols() == 1);
    unsigned int globalDimension = GetStructure()->GetDimension();

    const unsigned int numberOfNodes = mInterpolationType->Get(Node::DISPLACEMENTS).GetNumNodes();
    assert(rDerivativeShapeFunctions.rows() == numberOfNodes);

    Eigen::MatrixXd globalDerivativeShapeFunctions(globalDimension * numberOfNodes, 1);
    globalDerivativeShapeFunctions.setZero();

    for (unsigned int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        globalDerivativeShapeFunctions.row(globalDimension * iNode) = rDerivativeShapeFunctions.row(iNode);
    }

    Eigen::MatrixXd transformationMatrix = CalculateTransformationMatrix(globalDimension, numberOfNodes);

    Eigen::VectorXd result = rFactor * rEngineeringStress.GetData()[0] * globalDerivativeShapeFunctions;

    rResult.block(rRow, 0, globalDimension * numberOfNodes, 1) += transformationMatrix * result;

}

void NuTo::Element1DInXD::CheckElement()
{

    unsigned int numIntegrationPoints = GetNumIntegrationPoints();

    // check number of integration points
    assert(numIntegrationPoints > 0);

    Eigen::MatrixXd nodeCoordinates = ExtractNodeValues(0, Node::COORDINATES);

    double length = 0;
    for (unsigned int iIP = 0; iIP < numIntegrationPoints; ++iIP)
    {
        const Eigen::MatrixXd& derivativeShapeFunctions = mInterpolationType->Get(Node::COORDINATES).GetDerivativeShapeFunctionsNatural(iIP);
        double detJacobian = CalculateJacobian(derivativeShapeFunctions, nodeCoordinates);
        assert(detJacobian > 0);

        length += this->GetIntegrationPointWeight(iIP) * detJacobian;
    }

    // check element length
    assert(length > 1e-14 and "element with zero length (check nodes)");
}

