#include "nuto/mechanics/elements/Element1DInXD.h"
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixInt.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixDouble.h"
#include "nuto/mechanics/elements/ElementOutputVectorInt.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/structures/StructureBase.h"

#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveMatrix.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "nuto/mechanics/elements/EvaluateDataContinuum.h"

#include "nuto/math/FullMatrix.h"

#ifdef ENABLE_SERIALIZATION
#include "nuto/math/EigenBoostSerialization.h"
#endif


NuTo::Element1DInXD::Element1DInXD(const NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase*>& rNodes, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType, InterpolationType* rInterpolationType) :
        NuTo::ContinuumElement<1>(rStructure, rNodes, rElementDataType, rIpDataType, rInterpolationType)
{
    mRotationMatrix = CalculateRotationMatrix();
}

NuTo::Element::eElementType NuTo::Element1DInXD::GetEnumType() const
{
    return NuTo::Element::ELEMENT1DINXD;
}


Eigen::MatrixXd NuTo::Element1DInXD::CalculateRotationMatrix()
{
    const Eigen::VectorXd globalNodeCoordinates = ExtractGlobalNodeValues(0, Node::COORDINATES);

    const unsigned int globalDimension = GetStructure()->GetDimension();


    assert((globalDimension == 2 or globalDimension == 3) and "Need 2d or 3d coordinates");

    // the rotation matrix consists of three basis vectors that define the new coordinate system.
    Eigen::Vector3d basisVector00 = Eigen::Vector3d::Zero();
    Eigen::Vector3d basisVector01 = Eigen::Vector3d::Zero();
    Eigen::Vector3d basisVector02 = Eigen::Vector3d::Zero();

    // basisVector00 is aligned with the truss. Note: Trusses have to be straight!
    basisVector00.block(0, 0, globalDimension, 1) = globalNodeCoordinates.tail(globalDimension) - globalNodeCoordinates.head(globalDimension);

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

Eigen::VectorXd NuTo::Element1DInXD::ExtractNodeValues(int rTimeDerivative, Node::eDof rDofType) const
{
    Eigen::VectorXd globalNodeValues = ExtractGlobalNodeValues(rTimeDerivative, rDofType);
    const unsigned int structureDim = GetStructure()->GetDimension();
    const unsigned int numNodes     = GetInterpolationType()->Get(rDofType).GetNumNodes();

    Eigen::VectorXd nodeValues(numNodes);

    for (unsigned iNode = 0; iNode < numNodes; ++iNode)
    {
        Eigen::VectorXd tmp = (mRotationMatrix.transpose() * globalNodeValues.segment(iNode * structureDim, structureDim));
        nodeValues[iNode] = tmp.at(0,0);
    }

    return nodeValues;
}

const Eigen::VectorXd NuTo::Element1DInXD::ExtractGlobalNodeValues(int rTimeDerivative, Node::eDof rDofType) const
{

    const InterpolationBase& interpolationTypeDof = GetInterpolationType()->Get(rDofType);
    int numNodes = interpolationTypeDof.GetNumNodes();
    int numDofsPerNode = GetNumDofsPerNode(rDofType);
    const unsigned int structureDim = GetStructure()->GetDimension();


    Eigen::VectorXd globalNodeValues(numDofsPerNode * numNodes);

    for (int iNode = 0; iNode < numNodes; ++iNode)
    {
        const NodeBase& node = *GetNode(iNode, rDofType);

        switch (rDofType)
        {
        case Node::COORDINATES:
            globalNodeValues.segment(iNode * numDofsPerNode, structureDim) = node.GetCoordinates();
            break;
        case Node::DISPLACEMENTS:
            globalNodeValues.segment(iNode * numDofsPerNode, structureDim) = node.GetDisplacements(rTimeDerivative);
            break;
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for " + Node::DofToString(rDofType));
        }
    }

    return globalNodeValues;

}

void NuTo::Element1DInXD::CalculateElementOutputHessian0(BlockFullMatrix<double>& rHessian0, EvaluateDataContinuum<1>& rData, int rTheIP) const
{


    for (auto dofRow : mInterpolationType->GetActiveDofs())
    {
        for (auto dofCol : mInterpolationType->GetActiveDofs())
        {
            auto& hessian0 = rHessian0(dofRow, dofCol);
            switch (Node::CombineDofs(dofRow, dofCol))
            {
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS):
            {
                const unsigned int globalDimension = GetStructure()->GetDimension();
                const unsigned int numberOfNodes = mInterpolationType->Get(Node::DISPLACEMENTS).GetNumNodes();
                const unsigned int numberOfDofs = globalDimension * numberOfNodes;

                Eigen::MatrixXd transformationMatrix = CalculateTransformationMatrix(globalDimension, numberOfNodes);

                const Eigen::MatrixXd localHessian0 = rData.mDetJxWeightIPxSection * rData.mB.at(dofRow).transpose() * rData.mTangentStressStrain * rData.mB.at(dofRow);
                Eigen::MatrixXd globalHessian0(numberOfDofs, numberOfDofs);
                globalHessian0.setZero();

                for (unsigned iRow = 0; iRow < localHessian0.rows(); ++iRow)
                {
                    for (unsigned iCol = 0; iCol < localHessian0.cols(); ++iCol)
                    {
                        globalHessian0(globalDimension * iRow, globalDimension * iCol) = localHessian0(iRow, iCol);
                    }
                }
                hessian0 += transformationMatrix * globalHessian0 * transformationMatrix.transpose();


            }
                break;
            default:
                throw MechanicsException(__PRETTY_FUNCTION__, "Element output HESSIAN_0_TIME_DERIVATIVE for "
                        "(" + Node::DofToString(dofRow) + "," + Node::DofToString(dofCol) + ") not implemented.");
            }
        }
    }



}

void NuTo::Element1DInXD::CalculateElementOutputInternalGradient(BlockFullVector<double>& rInternalGradient, EvaluateDataContinuum<1>& rData, int rTheIP) const
{



    for (auto dofRow : mInterpolationType->GetActiveDofs())
    {
        switch (dofRow)
        {
        case Node::DISPLACEMENTS:
      {
            const unsigned int globalDimension = GetStructure()->GetDimension();
            const unsigned int numberOfNodes = mInterpolationType->Get(dofRow).GetNumNodes();

            Eigen::MatrixXd transformationMatrix = CalculateTransformationMatrix(globalDimension, numberOfNodes);

            const Eigen::VectorXd localInternalGradient = rData.mDetJxWeightIPxSection * rData.mB.at(dofRow).transpose() * rData.mEngineeringStress;
            Eigen::VectorXd globalInternalGradient(globalDimension * numberOfNodes);
            globalInternalGradient.setZero();

            for (unsigned iNode = 0; iNode < numberOfNodes; ++iNode)
                globalInternalGradient[iNode * globalDimension] = localInternalGradient[iNode];


            rInternalGradient[dofRow] += transformationMatrix * globalInternalGradient;
        }
            break;
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Element output INTERNAL_GRADIENT for " + Node::DofToString(dofRow) + " not implemented.");
        }
    }



}

int NuTo::Element1DInXD::GetNumDofsPerNode(Node::eDof rDofType) const
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

Eigen::VectorXd NuTo::Element1DInXD::InterpolateDofGlobal(const Eigen::VectorXd& rNaturalCoordinates, NuTo::Node::eDof rDofType) const
{
    return InterpolateDofGlobal(0, rNaturalCoordinates, rDofType);

}
Eigen::VectorXd NuTo::Element1DInXD::InterpolateDofGlobal(int rTimeDerivative, const Eigen::VectorXd& rNaturalCoordinates, NuTo::Node::eDof rDofType) const
{

    const InterpolationBase& interpolationType = mInterpolationType->Get(rDofType);
    const Eigen::VectorXd nodalValues = ExtractGlobalNodeValues(rTimeDerivative, rDofType);
    const Eigen::MatrixXd matrixNLocal = interpolationType.CalculateMatrixN(rNaturalCoordinates);

    int numNodes = GetNumNodes();
    int dimBlock = GetNumDofsPerNode(rDofType);



    Eigen::MatrixXd matrixNGlobal(dimBlock, numNodes * dimBlock);
    for (int iNode = 0, iBlock = 0; iNode < numNodes; ++iNode, iBlock += dimBlock)
    {
        matrixNGlobal.block(0, iBlock, dimBlock, dimBlock) = Eigen::MatrixXd::Identity(dimBlock, dimBlock) * matrixNLocal(iNode);
    }


    return matrixNGlobal * nodalValues;

}


void NuTo::Element1DInXD::CheckElement()
{
    unsigned int numIntegrationPoints = GetNumIntegrationPoints();

    // check number of integration points
    assert(numIntegrationPoints > 0);

    Eigen::MatrixXd nodeCoordinates = ExtractNodeValues(0, Node::COORDINATES);

    std::cout << "mRotationMatrix"<< mRotationMatrix << std::endl;
    double length = 0;
    for (unsigned int iIP = 0; iIP < numIntegrationPoints; ++iIP)
    {
        const Eigen::MatrixXd& derivativeShapeFunctions = mInterpolationType->Get(Node::COORDINATES).GetDerivativeShapeFunctionsNatural(iIP);
        Eigen::Matrix<double,1,1> detJacobian = CalculateJacobian(derivativeShapeFunctions, nodeCoordinates);
        assert(detJacobian(0,0) > 0 and "Jacobian needs to be greater than 0");

        length += this->GetIntegrationPointWeight(iIP) * detJacobian(0,0);
    }

    // check element length
    assert(length > 1e-14 and "element with zero length (check nodes)");

}

//void NuTo::Element1DInXD::CalculateGlobalRowDofs(BlockFullVector<int>& rGlobalRowDofs) const
//{
//    for (auto dof : mInterpolationType->GetActiveDofs())
//    {
//        const InterpolationBase& interpolationType = mInterpolationType->Get(dof);
//        int numNodes = interpolationType.GetNumNodes();
//        FullVector<int, Eigen::Dynamic>& dofWiseGlobalRowDofs = rGlobalRowDofs[dof];
//
//        dofWiseGlobalRowDofs.Resize(interpolationType.GetNumDofs());
//        dofWiseGlobalRowDofs.setZero();
//        const unsigned globalDimension = GetStructure()->GetDimension();
//        switch (dof)
//        {
//        case Node::DISPLACEMENTS:
//        {
//            for (int iNodeDof = 0; iNodeDof < numNodes; ++iNodeDof)
//            {
//                const NodeBase* nodePtr = mNodes[interpolationType.GetNodeIndex(iNodeDof)];
//                for (unsigned iDof = 0; iDof < globalDimension; ++iDof)
//                    dofWiseGlobalRowDofs[globalDimension * iNodeDof + iDof] = nodePtr->GetDofDisplacement(iDof);
//            }
//            break;
//        }
//        default:
//            throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for " + Node::DofToString(dof) + ".");
//        }
//    }
//
//
//}

#ifdef ENABLE_SERIALIZATION
template void NuTo::Element1DInXD::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Element1DInXD::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Element1DInXD::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Element1DInXD::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Element1DInXD::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void NuTo::Element1DInXD::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Element1DInXD::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Element1DInXD" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Element1D);
    ar & BOOST_SERIALIZATION_NVP(mRotationMatrix);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Element1DInXD" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Element1DInXD)
#endif
