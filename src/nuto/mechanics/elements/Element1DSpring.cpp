#include "nuto/mechanics/elements/Element1DSpring.h"
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixInt.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixDouble.h"
#include "nuto/mechanics/elements/ElementOutputVectorInt.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/math/FullMatrix.h"

NuTo::Element1DSpring::Element1DSpring(const NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase*>& rNodes, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType, InterpolationType* rInterpolationType) :
        NuTo::Element1D::Element1D(rStructure, rNodes, rElementDataType, rIpDataType, rInterpolationType)
{

}

NuTo::Element::eElementType NuTo::Element1DSpring::GetEnumType() const
{
    return NuTo::Element::ELEMENT1DSPRING;
}

const Eigen::MatrixXd NuTo::Element1DSpring::ExtractNodeValues(int rTimeDerivative, Node::eDof rDofType) const
{
    const InterpolationBase& interpolationTypeDof = GetInterpolationType()->Get(rDofType);
    int numNodes = interpolationTypeDof.GetNumNodes();
    int numDofsPerNode;

    switch (rDofType)
    {
    case NuTo::Node::COORDINATES:
        numDofsPerNode = GetStructure()->GetDimension();
        break;
    case NuTo::Node::DISPLACEMENTS:
        numDofsPerNode = GetStructure()->GetDimension();
        break;
    default:
        throw NuTo::MechanicsException("[NuTo::Element1DSpring::ExtractNodeValues] dof type not found.");
    }

    Eigen::MatrixXd globalNodeValues(numDofsPerNode, numNodes);

    for (int iNode = 0; iNode < numNodes; ++iNode)
    {
        const NodeBase* node = mNodes[interpolationTypeDof.GetNodeIndex(iNode)];

        switch (rDofType)
        {
        case Node::COORDINATES:
            globalNodeValues.block(0, iNode, numDofsPerNode, 1) = node->Get(Node::COORDINATES);
            break;

        case Node::DISPLACEMENTS:
            globalNodeValues.block(0, iNode, numDofsPerNode, 1) = node->GetDisplacements(rTimeDerivative);
            break;
        default:
            throw MechanicsException("[NuTo::Element1DSpring::ExtractNodeValues] Not implemented for " + Node::DofToString(rDofType) + ".");
        }
    }

    return globalNodeValues;

}

void NuTo::Element1DSpring::CheckElement()
{
    // Jacobian may be zero or negative.
}

NuTo::Error::eError NuTo::Element1DSpring::Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput)
{
    try
    {
        const SectionBase* section(GetSection());
        if (section == nullptr)
            throw MechanicsException("[NuTo::Element1DSpring::Evaluate] no section allocated for element.");

        const std::set<Node::eDof>& dofs = mInterpolationType->GetDofs();
        const std::set<Node::eDof>& activeDofs = mInterpolationType->GetActiveDofs();
        unsigned int numActiveDofs = mInterpolationType->GetNumActiveDofs();
        // extract all node values and store them
        std::map<Node::eDof, Eigen::MatrixXd> nodalValues;
        for (auto dof : dofs)
        {
            nodalValues[dof] = ExtractNodeValues(0, dof);
        }

        /*****************************************\
         *      EVALUATE CONSTITUTIVE LAW        *
         \*****************************************/

        int theIP = 0;
        ConstitutiveBase* constitutivePtr = GetConstitutiveLaw(theIP);

        double springStiffness = constitutivePtr->GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::SPRING_STIFFNESS);
        NuTo::FullVector<double, Eigen::Dynamic> springDirection = constitutivePtr->GetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter::SPRING_DIRECTION);

        /*****************************************\
        *           CALCULATE OUTPUT              *
         \*****************************************/

        for (auto it = rElementOutput.begin(); it != rElementOutput.end(); it++)
        {
            switch (it->first)
            {
            case Element::INTERNAL_GRADIENT:
            {
                it->second->GetFullVectorDouble().Resize(numActiveDofs);
            }
                break;
            case Element::HESSIAN_0_TIME_DERIVATIVE:
            {
                it->second->GetFullMatrixDouble().Resize(numActiveDofs, numActiveDofs);
                it->second->GetFullMatrixDouble().setZero();
                it->second->SetSymmetry(true);
                it->second->SetConstant(true);
                for (auto dof : activeDofs)
                {
                    unsigned int startIndex = mInterpolationType->Get(dof).GetLocalStartIndex();
                    switch (dof)
                    {
                    case Node::DISPLACEMENTS:
                    {
                        AddStiffnessMatrix(springStiffness, springDirection, startIndex, startIndex, it->second->GetFullMatrixDouble());
                    }
                        break;

                    default:
                        throw MechanicsException("[NuTo::Element1DSpring::Evaluate] Element output HESSIAN_0_TIME_DERIVATIVE for " + Node::DofToString(dof) + " not implemented.");

                    }
                }
            }
                break;

            case Element::HESSIAN_1_TIME_DERIVATIVE:

                break;

            case Element::HESSIAN_2_TIME_DERIVATIVE:

                break;
            case Element::LUMPED_HESSIAN_2_TIME_DERIVATIVE:
                break;
            case Element::UPDATE_STATIC_DATA:
            case Element::UPDATE_TMP_STATIC_DATA:
                break;
            case Element::IP_DATA:
                break;
            case Element::GLOBAL_ROW_DOF:
            {
                const Eigen::VectorXi& globalRowDofsEigen = CalculateGlobalRowDofs();
                std::vector<int> globalRowDofsStd(globalRowDofsEigen.data(), globalRowDofsEigen.data() + globalRowDofsEigen.rows());
                it->second->GetVectorInt() = globalRowDofsStd;
            }
                break;
            case Element::GLOBAL_COLUMN_DOF:
            {
                const Eigen::VectorXi& globalColumnDofsEigen = CalculateGlobalColumnDofs();
                std::vector<int> globalColumnDofsStd(globalColumnDofsEigen.data(), globalColumnDofsEigen.data() + globalColumnDofsEigen.rows());
                it->second->GetVectorInt() = globalColumnDofsStd;
            }
                break;

            default:
                throw MechanicsException("[NuTo::Element1DSpring::Evaluate] element output not implemented.");
            }
        }

    } catch (NuTo::MechanicsException& e)
    {
        std::stringstream ss;
        ss << mStructure->ElementGetId(this);
        e.AddMessage("[NuTo::Element1DSpring::Evaluate] Error evaluating element data of element " + ss.str() + ".");
        throw e;
    }
    return Error::SUCCESSFUL;

}

void NuTo::Element1DSpring::AddStiffnessMatrix(double rSpringStiffness, NuTo::FullVector<double, Eigen::Dynamic> rSpringDirection, unsigned int rRow, unsigned int rCol, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCoefficientMatrix) const
{
    unsigned int numNodes = mInterpolationType->Get(Node::DISPLACEMENTS).GetNumNodes();
    unsigned int globalDimension = GetStructure()->GetDimension();

    assert((globalDimension == 2 or globalDimension == 3) and "Need 2d or 3d coordinates");

    // the rotation matrix consists of three basis vectors that define the new coordinate system.
    Eigen::Vector3d basisVector00 = Eigen::Vector3d::Zero();
    Eigen::Vector3d basisVector01 = Eigen::Vector3d::Zero();
    Eigen::Vector3d basisVector02 = Eigen::Vector3d::Zero();

    rSpringDirection.normalize();
    basisVector00.block(0, 0, globalDimension, 1) = rSpringDirection;

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

    basisVector01.normalize();
    basisVector02.normalize();

    // implement roation
    Eigen::Matrix3d rotationMatrix3d = Eigen::Matrix3d::Zero();
    rotationMatrix3d.col(0) = basisVector00;
    rotationMatrix3d.col(1) = basisVector01;
    rotationMatrix3d.col(2) = basisVector02;

    Eigen::MatrixXd tmpStiffnessMatrix(globalDimension * numNodes, globalDimension * numNodes);
    tmpStiffnessMatrix.setZero();

    tmpStiffnessMatrix(0, 0) = 1.;
    tmpStiffnessMatrix(0, globalDimension) = -1.;
    tmpStiffnessMatrix(globalDimension, 0) = -1.;
    tmpStiffnessMatrix(globalDimension, globalDimension) = 1.;

    Eigen::MatrixXd transformationMatrix(globalDimension * numNodes, globalDimension * numNodes);
    transformationMatrix.setZero();

    for (unsigned int i = 0; i < numNodes; ++i)
    {
        transformationMatrix.block(globalDimension * i, globalDimension * i, globalDimension, globalDimension) = rotationMatrix3d.block(0, 0, globalDimension, globalDimension);
    }

    rCoefficientMatrix.block(rRow, rCol, globalDimension * numNodes, globalDimension * numNodes) += rSpringStiffness * transformationMatrix * tmpStiffnessMatrix * transformationMatrix.transpose();

}

const Eigen::VectorXi NuTo::Element1DSpring::CalculateGlobalRowDofs() const
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
                switch (GetStructure()->GetDimension())
                {
                case 1:
                    globalRowDofs[index++] = nodePtr->GetDof(Node::DISPLACEMENTS, 0);
                    break;
                case 2:
                    globalRowDofs[index++] = nodePtr->GetDof(Node::DISPLACEMENTS, 0);
                    globalRowDofs[index++] = nodePtr->GetDof(Node::DISPLACEMENTS, 1);
                    break;

                case 3:
                    globalRowDofs[index++] = nodePtr->GetDof(Node::DISPLACEMENTS, 0);
                    globalRowDofs[index++] = nodePtr->GetDof(Node::DISPLACEMENTS, 1);
                    globalRowDofs[index++] = nodePtr->GetDof(Node::DISPLACEMENTS, 2);
                    break;
                default:
                    throw MechanicsException("[NuTo::Element1DSpring::CalculateGlobalRowDofs] Global dimension must be either 1, 2, or 3");
                    break;
                }

            }
                break;
            default:
                throw MechanicsException("[NuTo::Element1DInXD::CalculateGlobalRowDofs] Not implemented for " + Node::DofToString(dof) + ".");

            }
        }
    }
    return globalRowDofs;
}

