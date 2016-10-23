
#include "nuto/base/ErrorEnum.h"

#include "nuto/math/FullMatrix.h"

#include "nuto/mechanics/elements/Element2DInterface.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/nodes/NodeEnum.h"

#include "nuto/mechanics/sections/SectionFibreMatrixBond.h"
#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/elements/ElementOutputBase.h"
#include "nuto/mechanics/elements/ElementOutputIpData.h"

#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/interpolationtypes/InterpolationBase.h"
#include "nuto/mechanics/interpolationtypes/InterpolationType.h"

#include "nuto/mechanics/elements/EvaluateDataContinuum.h"
#include "nuto/mechanics/elements/IpDataEnum.h"

#include "nuto/mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockFullMatrix.h"

#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveMatrixXd.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"

#include "nuto/mechanics/structures/StructureBase.h"
#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeComponent.h"
#include "nuto/visualize/VisualizeEnum.h"
#include "nuto/visualize/VisualizeUnstructuredGrid.h"
#endif // ENABLE_VISUALIZE


NuTo::Element2DInterface::Element2DInterface(const NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase*>& rNodes, const InterpolationType& rInterpolationType) :
        NuTo::ElementBase::ElementBase(rStructure, rInterpolationType),
        mNodes(rNodes)
{
    mTransformationMatrix = CalculateTransformationMatrix(GetStructure()->GetDimension(), mNodes.size());
}

NuTo::ConstitutiveOutputMap NuTo::Element2DInterface::GetConstitutiveOutputMap(std::map<Element::eOutput, std::shared_ptr<ElementOutputBase> >& rElementOutput)
{

//    const unsigned globalDimension = GetStructure()->GetDimension(); // --- unused so far
    ConstitutiveOutputMap constitutiveOutput;

    for (auto it : rElementOutput)
    {
        switch (it.first)
        {
        case Element::eOutput::INTERNAL_GRADIENT:
        {
            FillConstitutiveOutputMapInternalGradient(constitutiveOutput, it.second->GetBlockFullVectorDouble());
        }
            break;
        case Element::eOutput::HESSIAN_0_TIME_DERIVATIVE:
        {
            FillConstitutiveOutputMapHessian0(constitutiveOutput, it.second->GetBlockFullMatrixDouble());

        }
            break;
        case Element::eOutput::HESSIAN_1_TIME_DERIVATIVE:
        case Element::eOutput::HESSIAN_2_TIME_DERIVATIVE:
        case Element::eOutput::LUMPED_HESSIAN_2_TIME_DERIVATIVE:
            break;
        case Element::eOutput::UPDATE_STATIC_DATA:
            constitutiveOutput[NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA] = 0;
            break;
        case Element::eOutput::UPDATE_TMP_STATIC_DATA:
            constitutiveOutput[NuTo::Constitutive::eOutput::UPDATE_TMP_STATIC_DATA] = 0;
            break;
        case Element::eOutput::IP_DATA:
        {
            FillConstitutiveOutputMapIpData(constitutiveOutput, it.second->GetIpData());
        }
            break;
        case Element::eOutput::GLOBAL_ROW_DOF:
            CalculateGlobalRowDofs(it.second->GetBlockFullVectorInt());
            break;
        case Element::eOutput::GLOBAL_COLUMN_DOF:
            CalculateGlobalRowDofs(it.second->GetBlockFullVectorInt());
            break;
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "element output not implemented.");
        }
    }

    return constitutiveOutput;
}

NuTo::ConstitutiveInputMap NuTo::Element2DInterface::GetConstitutiveInputMap(const ConstitutiveOutputMap& rConstitutiveOutput) const
{
    ConstitutiveInputMap constitutiveInputMap = GetConstitutiveLaw(0).GetConstitutiveInputs(rConstitutiveOutput, GetInterpolationType());
    for (auto& itInput : constitutiveInputMap)
    {
        itInput.second = ConstitutiveIOBase::makeConstitutiveIO<2>(itInput.first);
    }
    return constitutiveInputMap;
}

NuTo::eError NuTo::Element2DInterface::Evaluate(const ConstitutiveInputMap& rInput, std::map<Element::eOutput, std::shared_ptr<ElementOutputBase> >& rElementOutput)
{

    if (mSection == nullptr)
        throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t no section allocated for element.");

    const unsigned globalDimension = GetStructure()->GetDimension();
    const std::set<Node::eDof>& dofs = mInterpolationType->GetDofs();

    // extract all node values and store them
    EvaluateData data;
    for (auto dof : dofs)
    {
        data.mNodalValues[dof] = ExtractNodeValues(0, dof);
    }


    auto constitutiveOutput = GetConstitutiveOutputMap(rElementOutput);
    auto constitutiveInput = GetConstitutiveInputMap(constitutiveOutput);
    constitutiveInput.Merge(rInput);

    for (int theIP = 0; theIP < GetNumIntegrationPoints(); ++theIP)
    {

        const Eigen::VectorXd nodeCoords = data.mNodalValues[Node::eDof::COORDINATES];

        // 0.5 * element length
        data.mDetJacobian = 0.5 * (nodeCoords.segment(0, globalDimension) - nodeCoords.segment(0.5 * nodeCoords.rows(), globalDimension)).norm();

        for (auto dof : mInterpolationType->GetDofs())
        {
            if (dof == Node::eDof::COORDINATES)
                continue;

            const InterpolationBase& interpolationType = mInterpolationType->Get(dof);
            data.mMatrixN[dof] = &interpolationType.GetMatrixN(theIP);
            const Eigen::MatrixXd shapeFunctions = *data.mMatrixN[dof];

            const unsigned numberOfNodes = interpolationType.GetNumNodes();
            Eigen::MatrixXd BMatrix(globalDimension, globalDimension * numberOfNodes);
            BMatrix.setZero();

            for (unsigned int iNode = 0; iNode < numberOfNodes; ++iNode)
            {
                const unsigned int index = globalDimension * iNode;
                for (unsigned int iDim = 0; iDim < globalDimension; ++iDim)
                {
                    BMatrix(iDim, index + iDim) = shapeFunctions(0, iNode);
                }

            }

            data.mMatrixB[dof] = BMatrix;

        }

        CalculateConstitutiveInputs(constitutiveInput, data);

        try
        {
            eError error = EvaluateConstitutiveLaw<2>(constitutiveInput, constitutiveOutput, theIP);
            if (error != eError::SUCCESSFUL)
                return error;
        } catch (NuTo::MechanicsException& e)
        {
            e.AddMessage(__PRETTY_FUNCTION__, "error evaluating the constitutive model.");
            throw e;
        }
        CalculateElementOutputs(rElementOutput, data, theIP, constitutiveOutput);
    }
    return eError::SUCCESSFUL;

}

NuTo::Element::eElementType NuTo::Element2DInterface::GetEnumType() const
{
    return Element::eElementType::ELEMENT2DINTERFACE;
}


NuTo::NodeBase* NuTo::Element2DInterface::GetNode(int rLocalNodeNumber)
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->GetNumNodes());
    return mNodes[rLocalNodeNumber];
}

const NuTo::NodeBase* NuTo::Element2DInterface::GetNode(int rLocalNodeNumber) const
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->GetNumNodes());
    return mNodes[rLocalNodeNumber];
}

NuTo::NodeBase* NuTo::Element2DInterface::GetNode(int rLocalNodeNumber, Node::eDof rDofType)
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->Get(rDofType).GetNumNodes());
    int nodeIndex = mInterpolationType->Get(rDofType).GetNodeIndex(rLocalNodeNumber);
    assert(nodeIndex < (int )mNodes.size());
    return mNodes[nodeIndex];
}

const NuTo::NodeBase* NuTo::Element2DInterface::GetNode(int rLocalNodeNumber, Node::eDof rDofType) const
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->Get(rDofType).GetNumNodes());
    int nodeIndex = mInterpolationType->Get(rDofType).GetNodeIndex(rLocalNodeNumber);
    assert(nodeIndex < (int )mNodes.size());
    assert(mNodes[nodeIndex] != nullptr);
    return mNodes[nodeIndex];
}

void NuTo::Element2DInterface::SetNode(int rLocalNodeNumber, NodeBase* rNode)
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->GetNumNodes());
    assert(rNode != nullptr);
    mNodes[rLocalNodeNumber] = rNode;
}

void NuTo::Element2DInterface::ResizeNodes(int rNewNumNodes)
{
    if (rNewNumNodes == (int) mNodes.size())
        return;

    if (rNewNumNodes > (int) mNodes.size())
    {
        // just resize (enlarge)
        mNodes.resize(rNewNumNodes);
    }
    else
    {
        throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Resize that reduces the number of nodes is not implemented yet.");
    }
}

void NuTo::Element2DInterface::ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
{
    throw MechanicsException(__PRETTY_FUNCTION__, "IMPLEMENT ME!");
}

const Eigen::VectorXd NuTo::Element2DInterface::GetIntegrationPointVolume() const
{
    throw MechanicsException(__PRETTY_FUNCTION__, "IMPLEMENT ME!");
}

void NuTo::Element2DInterface::CalculateGlobalRowDofs(BlockFullVector<int>& rGlobalRowDofs) const
{
    const unsigned globalDimension = GetStructure()->GetDimension();
    for (auto dof : mInterpolationType->GetActiveDofs())
     {
         const InterpolationBase& interpolationType = mInterpolationType->Get(dof);
         const int numNodes = interpolationType.GetNumNodes();
         FullVector<int, Eigen::Dynamic>& dofWiseGlobalRowDofs = rGlobalRowDofs[dof];

         dofWiseGlobalRowDofs.Resize(interpolationType.GetNumDofs());
         dofWiseGlobalRowDofs.setZero();
         switch (dof)
         {
         case Node::eDof::DISPLACEMENTS:
         {
             for (int iNodeDof = 0; iNodeDof < numNodes; ++iNodeDof)
             {
                 const NodeBase* nodePtr = mNodes[interpolationType.GetNodeIndex(iNodeDof)];
                 for (unsigned int iDof = 0; iDof < globalDimension; ++iDof)
                     dofWiseGlobalRowDofs[globalDimension * iNodeDof + iDof] = nodePtr->GetDof(Node::eDof::DISPLACEMENTS, iDof);
             }

         }
             break;
         default:
             throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for " + Node::DofToString(dof) + ".");
         }
     }
}

void NuTo::Element2DInterface::CalculateConstitutiveInputs(const ConstitutiveInputMap& rConstitutiveInput, EvaluateData& rData)
{
    for (auto& it : rConstitutiveInput)
     {
         switch (it.first)
        {
        case Constitutive::eInput::INTERFACE_SLIP:
        {
            Eigen::MatrixXd rotationMatrix = CalculateRotationMatrix();
            auto& slip = *static_cast<ConstitutiveMatrixXd*>(it.second.get());
            slip.AsMatrix() = rotationMatrix * (rData.mMatrixB[Node::eDof::DISPLACEMENTS] * rData.mNodalValues[Node::eDof::DISPLACEMENTS]);
        }
            break;
        case Constitutive::eInput::TIME_STEP:
        case Constitutive::eInput::CALCULATE_STATIC_DATA:

            break;
        default:
             throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive input for " + Constitutive::InputToString(it.first) + " not implemented.");
         }
     }
}

Eigen::MatrixXd NuTo::Element2DInterface::CalculateRotationMatrix()
{
    const Eigen::VectorXd globalNodeCoordinates = ExtractNodeValues(0, Node::eDof::COORDINATES);

    unsigned int globalDimension = GetStructure()->GetDimension();

    assert((globalDimension == 2 or globalDimension == 3) and "Need 2d or 3d coordinates");

    // the rotation matrix consists of three basis vectors that define the new coordinate system.
    Eigen::Vector3d basisVector00 = Eigen::Vector3d::Zero();
    Eigen::Vector3d basisVector01 = Eigen::Vector3d::Zero();
    Eigen::Vector3d basisVector02 = Eigen::Vector3d::Zero();

    // basisVector00 is aligned with the truss. Note: Trusses have to be straight!
    basisVector00.block(0, 0, globalDimension, 1) = globalNodeCoordinates.segment(globalDimension, globalDimension) - globalNodeCoordinates.head(globalDimension);

    // check if basisVector00 is linear independent to unitVectorZ and calculate the cross product to determine basisVector01
    if (basisVector00(0, 0) > 1e-8 or basisVector00(1, 0) > 1e-8)
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

Eigen::MatrixXd NuTo::Element2DInterface::CalculateTransformationMatrix(unsigned int rGlobalDimension, unsigned int rNumberOfNodes)
{

    Eigen::MatrixXd rotationMatrix;
    rotationMatrix = CalculateRotationMatrix();
    Eigen::MatrixXd transformationMatrix(rGlobalDimension * rNumberOfNodes, rGlobalDimension * rNumberOfNodes);
    transformationMatrix.setZero();

    for (unsigned int i = 0; i < rNumberOfNodes; ++i)
    {
        transformationMatrix.block(rGlobalDimension * i, rGlobalDimension * i, rGlobalDimension, rGlobalDimension) = rotationMatrix;
    }

    return transformationMatrix;

}

void NuTo::Element2DInterface::CalculateElementOutputs(std::map<Element::eOutput, std::shared_ptr<ElementOutputBase> >& rElementOutput, EvaluateData& rData, int rTheIP, const ConstitutiveOutputMap& constitutiveOutputMap) const
{
    rData.mDetJxWeightIPxSection = rData.mDetJacobian * GetIntegrationPointWeight(rTheIP) * mSection->GetCircumference();

        for (auto it : rElementOutput)
        {
            switch (it.first)
            {
            case Element::eOutput::INTERNAL_GRADIENT:
                CalculateElementOutputInternalGradient(it.second->GetBlockFullVectorDouble(), rData, rTheIP, constitutiveOutputMap);
                break;

            case Element::eOutput::HESSIAN_0_TIME_DERIVATIVE:
                CalculateElementOutputHessian0(it.second->GetBlockFullMatrixDouble(), rData, rTheIP, constitutiveOutputMap);
                break;
            case Element::eOutput::IP_DATA:
                //CalculateElementOutputIpData(it.second->GetIpData(), rData, rTheIP);
                break;
            case Element::eOutput::HESSIAN_1_TIME_DERIVATIVE:
            case Element::eOutput::HESSIAN_2_TIME_DERIVATIVE:
            case Element::eOutput::LUMPED_HESSIAN_2_TIME_DERIVATIVE:
            case Element::eOutput::UPDATE_STATIC_DATA:
            case Element::eOutput::UPDATE_TMP_STATIC_DATA:
            case Element::eOutput::GLOBAL_ROW_DOF:
            case Element::eOutput::GLOBAL_COLUMN_DOF:
                break;
            default:
                throw MechanicsException(__PRETTY_FUNCTION__, "element output not implemented.");
            }
        }

}

void NuTo::Element2DInterface::CalculateElementOutputInternalGradient(BlockFullVector<double>& rInternalGradient, EvaluateData& rData, int rTheIP, const ConstitutiveOutputMap& constitutiveOutputMap) const
{
    for (auto dofRow : mInterpolationType->GetActiveDofs())
    {
        switch (dofRow)
        {
        case Node::eDof::DISPLACEMENTS:
        {
            const auto& bondStress = *static_cast<ConstitutiveMatrixXd*>(constitutiveOutputMap.at(NuTo::Constitutive::eOutput::BOND_STRESS).get());
            rInternalGradient[dofRow] += mTransformationMatrix * rData.mDetJxWeightIPxSection *  rData.mMatrixB.at(dofRow).transpose() * bondStress;
            break;
        }
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Element output INTERNAL_GRADIENT for " + Node::DofToString(dofRow) + " not implemented.");
        }
    }
}

void NuTo::Element2DInterface::CalculateElementOutputHessian0(BlockFullMatrix<double>& rHessian0, EvaluateData& rData, int rTheIP, const ConstitutiveOutputMap& constitutiveOutputMap) const
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
                const auto& tangentBondStressSlip = *static_cast<ConstitutiveMatrixXd*>(constitutiveOutputMap.at(NuTo::Constitutive::eOutput::INTERFACE_CONSTITUTIVE_MATRIX).get());
                hessian0 += mTransformationMatrix * rData.mDetJxWeightIPxSection *  rData.mMatrixB.at(dofRow).transpose() * tangentBondStressSlip * rData.mMatrixB.at(dofRow) * mTransformationMatrix.transpose();
                break;
            }
            default:
                throw MechanicsException(__PRETTY_FUNCTION__, "Element output HESSIAN_0_TIME_DERIVATIVE for "
                        "(" + Node::DofToString(dofRow) + "," + Node::DofToString(dofCol) + ") not implemented.");
            }
        }
    }

}

void NuTo::Element2DInterface::CalculateElementOutputIpData(ElementOutputIpData& rIpData, EvaluateData& rData, int rTheIP) const
{
}

Eigen::VectorXd NuTo::Element2DInterface::ExtractNodeValues(int rTimeDerivative, Node::eDof rDofType) const
{
    const InterpolationBase& interpolationTypeDof = GetInterpolationType().Get(rDofType);

    const unsigned globalDimension = GetStructure()->GetDimension();
    int numNodes = interpolationTypeDof.GetNumNodes();
    int numDofs = interpolationTypeDof.GetNumDofs();
//    int numDofsPerNode = numDofs / numNodes; // --- unused so far

    Eigen::VectorXd nodeValues = Eigen::VectorXd::Constant(numDofs, -1337.0);


    for (int iNode = 0; iNode < numNodes; ++iNode)
    {
        const NodeBase* node = mNodes[interpolationTypeDof.GetNodeIndex(iNode)];

        switch (rDofType)
        {
        case Node::eDof::COORDINATES:
        {
            Eigen::MatrixXd tmp(node->Get(Node::eDof::COORDINATES));

            nodeValues.segment(iNode * globalDimension, globalDimension) = node->Get(Node::eDof::COORDINATES);
        }
            break;
        case Node::eDof::DISPLACEMENTS:
            nodeValues.segment(iNode * globalDimension, globalDimension) = node->Get(Node::eDof::DISPLACEMENTS);
            break;
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for " + Node::DofToString(rDofType));
        }
    }


    return nodeValues;

}


void NuTo::Element2DInterface::CheckElement()
{
}

void NuTo::Element2DInterface::FillConstitutiveOutputMapInternalGradient(ConstitutiveOutputMap& rConstitutiveOutput, BlockFullVector<double>& rInternalGradient) const
{
    for (auto dofRow : mInterpolationType->GetActiveDofs())
    {
        rInternalGradient[dofRow].setZero(mInterpolationType->Get(dofRow).GetNumDofs());
        switch (dofRow)
        {
        case Node::eDof::DISPLACEMENTS:
            rConstitutiveOutput[NuTo::Constitutive::eOutput::BOND_STRESS] = ConstitutiveIOBase::makeConstitutiveIO<2>(NuTo::Constitutive::eOutput::BOND_STRESS);
            break;
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive output INTERNAL_GRADIENT for " + Node::DofToString(dofRow) + " not implemented.");

        }
    }
}

void NuTo::Element2DInterface::FillConstitutiveOutputMapHessian0(ConstitutiveOutputMap& rConstitutiveOutput, BlockFullMatrix<double>& rHessian0) const
{


    for (auto dofRow : mInterpolationType->GetActiveDofs())
    {
        for (auto dofCol : mInterpolationType->GetActiveDofs())
        {
            NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& dofSubMatrix = rHessian0(dofRow, dofCol);
            dofSubMatrix.Resize(mInterpolationType->Get(dofRow).GetNumDofs(), mInterpolationType->Get(dofCol).GetNumDofs());
            dofSubMatrix.setZero();

            switch (Node::CombineDofs(dofRow, dofCol))
            {
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS):
                rConstitutiveOutput[NuTo::Constitutive::eOutput::INTERFACE_CONSTITUTIVE_MATRIX] = ConstitutiveIOBase::makeConstitutiveIO<2>(NuTo::Constitutive::eOutput::INTERFACE_CONSTITUTIVE_MATRIX);
                break;
            default:
                throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive output HESSIAN_0_TIME_DERIVATIVE for "
                        "(" + Node::DofToString(dofRow) + "," + Node::DofToString(dofCol) + ") not implemented.");
            }
        }
    }

}

#ifdef ENABLE_VISUALIZE
void NuTo::Element2DInterface::GetVisualizationCells(unsigned int& NumVisualizationPoints, std::vector<double>& VisualizationPointLocalCoordinates, unsigned int& NumVisualizationCells, std::vector<NuTo::eCellTypes>& VisualizationCellType,
        std::vector<unsigned int>& VisualizationCellsIncidence, std::vector<unsigned int>& VisualizationCellsIP) const
{

    const IntegrationTypeBase& integrationType = GetIntegrationType();
    const unsigned int globalDimension = GetStructure()->GetDimension();
    switch (integrationType.GetNumIntegrationPoints())
    {
    case 2:

        NumVisualizationPoints = 4;
        // Point 0
        VisualizationPointLocalCoordinates.push_back(-1);
        VisualizationPointLocalCoordinates.push_back(-1);
        if (globalDimension == 3)
            VisualizationPointLocalCoordinates.push_back(0);

        // Point 1
        VisualizationPointLocalCoordinates.push_back(+1);
        VisualizationPointLocalCoordinates.push_back(-1);
        if (globalDimension == 3)
                    VisualizationPointLocalCoordinates.push_back(0);
        // Point 2
        VisualizationPointLocalCoordinates.push_back(+1);
        VisualizationPointLocalCoordinates.push_back(+1);
        if (globalDimension == 3)
                    VisualizationPointLocalCoordinates.push_back(0);
        // Point 3
        VisualizationPointLocalCoordinates.push_back(-1);
        VisualizationPointLocalCoordinates.push_back(+1);
        if (globalDimension == 3)
                    VisualizationPointLocalCoordinates.push_back(0);

        NumVisualizationCells = 1;

        // cell 0
        VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
        VisualizationCellsIncidence.push_back(0);
        VisualizationCellsIncidence.push_back(1);
        VisualizationCellsIncidence.push_back(2);
        VisualizationCellsIncidence.push_back(3);
        VisualizationCellsIP.push_back(0);
        break;
    case 3:

        NumVisualizationPoints = 6;
        // Point 0
        VisualizationPointLocalCoordinates.push_back(-1);
        VisualizationPointLocalCoordinates.push_back(-1);

        // Point 1
        VisualizationPointLocalCoordinates.push_back(+0);
        VisualizationPointLocalCoordinates.push_back(-1);

        // Point 2
        VisualizationPointLocalCoordinates.push_back(+1);
        VisualizationPointLocalCoordinates.push_back(-1);

        // Point 3
        VisualizationPointLocalCoordinates.push_back(+1);
        VisualizationPointLocalCoordinates.push_back(+1);

        // Point 4
        VisualizationPointLocalCoordinates.push_back(+0);
        VisualizationPointLocalCoordinates.push_back(+1);

        // Point 5
        VisualizationPointLocalCoordinates.push_back(-1);
        VisualizationPointLocalCoordinates.push_back(+1);

        NumVisualizationCells = 2;

        // cell 0
        VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
        VisualizationCellsIncidence.push_back(0);
        VisualizationCellsIncidence.push_back(1);
        VisualizationCellsIncidence.push_back(4);
        VisualizationCellsIncidence.push_back(5);
        VisualizationCellsIP.push_back(0);

        // cell 0
        VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
        VisualizationCellsIncidence.push_back(1);
        VisualizationCellsIncidence.push_back(2);
        VisualizationCellsIncidence.push_back(3);
        VisualizationCellsIncidence.push_back(4);
        VisualizationCellsIP.push_back(1);

        break;

    case 4:// lobatto

        NumVisualizationPoints = 8;
        // Point 0
        VisualizationPointLocalCoordinates.push_back(-1);
        VisualizationPointLocalCoordinates.push_back(-1);

        // Point 1
        VisualizationPointLocalCoordinates.push_back(-1. + 2. / 3.);
        VisualizationPointLocalCoordinates.push_back(-1);

        // Point 2
        VisualizationPointLocalCoordinates.push_back(-1. + 4. / 3.);
        VisualizationPointLocalCoordinates.push_back(-1);

        // Point 3
        VisualizationPointLocalCoordinates.push_back(+1);
        VisualizationPointLocalCoordinates.push_back(-1);

        // Point 4
        VisualizationPointLocalCoordinates.push_back(+1);
        VisualizationPointLocalCoordinates.push_back(+1);

        // Point 5
        VisualizationPointLocalCoordinates.push_back(-1. + 4. / 3.);
        VisualizationPointLocalCoordinates.push_back(+1);

        // Point 6
        VisualizationPointLocalCoordinates.push_back(-1. + 2. / 3.);
        VisualizationPointLocalCoordinates.push_back(+1);

        // Point 7
        VisualizationPointLocalCoordinates.push_back(-1);
        VisualizationPointLocalCoordinates.push_back(+1);

        NumVisualizationCells = 3;

        // cell 0
        VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
        VisualizationCellsIncidence.push_back(0);
        VisualizationCellsIncidence.push_back(1);
        VisualizationCellsIncidence.push_back(6);
        VisualizationCellsIncidence.push_back(7);
        VisualizationCellsIP.push_back(0);

        // cell 1
        VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
        VisualizationCellsIncidence.push_back(1);
        VisualizationCellsIncidence.push_back(2);
        VisualizationCellsIncidence.push_back(5);
        VisualizationCellsIncidence.push_back(6);
        VisualizationCellsIP.push_back(1);

        // cell 2
        VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
        VisualizationCellsIncidence.push_back(2);
        VisualizationCellsIncidence.push_back(3);
        VisualizationCellsIncidence.push_back(4);
        VisualizationCellsIncidence.push_back(5);
        VisualizationCellsIP.push_back(1);
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Integration type not valid for this element.");
    }

}

void NuTo::Element2DInterface::Visualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList)
{

    // get visualization cells from integration type
    unsigned int NumVisualizationPoints;
    std::vector<double> VisualizationPointLocalCoordinates;
    unsigned int NumVisualizationCells;
    std::vector<NuTo::eCellTypes> VisualizationCellType;
    std::vector<unsigned int> VisualizationCellsIncidence;
    std::vector<unsigned int> VisualizationCellsIP;
    //get the visualization cells either from the integration type (standard)
    // or (if the routine is rewritten for e.g. XFEM or lattice elements, from other element data
    GetVisualizationCells(NumVisualizationPoints, VisualizationPointLocalCoordinates, NumVisualizationCells, VisualizationCellType, VisualizationCellsIncidence, VisualizationCellsIP);

    // calculate global point coordinates and store point at the visualize object
    int dimension(VisualizationPointLocalCoordinates.size() / NumVisualizationPoints);
    assert(VisualizationPointLocalCoordinates.size() == NumVisualizationPoints * dimension);

    // TODO: fix that by proper visualization point (natural!) coordinates, local might be misleading here
    Eigen::MatrixXd visualizationPointNaturalCoordinates = Eigen::MatrixXd::Map(VisualizationPointLocalCoordinates.data(), dimension, NumVisualizationPoints);

    std::vector<unsigned int> PointIdVec;
    for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
    {
        if (dimension != 1 and dimension != 2 and dimension != 3)
            throw NuTo::MechanicsException("[NuTo::ElementBase::Visualize] invalid dimension of local coordinates");

        Eigen::Vector3d GlobalPointCoor = Eigen::Vector3d::Zero();

        GlobalPointCoor.head(dimension) = ExtractNodeValues(0, Node::eDof::COORDINATES).segment(dimension*PointCount, dimension);
        unsigned int PointId = rVisualize.AddPoint(GlobalPointCoor.data());
        PointIdVec.push_back(PointId);
    }

    // store cells at the visualize object
    assert(VisualizationCellType.size() == NumVisualizationCells);
    assert(VisualizationCellsIP.size() == NumVisualizationCells);
    std::vector<unsigned int> CellIdVec;
    unsigned int Pos = 0;
    for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
    {
        switch (VisualizationCellType[CellCount])
        {
        case NuTo::eCellTypes::QUAD:
        {
            assert(Pos + 4 <= VisualizationCellsIncidence.size());
            unsigned int Points[4];
            for (unsigned int PointCount = 0; PointCount < 4; PointCount++)
            {
                Points[PointCount] = PointIdVec[VisualizationCellsIncidence[Pos + PointCount]];
            }
            unsigned int CellId = rVisualize.AddQuadCell(Points);
            CellIdVec.push_back(CellId);
            Pos += 4;
        }
            break;
        default:
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "unsupported visualization cell type");
        }
    }

    //determine the ipdata and determine the map
    std::map<NuTo::Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutput;
    elementOutput[Element::eOutput::IP_DATA] = std::make_shared<ElementOutputIpData>();
    auto& elementIpDataMap = elementOutput.at(Element::eOutput::IP_DATA)->GetIpData().GetIpDataMap();

    for (auto const &it : rVisualizationList)
    {
        switch (it.get()->GetComponentEnum())
        {
        case NuTo::eVisualizeWhat::BOND_STRESS:
            elementIpDataMap[IpData::eIpStaticDataType::BOND_STRESS];
            break;
        case NuTo::eVisualizeWhat::ENGINEERING_STRESS:
        case NuTo::eVisualizeWhat::DISPLACEMENTS:
        case NuTo::eVisualizeWhat::SECTION:
        case NuTo::eVisualizeWhat::ELEMENT:
        default:
            //do nothing
            break;
        }
    }

    //calculate the element solution
    ConstitutiveInputMap input;
    input[Constitutive::eInput::CALCULATE_STATIC_DATA] = std::make_unique<ConstitutiveCalculateStaticData>(eCalculateStaticData::EULER_BACKWARD);


    Evaluate(input, elementOutput);


    // store data
    for (auto const &it : rVisualizationList)
    {
        switch (it.get()->GetComponentEnum())
        {
        case NuTo::eVisualizeWhat::DISPLACEMENTS:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                Eigen::Vector3d GlobalDisplacements = Eigen::Vector3d::Zero();

                GlobalDisplacements.head(dimension) = ExtractNodeValues(0, Node::eDof::DISPLACEMENTS).segment(dimension*PointCount, dimension);
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(), GlobalDisplacements.data());
            }
            break;
        case NuTo::eVisualizeWhat::ENGINEERING_STRESS:
            break;
        case NuTo::eVisualizeWhat::CONSTITUTIVE:
        {
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                unsigned int CellId = CellIdVec[CellCount];
                int constitutiveId = mStructure->ConstitutiveLawGetId(&GetConstitutiveLaw(theIp));

                rVisualize.SetCellDataScalar(CellId, it.get()->GetComponentName(), constitutiveId);
            }
        }
            break;
        case NuTo::eVisualizeWhat::SECTION:
        {
            int sectionId = mStructure->SectionGetId(&GetSection());
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataScalar(CellId, it.get()->GetComponentName(), sectionId);
            }
        }
            break;
        case NuTo::eVisualizeWhat::ELEMENT:
        {
            int elementId = this->ElementGetId();
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataScalar(CellId, it.get()->GetComponentName(), elementId);
            }
        }
            break;
        case NuTo::eVisualizeWhat::BOND_STRESS:
        {
            const auto& bondStress = elementIpDataMap.at(IpData::eIpStaticDataType::BOND_STRESS);
            assert(bondStress.size() != 0);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                double bondStressTensor[9];
                bondStressTensor[0] = bondStress(0, theIp);
                bondStressTensor[1] = bondStress(1, theIp);
                bondStressTensor[2] = 0.0;
                bondStressTensor[3] = 0.0;
                bondStressTensor[4] = 0.0;
                bondStressTensor[5] = 0.0;
                bondStressTensor[6] = 0.0;
                bondStressTensor[7] = 0.0;
                bondStressTensor[8] = 0.0;

                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataTensor(CellId, it.get()->GetComponentName(), bondStressTensor);
            }
        }
            break;
        default:
            throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t unsupported datatype for visualization.");
        }
    }

}
#endif // ENABLE_VISUALIZE
void NuTo::Element2DInterface::FillConstitutiveOutputMapIpData(ConstitutiveOutputMap& rConstitutiveOutput, ElementOutputIpData& rIpData) const
{

    for (auto& it : rIpData.GetIpDataMap()) // this reference here is _EXTREMLY_ important, since the GetIpDataMap() contains a
    {                                       // FullMatrix VALUE and you want to access this value by reference. Without the &, a tmp copy would be made.
        switch (it.first)
        {
        case NuTo::IpData::eIpStaticDataType::BOND_STRESS:
            it.second.Resize(6, GetNumIntegrationPoints());
            rConstitutiveOutput[NuTo::Constitutive::eOutput::BOND_STRESS] = ConstitutiveIOBase::makeConstitutiveIO<2>(NuTo::Constitutive::eOutput::BOND_STRESS);
            break;
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "this ip data type is not implemented.");
        }
    }

}
void NuTo::Element2DInterface::SetSection(const NuTo::SectionBase& rSection)
{
    mSection = &rSection;
}
const NuTo::SectionBase& NuTo::Element2DInterface::GetSection() const
{
    if (mSection != nullptr)
        return *mSection;

    Info();
    throw MechanicsException(__PRETTY_FUNCTION__, "This element has no section assigned yet.");
}




