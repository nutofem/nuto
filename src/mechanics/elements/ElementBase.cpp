#include "mechanics/elements/ElementBase.h"

#include <iostream>

#include "base/Exception.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "mechanics/constitutive/inputoutput/ConstitutivePlaneState.h"
#include "mechanics/elements/ElementOutputIpData.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/elements/IpDataEnum.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "mechanics/interpolationtypes/InterpolationBase.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/sections/Section.h"
#include "visualize/ComponentName.h"

#include "math/EigenCompanion.h"

#include <eigen3/Eigen/Core>

#ifdef ENABLE_VISUALIZE
#include "visualize/Point.h"
#include "visualize/Cell.h"
#include "visualize/UnstructuredGrid.h"
#endif

using namespace NuTo;

namespace NuTo
{
std::ostream& operator<<(std::ostream& out, const ElementBase& element)
{
    element.Info(out);
    return out;
}
} /* NuTo */

NuTo::ElementBase::ElementBase(const InterpolationType& interpolationType, const IntegrationTypeBase& integrationType)
    : mInterpolationType(&interpolationType)
    , mIPData(integrationType)
{
}


void NuTo::ElementBase::Evaluate(std::map<ElementEnum::eOutput, std::shared_ptr<ElementOutputBase>>& rOutput)
{
    ConstitutiveInputMap input;
    input[Constitutive::eInput::CALCULATE_STATIC_DATA] =
            std::make_unique<ConstitutiveCalculateStaticData>(eCalculateStaticData::EULER_BACKWARD);

    return this->Evaluate(input, rOutput);
}

const NuTo::ConstitutiveBase& NuTo::ElementBase::GetConstitutiveLaw(unsigned int rIP) const
{
    return mIPData.GetIPConstitutiveLaw(rIP).GetConstitutiveLaw();
}

NuTo::ConstitutiveBase& NuTo::ElementBase::GetConstitutiveLaw(unsigned int rIP)
{
    return mIPData.GetIPConstitutiveLaw(rIP).GetConstitutiveLaw();
}

NuTo::IPData& NuTo::ElementBase::GetIPData()
{
    return mIPData;
}

void NuTo::ElementBase::SetConstitutiveLaw(ConstitutiveBase& rConstitutiveLaw)
{
    mIPData.SetConstitutiveLaw(rConstitutiveLaw);
}

bool NuTo::ElementBase::HasConstitutiveLawAssigned(unsigned int rIP) const
{
    return mIPData.HasConstitutiveLawAssigned(rIP);
}

void NuTo::ElementBase::SetSection(std::shared_ptr<const Section>)
{
    throw Exception(__PRETTY_FUNCTION__, "This element type has so section.");
}

std::shared_ptr<const Section> NuTo::ElementBase::GetSection() const
{
    throw Exception(__PRETTY_FUNCTION__, "This element type has so section.");
}


int NuTo::ElementBase::GetNumNodes() const
{
    return mInterpolationType->GetNumNodes();
}

int NuTo::ElementBase::GetNumNodes(Node::eDof rDofType) const
{
    return mInterpolationType->Get(rDofType).GetNumNodes();
}

Eigen::VectorXd NuTo::ElementBase::ExtractNodeValues(int, Node::eDof) const
{
    throw NuTo::Exception("[NuTo::ElementBase::ExtractNodeValues] not implemented.");
}


Eigen::VectorXd NuTo::ElementBase::InterpolateDofGlobal(const Eigen::VectorXd& rNaturalCoordinates,
                                                        Node::eDof rDofType) const
{
    return InterpolateDofGlobal(0, rNaturalCoordinates, rDofType);
}

Eigen::VectorXd NuTo::ElementBase::InterpolateDofGlobal(int rTimeDerivative, const Eigen::VectorXd& rNaturalCoordinates,
                                                        Node::eDof rDofType) const
{

    const InterpolationBase& interpolationType = mInterpolationType->Get(rDofType);
    Eigen::MatrixXd nodalValues = ExtractNodeValues(rTimeDerivative, rDofType);
    Eigen::MatrixXd matrixN = interpolationType.MatrixN(rNaturalCoordinates);

    return matrixN * nodalValues;
}

Eigen::Vector3d NuTo::ElementBase::InterpolateDof3D(const Eigen::VectorXd& rNaturalCoordinates,
                                                    Node::eDof rDofType) const
{
    return InterpolateDof3D(0, rNaturalCoordinates, rDofType);
}

Eigen::Vector3d NuTo::ElementBase::InterpolateDof3D(int rTimeDerivative, const Eigen::VectorXd& rNaturalCoordinates,
                                                    Node::eDof rDofType) const
{
    return EigenCompanion::To3D(InterpolateDofGlobal(rTimeDerivative, rNaturalCoordinates, rDofType));
}


void NuTo::ElementBase::SetIntegrationType(const NuTo::IntegrationTypeBase& rIntegrationType)
{
    // check compatibility between element type and constitutive law
    if (GetLocalDimension() == rIntegrationType.GetDimension())
    {
        mIPData.SetIntegrationType(rIntegrationType);
    }
    else
    {
        throw Exception(__PRETTY_FUNCTION__, "Integration Type does not match element type of element.");
    }
}

const NuTo::IntegrationTypeBase& NuTo::ElementBase::GetIntegrationType() const
{
    return mIPData.GetIntegrationType();
}

void NuTo::ElementBase::SetInterpolationType(const InterpolationType& rInterpolationType)
{
    mInterpolationType = &rInterpolationType;
}

const NuTo::InterpolationType& NuTo::ElementBase::GetInterpolationType() const
{
    // mInterpolationType only set via references. No nullptr check required.
    return *mInterpolationType;
}


int NuTo::ElementBase::GetNumIntegrationPoints() const
{
    return mIPData.GetIntegrationType().GetNumIntegrationPoints();
}

double NuTo::ElementBase::GetIntegrationPointWeight(unsigned int rIP) const
{
    return mIPData.GetIntegrationType().GetIntegrationPointWeight(rIP);
}


template <int TDim>
void NuTo::ElementBase::EvaluateConstitutiveLaw(const NuTo::ConstitutiveInputMap& rConstitutiveInput,
                                                NuTo::ConstitutiveOutputMap& rConstitutiveOutput, unsigned int IP)
{
    Constitutive::IPConstitutiveLawBase& ipConstitutiveLaw = mIPData.GetIPConstitutiveLaw(IP);

    for (auto& itOutput : rConstitutiveOutput)
        if (itOutput.second != nullptr) // check nullptr because of static data
            itOutput.second->SetIsCalculated(false);
    ipConstitutiveLaw.Evaluate<TDim>(rConstitutiveInput, rConstitutiveOutput);

    for (auto& itOutput : rConstitutiveOutput)
        if (itOutput.second != nullptr && !itOutput.second->GetIsCalculated()) // check nullptr because of static data
            throw Exception(__PRETTY_FUNCTION__, "Output " + Constitutive::OutputToString(itOutput.first) +
                                                                  " not calculated by constitutive law");
}


template void NuTo::ElementBase::EvaluateConstitutiveLaw<1>(const NuTo::ConstitutiveInputMap&,
                                                            NuTo::ConstitutiveOutputMap&, unsigned int);
template void NuTo::ElementBase::EvaluateConstitutiveLaw<2>(const NuTo::ConstitutiveInputMap&,
                                                            NuTo::ConstitutiveOutputMap&, unsigned int);
template void NuTo::ElementBase::EvaluateConstitutiveLaw<3>(const NuTo::ConstitutiveInputMap&,
                                                            NuTo::ConstitutiveOutputMap&, unsigned int);


const Eigen::Vector3d NuTo::ElementBase::GetGlobalIntegrationPointCoordinates(int rIpNum) const
{
    const auto naturalCoords = GetIntegrationType().GetLocalIntegrationPointCoordinates(rIpNum);
    const Eigen::MatrixXd& matrixN = mInterpolationType->Get(Node::eDof::COORDINATES).MatrixN(naturalCoords);
    Eigen::VectorXd nodeCoordinates = ExtractNodeValues(0, Node::eDof::COORDINATES);

    Eigen::Vector3d globalIntegrationPointCoordinates = Eigen::Vector3d::Zero();

    globalIntegrationPointCoordinates.segment(0, GetLocalDimension()) = matrixN * nodeCoordinates;

    return globalIntegrationPointCoordinates;
}

bool NuTo::ElementBase::GetLocalPointCoordinates(const double*, double*) const
{
    throw NuTo::Exception(
            "[NuTo::ElementBase::GetLocalPointCoordinates] not implemented for this element type.");
}

NuTo::NodeBase* NuTo::ElementBase::GetBoundaryControlNode() const
{
    throw NuTo::Exception(__PRETTY_FUNCTION__, "Not implemented for this element type.");
}


#ifdef ENABLE_VISUALIZE
void NuTo::ElementBase::GetVisualizationCells(unsigned int& NumVisualizationPoints,
                                              std::vector<double>& VisualizationPointLocalCoordinates,
                                              unsigned int& NumVisualizationCells,
                                              std::vector<NuTo::eCellTypes>& VisualizationCellType,
                                              std::vector<unsigned int>& VisualizationCellsIncidence,
                                              std::vector<unsigned int>& VisualizationCellsIP) const
{
    GetIntegrationType().GetVisualizationCells(NumVisualizationPoints, VisualizationPointLocalCoordinates,
                                               NumVisualizationCells, VisualizationCellType,
                                               VisualizationCellsIncidence, VisualizationCellsIP);
}


NuTo::IpData::eIpStaticDataType ToIpDataEnum(NuTo::eVisualizeWhat what)
{
    switch (what)
    {
    case NuTo::eVisualizeWhat::BOND_STRESS:
        return IpData::eIpStaticDataType::BOND_STRESS;
    case NuTo::eVisualizeWhat::DAMAGE:
        return IpData::eIpStaticDataType::DAMAGE;
    case NuTo::eVisualizeWhat::ENGINEERING_PLASTIC_STRAIN:
        return IpData::eIpStaticDataType::ENGINEERING_PLASTIC_STRAIN;
    case NuTo::eVisualizeWhat::ENGINEERING_STRAIN:
        return IpData::eIpStaticDataType::ENGINEERING_STRAIN;
    case NuTo::eVisualizeWhat::SHRINKAGE_STRAIN:
        return IpData::eIpStaticDataType::SHRINKAGE_STRAIN;
    case NuTo::eVisualizeWhat::THERMAL_STRAIN:
        return IpData::eIpStaticDataType::THERMAL_STRAIN;
    case NuTo::eVisualizeWhat::ENGINEERING_STRESS:
        return IpData::eIpStaticDataType::ENGINEERING_STRESS;
    case NuTo::eVisualizeWhat::HEAT_FLUX:
        return IpData::eIpStaticDataType::HEAT_FLUX;
    case NuTo::eVisualizeWhat::ELECTRIC_FIELD:
        return IpData::eIpStaticDataType::ELECTRIC_FIELD;
    case NuTo::eVisualizeWhat::ELECTRIC_DISPLACEMENT:
        return IpData::eIpStaticDataType::ELECTRIC_DISPLACEMENT;
    case NuTo::eVisualizeWhat::LOCAL_EQ_STRAIN:
        return IpData::eIpStaticDataType::LOCAL_EQ_STRAIN;
    case NuTo::eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS:
        return IpData::eIpStaticDataType::ENGINEERING_STRESS;
    case NuTo::eVisualizeWhat::TOTAL_INELASTIC_EQ_STRAIN:
        return IpData::eIpStaticDataType::TOTAL_INELASTIC_EQ_STRAIN;
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "No conversion from eVisualizeWhat to eIpStaticDataType");
    }
}

NuTo::Node::eDof ToNodeEnum(NuTo::eVisualizeWhat what)
{
    switch (what)
    {
    case NuTo::eVisualizeWhat::TEMPERATURE:
        return NuTo::Node::eDof::TEMPERATURE;
    case NuTo::eVisualizeWhat::NONLOCAL_EQ_STRAIN:
        return NuTo::Node::eDof::NONLOCALEQSTRAIN;
    case NuTo::eVisualizeWhat::DISPLACEMENTS:
        return NuTo::Node::eDof::DISPLACEMENTS;
    case NuTo::eVisualizeWhat::CRACK_PHASE_FIELD:
        return NuTo::Node::eDof::CRACKPHASEFIELD;
    case NuTo::eVisualizeWhat::RELATIVE_HUMIDITY:
        return NuTo::Node::eDof::RELATIVEHUMIDITY;
    case NuTo::eVisualizeWhat::ELECTRIC_POTENTIAL:
        return NuTo::Node::eDof::ELECTRICPOTENTIAL;
    case NuTo::eVisualizeWhat::WATER_VOLUME_FRACTION:
        return NuTo::Node::eDof::WATERVOLUMEFRACTION;
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "No conversion from eVisualizeWhat to eDof");
    }
}

void NuTo::ElementBase::Visualize(Visualize::UnstructuredGrid& visualizer,
                                  const std::vector<eVisualizeWhat>& visualizeComponents)
{
    IntegrationTypeBase::IpCellInfo ipCellInfo = GetIntegrationType().GetVisualizationCells();
    auto& cells = ipCellInfo.cells;
    auto& points = ipCellInfo.vertices;
    if (ipCellInfo.cells.size() == 0)
        return; // nothing to visualize

    // add visualization points and store their id together with their local coordinates
    for (auto& pointInfo : points)
    {
        Eigen::Vector3d globalCoords = InterpolateDof3D(pointInfo.localCoords, Node::eDof::COORDINATES);
        pointInfo.visualizePointId = visualizer.AddPoint(globalCoords);
    }

    for (auto& cellInfo : cells)
    {
        // transform the ids of the cell points
        std::vector<int> globalPointIds;
        globalPointIds.reserve(cellInfo.pointIds.size());
        for (int localPointId : cellInfo.pointIds)
            globalPointIds.push_back(points[localPointId].visualizePointId);
        cellInfo.visualizeCellId = visualizer.AddCell(globalPointIds, cellInfo.cellType);
    }

    // determine the ipdata and determine the map
    std::map<NuTo::ElementEnum::eOutput, std::shared_ptr<ElementOutputBase>> elementOutput;
    elementOutput[ElementEnum::eOutput::IP_DATA] = std::make_shared<ElementOutputIpData>();

    auto& elementIpDataMap = elementOutput.at(ElementEnum::eOutput::IP_DATA)->GetIpData().GetIpDataMap();

    bool evaluateStress(false);

    for (auto component : visualizeComponents)
    {
        switch (component)
        {
        case NuTo::eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS:
        case NuTo::eVisualizeWhat::ENGINEERING_STRESS:
            if (evaluateStress == false)
                evaluateStress = true;
        // no break here...
        case NuTo::eVisualizeWhat::BOND_STRESS:
        case NuTo::eVisualizeWhat::DAMAGE:
        case NuTo::eVisualizeWhat::ENGINEERING_PLASTIC_STRAIN:
        case NuTo::eVisualizeWhat::ENGINEERING_STRAIN:
        case NuTo::eVisualizeWhat::SHRINKAGE_STRAIN:
        case NuTo::eVisualizeWhat::THERMAL_STRAIN:
        case NuTo::eVisualizeWhat::HEAT_FLUX:
        case NuTo::eVisualizeWhat::ELECTRIC_FIELD:
        case NuTo::eVisualizeWhat::ELECTRIC_DISPLACEMENT:
        case NuTo::eVisualizeWhat::LOCAL_EQ_STRAIN:
        case NuTo::eVisualizeWhat::TOTAL_INELASTIC_EQ_STRAIN:
            elementIpDataMap[ToIpDataEnum(component)];
            break;


        case NuTo::eVisualizeWhat::ACCELERATION:
        case NuTo::eVisualizeWhat::ANGULAR_VELOCITY:
        case NuTo::eVisualizeWhat::ANGULAR_ACCELERATION:
        case NuTo::eVisualizeWhat::DISPLACEMENTS:
        case NuTo::eVisualizeWhat::NONLOCAL_EQ_STRAIN:
        case NuTo::eVisualizeWhat::PARTICLE_RADIUS:
        case NuTo::eVisualizeWhat::RELATIVE_HUMIDITY:
        case NuTo::eVisualizeWhat::ROTATION:
        case NuTo::eVisualizeWhat::TEMPERATURE:
        case NuTo::eVisualizeWhat::ELECTRIC_POTENTIAL:
        case NuTo::eVisualizeWhat::VELOCITY:
        case NuTo::eVisualizeWhat::WATER_VOLUME_FRACTION:
        default:
            break;
        }
    }

    // calculate the element solution
    ConstitutiveInputMap input;
    input[Constitutive::eInput::CALCULATE_STATIC_DATA] =
            std::make_unique<ConstitutiveCalculateStaticData>(eCalculateStaticData::USE_PREVIOUS);
    Evaluate(input, elementOutput);
    //    Evaluate(elementOutput);

    // assign the outputs


    // store data
    for (auto component : visualizeComponents)
    {
        auto componentName = GetComponentName(component);
        switch (component)
        {
        case NuTo::eVisualizeWhat::LOCAL_EQ_STRAIN:
        case NuTo::eVisualizeWhat::DAMAGE:
        case NuTo::eVisualizeWhat::SHRINKAGE_STRAIN:
        case NuTo::eVisualizeWhat::TOTAL_INELASTIC_EQ_STRAIN:
        case NuTo::eVisualizeWhat::ENGINEERING_STRESS:
        case NuTo::eVisualizeWhat::BOND_STRESS:
        case NuTo::eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS:
        case NuTo::eVisualizeWhat::HEAT_FLUX:
        case NuTo::eVisualizeWhat::ELECTRIC_FIELD:
        case NuTo::eVisualizeWhat::ELECTRIC_DISPLACEMENT:
        {
            const Eigen::MatrixXd& data = elementIpDataMap.at(ToIpDataEnum(component));
            assert(data.size() != 0);
            for (auto cell : cells)
                visualizer.SetCellData(cell.visualizeCellId, componentName, data.col(cell.ipId));
        }
        break;

        case NuTo::eVisualizeWhat::ENGINEERING_PLASTIC_STRAIN:
        case NuTo::eVisualizeWhat::THERMAL_STRAIN:
        case NuTo::eVisualizeWhat::ENGINEERING_STRAIN:
        {
            Eigen::MatrixXd data = elementIpDataMap.at(ToIpDataEnum(component));
            assert(data.rows() == 6);
            data.bottomRows(3) /= 2.; // transform engineering gamma to epsilon
            for (auto cell : cells)
                visualizer.SetCellData(cell.visualizeCellId, componentName, data.col(cell.ipId));
        }
        break;


        case NuTo::eVisualizeWhat::NONLOCAL_EQ_STRAIN:
        case NuTo::eVisualizeWhat::TEMPERATURE:
        case NuTo::eVisualizeWhat::CRACK_PHASE_FIELD:
        case NuTo::eVisualizeWhat::RELATIVE_HUMIDITY:
        case NuTo::eVisualizeWhat::ELECTRIC_POTENTIAL:
        case NuTo::eVisualizeWhat::WATER_VOLUME_FRACTION:
        {
            auto nodeDof = ToNodeEnum(component);
            for (auto point : points)
                visualizer.SetPointData(point.visualizePointId, componentName,
                                        InterpolateDofGlobal(point.localCoords, nodeDof));
        }
        break;

        case NuTo::eVisualizeWhat::DISPLACEMENTS:
        {
            for (auto point : points)
                visualizer.SetPointData(point.visualizePointId, componentName,
                                        InterpolateDof3D(point.localCoords, Node::eDof::DISPLACEMENTS));
        }
        break;

        case NuTo::eVisualizeWhat::VELOCITY:
        {
            for (auto point : points)
                visualizer.SetPointData(point.visualizePointId, componentName,
                                        InterpolateDof3D(1, point.localCoords, Node::eDof::DISPLACEMENTS));
        }
        break;

        case NuTo::eVisualizeWhat::ACCELERATION:
        {
            for (auto point : points)
                visualizer.SetPointData(point.visualizePointId, componentName,
                                        InterpolateDof3D(2, point.localCoords, Node::eDof::DISPLACEMENTS));
        }
        break;

        case NuTo::eVisualizeWhat::PARTICLE_RADIUS:
        case NuTo::eVisualizeWhat::ROTATION:
        case NuTo::eVisualizeWhat::ANGULAR_VELOCITY:
        case NuTo::eVisualizeWhat::ANGULAR_ACCELERATION:
            // do nothing
            break;
        default:
            throw NuTo::Exception(__PRETTY_FUNCTION__,
                                           "visualization of " + componentName + " not implemented.");
        }
    }
}

void NuTo::ElementBase::VisualizeExtrapolateToNodes(Visualize::UnstructuredGrid&,
                                                    const std::vector<eVisualizeWhat>&)
{
    throw NuTo::Exception(
            std::string(__PRETTY_FUNCTION__) +
            ": \t This function is not ready to be used yet. Choose a different visualization type!");
}

void NuTo::ElementBase::VisualizeIntegrationPointData(Visualize::UnstructuredGrid& visualizer,
                                                      const std::vector<eVisualizeWhat>& visualizeComponents)
{

    struct IpInfo
    {
        int pointId;
        int cellId;
        Eigen::Vector3d localCoords;
    };


    const int numIp = GetNumIntegrationPoints();
    std::vector<IpInfo> ipInfo(numIp);

    for (int i = 0; i < numIp; ++i)
    {
        IpInfo info;
        info.localCoords = GetIntegrationType().GetLocalIntegrationPointCoordinates(i);
        info.pointId = visualizer.AddPoint(info.localCoords);
        info.cellId = visualizer.AddCell({info.pointId}, eCellTypes::VERTEX);
    }

    // determine the ipdata and determine the map
    std::map<NuTo::ElementEnum::eOutput, std::shared_ptr<ElementOutputBase>> elementOutput;
    elementOutput[ElementEnum::eOutput::IP_DATA] = std::make_shared<ElementOutputIpData>();
    auto& elementIpDataMap = elementOutput[ElementEnum::eOutput::IP_DATA]->GetIpData().GetIpDataMap();

    for (auto component : visualizeComponents)
    {
        switch (component)
        {
        case NuTo::eVisualizeWhat::ENGINEERING_STRAIN:
        case NuTo::eVisualizeWhat::DAMAGE:
        case NuTo::eVisualizeWhat::SHRINKAGE_STRAIN:
            elementIpDataMap[ToIpDataEnum(component)];
            break;
        default:
            throw NuTo::Exception(__PRETTY_FUNCTION__,
                                           "Visualization component " + GetComponentName(component) +
                                                   " is not implemented or not known at the integration points.");
            break;
        }
    }


    // calculate the element solution
    Evaluate(elementOutput);

    // store data
    for (auto component : visualizeComponents)
    {
        auto componentName = GetComponentName(component);
        switch (component)
        {
        case NuTo::eVisualizeWhat::ENGINEERING_STRAIN:
        case NuTo::eVisualizeWhat::DAMAGE:
        case NuTo::eVisualizeWhat::SHRINKAGE_STRAIN:
        {
            const Eigen::MatrixXd& data = elementIpDataMap.at(ToIpDataEnum(component));
            assert(data.size() != 0);
            for (int i = 0; i < numIp; ++i)
                visualizer.SetCellData(ipInfo[i].cellId, componentName, data.col(i));
        }
        break;
        case NuTo::eVisualizeWhat::DISPLACEMENTS:
        case NuTo::eVisualizeWhat::TEMPERATURE:
        case NuTo::eVisualizeWhat::NONLOCAL_EQ_STRAIN:
        {
            for (int i = 0; i < numIp; ++i)
                visualizer.SetPointData(ipInfo[i].pointId, componentName,
                                        InterpolateDof3D(ipInfo[i].localCoords, Node::eDof::DISPLACEMENTS));
        }
        break;
        default:
            throw NuTo::Exception(__PRETTY_FUNCTION__,
                                           "Visualization component " + componentName +
                                                   " is not implemented or not known at the integration points.");
        }
    }
}
#endif // ENABLE_VISUALIZE


void NuTo::ElementBase::GetIntegratedStress(Eigen::MatrixXd& rStress)
{
    std::map<ElementEnum::eOutput, std::shared_ptr<ElementOutputBase>> elementOutput;
    elementOutput[ElementEnum::eOutput::IP_DATA] =
            std::make_shared<ElementOutputIpData>(IpData::eIpStaticDataType::ENGINEERING_STRESS);

    this->Evaluate(elementOutput);

    const auto& ipStress = elementOutput.at(ElementEnum::eOutput::IP_DATA)
                                   ->GetIpData()
                                   .GetIpDataMap()[IpData::eIpStaticDataType::ENGINEERING_STRESS];
    Eigen::VectorXd ipVolume = this->GetIntegrationPointVolume();

    rStress.resize(ipStress.rows(), 1);
    rStress.setZero();
    for (int countIP = 0; countIP < ipStress.cols(); countIP++)
    {
        rStress += (ipStress.col(countIP) * (ipVolume[countIP]));
    }
}

void NuTo::ElementBase::GetIntegratedStrain(Eigen::MatrixXd& rStrain)
{
    std::map<ElementEnum::eOutput, std::shared_ptr<ElementOutputBase>> elementOutput;
    elementOutput[ElementEnum::eOutput::IP_DATA] =
            std::make_shared<ElementOutputIpData>(IpData::eIpStaticDataType::ENGINEERING_STRAIN);

    this->Evaluate(elementOutput);

    const auto& ipStress = elementOutput.at(ElementEnum::eOutput::IP_DATA)
                                   ->GetIpData()
                                   .GetIpDataMap()[IpData::eIpStaticDataType::ENGINEERING_STRAIN];
    Eigen::VectorXd ipVolume = this->GetIntegrationPointVolume();

    rStrain.resize(ipStress.rows(), 1);
    rStrain.setZero();
    for (int countIP = 0; countIP < ipStress.cols(); countIP++)
    {
        rStrain += (ipStress.col(countIP) * (ipVolume[countIP]));
    }
}


void NuTo::ElementBase::Info() const
{
    std::cout << "[" << __PRETTY_FUNCTION__ << "] \n";
    std::cout << "InterpolationTypeInfo:\n" << GetInterpolationType().Info() << "\n";

    for (int iNode = 0; iNode < GetNumNodes(); ++iNode)
    {
        const NodeBase* node = GetNode(iNode);
        std::cout << "NodeInfo of local node " << iNode << ": \n";
        std::cout << *node << "\n";
    }
}

void NuTo::ElementBase::Info(std::ostream& out) const
{
    out << "InterpolationTypeInfo:\n" << GetInterpolationType().Info() << "\n";

    for (int iNode = 0; iNode < GetNumNodes(); ++iNode)
    {
        const NodeBase* node = GetNode(iNode);
        out << "NodeInfo of local node " << iNode << ": \n";
        out << *node << "\n";
    }
}

void NuTo::ElementBase::ReorderNodes()
{
    const Eigen::MatrixX2i& reorderIndices = mInterpolationType->GetNodeRenumberingIndices();
    for (int i = 0; i < reorderIndices.rows(); ++i)
    {
        int i0 = reorderIndices(i, 0);
        int i1 = reorderIndices(i, 1);
        // swap nodes i0 and i1
        NodeBase* tmpNode0 = GetNode(i0);
        NodeBase* tmpNode1 = GetNode(i1);
        SetNode(i0, tmpNode1);
        SetNode(i1, tmpNode0);
    }
}


void NuTo::ElementBase::AddPlaneStateToInput(ConstitutiveInputMap& constitutiveInput) const
{
    auto planeState = NuTo::Constitutive::eInput::PLANE_STATE;
    if (!GetSection()->IsPlaneStrain())
    {
        // plane stress is default
        constitutiveInput[planeState] = ConstitutiveIOBase::makeConstitutiveIO<2>(planeState);
    }
    else
    {
        constitutiveInput[planeState] = ConstitutiveIOBase::makeConstitutiveIO<2>(planeState);
        auto& value = *static_cast<ConstitutivePlaneState*>(constitutiveInput[planeState].get());
        value.SetPlaneState(NuTo::ePlaneState::PLANE_STRAIN);
    }
}
