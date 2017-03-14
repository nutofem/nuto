#include "mechanics/elements/ElementBase.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include <boost/foreach.hpp>
#include <boost/assign/ptr_map_inserter.hpp>

#include "mechanics/MechanicsException.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "mechanics/constitutive/inputoutput/ConstitutivePlaneState.h"
#include "mechanics/constraints/ConstraintBase.h"
#include "mechanics/elements/ElementOutputIpData.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/elements/IpDataEnum.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "mechanics/interpolationtypes/InterpolationBase.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/groups/GroupBase.h"
#include "mechanics/loads/LoadBase.h"
#include "mechanics/sections/SectionBase.h"
#include "mechanics/sections/SectionEnum.h"
#include "mechanics/structures/StructureBase.h"
#include "visualize/VisualizeEnum.h"


#include <eigen3/Eigen/QR>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Dense>

#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeUnstructuredGrid.h"
#include "visualize/VisualizeComponent.h"
#include "visualize/VisualizeException.h"
#endif

NuTo::ElementBase::ElementBase(const StructureBase* rStructure, const InterpolationType& rInterpolationType) :
        mStructure(rStructure),
        mInterpolationType(&rInterpolationType),
        mIPData(rInterpolationType.GetCurrentIntegrationType())
{}


#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ElementBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ElementBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ElementBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ElementBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ElementBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ElementBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ElementBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ElementBase " << std::endl;
#endif
    ar & boost::serialization::make_nvp ("mStructure", const_cast<StructureBase*&>(mStructure));
    ar & boost::serialization::make_nvp ("mInterpolationType", const_cast<InterpolationType*&>(mInterpolationType));
    ar & boost::serialization::make_nvp ("mElementData", mElementData);

    // the element data has to be saved on the main structure due to problems with a recursion on the stack (nonlocal data contains ptr to elements)
    // the idea is to first serialize all the elements in the table, and afterwards update the pointers of the element data in the element data routine
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ElementBase" << std::endl;
#endif
}
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ElementBase)
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ElementBase)
#endif // ENABLE_SERIALIZATION


void NuTo::ElementBase::Evaluate(std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rOutput)
{
    ConstitutiveInputMap input;
    input[Constitutive::eInput::CALCULATE_STATIC_DATA] = std::make_unique<ConstitutiveCalculateStaticData>(eCalculateStaticData::EULER_BACKWARD);

    return this->Evaluate(input, rOutput);
}

int NuTo::ElementBase::ElementGetId() const
{
    return mStructure->ElementGetId(this);
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

void NuTo::ElementBase::SetSection(const SectionBase& rSection)
{
    throw MechanicsException(__PRETTY_FUNCTION__, "This element type has so section.");
}

const NuTo::SectionBase& NuTo::ElementBase::GetSection() const
{
    throw MechanicsException(__PRETTY_FUNCTION__, "This element type has so section.");
}


int NuTo::ElementBase::GetNumNodes() const
{
    return mInterpolationType->GetNumNodes();
}

int NuTo::ElementBase::GetNumNodes(Node::eDof rDofType) const
{
    return mInterpolationType->Get(rDofType).GetNumNodes();
}

Eigen::VectorXd NuTo::ElementBase::ExtractNodeValues(int rTimeDerivative, Node::eDof rDofType) const
{
    throw NuTo::MechanicsException("[NuTo::ElementBase::ExtractNodeValues] not implemented.");
}


Eigen::VectorXd NuTo::ElementBase::InterpolateDofGlobal(const Eigen::VectorXd& rNaturalCoordinates, Node::eDof rDofType) const
{
    //assert(GetLocalDimension() == GetStructure()->GetDimension() && "Global and local dimensions do not agree.");

    return InterpolateDofGlobal(0, rNaturalCoordinates, rDofType);
}

Eigen::VectorXd NuTo::ElementBase::InterpolateDofGlobal(int rTimeDerivative, const Eigen::VectorXd& rNaturalCoordinates, Node::eDof rDofType) const
{

    const InterpolationBase& interpolationType = mInterpolationType->Get(rDofType);
    Eigen::MatrixXd nodalValues = ExtractNodeValues(rTimeDerivative, rDofType);
    Eigen::MatrixXd matrixN = interpolationType.CalculateMatrixN(rNaturalCoordinates);

    return matrixN * nodalValues;
}

Eigen::Vector3d NuTo::ElementBase::InterpolateDof3D(const Eigen::VectorXd& rNaturalCoordinates, Node::eDof rDofType) const
{
    return InterpolateDof3D(0, rNaturalCoordinates, rDofType);
}

Eigen::Vector3d NuTo::ElementBase::InterpolateDof3D(int rTimeDerivative, const Eigen::VectorXd& rNaturalCoordinates, Node::eDof rDofType) const
{

    Eigen::VectorXd interpolatedDofs = InterpolateDofGlobal(rTimeDerivative, rNaturalCoordinates, rDofType);
    Eigen::Vector3d interpolation3D = Eigen::Vector3d::Zero();

    interpolation3D.block(0,0, interpolatedDofs.rows(), 1) = interpolatedDofs;

    return interpolation3D;
}


void NuTo::ElementBase::SetIntegrationType(const NuTo::IntegrationTypeBase& rIntegrationType)
{
    //check compatibility between element type and constitutive law
    if (GetLocalDimension() ==  rIntegrationType.GetDimension())
    {
        mIPData.SetIntegrationType(rIntegrationType);
    } else
    {
        std::stringstream message;
        message << "[NuTo::ElementBase::SetIntegrationType] Integration Type does not match element type of element " << mStructure->ElementGetId(this) << "." << std::endl;
        throw MechanicsException(message.str());
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


template<int TDim>
void NuTo::ElementBase::EvaluateConstitutiveLaw(
        const NuTo::ConstitutiveInputMap& rConstitutiveInput,
        NuTo::ConstitutiveOutputMap& rConstitutiveOutput, unsigned int IP)
{
    Constitutive::IPConstitutiveLawBase& ipConstitutiveLaw = mIPData.GetIPConstitutiveLaw(IP);

    for (auto& itOutput : rConstitutiveOutput)
        if(itOutput.second!=nullptr) //check nullptr because of static data
            itOutput.second->SetIsCalculated(false);
    ipConstitutiveLaw.Evaluate<TDim>(rConstitutiveInput, rConstitutiveOutput);

    for(auto& itOutput : rConstitutiveOutput)
        if(itOutput.second!=nullptr && !itOutput.second->GetIsCalculated()) //check nullptr because of static data
            throw MechanicsException(__PRETTY_FUNCTION__,
                    "Output "+Constitutive::OutputToString(itOutput.first)+" not calculated by constitutive law");
}


template void NuTo::ElementBase::EvaluateConstitutiveLaw<1>(const NuTo::ConstitutiveInputMap&,
        NuTo::ConstitutiveOutputMap&, unsigned int);
template void NuTo::ElementBase::EvaluateConstitutiveLaw<2>(const NuTo::ConstitutiveInputMap&,
        NuTo::ConstitutiveOutputMap&, unsigned int);
template void NuTo::ElementBase::EvaluateConstitutiveLaw<3>(const NuTo::ConstitutiveInputMap&,
        NuTo::ConstitutiveOutputMap&, unsigned int);


const Eigen::Vector3d NuTo::ElementBase::GetGlobalIntegrationPointCoordinates(int rIpNum) const
{
    const Eigen::MatrixXd& matrixN = mInterpolationType->Get(Node::eDof::COORDINATES).GetMatrixN(rIpNum);
    Eigen::VectorXd nodeCoordinates = ExtractNodeValues(0, Node::eDof::COORDINATES);

    Eigen::Vector3d globalIntegrationPointCoordinates = Eigen::Vector3d::Zero();

    globalIntegrationPointCoordinates.segment(0, GetLocalDimension()) = matrixN * nodeCoordinates;

    return globalIntegrationPointCoordinates;
}

bool NuTo::ElementBase::GetLocalPointCoordinates(const double* rGlobCoords, double* rLocCoords) const
{
    throw NuTo::MechanicsException("[NuTo::ElementBase::GetLocalPointCoordinates] not implemented for this element type.");
}

NuTo::NodeBase *NuTo::ElementBase::GetBoundaryControlNode() const
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Not implemented for this element type.");
}


#ifdef ENABLE_VISUALIZE
void NuTo::ElementBase::GetVisualizationCells(unsigned int& NumVisualizationPoints, std::vector<double>& VisualizationPointLocalCoordinates, unsigned int& NumVisualizationCells, std::vector<NuTo::eCellTypes>& VisualizationCellType, std::vector<unsigned int>& VisualizationCellsIncidence,
        std::vector<unsigned int>& VisualizationCellsIP) const
{
    GetIntegrationType().GetVisualizationCells(NumVisualizationPoints, VisualizationPointLocalCoordinates, NumVisualizationCells, VisualizationCellType, VisualizationCellsIncidence, VisualizationCellsIP);
}

void NuTo::ElementBase::Visualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList)
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
    if (NumVisualizationPoints == 0) return; // nothing to visualize
    int dimension(VisualizationPointLocalCoordinates.size() / NumVisualizationPoints);
    assert(VisualizationPointLocalCoordinates.size() == NumVisualizationPoints * dimension);

    // TODO: fix that by proper visualization point (natural!) coordinates, local might be misleading here
    Eigen::MatrixXd visualizationPointNaturalCoordinates = Eigen::MatrixXd::Map(VisualizationPointLocalCoordinates.data(), dimension, NumVisualizationPoints);

    std::vector<unsigned int> PointIdVec;
    for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
    {
        if (dimension != 1 and dimension != 2 and dimension != 3)
            throw NuTo::MechanicsException("[NuTo::ElementBase::Visualize] invalid dimension of local coordinates");

        const Eigen::VectorXd& coords = visualizationPointNaturalCoordinates.col(PointCount);
        Eigen::Matrix<double, 3, 1> GlobalPointCoor = this->InterpolateDof3D(coords, Node::eDof::COORDINATES);
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
        case NuTo::eCellTypes::LINE:
        {
            assert(Pos + 2 <= VisualizationCellsIncidence.size());
            unsigned int Points[2];
            Points[0] = PointIdVec[VisualizationCellsIncidence[Pos]];
            Points[1] = PointIdVec[VisualizationCellsIncidence[Pos + 1]];
            unsigned int CellId = rVisualize.AddLineCell(Points);
            CellIdVec.push_back(CellId);
            Pos += 2;
        }
            break;
        case NuTo::eCellTypes::TRIANGLE:
        {
            assert(Pos + 3 <= VisualizationCellsIncidence.size());
            unsigned int Points[3];
            for (unsigned int PointCount = 0; PointCount < 3; PointCount++)
            {
                Points[PointCount] = PointIdVec[VisualizationCellsIncidence[Pos + PointCount]];
            }
            unsigned int CellId = rVisualize.AddTriangleCell(Points);
            CellIdVec.push_back(CellId);
            Pos += 3;
        }
            break;
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
        case NuTo::eCellTypes::TETRAEDER:
        {
            assert(Pos + 4 <= VisualizationCellsIncidence.size());
            unsigned int Points[4];
            for (unsigned int PointCount = 0; PointCount < 4; PointCount++)
            {
                Points[PointCount] = PointIdVec[VisualizationCellsIncidence[Pos + PointCount]];
            }
            unsigned int CellId = rVisualize.AddTetraCell(Points);
            CellIdVec.push_back(CellId);
            Pos += 4;
        }
            break;
        case NuTo::eCellTypes::HEXAHEDRON:
        {
            assert(Pos + 8 <= VisualizationCellsIncidence.size());
            unsigned int Points[8];
            for (unsigned int PointCount = 0; PointCount < 8; PointCount++)
            {
                Points[PointCount] = PointIdVec[VisualizationCellsIncidence[Pos + PointCount]];
            }
            unsigned int CellId = rVisualize.AddHexahedronCell(Points);
            CellIdVec.push_back(CellId);
            Pos += 8;
        }
            break;
        default:
            throw NuTo::MechanicsException("[NuTo::ElementBase::Visualize] unsupported visualization cell type");
        }
    }

    //determine the ipdata and determine the map
    std::map<NuTo::Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutput;
    elementOutput[Element::eOutput::IP_DATA] = std::make_shared<ElementOutputIpData>();

    auto& elementIpDataMap = elementOutput.at(Element::eOutput::IP_DATA)->GetIpData().GetIpDataMap();

    bool evaluateStress(false);


    for (auto const &it : rVisualizationList)
    {
        switch (it.get()->GetComponentEnum())
        {
        case NuTo::eVisualizeWhat::BOND_STRESS:
            elementIpDataMap[IpData::eIpStaticDataType::BOND_STRESS];
            break;
        case NuTo::eVisualizeWhat::DAMAGE:
            elementIpDataMap[IpData::eIpStaticDataType::DAMAGE];
            break;
        case NuTo::eVisualizeWhat::ENGINEERING_PLASTIC_STRAIN:
            elementIpDataMap[IpData::eIpStaticDataType::ENGINEERING_PLASTIC_STRAIN];
            break;
        case NuTo::eVisualizeWhat::ENGINEERING_STRAIN:
            elementIpDataMap[IpData::eIpStaticDataType::ENGINEERING_STRAIN];
            break;
        case NuTo::eVisualizeWhat::SHRINKAGE_STRAIN:
            elementIpDataMap[IpData::eIpStaticDataType::SHRINKAGE_STRAIN];
            break;
        case NuTo::eVisualizeWhat::THERMAL_STRAIN:
            elementIpDataMap[IpData::eIpStaticDataType::THERMAL_STRAIN];
            break;
        case NuTo::eVisualizeWhat::ENGINEERING_STRESS:
            if (evaluateStress == false)
            {
                elementIpDataMap[IpData::eIpStaticDataType::ENGINEERING_STRESS];
                evaluateStress = true;
            }
            break;
        case NuTo::eVisualizeWhat::HEAT_FLUX:
            elementIpDataMap[IpData::eIpStaticDataType::HEAT_FLUX];
            break;
        case NuTo::eVisualizeWhat::LOCAL_EQ_STRAIN:
            elementIpDataMap[IpData::eIpStaticDataType::LOCAL_EQ_STRAIN];
            break;
        case NuTo::eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS:
            if (evaluateStress == false)
            {
                elementIpDataMap[IpData::eIpStaticDataType::ENGINEERING_STRESS];
                evaluateStress = true;
            }
            break;
        case NuTo::eVisualizeWhat::TOTAL_INELASTIC_EQ_STRAIN:
            elementIpDataMap[IpData::eIpStaticDataType::TOTAL_INELASTIC_EQ_STRAIN];
            break;

        case NuTo::eVisualizeWhat::ACCELERATION:
        case NuTo::eVisualizeWhat::ANGULAR_VELOCITY:
        case NuTo::eVisualizeWhat::ANGULAR_ACCELERATION:
        case NuTo::eVisualizeWhat::CONSTITUTIVE:
        case NuTo::eVisualizeWhat::DISPLACEMENTS:
        case NuTo::eVisualizeWhat::ELEMENT:
        case NuTo::eVisualizeWhat::NONLOCAL_EQ_STRAIN:
        case NuTo::eVisualizeWhat::PARTICLE_RADIUS:
        case NuTo::eVisualizeWhat::RELATIVE_HUMIDITY:
        case NuTo::eVisualizeWhat::ROTATION:
        case NuTo::eVisualizeWhat::SECTION:
        case NuTo::eVisualizeWhat::TEMPERATURE:
        case NuTo::eVisualizeWhat::VELOCITY:
        case NuTo::eVisualizeWhat::WATER_VOLUME_FRACTION:
        default:
            //do nothing
            ;
            break;
        }
    }

    //calculate the element solution
    ConstitutiveInputMap input;
    input[Constitutive::eInput::CALCULATE_STATIC_DATA] = std::make_unique<ConstitutiveCalculateStaticData>(
            eCalculateStaticData::USE_PREVIOUS);
    Evaluate(input, elementOutput);
//    Evaluate(elementOutput);

    //assign the outputs

    // store data
    for (auto const &it : rVisualizationList)
    {
        switch (it.get()->GetComponentEnum())
        {
        case NuTo::eVisualizeWhat::DAMAGE:
        {
            const auto& damage = elementIpDataMap.at(IpData::eIpStaticDataType::DAMAGE);
            assert(damage.size() != 0);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataScalar(CellId, it.get()->GetComponentName(), damage.data()[theIp]);
            }
        }
        break;
        case NuTo::eVisualizeWhat::LOCAL_EQ_STRAIN:
        {
            const auto& localEqStrain = elementIpDataMap.at(IpData::eIpStaticDataType::LOCAL_EQ_STRAIN);
            assert(localEqStrain.size() != 0);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataScalar(CellId, it.get()->GetComponentName(), localEqStrain.data()[theIp]);
            }
        }
        break;
        case NuTo::eVisualizeWhat::DISPLACEMENTS:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                const Eigen::VectorXd& coords = visualizationPointNaturalCoordinates.col(PointCount);
                Eigen::Vector3d GlobalDisplacements = this->InterpolateDof3D(coords, Node::eDof::DISPLACEMENTS);
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(), GlobalDisplacements.data());
            }
            break;
        case NuTo::eVisualizeWhat::VELOCITY:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                const Eigen::VectorXd& coords = visualizationPointNaturalCoordinates.col(PointCount);
                Eigen::Vector3d GlobalVelocity = this->InterpolateDof3D(1, coords, Node::eDof::DISPLACEMENTS);
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(), GlobalVelocity.data());
            }
            break;
        case NuTo::eVisualizeWhat::ACCELERATION:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                const Eigen::VectorXd& coords = visualizationPointNaturalCoordinates.col(PointCount);
                Eigen::Vector3d GlobalAcceleration = this->InterpolateDof3D(2, coords, Node::eDof::DISPLACEMENTS);
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(), GlobalAcceleration.data());
            }
            break;
        case NuTo::eVisualizeWhat::ENGINEERING_STRAIN:
        {
            const auto& engineeringStrain = elementIpDataMap.at(IpData::eIpStaticDataType::ENGINEERING_STRAIN);
            assert(engineeringStrain.size() != 0);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                 double EngineeringStrainTensor[9];
                EngineeringStrainTensor[0] =       engineeringStrain(0, theIp);
                EngineeringStrainTensor[1] = 0.5 * engineeringStrain(3, theIp);
                EngineeringStrainTensor[2] = 0.5 * engineeringStrain(5, theIp);
                EngineeringStrainTensor[3] = 0.5 * engineeringStrain(3, theIp);
                EngineeringStrainTensor[4] =       engineeringStrain(1, theIp);
                EngineeringStrainTensor[5] = 0.5 * engineeringStrain(4, theIp);
                EngineeringStrainTensor[6] = 0.5 * engineeringStrain(5, theIp);
                EngineeringStrainTensor[7] = 0.5 * engineeringStrain(4, theIp);
                EngineeringStrainTensor[8] =       engineeringStrain(2, theIp);

                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataTensor(CellId, it.get()->GetComponentName(), EngineeringStrainTensor);
            }
        }
            break;
        case NuTo::eVisualizeWhat::SHRINKAGE_STRAIN:
        {
            const auto& shrinkageStrain = elementIpDataMap.at(IpData::eIpStaticDataType::SHRINKAGE_STRAIN);
            assert(shrinkageStrain.size() != 0);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                 double shrinkageStrainTensor[9];
                shrinkageStrainTensor[0] =       shrinkageStrain(0, theIp);
                shrinkageStrainTensor[1] = 0.5 * shrinkageStrain(3, theIp);
                shrinkageStrainTensor[2] = 0.5 * shrinkageStrain(5, theIp);
                shrinkageStrainTensor[3] = 0.5 * shrinkageStrain(3, theIp);
                shrinkageStrainTensor[4] =       shrinkageStrain(1, theIp);
                shrinkageStrainTensor[5] = 0.5 * shrinkageStrain(4, theIp);
                shrinkageStrainTensor[6] = 0.5 * shrinkageStrain(5, theIp);
                shrinkageStrainTensor[7] = 0.5 * shrinkageStrain(4, theIp);
                shrinkageStrainTensor[8] =       shrinkageStrain(2, theIp);

                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataTensor(CellId, it.get()->GetComponentName(), shrinkageStrainTensor);
            }
        }
            break;
        case NuTo::eVisualizeWhat::THERMAL_STRAIN:
        {
            const auto& thermalStrain = elementIpDataMap.at(IpData::eIpStaticDataType::THERMAL_STRAIN);
            assert(thermalStrain.size() != 0);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                double thermalStrainTensor[9];
                thermalStrainTensor[0] =       thermalStrain(0, theIp);
                thermalStrainTensor[1] = 0.5 * thermalStrain(3, theIp);
                thermalStrainTensor[2] = 0.5 * thermalStrain(5, theIp);
                thermalStrainTensor[3] = 0.5 * thermalStrain(3, theIp);
                thermalStrainTensor[4] =       thermalStrain(1, theIp);
                thermalStrainTensor[5] = 0.5 * thermalStrain(4, theIp);
                thermalStrainTensor[6] = 0.5 * thermalStrain(5, theIp);
                thermalStrainTensor[7] = 0.5 * thermalStrain(4, theIp);
                thermalStrainTensor[8] =       thermalStrain(2, theIp);

                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataTensor(CellId, it.get()->GetComponentName(), thermalStrainTensor);
            }
        }
            break;
        case NuTo::eVisualizeWhat::ENGINEERING_PLASTIC_STRAIN:
        {
            const auto& engineeringPlasticStrain = elementIpDataMap.at(IpData::eIpStaticDataType::ENGINEERING_PLASTIC_STRAIN);
            assert(engineeringPlasticStrain.size() != 0);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                double EngineeringStrainTensor[9];
                EngineeringStrainTensor[0] =       engineeringPlasticStrain(0, theIp);
                EngineeringStrainTensor[1] = 0.5 * engineeringPlasticStrain(3, theIp);
                EngineeringStrainTensor[2] = 0.5 * engineeringPlasticStrain(5, theIp);
                EngineeringStrainTensor[3] = 0.5 * engineeringPlasticStrain(3, theIp);
                EngineeringStrainTensor[4] =       engineeringPlasticStrain(1, theIp);
                EngineeringStrainTensor[5] = 0.5 * engineeringPlasticStrain(4, theIp);
                EngineeringStrainTensor[6] = 0.5 * engineeringPlasticStrain(5, theIp);
                EngineeringStrainTensor[7] = 0.5 * engineeringPlasticStrain(4, theIp);
                EngineeringStrainTensor[8] =       engineeringPlasticStrain(2, theIp);

                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataTensor(CellId, it.get()->GetComponentName(), EngineeringStrainTensor);
            }
        }
            break;
        case NuTo::eVisualizeWhat::TOTAL_INELASTIC_EQ_STRAIN:
        {
            const auto& totalInelasticEqStrain = elementIpDataMap.at(IpData::eIpStaticDataType::TOTAL_INELASTIC_EQ_STRAIN);
            assert(totalInelasticEqStrain.size() != 0);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataScalar(CellId, it.get()->GetComponentName(), totalInelasticEqStrain.data()[theIp]);
            }
        }
            break;
        case NuTo::eVisualizeWhat::ENGINEERING_STRESS:
        {
            const auto& engineeringStress = elementIpDataMap.at(IpData::eIpStaticDataType::ENGINEERING_STRESS);
            assert(engineeringStress.size() != 0);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                double EngineeringStressTensor[9];
                EngineeringStressTensor[0] = engineeringStress(0, theIp);
                EngineeringStressTensor[1] = engineeringStress(3, theIp);
                EngineeringStressTensor[2] = engineeringStress(5, theIp);
                EngineeringStressTensor[3] = engineeringStress(3, theIp);
                EngineeringStressTensor[4] = engineeringStress(1, theIp);
                EngineeringStressTensor[5] = engineeringStress(4, theIp);
                EngineeringStressTensor[6] = engineeringStress(5, theIp);
                EngineeringStressTensor[7] = engineeringStress(4, theIp);
                EngineeringStressTensor[8] = engineeringStress(2, theIp);

                //std::cout<<"[NuTo::ElementBase::VisualizeEngineeringStressTensor]" << EngineeringStressTensor[0] << EngineeringStressTensor[1] << std::endl;
                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataTensor(CellId, it.get()->GetComponentName(), EngineeringStressTensor);
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
        case NuTo::eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS:
        {
            const auto& engineeringStress = elementIpDataMap.at(IpData::eIpStaticDataType::ENGINEERING_STRESS);
            assert(engineeringStress.size() != 0);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                Eigen::Matrix<double, 3, 3> EngineeringStressTensor;

                EngineeringStressTensor(0, 0) = engineeringStress(0, theIp);
                EngineeringStressTensor(1, 0) = engineeringStress(3, theIp);
                EngineeringStressTensor(2, 0) = engineeringStress(5, theIp);
                EngineeringStressTensor(0, 1) = engineeringStress(3, theIp);
                EngineeringStressTensor(1, 1) = engineeringStress(1, theIp);
                EngineeringStressTensor(2, 1) = engineeringStress(4, theIp);
                EngineeringStressTensor(0, 2) = engineeringStress(5, theIp);
                EngineeringStressTensor(1, 2) = engineeringStress(4, theIp);
                EngineeringStressTensor(2, 2) = engineeringStress(2, theIp);

                //std::cout<<"[NuTo::ElementBase::VisualizeEngineeringStressTensor]" << EngineeringStressTensor[0] << EngineeringStressTensor[1] << std::endl;
                unsigned int CellId = CellIdVec[CellCount];
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> mySolver(EngineeringStressTensor, false);
                Eigen::Matrix<double, 3, 1> eigenValues(mySolver.eigenvalues());
                rVisualize.SetCellDataVector(CellId, it.get()->GetComponentName(), eigenValues.data());
            }
        }
            break;
        case NuTo::eVisualizeWhat::NONLOCAL_EQ_STRAIN:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                Eigen::VectorXd nonlocalEqStrain(1);
                if (mInterpolationType->IsDof(Node::eDof::NONLOCALEQSTRAIN))
                {
                    // calculate only if element has nonlocal eq strain dofs
                    const Eigen::VectorXd& coords = visualizationPointNaturalCoordinates.col(PointCount);
                    nonlocalEqStrain = InterpolateDofGlobal(coords, Node::eDof::NONLOCALEQSTRAIN);
                    assert(nonlocalEqStrain.rows() == 1);
                }
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataScalar(PointId, it.get()->GetComponentName(), nonlocalEqStrain[0]);
            }
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
        case NuTo::eVisualizeWhat::PARTICLE_RADIUS:
            //do nothing
            break;
        case NuTo::eVisualizeWhat::ROTATION:
        case NuTo::eVisualizeWhat::ANGULAR_VELOCITY:
        case NuTo::eVisualizeWhat::ANGULAR_ACCELERATION:
            //do nothing
            break;
        case NuTo::eVisualizeWhat::HEAT_FLUX:
        {
            const auto& heatFlux = elementIpDataMap.at(IpData::eIpStaticDataType::HEAT_FLUX);
            assert(heatFlux.size() != 0);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                unsigned int CellId = CellIdVec[CellCount];
                double heatFlux_VISUALIZE_I_HATE_YOU_SO_MUCH_OMG[3];
                heatFlux_VISUALIZE_I_HATE_YOU_SO_MUCH_OMG[0] = heatFlux(0, theIp);
                heatFlux_VISUALIZE_I_HATE_YOU_SO_MUCH_OMG[1] = heatFlux(1, theIp);
                heatFlux_VISUALIZE_I_HATE_YOU_SO_MUCH_OMG[2] = heatFlux(2, theIp);

                rVisualize.SetCellDataVector(CellId, it.get()->GetComponentName(),heatFlux_VISUALIZE_I_HATE_YOU_SO_MUCH_OMG);
            }
        }
            break;
        case NuTo::eVisualizeWhat::TEMPERATURE:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                const Eigen::VectorXd& coords = visualizationPointNaturalCoordinates.col(PointCount);
                Eigen::VectorXd temperature = InterpolateDofGlobal(coords, Node::eDof::TEMPERATURE);
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataScalar(PointId, it.get()->GetComponentName(), temperature[0]);
            }
            break;

        case NuTo::eVisualizeWhat::CRACK_PHASE_FIELD:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                const Eigen::VectorXd& coords = visualizationPointNaturalCoordinates.col(PointCount);
                Eigen::VectorXd damage = InterpolateDofGlobal(coords, Node::eDof::CRACKPHASEFIELD);
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataScalar(PointId, it.get()->GetComponentName(), damage[0]);
            }
            break;

        case NuTo::eVisualizeWhat::RELATIVE_HUMIDITY:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                Eigen::Matrix<double, Eigen::Dynamic, 1> relativeHumidity;
                if (mInterpolationType->IsDof(Node::eDof::RELATIVEHUMIDITY))
                {
                    // calculate only if element has nonlocal eq strain dofs
                    const Eigen::VectorXd& coords = visualizationPointNaturalCoordinates.col(PointCount);
                    relativeHumidity = InterpolateDofGlobal(coords, Node::eDof::RELATIVEHUMIDITY);
                    assert(relativeHumidity.rows() == 1);
                }
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataScalar(PointId, it.get()->GetComponentName(), relativeHumidity[0]);
            }
            break;
        case NuTo::eVisualizeWhat::WATER_VOLUME_FRACTION:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                Eigen::Matrix<double, Eigen::Dynamic, 1> waterVolumeFraction;
                if (mInterpolationType->IsDof(Node::eDof::WATERVOLUMEFRACTION))
                {
                    // calculate only if element has nonlocal eq strain dofs
                    const Eigen::VectorXd& coords = visualizationPointNaturalCoordinates.col(PointCount);
                    waterVolumeFraction = InterpolateDofGlobal(coords, Node::eDof::WATERVOLUMEFRACTION);
                    assert(waterVolumeFraction.rows() == 1);
                }
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataScalar(PointId, it.get()->GetComponentName(), waterVolumeFraction[0]);
            }
            break;
        default:
            //VHIRTHAMTODO: Create enum to string function and replace static cast!
            std::cout << static_cast<int>(it.get()->GetComponentEnum()) << "\n";
            throw NuTo::MechanicsException("[NuTo::ElementBase::Visualize] unsupported datatype for visualization.");
        }
    }
}

void NuTo::ElementBase::VisualizeExtrapolateToNodes(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList)
{
    throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) +": \t This function is not ready to be used yet. Choose a different visualization type!");
    /*
    // get visualization cells from integration type
    unsigned int NumVisualizationPoints;
    std::vector<double> VisualizationPointLocalCoordinates;
    unsigned int NumVisualizationCells;
    std::vector<NuTo::eCellTypes> VisualizationCellType;
    std::vector<unsigned int> VisualizationCellsIncidence;
    std::vector<unsigned int> VisualizationCellsIP;
    int dimension = GetLocalDimension();

    NumVisualizationPoints = 3;

    for (unsigned int iNode = 0; iNode < NumVisualizationPoints; ++iNode)
    {
        for (int iDim = 0; iDim < dimension; ++iDim)
        {
            VisualizationPointLocalCoordinates.push_back(GetInterpolationType()->GetNaturalNodeCoordinates(iNode).at(iDim,0));
        }

    }

    NumVisualizationCells = 1;

    // cell 0
    VisualizationCellType.push_back(NuTo::eCellTypes::TRIANGLE);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIP.push_back(0);


    // TODO: fix that by proper visualization point (natural!) coordinates, local might be misleading here
    Eigen::MatrixXd visualizationPointNaturalCoordinates = Eigen::MatrixXd::Map(VisualizationPointLocalCoordinates.data(), dimension, NumVisualizationPoints);

    std::vector<unsigned int> PointIdVec;
    for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
    {

        const Eigen::VectorXd& coords = visualizationPointNaturalCoordinates.col(PointCount);
        Eigen::Matrix<double, 3, 1> GlobalPointCoor = this->InterpolateDof3D(coords, Node::eDof::COORDINATES);
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
        case NuTo::eCellTypes::TRIANGLE:
        {
            assert(Pos + 3 <= VisualizationCellsIncidence.size());
            unsigned int Points[3];
            for (unsigned int PointCount = 0; PointCount < 3; PointCount++)
            {
                Points[PointCount] = PointIdVec[VisualizationCellsIncidence[Pos + PointCount]];
            }
            unsigned int CellId = rVisualize.AddTriangleCell(Points);
            CellIdVec.push_back(CellId);
            Pos += 3;
        }
            break;
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
            throw NuTo::MechanicsException("[NuTo::ElementBase::Visualize] unsupported visualization cell type");
        }
    }

    //determine the ipdata and determine the map
    boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;


    for (auto const &it : rVisualizationList)
    {
        switch (it.get()->GetComponentEnum())
        {
        case NuTo::eVisualizeWhat::ENGINEERING_STRAIN:
            boost::assign::ptr_map_insert<ElementOutputIpData>( elementOutput )( Element::IP_DATA ,IpData::ENGINEERING_STRAIN);
        break;
        case NuTo::eVisualizeWhat::DISPLACEMENTS:
        case NuTo::eVisualizeWhat::SECTION:
        case NuTo::eVisualizeWhat::ELEMENT:
        default:
            //do nothing
            break;
        }
    }

    //calculate the element solution
    Evaluate(elementOutput);

    //assign the outputs
    Eigen::MatrixXd* engineeringStrain(nullptr);

    for (auto itElementOutput=elementOutput.begin(); itElementOutput!=elementOutput.end(); itElementOutput++)
    {
        switch (itElementOutput->second->GetIpDataType())
        {
        case NuTo::IpData::ENGINEERING_STRAIN:
            engineeringStrain = &(itElementOutput->second->GetFullMatrixDouble());
        break;
        default:
            throw MechanicsException("[NuTo::ElementBase::Visualize] other ipdatatypes not supported.");
        }
    }

    // store data
    for (auto const &it : rVisualizationList)
    {
        switch (it.get()->GetComponentEnum())
        {
        case NuTo::eVisualizeWhat::DISPLACEMENTS:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                const Eigen::VectorXd& coords = visualizationPointNaturalCoordinates.col(PointCount);
                Eigen::Vector3d GlobalDisplacements = this->InterpolateDof3D(coords, Node::eDof::DISPLACEMENTS);
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(), GlobalDisplacements.data());
            }
            break;
        case NuTo::eVisualizeWhat::ENGINEERING_STRAIN:
        {

            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                const Eigen::VectorXd& coords = visualizationPointNaturalCoordinates.col(PointCount);
                const double* EngineeringStrainVector = &(engineeringStrain->data()[PointCount * 6]);
                Eigen::Vector3d normalStrains;

                normalStrains(0,0) = EngineeringStrainVector[0];
                normalStrains(1,0) = EngineeringStrainVector[1];
                normalStrains(2,0) = EngineeringStrainVector[2];

                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(), normalStrains.data());
            }

        }
            break;
            case NuTo::eVisualizeWhat::SECTION:
        {
            int sectionId = mStructure->SectionGetId(GetSection());
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
            default:
            std::cout << it.get()->GetComponentEnum() << "\n";
            throw NuTo::MechanicsException("[NuTo::ElementBase::Visualize] unsupported datatype for visualization.");
        }
    }
    */
}


void NuTo::ElementBase::VisualizeIntegrationPointData(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList)
{
    //
    //  This function is still in beta and only works for engineering strain. Implementation is still in progress...
    //

    // get visualization cells from integration type
    unsigned int NumVisualizationPoints = GetInterpolationType().GetCurrentIntegrationType().GetNumIntegrationPoints();
    unsigned int NumVisualizationCells = NumVisualizationPoints;
    std::vector<NuTo::eCellTypes> VisualizationCellType;


    std::vector<double> VisualizationPointLocalCoordinates;

    std::vector<unsigned int> VisualizationCellsIncidence;
    std::vector<unsigned int> VisualizationCellsIP;
    int dimension = GetLocalDimension();



    for (unsigned int iIp = 0; iIp < NumVisualizationPoints; ++iIp)
    {
        auto naturalIpCoords = GetInterpolationType().GetCurrentIntegrationType().GetLocalIntegrationPointCoordinates(iIp);

        for (int iDim = 0; iDim < dimension; ++iDim)
        {
            VisualizationPointLocalCoordinates.push_back(naturalIpCoords[iDim]);
        }

        VisualizationCellType.push_back(NuTo::eCellTypes::VERTEX);
        VisualizationCellsIncidence.push_back(iIp);
        VisualizationCellsIP.push_back(iIp);
    }


    // TODO: fix that by proper visualization point (natural!) coordinates, local might be misleading here
        Eigen::MatrixXd visualizationPointNaturalCoordinates = Eigen::MatrixXd::Map(VisualizationPointLocalCoordinates.data(), dimension, NumVisualizationPoints);

        std::vector<unsigned int> PointIdVec;
        for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
        {

            const Eigen::VectorXd& coords = visualizationPointNaturalCoordinates.col(PointCount);
            Eigen::Matrix<double, 3, 1> GlobalPointCoor = this->InterpolateDof3D(coords, Node::eDof::COORDINATES);
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
            case NuTo::eCellTypes::VERTEX:
            {
                assert(Pos + 1 <= VisualizationCellsIncidence.size());
                unsigned int Points[1];
                for (unsigned int PointCount = 0; PointCount < 1; PointCount++)
                {
                    Points[PointCount] = PointIdVec[VisualizationCellsIncidence[Pos + PointCount]];
                }
                unsigned int CellId = rVisualize.AddVertexCell(Points);
                CellIdVec.push_back(CellId);
                Pos += 1;
            }
                break;
            default:
                throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) + ": \t unsupported visualization cell type");
            }
        }

        //determine the ipdata and determine the map
        std::map<NuTo::Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutput;
        elementOutput[Element::eOutput::IP_DATA] = std::make_shared<ElementOutputIpData>();
        auto& elementIpDataMap = elementOutput[Element::eOutput::IP_DATA]->GetIpData().GetIpDataMap();

        for (auto const &it : rVisualizationList)
        {
            switch (it.get()->GetComponentEnum())
            {
                case NuTo::eVisualizeWhat::ENGINEERING_STRAIN:
                    elementIpDataMap[IpData::eIpStaticDataType::ENGINEERING_STRAIN];
                    break;
                case NuTo::eVisualizeWhat::DAMAGE:
                    elementIpDataMap[IpData::eIpStaticDataType::DAMAGE];
                    break;
                case NuTo::eVisualizeWhat::SHRINKAGE_STRAIN:
                    elementIpDataMap[IpData::eIpStaticDataType::SHRINKAGE_STRAIN];
                    break;
                default:
                    throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) + ": \t Visualization component " + it.get()->GetComponentName() + " is not implemented or not known at the integration points.");
                    break;
            }
        }

        //calculate the element solution
        Evaluate(elementOutput);

    // store data
    for (auto const &it : rVisualizationList)
    {
        switch (it.get()->GetComponentEnum())
        {
            case NuTo::eVisualizeWhat::ENGINEERING_STRAIN:
            {
                const auto& engineeringStrain = elementIpDataMap.at(IpData::eIpStaticDataType::ENGINEERING_STRAIN);
                assert(engineeringStrain.size() != 0);
                for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
                {
                    unsigned int theIp = VisualizationCellsIP[CellCount];
                    double EngineeringStrainTensor[9];
                    EngineeringStrainTensor[0] =       engineeringStrain(0, theIp);
                    EngineeringStrainTensor[1] = 0.5 * engineeringStrain(3, theIp);
                    EngineeringStrainTensor[2] = 0.5 * engineeringStrain(5, theIp);
                    EngineeringStrainTensor[3] = 0.5 * engineeringStrain(3, theIp);
                    EngineeringStrainTensor[4] =       engineeringStrain(1, theIp);
                    EngineeringStrainTensor[5] = 0.5 * engineeringStrain(4, theIp);
                    EngineeringStrainTensor[6] = 0.5 * engineeringStrain(5, theIp);
                    EngineeringStrainTensor[7] = 0.5 * engineeringStrain(4, theIp);
                    EngineeringStrainTensor[8] =       engineeringStrain(2, theIp);

                    unsigned int CellId = CellIdVec[CellCount];
                    rVisualize.SetCellDataTensor(CellId, it.get()->GetComponentName(), EngineeringStrainTensor);
                }
            }
                break;

            case NuTo::eVisualizeWhat::DAMAGE:
            {
                const auto& damage = elementIpDataMap.at(IpData::eIpStaticDataType::DAMAGE);
                for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
                {
                    unsigned int theIp = VisualizationCellsIP[CellCount];
                    double Damage =       damage(0, theIp);
                    unsigned int CellId = CellIdVec[CellCount];

                    rVisualize.SetCellDataScalar(CellId, it.get()->GetComponentName(), Damage);
                }
            }
                break;
            case NuTo::eVisualizeWhat::SHRINKAGE_STRAIN:
            {
                const auto& shrinkageStrain = elementIpDataMap.at(IpData::eIpStaticDataType::SHRINKAGE_STRAIN);
                assert(shrinkageStrain.size() != 0);
                for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
                {
                    unsigned int theIp = VisualizationCellsIP[CellCount];
                    double shrinkageStrainTensor[9];
                    shrinkageStrainTensor[0] =       shrinkageStrain(0, theIp);
                    shrinkageStrainTensor[1] = 0.5 * shrinkageStrain(3, theIp);
                    shrinkageStrainTensor[2] = 0.5 * shrinkageStrain(5, theIp);
                    shrinkageStrainTensor[3] = 0.5 * shrinkageStrain(3, theIp);
                    shrinkageStrainTensor[4] =       shrinkageStrain(1, theIp);
                    shrinkageStrainTensor[5] = 0.5 * shrinkageStrain(4, theIp);
                    shrinkageStrainTensor[6] = 0.5 * shrinkageStrain(5, theIp);
                    shrinkageStrainTensor[7] = 0.5 * shrinkageStrain(4, theIp);
                    shrinkageStrainTensor[8] =       shrinkageStrain(2, theIp);

                    unsigned int CellId = CellIdVec[CellCount];
                    rVisualize.SetCellDataTensor(CellId, it.get()->GetComponentName(), shrinkageStrainTensor);
                }
            }
                break;
            default:
                throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) + ": \t Visualization component " + it.get()->GetComponentName() + " is not implemented or not known at the integration points.");
        }
    }

}
#endif // ENABLE_VISUALIZE


void NuTo::ElementBase::GetIntegratedStress(Eigen::MatrixXd& rStress)
{
    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutput;
    elementOutput[Element::eOutput::IP_DATA] = std::make_shared<ElementOutputIpData>(IpData::eIpStaticDataType::ENGINEERING_STRESS);

    this->Evaluate(elementOutput);

    const auto& ipStress = elementOutput.at(Element::eOutput::IP_DATA)->GetIpData().GetIpDataMap()[IpData::eIpStaticDataType::ENGINEERING_STRESS];
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
    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutput;
    elementOutput[Element::eOutput::IP_DATA] = std::make_shared<ElementOutputIpData>(IpData::eIpStaticDataType::ENGINEERING_STRAIN);

    this->Evaluate(elementOutput);

    const auto& ipStress = elementOutput.at(Element::eOutput::IP_DATA)->GetIpData().GetIpDataMap()[IpData::eIpStaticDataType::ENGINEERING_STRAIN];
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
    mStructure->GetLogger() << "[" << __PRETTY_FUNCTION__ << "] \n";
    mStructure->GetLogger() << "InterpolationTypeInfo:\n" << GetInterpolationType().Info() << "\n";
    mStructure->GetLogger() << Element::ElementTypeToString(GetEnumType()) << "\n";

    for (int iNode = 0; iNode < GetNumNodes(); ++iNode)
    {
        const NodeBase* node = GetNode(iNode);
        mStructure->GetLogger() << "NodeInfo of local node " << iNode << ": \n";
        mStructure->GetLogger() << node->GetNodeTypeStr() << "\n";
    }
//    mStructure->GetLogger() << "InterpolationTypeInfo: \n " << GetInterpolationType()->Info();

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

        if (mStructure->GetVerboseLevel() > 5)
            mStructure->GetLogger() << "[NuTo::ElementBase::ReorderNodes] Swapped nodes " << i0 << " and " << i1 << ".\n";
    }
    if (mStructure->GetVerboseLevel() > 5)
    	mStructure->GetLogger() << "\n";
}


void NuTo::ElementBase::AddPlaneStateToInput(ConstitutiveInputMap& constitutiveInput) const
{
    auto planeState = NuTo::Constitutive::eInput::PLANE_STATE;
    if (GetSection().GetType() == NuTo::eSectionType::PLANE_STRESS)
    {
        // plane stress is default
        constitutiveInput[planeState] = ConstitutiveIOBase::makeConstitutiveIO<2>(planeState);
    }
    if (GetSection().GetType() == NuTo::eSectionType::PLANE_STRAIN)
    {
        constitutiveInput[planeState] = ConstitutiveIOBase::makeConstitutiveIO<2>(planeState);
        auto& value = *static_cast<ConstitutivePlaneState*>(constitutiveInput[planeState].get());
        value.SetPlaneState(NuTo::ePlaneState::PLANE_STRAIN);
    }
}
