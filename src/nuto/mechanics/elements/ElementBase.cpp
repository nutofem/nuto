// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include <boost/assign/ptr_map_inserter.hpp>

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/constraints/ConstraintBase.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/ElementDataConstitutiveIp.h"
#include "nuto/mechanics/elements/ElementDataVariableConstitutiveIp.h"
#include "nuto/mechanics/elements/ElementDataConstitutiveIpCrack.h"
#include "nuto/mechanics/elements/ElementDataConstitutiveIpNonlocal.h"
#include "nuto/mechanics/elements/ElementOutputIpData.h"
#include "nuto/mechanics/groups/GroupBase.h"
#include "nuto/mechanics/loads/LoadBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/structures/StructureBase.h"



#include <eigen3/Eigen/QR>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Dense>

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeException.h"
#include "nuto/visualize/VisualizeComponent.h"
#endif

NuTo::ElementBase::ElementBase(const StructureBase* rStructure, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType, const InterpolationType* rInterpolationType) :
        mStructure(rStructure), mInterpolationType(rInterpolationType)
{

    const IntegrationTypeBase* integrationType = mInterpolationType->GetCurrentIntegrationType();

    //allocate element data
    ElementDataBase* ptrElementData;
    switch (rElementDataType)
    {
    case NuTo::ElementData::NOELEMENTDATA:
        throw MechanicsException("[NuTo::ElementBase::ElementBase] no elements without element data implemented.");
        break;
    case NuTo::ElementData::CONSTITUTIVELAWIP:
        ptrElementData = new NuTo::ElementDataConstitutiveIp(this, integrationType, rIpDataType);
        break;
    case NuTo::ElementData::CONSTITUTIVELAWIPNONLOCAL:
        ptrElementData = new NuTo::ElementDataConstitutiveIpNonlocal(this, integrationType, rIpDataType);
        break;
    case NuTo::ElementData::CONSTITUTIVELAWIPCRACK:
        ptrElementData = new NuTo::ElementDataConstitutiveIpCrack(this, integrationType, rIpDataType);
        break;
    case NuTo::ElementData::VARIABLECONSTITUTIVELAWIP:
        ptrElementData = new NuTo::ElementDataVariableConstitutiveIp(this, integrationType, rIpDataType);
        break;

    default:
        throw MechanicsException("[NuTo::ElementBase::ElementWithDataBase] unsupported element data type.");
    }
    mElementData = ptrElementData;
}

NuTo::ElementBase::ElementBase(const StructureBase* rStructure, ElementData::eElementDataType rElementDataType, int rNumIp, IpData::eIpDataType rIpDataType, const InterpolationType* rInterpolationType) :
        mStructure(rStructure), mInterpolationType(rInterpolationType)
{

    //allocate element data
    ElementDataBase* ptrElementData;
    switch (rElementDataType)
    {
    case NuTo::ElementData::NOELEMENTDATA:
        throw MechanicsException("[NuTo::ElementBase::ElementBase] no elements without element data implemented.");
        break;
    case NuTo::ElementData::CONSTITUTIVELAWIP:
        ptrElementData = new NuTo::ElementDataConstitutiveIp(this, rNumIp, rIpDataType);
        break;
    case NuTo::ElementData::CONSTITUTIVELAWIPNONLOCAL:
        ptrElementData = new NuTo::ElementDataConstitutiveIpNonlocal(this, rNumIp, rIpDataType);
        break;
    case NuTo::ElementData::CONSTITUTIVELAWIPCRACK:
        ptrElementData = new NuTo::ElementDataConstitutiveIpCrack(this, rNumIp, rIpDataType);
        break;
    case NuTo::ElementData::VARIABLECONSTITUTIVELAWIP:
        ptrElementData = new NuTo::ElementDataVariableConstitutiveIp(this, rNumIp, rIpDataType);
        break;

    default:
        throw MechanicsException("[NuTo::ElementBase::ElementBase] unsupported element data type.");
    }
    mElementData = ptrElementData;
}

// destructor
NuTo::ElementBase::~ElementBase()
{
    if (this->mElementData != 0)
    {
//    	printf("Delete %p",static_cast<void*>(mElementData));
        delete this->mElementData;
    }
}

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


//! @brief calculates output data for the element with a standard input (EULER_BACKWARD static data)
//! @param rOutput ...  coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
NuTo::Error::eError NuTo::ElementBase::Evaluate(std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rOutput)
{
    ConstitutiveInputMap input;
    ConstitutiveCalculateStaticData calculateImplicitly(CalculateStaticData::EULER_BACKWARD);
    input[Constitutive::Input::CALCULATE_STATIC_DATA] = &calculateImplicitly;

    return this->Evaluate(input, rOutput);
}

//! @brief returns the id number of the element
//! @return id
int NuTo::ElementBase::ElementGetId() const
{
    return mStructure->ElementGetId(this);
}

//! @brief returns a pointer to the constitutive law for an integration point
//! @param integration point number (counting from zero)
//! @return pointer to constitutive law
const NuTo::ConstitutiveBase* NuTo::ElementBase::GetConstitutiveLaw(int rIp) const
{
    return mElementData->GetConstitutiveLaw(rIp);
}

//! @brief returns a pointer to the constitutive law for an integration point
//! @param integration point number (counting from zero)
//! @return pointer to constitutive law
NuTo::ConstitutiveBase* NuTo::ElementBase::GetConstitutiveLaw(int rIp)
{
    return mElementData->GetConstitutiveLaw(rIp);
}

//! @brief sets the constitutive law for an element
//! @param rConstitutiveLaw Pointer to constitutive law entry
void NuTo::ElementBase::SetConstitutiveLaw(ConstitutiveBase* rConstitutiveLaw)
{
    //check compatibility between element type and constitutive law
    if (rConstitutiveLaw->CheckElementCompatibility(this->GetEnumType()))
    {
        //printf("check element constitutive was positive.\n");

        mElementData->SetConstitutiveLaw(this, rConstitutiveLaw);
    } else
    {
        std::stringstream message;
        message << "[NuTo::ElementBase::SetConstitutiveLaw] Constitutive Law " << mStructure->ConstitutiveLawGetId(rConstitutiveLaw) << " does not match element type of element " << mStructure->ElementGetId(this) << "." << std::endl;
        throw MechanicsException(message.str());
    }
}

//! @brief sets the constitutive law for an element
//! @param rIp id of integration point
//! @param rConstitutiveLaw Pointer to constitutive law entry
void NuTo::ElementBase::SetConstitutiveLaw(int rIp, ConstitutiveBase* rConstitutiveLaw)
{
    //check compatibility between element type and constitutive law
    if (rConstitutiveLaw->CheckElementCompatibility(this->GetEnumType()))
    {
        //printf("check element constitutive was positive.\n");

        mElementData->SetConstitutiveLaw(this, rIp, rConstitutiveLaw);
    } else
    {
        std::stringstream message;
        message << "[NuTo::ElementBase::SetConstitutiveLaw] Constitutive Law " << mStructure->ConstitutiveLawGetId(rConstitutiveLaw) << " does not match element type of element " << mStructure->ElementGetId(this) << "." << std::endl;
        throw MechanicsException(message.str());
    }
}

//! @brief returns true, if the constitutive law has been assigned
bool NuTo::ElementBase::HasConstitutiveLawAssigned(int rIp) const
{
    return mElementData->HasConstitutiveLawAssigned(rIp);
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

//    std::cout << matrixN * nodalValues << std::endl;

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

//! @brief calculates the area of a plane element via the nodes (probably faster than sum over integration points)
//! @return Area
double NuTo::ElementBase::CalculateArea() const
{
    throw MechanicsException("[NuTo::ElementBase::CalculateArea] The area can only be calculated for 2D elements.");
}

//! @brief sets the integration type of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need an integration type
//! @param rIntegrationType pointer to integration type
void NuTo::ElementBase::SetIntegrationType(const NuTo::IntegrationTypeBase* rIntegrationType, NuTo::IpData::eIpDataType rIpDataType)
{
    //check compatibility between element type and constitutive law
    if (rIntegrationType->CheckElementCompatibility(this->GetEnumType()))
    {
        mElementData->SetIntegrationType(this, rIntegrationType, rIpDataType);
    } else
    {
        std::stringstream message;
        message << "[NuTo::ElementBase::SetIntegrationType] Integration Type does not match element type of element " << mStructure->ElementGetId(this) << "." << std::endl;
        throw MechanicsException(message.str());
    }
}

//! @brief returns a pointer to the integration type of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need an integration type
//! @return pointer to integration type
const NuTo::IntegrationTypeBase* NuTo::ElementBase::GetIntegrationType() const
{
    return mElementData->GetIntegrationType();
}

//! @brief sets the interpolation type of an element
//! @param rInterpolationType interpolation type
void NuTo::ElementBase::SetInterpolationType(const InterpolationType* rInterpolationType)
{
    mInterpolationType = rInterpolationType;
}

//! @brief returns a pointer to the interpolation type of an element
//! @return pointer to interpolation type
const NuTo::InterpolationType* NuTo::ElementBase::GetInterpolationType() const
{
    return mInterpolationType;
}

//! @brief returns ip data type of the element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need an integration type
//! @return pointer to integration type
NuTo::IpData::eIpDataType NuTo::ElementBase::GetIpDataType(int rIp) const
{
    return mElementData->GetIpDataType(rIp);
}

//! @brief returns the enum of element data type
//! @return enum of ElementDataType
const NuTo::ElementData::eElementDataType NuTo::ElementBase::GetElementDataType() const
{
    return mElementData->GetElementDataType();
}

//! @brief returns the number of integration points
//! @return number of integration points
int NuTo::ElementBase::GetNumIntegrationPoints() const
{
    return this->mElementData->GetNumIntegrationPoints();
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point
//! @return weight
double NuTo::ElementBase::GetIntegrationPointWeight(int rIpNum) const
{
    return this->mElementData->GetIntegrationPointWeight(rIpNum);
}

//! @brief returns the static data of an integration point
//! @param rIp integration point
//! @return static data pointer
const NuTo::ConstitutiveStaticDataBase* NuTo::ElementBase::GetStaticData(int rIp) const
{
    return this->mElementData->GetStaticData(rIp);
}

//! @brief returns the static data of an integration point
//! @param rIp integration point
//! @return static data pointer
NuTo::ConstitutiveStaticDataBase* NuTo::ElementBase::GetStaticData(int rIp)
{
    return this->mElementData->GetStaticData(rIp);
}

NuTo::IpDataStaticDataBase& NuTo::ElementBase::GetStaticDataBase(int rIp)
{
    return mElementData->GetStaticDataBase(rIp);
}

const NuTo::IpDataStaticDataBase& NuTo::ElementBase::GetStaticDataBase(int rIp) const
{
    return mElementData->GetStaticDataBase(rIp);
}

//! @brief sets the static data for an integration point of an element
//! @param rIp integration point
//! @param rStaticData static data
void NuTo::ElementBase::SetStaticData(int rIp, ConstitutiveStaticDataBase* rStaticData)
{
    if (rStaticData->CheckConstitutiveCompatibility(mElementData->GetConstitutiveLaw(rIp)->GetType(), this->GetEnumType()))
        return this->mElementData->SetStaticData(rIp, rStaticData);
    else
        throw MechanicsException("[NuTo::ElementBase::SetStaticData] Static data is not compatible with the element and or constitutive model");
}

const Eigen::Vector3d NuTo::ElementBase::GetGlobalIntegrationPointCoordinates(int rIpNum) const
{
    const Eigen::MatrixXd& matrixN = mInterpolationType->Get(Node::COORDINATES).GetMatrixN(rIpNum);
    Eigen::VectorXd nodeCoordinates = ExtractNodeValues(0, Node::COORDINATES);

    Eigen::Vector3d globalIntegrationPointCoordinates = Eigen::Vector3d::Zero();

    globalIntegrationPointCoordinates.segment(0, GetLocalDimension()) = matrixN * nodeCoordinates;

    return globalIntegrationPointCoordinates;
}

//! @brief returns the natural coordinates of an given point
//! implemented with an exception for all elements, reimplementation required for those elements
//! @param rGlobCoords (input) ... pointer to the array of coordinates
//! @param rLocCoords (output) ... coordinates to be returned
//! @return True if coordinates are within the element, False otherwise
bool NuTo::ElementBase::GetLocalPointCoordinates(const double* rGlobCoords, double* rLocCoords) const
{
    throw NuTo::MechanicsException("[NuTo::ElementBase::GetLocalPointCoordinates] not implemented for this element type.");
}

//! @brief Gets the control node of an boundary element, if it has one
//! @return boundary control node
NuTo::NodeBase *NuTo::ElementBase::GetBoundaryControlNode() const
{
    throw NuTo::MechanicsException(__PRETTY_FUNCTION__,"Not implemented for this element type.");
}

//VHIRTHAMTODO Remove???
////! @brief sets the water volume fraction at the boundary surface
////! @return water volume fraction at the boundary surface
//double NuTo::ElementBase::GetBoundaryWaterVolumeFraction() const
//{
//    throw NuTo::MechanicsException("[NuTo::ElementBase::GetBoundaryWaterVolumeFraction] not implemented for this element type.");
//}

////! @brief sets the water volume fraction at the boundary surface
////! @param water volume fraction at the boundary surface
//void NuTo::ElementBase::SetBoundaryWaterVolumeFraction(double rBoundaryWaterVolumeFraction)
//{
//    throw NuTo::MechanicsException("[NuTo::ElementBase::SetBoundaryWaterVolumeFraction] not implemented for this element type.");
//}

////! @brief sets the relative humidity at the boundary surface
////! @param relative humidity at the boundary surface
//double NuTo::ElementBase::GetBoundaryRelativeHumidity() const
//{
//    throw NuTo::MechanicsException("[NuTo::ElementBase::GetBoundaryRelativeHumidity] not implemented for this element type.");
//}

////! @brief sets the relative humidity at the boundary surface
////! @param relative humidity at the boundary surface
//void NuTo::ElementBase::SetBoundaryRelativeHumidity(double rBoundaryRelativeHumidity)
//{
//    throw NuTo::MechanicsException("[NuTo::ElementBase::SetBoundaryRelativeHumidity] not implemented for this element type.");
//}

//! @brief checks if a node is inside of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! @param rGlobCoords (input) ... pointer to the array of coordinates
//! @return True if coordinates are within the element, False otherwise
bool NuTo::ElementBase::CheckPointInside(const double* rGlobCoords) const
{
    throw NuTo::MechanicsException("[NuTo::ElementBase::CheckPointInside] not implemented for this element type.");
}

#ifdef ENABLE_VISUALIZE
void NuTo::ElementBase::GetVisualizationCells(unsigned int& NumVisualizationPoints, std::vector<double>& VisualizationPointLocalCoordinates, unsigned int& NumVisualizationCells, std::vector<NuTo::CellBase::eCellTypes>& VisualizationCellType, std::vector<unsigned int>& VisualizationCellsIncidence,
        std::vector<unsigned int>& VisualizationCellsIP) const
{
    const IntegrationTypeBase* integrationType = this->mElementData->GetIntegrationType();
    if (integrationType == 0 || integrationType->GetNumIntegrationPoints() == 0)
        throw MechanicsException("[NuTo::ElementBase::GetVisualizationCells] integration type for this element is not set.");

    //get integration cells from the integration type
    integrationType->GetVisualizationCells(NumVisualizationPoints, VisualizationPointLocalCoordinates, NumVisualizationCells, VisualizationCellType, VisualizationCellsIncidence, VisualizationCellsIP);
}

void NuTo::ElementBase::Visualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList)
{
    // get visualization cells from integration type
    unsigned int NumVisualizationPoints;
    std::vector<double> VisualizationPointLocalCoordinates;
    unsigned int NumVisualizationCells;
    std::vector<NuTo::CellBase::eCellTypes> VisualizationCellType;
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
        case NuTo::CellBase::LINE:
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
        case NuTo::CellBase::TRIANGLE:
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
        case NuTo::CellBase::QUAD:
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
        case NuTo::CellBase::TETRAEDER:
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
        case NuTo::CellBase::HEXAHEDRON:
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
    elementOutput[Element::IP_DATA] = std::make_shared<ElementOutputIpData>();

    auto& elementIpDataMap = elementOutput.at(Element::IP_DATA)->GetIpData().GetIpDataMap();

    bool evaluateStress(false);


    for (auto const &it : rVisualizationList)
    {
        switch (it.get()->GetComponentEnum())
        {
        case NuTo::VisualizeBase::BOND_STRESS:
            elementIpDataMap[IpData::BOND_STRESS];
            break;
        case NuTo::VisualizeBase::DAMAGE:
            elementIpDataMap[IpData::DAMAGE];
            break;
        case NuTo::VisualizeBase::ENGINEERING_STRAIN:
            elementIpDataMap[IpData::ENGINEERING_STRAIN];
            break;
        case NuTo::VisualizeBase::LOCAL_EQ_STRAIN:
            elementIpDataMap[IpData::LOCAL_EQ_STRAIN];
            break;
        case NuTo::VisualizeBase::ENGINEERING_PLASTIC_STRAIN:
            elementIpDataMap[IpData::ENGINEERING_PLASTIC_STRAIN];
            break;
        case NuTo::VisualizeBase::TOTAL_INELASTIC_EQ_STRAIN:
            elementIpDataMap[IpData::TOTAL_INELASTIC_EQ_STRAIN];
            break;
        case NuTo::VisualizeBase::ENGINEERING_STRESS:
            if (evaluateStress == false)
            {
                elementIpDataMap[IpData::ENGINEERING_STRESS];
                evaluateStress = true;
            }
            break;
        case NuTo::VisualizeBase::PRINCIPAL_ENGINEERING_STRESS:
            if (evaluateStress == false)
            {
                elementIpDataMap[IpData::ENGINEERING_STRESS];
                evaluateStress = true;
            }
            break;
        case NuTo::VisualizeBase::HEAT_FLUX:
            elementIpDataMap[IpData::HEAT_FLUX];
            break;
        case NuTo::VisualizeBase::DISPLACEMENTS:
        case NuTo::VisualizeBase::NONLOCAL_WEIGHT:
        case NuTo::VisualizeBase::NONLOCAL_EQ_STRAIN:
        case NuTo::VisualizeBase::CONSTITUTIVE:
        case NuTo::VisualizeBase::SECTION:
        case NuTo::VisualizeBase::ELEMENT:
        case NuTo::VisualizeBase::CRACK:
        case NuTo::VisualizeBase::TEMPERATURE:
        case NuTo::VisualizeBase::ROTATION:
        case NuTo::VisualizeBase::VELOCITY:
        case NuTo::VisualizeBase::ACCELERATION:
        case NuTo::VisualizeBase::ANGULAR_VELOCITY:
        case NuTo::VisualizeBase::ANGULAR_ACCELERATION:
        case NuTo::VisualizeBase::PARTICLE_RADIUS:
        case NuTo::VisualizeBase::RELATIVE_HUMIDITY:
        case NuTo::VisualizeBase::WATER_VOLUME_FRACTION:
        default:
            //do nothing
            ;
            break;
        }
    }

    //calculate the element solution
//    ConstitutiveInputMap input;
//    ConstitutiveCalculateStaticData calc(CalculateStaticData::USE_PREVIOUS);
//    input[Constitutive::Input::CALCULATE_STATIC_DATA] = &calc;
//    Evaluate(input, elementOutput);
    Evaluate(elementOutput);

    //assign the outputs

    // store data
    for (auto const &it : rVisualizationList)
    {
        switch (it.get()->GetComponentEnum())
        {
        case NuTo::VisualizeBase::DAMAGE:
        {
            const auto& damage = elementIpDataMap.at(IpData::DAMAGE);
            assert(damage.size() != 0);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataScalar(CellId, it.get()->GetComponentName(), damage.data()[theIp]);
            }
        }
        break;
        case NuTo::VisualizeBase::LOCAL_EQ_STRAIN:
        {
            const auto& localEqStrain = elementIpDataMap.at(IpData::LOCAL_EQ_STRAIN);
            assert(localEqStrain.size() != 0);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataScalar(CellId, it.get()->GetComponentName(), localEqStrain.data()[theIp]);
            }
        }
        break;
        case NuTo::VisualizeBase::DISPLACEMENTS:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                const Eigen::VectorXd& coords = visualizationPointNaturalCoordinates.col(PointCount);
                Eigen::Vector3d GlobalDisplacements = this->InterpolateDof3D(coords, Node::eDof::DISPLACEMENTS);
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(), GlobalDisplacements.data());
            }
            break;
        case NuTo::VisualizeBase::VELOCITY:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                const Eigen::VectorXd& coords = visualizationPointNaturalCoordinates.col(PointCount);
                Eigen::Vector3d GlobalVelocity = this->InterpolateDof3D(1, coords, Node::eDof::DISPLACEMENTS);
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(), GlobalVelocity.data());
            }
            break;
        case NuTo::VisualizeBase::ACCELERATION:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                const Eigen::VectorXd& coords = visualizationPointNaturalCoordinates.col(PointCount);
                Eigen::Vector3d GlobalAcceleration = this->InterpolateDof3D(2, coords, Node::eDof::DISPLACEMENTS);
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(), GlobalAcceleration.data());
            }
            break;
        case NuTo::VisualizeBase::ENGINEERING_STRAIN:
        {
            const auto& engineeringStrain = elementIpDataMap.at(IpData::ENGINEERING_STRAIN);
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
        case NuTo::VisualizeBase::ENGINEERING_PLASTIC_STRAIN:
        {
            const auto& engineeringPlasticStrain = elementIpDataMap.at(IpData::ENGINEERING_PLASTIC_STRAIN);
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
        case NuTo::VisualizeBase::TOTAL_INELASTIC_EQ_STRAIN:
        {
            const auto& totalInelasticEqStrain = elementIpDataMap.at(IpData::TOTAL_INELASTIC_EQ_STRAIN);
            assert(totalInelasticEqStrain.size() != 0);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataScalar(CellId, it.get()->GetComponentName(), totalInelasticEqStrain.data()[theIp]);
            }
        }
            break;
        case NuTo::VisualizeBase::ENGINEERING_STRESS:
        {
            const auto& engineeringStress = elementIpDataMap.at(IpData::ENGINEERING_STRESS);
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
        case NuTo::VisualizeBase::BOND_STRESS:
        {
            const auto& bondStress = elementIpDataMap.at(IpData::BOND_STRESS);
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
        case NuTo::VisualizeBase::PRINCIPAL_ENGINEERING_STRESS:
        {
            const auto& engineeringStress = elementIpDataMap.at(IpData::ENGINEERING_STRESS);
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
        case NuTo::VisualizeBase::NONLOCAL_WEIGHT:
        {
            // get constitutive model and nonlocal element
            const ElementBase* visualizeElement(it.get()->GetElement());

            // get local number within nonlocal elements for current element
            try
            {
                const std::vector<const NuTo::ElementBase*>& nonlocalElements(visualizeElement->GetNonlocalElements());
                int countNonlocalElement;
                for (countNonlocalElement = 0; countNonlocalElement < (int) nonlocalElements.size(); countNonlocalElement++)
                {
                    if (this == nonlocalElements[countNonlocalElement])
                    {
                        //visualize the weights
                        break;
                    }
                }
                if (countNonlocalElement < (int) nonlocalElements.size())
                {
                    // get the weights
                    for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
                    {
                        unsigned int theIp = VisualizationCellsIP[CellCount];
                        unsigned int CellId = CellIdVec[CellCount];
                        double weight(visualizeElement->GetNonlocalWeights(it.get()->GetIp(), countNonlocalElement)[theIp]);
                        rVisualize.SetCellDataScalar(CellId, it.get()->GetComponentName(), weight);
                    }
                } else
                {
                    // otherwise they just remain zero
                    for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
                    {
                        unsigned int CellId = CellIdVec[CellCount];
                        rVisualize.SetCellDataScalar(CellId, it.get()->GetComponentName(), 0);
                    }
                }
            } catch (NuTo::MechanicsException &e)
            {
                e.AddMessage("[NuTo::ElementBase::Visualize] error getting visualization data for nonlocal weights.");
                throw e;
            }
#ifdef ENABLE_VISUALIZE
            catch (NuTo::VisualizeException &e)
            {
                e.AddMessage("[NuTo::ElementBase::Visualize] error getting visualization data for nonlocal weights.");
                throw e;
            }
#endif
            catch (...)
            {
                throw NuTo::MechanicsException("[NuTo::ElementBase::Visualize] error getting visualization data for nonlocal weights.");
            }
        }
            break;
        case NuTo::VisualizeBase::NONLOCAL_EQ_STRAIN:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                Eigen::VectorXd nonlocalEqStrain(1);
                if (mInterpolationType->IsDof(Node::NONLOCALEQSTRAIN))
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
        case NuTo::VisualizeBase::CONSTITUTIVE:
        {
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                unsigned int CellId = CellIdVec[CellCount];
                int constitutiveId = mStructure->ConstitutiveLawGetId(GetConstitutiveLaw(theIp));

                rVisualize.SetCellDataScalar(CellId, it.get()->GetComponentName(), constitutiveId);
            }
        }
            break;
        case NuTo::VisualizeBase::SECTION:
        {
            int sectionId = mStructure->SectionGetId(GetSection());
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataScalar(CellId, it.get()->GetComponentName(), sectionId);
            }
        }
            break;
        case NuTo::VisualizeBase::ELEMENT:
        {
            int elementId = this->ElementGetId();
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataScalar(CellId, it.get()->GetComponentName(), elementId);
            }
        }
            break;
        case NuTo::VisualizeBase::CRACK:
        {
            std::vector<NuTo::CrackBase*> elementCracks = this->GetDataPtr()->GetCracks();
            //! @todo [DA] change SetCellDataVector to variable data size
            if (2 < elementCracks.size())
                throw NuTo::VisualizeException("[NuTo::ElementBase::Visualize::CRACK] cannot visualize more than 3 cracks per elements");
            double elementCrackIds[] =
            { -1.0, -1.0, -1.0 };
            unsigned int crackCount = 0;
            BOOST_FOREACH(const NuTo::CrackBase* thisCrack, elementCracks)
                elementCrackIds[crackCount++] = this->mStructure->CrackGetId(thisCrack);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
                rVisualize.SetCellDataVector(CellIdVec[CellCount], it.get()->GetComponentName(), elementCrackIds);
        }
            break;
        case NuTo::VisualizeBase::PARTICLE_RADIUS:
            //do nothing
            break;
        case NuTo::VisualizeBase::ROTATION:
        case NuTo::VisualizeBase::ANGULAR_VELOCITY:
        case NuTo::VisualizeBase::ANGULAR_ACCELERATION:
            //do nothing
            break;
        case NuTo::VisualizeBase::HEAT_FLUX:
        {
            const auto& heatFlux = elementIpDataMap.at(IpData::HEAT_FLUX);
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
        case NuTo::VisualizeBase::TEMPERATURE:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                const Eigen::VectorXd& coords = visualizationPointNaturalCoordinates.col(PointCount);
                Eigen::VectorXd temperature = InterpolateDofGlobal(coords, Node::eDof::TEMPERATURE);
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataScalar(PointId, it.get()->GetComponentName(), temperature[0]);
            }
            break;

        case NuTo::VisualizeBase::RELATIVE_HUMIDITY:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                Eigen::Matrix<double, Eigen::Dynamic, 1> relativeHumidity;
                if (mInterpolationType->IsDof(Node::RELATIVEHUMIDITY))
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
        case NuTo::VisualizeBase::WATER_VOLUME_FRACTION:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                Eigen::Matrix<double, Eigen::Dynamic, 1> waterVolumeFraction;
                if (mInterpolationType->IsDof(Node::WATERVOLUMEFRACTION))
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
            std::cout << it.get()->GetComponentEnum() << "\n";
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
    std::vector<NuTo::CellBase::eCellTypes> VisualizationCellType;
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
    VisualizationCellType.push_back(NuTo::CellBase::TRIANGLE);
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
        case NuTo::CellBase::TRIANGLE:
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
        case NuTo::CellBase::QUAD:
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
        case NuTo::VisualizeBase::ENGINEERING_STRAIN:
            boost::assign::ptr_map_insert<ElementOutputIpData>( elementOutput )( Element::IP_DATA ,IpData::ENGINEERING_STRAIN);
        break;
        case NuTo::VisualizeBase::DISPLACEMENTS:
        case NuTo::VisualizeBase::SECTION:
        case NuTo::VisualizeBase::ELEMENT:
        default:
            //do nothing
            break;
        }
    }

    //calculate the element solution
    Evaluate(elementOutput);

    //assign the outputs
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>* engineeringStrain(nullptr);

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
        case NuTo::VisualizeBase::DISPLACEMENTS:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                const Eigen::VectorXd& coords = visualizationPointNaturalCoordinates.col(PointCount);
                Eigen::Vector3d GlobalDisplacements = this->InterpolateDof3D(coords, Node::eDof::DISPLACEMENTS);
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataVector(PointId, it.get()->GetComponentName(), GlobalDisplacements.data());
            }
            break;
        case NuTo::VisualizeBase::ENGINEERING_STRAIN:
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
            case NuTo::VisualizeBase::SECTION:
        {
            int sectionId = mStructure->SectionGetId(GetSection());
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataScalar(CellId, it.get()->GetComponentName(), sectionId);
            }
        }
            break;
        case NuTo::VisualizeBase::ELEMENT:
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
    unsigned int NumVisualizationPoints = GetInterpolationType()->GetCurrentIntegrationType()->GetNumIntegrationPoints();
    unsigned int NumVisualizationCells = NumVisualizationPoints;
    std::vector<NuTo::CellBase::eCellTypes> VisualizationCellType;


    std::vector<double> VisualizationPointLocalCoordinates;

    std::vector<unsigned int> VisualizationCellsIncidence;
    std::vector<unsigned int> VisualizationCellsIP;
    int dimension = GetLocalDimension();



    for (unsigned int iIp = 0; iIp < NumVisualizationPoints; ++iIp)
    {
        double naturalIpCoords[2];
        GetInterpolationType()->GetCurrentIntegrationType()->GetLocalIntegrationPointCoordinates2D(iIp, naturalIpCoords);

        for (int iDim = 0; iDim < dimension; ++iDim)
        {
            VisualizationPointLocalCoordinates.push_back(naturalIpCoords[iDim]);
        }

        VisualizationCellType.push_back(NuTo::CellBase::VERTEX);
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
            case NuTo::CellBase::VERTEX:
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
        elementOutput[Element::IP_DATA] = std::make_shared<ElementOutputIpData>();
        auto& elementIpDataMap = elementOutput[Element::IP_DATA]->GetIpData().GetIpDataMap();

        for (auto const &it : rVisualizationList)
        {
            switch (it.get()->GetComponentEnum())
            {
            case NuTo::VisualizeBase::ENGINEERING_STRAIN:
                elementIpDataMap[IpData::ENGINEERING_STRAIN];
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
            case NuTo::VisualizeBase::ENGINEERING_STRAIN:
            {
                const auto& engineeringStrain = elementIpDataMap.at(IpData::ENGINEERING_STRAIN);
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
                default:
                throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) + ": \t Visualization component " + it.get()->GetComponentName() + " is not implemented or not known at the integration points.");
            }
        }

}
#endif // ENABLE_VISUALIZE

//! @brief adds the nonlocal IP to an integration point
//! @param rLocalIpNumber local Ip
//! @param rConstitutive constitutive model for which nonlocal data is to be calculated
//! @param rNonlocalElement element of the nonlocal ip
//! @param rNonlocalIp local ip number of the nonlocal ip
//! @param rWeight weight
void NuTo::ElementBase::SetNonlocalWeight(int rLocalIpNumber, const ElementBase* rNonlocalElement, int rNonlocalIp, double rWeight)
{
    this->mElementData->SetNonlocalWeight(rLocalIpNumber, rNonlocalElement, rNonlocalIp, rWeight);
}

//! @brief returns a vector of weights for an ip
//! @param rIp local Ip
//! @param rNonlocalElement nonlocal element (must be in the range of the nonlocal element size stored at the element data level)
//! @return weights for each integration point of the nonlocal element
const std::vector<double>& NuTo::ElementBase::GetNonlocalWeights(int rIp, int rNonlocalElement) const
{
    return mElementData->GetNonlocalWeights(rIp, rNonlocalElement);
}

//! @brief returns a vector of the nonlocal elements
//! @return nonlocal elements
const std::vector<const NuTo::ElementBase*>& NuTo::ElementBase::GetNonlocalElements() const
{
    return this->mElementData->GetNonlocalElements();
}

//! @brief delete the nonlocal elements
void NuTo::ElementBase::DeleteNonlocalElements()
{
    this->mElementData->DeleteNonlocalElements();
}

//! @brief returns a vector of the nonlocal elements
//! @param rConstitutive constitutive model for the nonlocale elements
//! @return nonlocal elements
int NuTo::ElementBase::GetNumNonlocalElements() const
{
    return this->mElementData->GetNumNonlocalElements();
}

//! @brief integrates the stress over the element
//! @param rStress integrated stress
void NuTo::ElementBase::GetIntegratedStress(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rStress)
{
    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutput;
    elementOutput[Element::IP_DATA] = std::make_shared<ElementOutputIpData>(IpData::ENGINEERING_STRESS);

    this->Evaluate(elementOutput);

    const auto& ipStress = elementOutput.at(Element::IP_DATA)->GetIpData().GetIpDataMap()[IpData::ENGINEERING_STRESS];
    Eigen::VectorXd ipVolume = this->GetIntegrationPointVolume();

    rStress.Resize(ipStress.GetNumRows(), 1);
    for (int countIP = 0; countIP < ipStress.GetNumColumns(); countIP++)
    {
        rStress += (ipStress.col(countIP) * (ipVolume[countIP]));
    }
}

//! @brief integrates the strain over the element
//! @param rStrain integrated strain
void NuTo::ElementBase::GetIntegratedStrain(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rStrain)
{
    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutput;
    elementOutput[Element::IP_DATA] = std::make_shared<ElementOutputIpData>(IpData::ENGINEERING_STRAIN);

    this->Evaluate(elementOutput);

    const auto& ipStress = elementOutput.at(Element::IP_DATA)->GetIpData().GetIpDataMap()[IpData::ENGINEERING_STRAIN];
    Eigen::VectorXd ipVolume = this->GetIntegrationPointVolume();

    rStrain.Resize(ipStress.GetNumRows(), 1);
    for (int countIP = 0; countIP < ipStress.GetNumColumns(); countIP++)
    {
        rStrain += (ipStress.col(countIP) * (ipVolume[countIP]));
    }
}

//! @brief Returns the vector of crack pointers of an element
//! @return crack pointer vector
const std::vector<NuTo::CrackBase*> NuTo::ElementBase::GetCracks() const
{
    return this->mElementData->GetCracks();
}

//! @brief Set the information that the element is already cracked or not
//! @param bool cracked or not
void NuTo::ElementBase::IsCracked(const bool rIsCracked)
{
    this->mElementData->IsCracked(rIsCracked);
}

//! @brief Give the information if the element is already cracked or not
//! @return bool cracked or not
const bool NuTo::ElementBase::IsCracked() const
{
    return this->mElementData->IsCracked();
}


void NuTo::ElementBase::Info() const
{
    mStructure->GetLogger() << "[" << __PRETTY_FUNCTION__ << "] \n";
    mStructure->GetLogger() << "InterpolationTypeInfo:\n" << GetInterpolationType()->Info() << "\n";

    for (int iNode = 0; iNode < GetNumNodes(); ++iNode)
    {
        const NodeBase* node = GetNode(iNode);
        mStructure->GetLogger() << "NodeInfo of local node " << iNode << ": \n";
        mStructure->GetLogger() << node->GetNodeTypeStr() << "\n";
    }
//    mStructure->GetLogger() << "InterpolationTypeInfo: \n " << GetInterpolationType()->Info();

}

//! @brief returns the Element Data Vector
//! this was necessary due to recursive problems for serialization (nonlocal data)
//! this method should only be called from the serialization routine of the structure
NuTo::ElementDataBase* NuTo::ElementBase::GetDataPtr() const
{
    return mElementData;
}

//! @brief sets the Element Data Vector
//! this was necessary due to recursive problems for serialization (nonlocal data)
//! this method should only be called from the serialization routine of the structure
void NuTo::ElementBase::SetDataPtr(NuTo::ElementDataBase* rElementData)
{
    if (mElementData != 0)
        delete mElementData;
    mElementData = rElementData;
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
