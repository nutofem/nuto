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
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
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
#endif

NuTo::ElementBase::ElementBase(const StructureBase* rStructure, ElementData::eElementDataType rElementDataType,
        IpData::eIpDataType rIpDataType, const InterpolationType* rInterpolationType)
        : mStructure(rStructure), mInterpolationType(rInterpolationType)
{

    const IntegrationTypeBase* integrationType = mInterpolationType->GetCurrentIntegrationType();

    //allocate element data
	ElementDataBase* ptrElementData;
	switch(rElementDataType)
	{
	case NuTo::ElementData::NOELEMENTDATA:
		throw MechanicsException("[NuTo::ElementBase::ElementBase] no elements without element data implemented.");
    break;
	case NuTo::ElementData::CONSTITUTIVELAWIP:
		ptrElementData = new NuTo::ElementDataConstitutiveIp(this,integrationType,rIpDataType);
    break;
	case NuTo::ElementData::CONSTITUTIVELAWIPNONLOCAL:
		ptrElementData = new NuTo::ElementDataConstitutiveIpNonlocal(this,integrationType,rIpDataType);
    break;
	case NuTo::ElementData::CONSTITUTIVELAWIPCRACK:
		ptrElementData = new NuTo::ElementDataConstitutiveIpCrack(this,integrationType,rIpDataType);
    break;
	case NuTo::ElementData::VARIABLECONSTITUTIVELAWIP:
		ptrElementData = new NuTo::ElementDataVariableConstitutiveIp(this,integrationType,rIpDataType);
    break;

    default:
		throw MechanicsException("[NuTo::ElementBase::ElementWithDataBase] unsupported element data type.");
	}
	mElementData = ptrElementData;
}

NuTo::ElementBase::ElementBase(const StructureBase* rStructure, ElementData::eElementDataType rElementDataType,
		int rNumIp, IpData::eIpDataType rIpDataType, const InterpolationType* rInterpolationType)
        : mStructure(rStructure), mInterpolationType(rInterpolationType)
{

    //allocate element data
	ElementDataBase* ptrElementData;
	switch(rElementDataType)
	{
	case NuTo::ElementData::NOELEMENTDATA:
		throw MechanicsException("[NuTo::ElementBase::ElementBase] no elements without element data implemented.");
    break;
	case NuTo::ElementData::CONSTITUTIVELAWIP:
		ptrElementData = new NuTo::ElementDataConstitutiveIp(this,rNumIp,rIpDataType);
    break;
	case NuTo::ElementData::CONSTITUTIVELAWIPNONLOCAL:
		ptrElementData = new NuTo::ElementDataConstitutiveIpNonlocal(this,rNumIp,rIpDataType);
    break;
	case NuTo::ElementData::CONSTITUTIVELAWIPCRACK:
		ptrElementData = new NuTo::ElementDataConstitutiveIpCrack(this,rNumIp,rIpDataType);
    break;
	case NuTo::ElementData::VARIABLECONSTITUTIVELAWIP:
		ptrElementData = new NuTo::ElementDataVariableConstitutiveIp(this,rNumIp,rIpDataType);
    break;

    default:
		throw MechanicsException("[NuTo::ElementBase::ElementBase] unsupported element data type.");
	}
	mElementData = ptrElementData;
}

// destructor
NuTo::ElementBase::~ElementBase()
{
    if(this->mElementData != 0)
    {
//    	printf("Delete %p",static_cast<void*>(mElementData));
        delete this->mElementData;
    }
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ElementBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ElementBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ElementBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ElementBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ElementBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ElementBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ElementBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ElementBase " << std::endl;
#endif
    ar & BOOST_SERIALIZATION_NVP(mStructure);
    // the element data has to be saved on the main structure due to problems with a recursion on the stack (nonlocal data contains ptr to elements)
    // the idea is to first serialize all the elements in the table, and afterwards update the pointers of the element data in the element data routine
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ElementBase" << std::endl;
#endif
}
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ElementBase)
#endif // ENABLE_SERIALIZATION

//! @brief returns the id number of the element
//! @return id
int NuTo::ElementBase::ElementGetId()const
{
	return mStructure->ElementGetId(this);
}

//! @brief returns a pointer to the constitutive law for an integration point
//! @param integration point number (counting from zero)
//! @return pointer to constitutive law
const NuTo::ConstitutiveBase* NuTo::ElementBase::GetConstitutiveLaw(int rIp)const
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
	}
	else
	{
		std::stringstream message;
		message << "[NuTo::ElementBase::SetConstitutiveLaw] Constitutive Law " << mStructure->ConstitutiveLawGetId(rConstitutiveLaw)
				<<" does not match element type of element "<< mStructure->ElementGetId(this) <<"." <<std::endl;
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

		mElementData->SetConstitutiveLaw(this,rIp, rConstitutiveLaw);
	}
	else
	{
		std::stringstream message;
		message << "[NuTo::ElementBase::SetConstitutiveLaw] Constitutive Law " << mStructure->ConstitutiveLawGetId(rConstitutiveLaw)
				<<" does not match element type of element "<< mStructure->ElementGetId(this) <<"." <<std::endl;
	    throw MechanicsException(message.str());
	}
}

//! @brief returns true, if the constitutive law has been assigned
bool NuTo::ElementBase::HasConstitutiveLawAssigned(int rIp)const
{
	return mElementData->HasConstitutiveLawAssigned(rIp);
}

int NuTo::ElementBase::GetNumNodes() const
{
    return mInterpolationType->GetNumNodes();
}

int NuTo::ElementBase::GetNumNodes(Node::eAttributes rDofType) const
{
    return mInterpolationType->Get(rDofType).GetNumNodes();
}

const Eigen::MatrixXd NuTo::ElementBase::ExtractNodeValues(Node::eAttributes rDofType) const
{
    return this->ExtractNodeValues(0, rDofType);
}

const Eigen::MatrixXd NuTo::ElementBase::ExtractNodeValues(int rTimeDerivative, Node::eAttributes rDofType) const
{
    throw NuTo::MechanicsException("[NuTo::ElementBase::ExtractNodeValues] not implemented.");
}

const Eigen::VectorXd NuTo::ElementBase::InterpolateDof(const Eigen::VectorXd& rNaturalCoordinates, Node::eAttributes rDofType) const
{
    return InterpolateDof(0, rNaturalCoordinates, rDofType);
}

const Eigen::VectorXd NuTo::ElementBase::InterpolateDof(int rTimeDerivative, const Eigen::VectorXd& rNaturalCoordinates, Node::eAttributes rDofType) const
{
    const InterpolationBase& interpolationType = mInterpolationType->Get(rDofType);
    const Eigen::MatrixXd nodalValues = ExtractNodeValues(rTimeDerivative, rDofType);
    const Eigen::VectorXd shapeFunctions = interpolationType.CalculateShapeFunctions(rNaturalCoordinates);

    return nodalValues * shapeFunctions;
}

const Eigen::Vector3d NuTo::ElementBase::InterpolateDof3D(const Eigen::VectorXd& rNaturalCoordinates, Node::eAttributes rDofType) const
{
    return InterpolateDof3D(0, rNaturalCoordinates, rDofType);
}

const Eigen::Vector3d NuTo::ElementBase::InterpolateDof3D(int rTimeDerivative, const Eigen::VectorXd& rNaturalCoordinates, Node::eAttributes rDofType) const
{
    const InterpolationBase& interpolationType = mInterpolationType->Get(rDofType);
    const Eigen::MatrixXd nodalValues = ExtractNodeValues(rTimeDerivative, rDofType);
    const Eigen::VectorXd shapeFunctions = interpolationType.CalculateShapeFunctions(rNaturalCoordinates);

    Eigen::VectorXd result = nodalValues * shapeFunctions;
    Eigen::Vector3d interpolation3D = Eigen::Vector3d::Zero();


    for (int iRow = 0; iRow < result.rows(); ++iRow)
        interpolation3D[iRow] = result[iRow];

    return interpolation3D;
}

//! @brief calculates the area of a plane element via the nodes (probably faster than sum over integration points)
//! @return Area
double NuTo::ElementBase::CalculateArea()const
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
	}
	else
	{
		std::stringstream message;
		message << "[NuTo::ElementBase::SetIntegrationType] Integration Type does not match element type of element "<< mStructure->ElementGetId(this) <<"." <<std::endl;
	    throw MechanicsException(message.str());
	}
}

//! @brief returns a pointer to the integration type of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need an integration type
//! @return pointer to integration type
const NuTo::IntegrationTypeBase* NuTo::ElementBase::GetIntegrationType()const
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
const NuTo::InterpolationType* NuTo::ElementBase::GetInterpolationType () const
{
    return mInterpolationType;
}

//! @brief returns ip data type of the element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need an integration type
//! @return pointer to integration type
NuTo::IpData::eIpDataType NuTo::ElementBase::GetIpDataType(int  rIp)const
{
	return mElementData->GetIpDataType(rIp);
}

//! @brief returns the enum of element data type
//! @return enum of ElementDataType
const NuTo::ElementData::eElementDataType NuTo::ElementBase::GetElementDataType()const
{
    return mElementData->GetElementDataType();
}

//! @brief returns the number of integration points
//! @return number of integration points
int NuTo::ElementBase::GetNumIntegrationPoints()const
{
    return this->mElementData->GetNumIntegrationPoints();
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point
//! @return weight
double NuTo::ElementBase::GetIntegrationPointWeight(int rIpNum)const
{
    return this->mElementData->GetIntegrationPointWeight(rIpNum);
}

//! @brief calculate the length of an edge (belonging to an integration point for lattice elements)
//! @param rIp integration point
//! @return edge length
double NuTo::ElementBase::GetIpEdgeLength(int rIp)const
{
	throw MechanicsException("[NuTo::ElementBase::GetIpEdgeLength] not implement for this type of element. (implemented for lattice elements).");
}

//! @brief returns the static data of an integration point
//! @param rIp integration point
//! @return static data pointer
const NuTo::ConstitutiveStaticDataBase* NuTo::ElementBase::GetStaticData(int rIp)const
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

//! @brief sets the static data for an integration point of an element
//! @param rIp integration point
//! @param rStaticData static data
void NuTo::ElementBase::SetStaticData(int rIp, ConstitutiveStaticDataBase* rStaticData)
{
	if (rStaticData->CheckConstitutiveCompatibility(mElementData->GetConstitutiveLaw(rIp)->GetType(), this->GetEnumType()))
	    return this->mElementData->SetStaticData(rIp,rStaticData);
	else
		throw MechanicsException("[NuTo::ElementBase::SetStaticData] Static data is not compatible with the element and or constitutive model");
}

const Eigen::Vector3d NuTo::ElementBase::GetGlobalIntegrationPointCoordinates(int rIpNum) const
{
    Eigen::VectorXd shapeFunctions = mInterpolationType->Get(Node::COORDINATES).GetShapeFunctions(rIpNum);
    Eigen::MatrixXd nodeCoordinates = ExtractNodeValues(0, Node::COORDINATES);
    assert(shapeFunctions.rows() == nodeCoordinates.cols());

    Eigen::Vector3d globalIntegrationPointCoordinates = Eigen::Vector3d::Zero();
    for (int iRow = 0; iRow < nodeCoordinates.rows(); iRow++)
        globalIntegrationPointCoordinates[iRow] = nodeCoordinates.row(iRow) * shapeFunctions;

    return globalIntegrationPointCoordinates;
}

//! @brief returns the natural coordinates of an given point
//! implemented with an exception for all elements, reimplementation required for those elements
//! @param rGlobCoords (input) ... pointer to the array of coordinates
//! @param rLocCoords (output) ... coordinates to be returned
//! @return True if coordinates are within the element, False otherwise
bool NuTo::ElementBase::GetLocalPointCoordinates(const double* rGlobCoords,  double* rLocCoords)const
{
    throw NuTo::MechanicsException("[NuTo::ElementBase::GetLocalPointCoordinates] not implemented for this element type.");
}

//! @brief Gets the additional node of an boundary element, if it has one
//! @return Additional boundary node
NuTo::NodeBase *NuTo::ElementBase::GetAdditionalBoundaryNode() const
{
    throw NuTo::MechanicsException("[NuTo::ElementBase::GetAdditionalBoundaryNode] not implemented for this element type.");
}

//! @brief sets the water volume fraction at the boundary surface
//! @return water volume fraction at the boundary surface
double NuTo::ElementBase::GetBoundaryWaterVolumeFraction() const
{
    throw NuTo::MechanicsException("[NuTo::ElementBase::GetBoundaryWaterVolumeFraction] not implemented for this element type.");
}

//! @brief sets the water volume fraction at the boundary surface
//! @param water volume fraction at the boundary surface
void NuTo::ElementBase::SetBoundaryWaterVolumeFraction(double rBoundaryWaterVolumeFraction)
{
    throw NuTo::MechanicsException("[NuTo::ElementBase::SetBoundaryWaterVolumeFraction] not implemented for this element type.");
}

//! @brief sets the relative humidity at the boundary surface
//! @param relative humidity at the boundary surface
double NuTo::ElementBase::GetBoundaryRelativeHumidity() const
{
    throw NuTo::MechanicsException("[NuTo::ElementBase::GetBoundaryRelativeHumidity] not implemented for this element type.");
}

//! @brief sets the relative humidity at the boundary surface
//! @param relative humidity at the boundary surface
void NuTo::ElementBase::SetBoundaryRelativeHumidity(double rBoundaryRelativeHumidity)
{
    throw NuTo::MechanicsException("[NuTo::ElementBase::SetBoundaryRelativeHumidity] not implemented for this element type.");
}

//! @brief checks if a node is inside of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! @param rGlobCoords (input) ... pointer to the array of coordinates
//! @return True if coordinates are within the element, False otherwise
bool NuTo::ElementBase::CheckPointInside(const double* rGlobCoords)const
{
	throw NuTo::MechanicsException("[NuTo::ElementBase::CheckPointInside] not implemented for this element type.");
}

#ifdef ENABLE_VISUALIZE
void NuTo::ElementBase::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::CellBase::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    const IntegrationTypeBase* integrationType = this->mElementData->GetIntegrationType();
    if (integrationType==0 || integrationType->GetNumIntegrationPoints()==0)
    	throw MechanicsException("[NuTo::ElementBase::GetVisualizationCells] integration type for this element is not set.");

    //get integration cells from the integration type
    integrationType->GetVisualizationCells(
            NumVisualizationPoints,
            VisualizationPointLocalCoordinates,
            NumVisualizationCells,
            VisualizationCellType,
            VisualizationCellsIncidence,
            VisualizationCellsIP);
}


void NuTo::ElementBase::Visualize(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat)
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
    GetVisualizationCells(
        NumVisualizationPoints,
        VisualizationPointLocalCoordinates,
        NumVisualizationCells,
        VisualizationCellType,
        VisualizationCellsIncidence,
        VisualizationCellsIP);


    // calculate global point coordinates and store point at the visualize object
    int dimension (VisualizationPointLocalCoordinates.size()/NumVisualizationPoints);
    assert(VisualizationPointLocalCoordinates.size() == NumVisualizationPoints * dimension);

    // TODO: fix that by proper visualization point (natural!) coordinates, local might be misleading here
    Eigen::MatrixXd visualizationPointNaturalCoordinates = Eigen::MatrixXd::Map(VisualizationPointLocalCoordinates.data(), dimension, NumVisualizationPoints);

    std::vector<unsigned int> PointIdVec;
    for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
    {
        if (dimension != 1 and dimension != 2 and dimension != 3)
            throw NuTo::MechanicsException("[NuTo::ElementBase::Visualize] invalid dimension of local coordinates");

        const Eigen::VectorXd& coords = visualizationPointNaturalCoordinates.col(PointCount);
        Eigen::Matrix<double, 3, 1> GlobalPointCoor = this->InterpolateDof3D(coords, Node::eAttributes::COORDINATES);
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
            Points[0] = PointIdVec[VisualizationCellsIncidence[Pos    ]];
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
	boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;
	bool evaluateStress(false);
    boost::ptr_list<VisualizeComponentBase>::const_iterator WhatIter = rWhat.begin();
    for (auto WhatIter = rWhat.begin();WhatIter != rWhat.end(); WhatIter++)
    {
        switch (WhatIter->GetComponentEnum())
        {
		case NuTo::VisualizeBase::DAMAGE:
			boost::assign::ptr_map_insert<ElementOutputIpData>( elementOutput )( Element::IP_DATA ,IpData::DAMAGE);
		break;
		case NuTo::VisualizeBase::ENGINEERING_STRAIN:
			boost::assign::ptr_map_insert<ElementOutputIpData>( elementOutput )( Element::IP_DATA ,IpData::ENGINEERING_STRAIN);
		break;
		case NuTo::VisualizeBase::ENGINEERING_PLASTIC_STRAIN:
			boost::assign::ptr_map_insert<ElementOutputIpData>( elementOutput )( Element::IP_DATA ,IpData::ENGINEERING_PLASTIC_STRAIN);
		break;
		case NuTo::VisualizeBase::ENGINEERING_STRESS:
			if (evaluateStress==false)
			{
				boost::assign::ptr_map_insert<ElementOutputIpData>( elementOutput )( Element::IP_DATA ,IpData::ENGINEERING_STRESS);
				evaluateStress=true;
			}
		break;
		case NuTo::VisualizeBase::PRINCIPAL_ENGINEERING_STRESS:
			if (evaluateStress==false)
			{
				boost::assign::ptr_map_insert<ElementOutputIpData>( elementOutput )( Element::IP_DATA ,IpData::ENGINEERING_STRESS);
				evaluateStress=true;
			}
		break;
		case NuTo::VisualizeBase::HEAT_FLUX:
			boost::assign::ptr_map_insert<ElementOutputIpData>( elementOutput )( Element::IP_DATA ,IpData::HEAT_FLUX);
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
    Evaluate(elementOutput);













    //assign the outputs
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>* damage(0);
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>* engineeringStrain(0);
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>* engineeringPlasticStrain(0);
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>* engineeringStress(0);
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>* heatFlux(0);
    for (auto itElementOutput=elementOutput.begin(); itElementOutput!=elementOutput.end(); itElementOutput++)
    {
        switch (itElementOutput->second->GetIpDataType())
        {
		case NuTo::IpData::DAMAGE:
			damage = &(itElementOutput->second->GetFullMatrixDouble());
		break;
		case NuTo::IpData::ENGINEERING_STRAIN:
			engineeringStrain = &(itElementOutput->second->GetFullMatrixDouble());
		break;
		case NuTo::IpData::ENGINEERING_PLASTIC_STRAIN:
			engineeringPlasticStrain = &(itElementOutput->second->GetFullMatrixDouble());
		break;
		case NuTo::IpData::ENGINEERING_STRESS:
			engineeringStress = &(itElementOutput->second->GetFullMatrixDouble());
		break;
		case NuTo::IpData::HEAT_FLUX:
			heatFlux = &(itElementOutput->second->GetFullMatrixDouble());
		break;
		default:
			throw MechanicsException("[NuTo::ElementBase::Visualize] other ipdatatypes not supported.");
        }
    }

    // store data
    for (auto WhatIter = rWhat.begin();WhatIter != rWhat.end(); WhatIter++)
    {
        switch (WhatIter->GetComponentEnum())
        {
        case NuTo::VisualizeBase::DAMAGE:
        {
        	assert(damage!=0);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataScalar(CellId, WhatIter->GetComponentName(), damage->data()[theIp]);
            }
        }
        break;
        case NuTo::VisualizeBase::DISPLACEMENTS:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                const Eigen::VectorXd& coords = visualizationPointNaturalCoordinates.col(PointCount);
                Eigen::Vector3d GlobalDisplacements = this->InterpolateDof3D(coords, Node::eAttributes::DISPLACEMENTS);
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataVector(PointId, WhatIter->GetComponentName(), GlobalDisplacements.data());
            }
            break;
        case NuTo::VisualizeBase::VELOCITY:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                const Eigen::VectorXd& coords = visualizationPointNaturalCoordinates.col(PointCount);
                Eigen::Vector3d GlobalVelocity = this->InterpolateDof3D(1,coords, Node::eAttributes::DISPLACEMENTS);
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataVector(PointId, WhatIter->GetComponentName(), GlobalVelocity.data());
            }
            break;
        case NuTo::VisualizeBase::ACCELERATION:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                const Eigen::VectorXd& coords = visualizationPointNaturalCoordinates.col(PointCount);
                Eigen::Vector3d GlobalAcceleration = this->InterpolateDof3D(2,coords, Node::eAttributes::DISPLACEMENTS);
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataVector(PointId, WhatIter->GetComponentName(), GlobalAcceleration.data());
            }
            break;
        case NuTo::VisualizeBase::ENGINEERING_STRAIN:
        {
        	assert(engineeringStrain!=0);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                const double* EngineeringStrainVector = &(engineeringStrain->data()[theIp*6]);
                double EngineeringStrainTensor[9];
                EngineeringStrainTensor[0] = EngineeringStrainVector[0];
                EngineeringStrainTensor[1] = 0.5 * EngineeringStrainVector[3];
                EngineeringStrainTensor[2] = 0.5 * EngineeringStrainVector[5];
                EngineeringStrainTensor[3] = 0.5 * EngineeringStrainVector[3];
                EngineeringStrainTensor[4] = EngineeringStrainVector[1];
                EngineeringStrainTensor[5] = 0.5 * EngineeringStrainVector[4];
                EngineeringStrainTensor[6] = 0.5 * EngineeringStrainVector[5];
                EngineeringStrainTensor[7] = 0.5 * EngineeringStrainVector[4];
                EngineeringStrainTensor[8] = EngineeringStrainVector[2];

                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataTensor(CellId, WhatIter->GetComponentName(), EngineeringStrainTensor);
            }
        }
        break;
        case NuTo::VisualizeBase::ENGINEERING_PLASTIC_STRAIN:
        {
        	assert(engineeringPlasticStrain!=0);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                const double* EngineeringStrainVector = &(engineeringPlasticStrain->data()[theIp*6]);
                double EngineeringStrainTensor[9];
                EngineeringStrainTensor[0] = EngineeringStrainVector[0];
                EngineeringStrainTensor[1] = 0.5 * EngineeringStrainVector[3];
                EngineeringStrainTensor[2] = 0.5 * EngineeringStrainVector[5];
                EngineeringStrainTensor[3] = 0.5 * EngineeringStrainVector[3];
                EngineeringStrainTensor[4] = EngineeringStrainVector[1];
                EngineeringStrainTensor[5] = 0.5 * EngineeringStrainVector[4];
                EngineeringStrainTensor[6] = 0.5 * EngineeringStrainVector[5];
                EngineeringStrainTensor[7] = 0.5 * EngineeringStrainVector[4];
                EngineeringStrainTensor[8] = EngineeringStrainVector[2];

                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataTensor(CellId, WhatIter->GetComponentName(), EngineeringStrainTensor);
            }
        }
        break;
        case NuTo::VisualizeBase::ENGINEERING_STRESS:
        {
        	assert(engineeringStress!=0);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                const double* EngineeringStressVector = &(engineeringStress->data()[theIp*6]);
                double EngineeringStressTensor[9];
                EngineeringStressTensor[0] = EngineeringStressVector[0];
                EngineeringStressTensor[1] = EngineeringStressVector[3];
                EngineeringStressTensor[2] = EngineeringStressVector[5];
                EngineeringStressTensor[3] = EngineeringStressVector[3];
                EngineeringStressTensor[4] = EngineeringStressVector[1];
                EngineeringStressTensor[5] = EngineeringStressVector[4];
                EngineeringStressTensor[6] = EngineeringStressVector[5];
                EngineeringStressTensor[7] = EngineeringStressVector[4];
                EngineeringStressTensor[8] = EngineeringStressVector[2];

                //std::cout<<"[NuTo::ElementBase::VisualizeEngineeringStressTensor]" << EngineeringStressTensor[0] << EngineeringStressTensor[1] << std::endl;
                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataTensor(CellId, WhatIter->GetComponentName(), EngineeringStressTensor);
            }
        }
        break;
        case NuTo::VisualizeBase::PRINCIPAL_ENGINEERING_STRESS:
        {
        	assert(engineeringStress!=0);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                const double* EngineeringStressVector = &(engineeringStress->data()[theIp*6]);
                Eigen::Matrix<double,3,3>  EngineeringStressTensor;

                EngineeringStressTensor(0,0) = EngineeringStressVector[0];
                EngineeringStressTensor(1,0) = EngineeringStressVector[3];
                EngineeringStressTensor(2,0) = EngineeringStressVector[5];
                EngineeringStressTensor(0,1) = EngineeringStressVector[3];
                EngineeringStressTensor(1,1) = EngineeringStressVector[1];
                EngineeringStressTensor(2,1) = EngineeringStressVector[4];
                EngineeringStressTensor(0,2) = EngineeringStressVector[5];
                EngineeringStressTensor(1,2) = EngineeringStressVector[4];
                EngineeringStressTensor(2,2) = EngineeringStressVector[2];

                //std::cout<<"[NuTo::ElementBase::VisualizeEngineeringStressTensor]" << EngineeringStressTensor[0] << EngineeringStressTensor[1] << std::endl;
                unsigned int CellId = CellIdVec[CellCount];
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> mySolver(EngineeringStressTensor,false);
                Eigen::Matrix<double,3,1>  eigenValues(mySolver.eigenvalues());
                rVisualize.SetCellDataVector(CellId, WhatIter->GetComponentName(), eigenValues.data());
            }
        }
        break;
        case NuTo::VisualizeBase::NONLOCAL_WEIGHT:
        {
        	// get constitutive model and nonlocal element
        	const ElementBase* visualizeElement(WhatIter->GetElement());

        	// get local number within nonlocal elements for current element
        	try
        	{
            	const std::vector<const NuTo::ElementBase*>& nonlocalElements(visualizeElement->GetNonlocalElements());
            	int countNonlocalElement;
            	for (countNonlocalElement=0; countNonlocalElement<(int)nonlocalElements.size(); countNonlocalElement++)
            	{
            		if (this==nonlocalElements[countNonlocalElement])
            		{
            			//visualize the weights
            			break;
            		}
            	}
            	if (countNonlocalElement<(int)nonlocalElements.size())
            	{
            		// get the weights
                    for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
                     {
                         unsigned int theIp = VisualizationCellsIP[CellCount];
                         unsigned int CellId = CellIdVec[CellCount];
                         double weight(visualizeElement->GetNonlocalWeights(WhatIter->GetIp(),countNonlocalElement)[theIp]);
                         rVisualize.SetCellDataScalar(CellId, WhatIter->GetComponentName(),weight);
                     }
            	}
            	else
            	{
            		// otherwise they just remain zero
                    for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
                     {
                         unsigned int CellId = CellIdVec[CellCount];
                         rVisualize.SetCellDataScalar(CellId, WhatIter->GetComponentName(), 0);
                     }
            	}
        	}
            catch (NuTo::MechanicsException &e)
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
            catch(...)
            {
            	throw NuTo::MechanicsException("[NuTo::ElementBase::Visualize] error getting visualization data for nonlocal weights.");
            }
        }
        break;
        case NuTo::VisualizeBase::NONLOCAL_EQ_STRAIN:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                Eigen::Matrix<double, Eigen::Dynamic, 1> nonlocalEqStrain(1);
                if (mInterpolationType->IsDof(Node::NONLOCALEQSTRAIN))
                {
                    // calculate only if element has nonlocal eq strain dofs
                    const Eigen::VectorXd& coords = visualizationPointNaturalCoordinates.col(PointCount);
                    nonlocalEqStrain =  InterpolateDof(coords, Node::eAttributes::NONLOCALEQSTRAIN);
                    assert(nonlocalEqStrain.rows() == 1);
                }
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataScalar(PointId, WhatIter->GetComponentName(), nonlocalEqStrain.at(0,0));
            }
            break;
        case NuTo::VisualizeBase::CONSTITUTIVE:
        {
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                unsigned int CellId = CellIdVec[CellCount];
                int constitutiveId = mStructure->ConstitutiveLawGetId(GetConstitutiveLaw(theIp));

                rVisualize.SetCellDataScalar(CellId, WhatIter->GetComponentName(), constitutiveId);
            }
        }
        break;
        case NuTo::VisualizeBase::SECTION:
        {
            int sectionId = mStructure->SectionGetId(GetSection());
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataScalar(CellId, WhatIter->GetComponentName(), sectionId);
            }
        }
        break;
        case NuTo::VisualizeBase::ELEMENT:
        {
            int elementId = this->ElementGetId();
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataScalar(CellId, WhatIter->GetComponentName(), elementId);
            }
        }
        break;
        case NuTo::VisualizeBase::CRACK:
        {
        	std::vector<NuTo::CrackBase*> elementCracks = this->GetDataPtr()->GetCracks();
        	//! @todo [DA] change SetCellDataVector to variable data size
        	if(2 < elementCracks.size())
            	throw NuTo::VisualizeException("[NuTo::ElementBase::Visualize::CRACK] cannot visualize more than 3 cracks per elements");
        	double elementCrackIds[]={-1.0,-1.0,-1.0};
        	unsigned int crackCount=0;
        	BOOST_FOREACH(const NuTo::CrackBase* thisCrack, elementCracks)
            	elementCrackIds[crackCount++] = this->mStructure->CrackGetId(thisCrack);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
                rVisualize.SetCellDataVector(CellIdVec[CellCount], WhatIter->GetComponentName(), elementCrackIds );
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
        	assert(heatFlux!=0);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataVector(CellId, WhatIter->GetComponentName(), &(heatFlux->data()[theIp*3]));
            }
        }
        break;
        case NuTo::VisualizeBase::TEMPERATURE:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                const Eigen::VectorXd& coords = visualizationPointNaturalCoordinates.col(PointCount);
                Eigen::VectorXd temperature = InterpolateDof(coords, Node::eAttributes::TEMPERATURES);
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataScalar(PointId, WhatIter->GetComponentName(), temperature.at(0,0));
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
                    relativeHumidity =  InterpolateDof(coords, Node::eAttributes::RELATIVEHUMIDITY);
                    assert(relativeHumidity.rows() == 1);
                }
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataScalar(PointId, WhatIter->GetComponentName(), relativeHumidity.at(0,0));
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
                    waterVolumeFraction =  InterpolateDof(coords, Node::eAttributes::WATERVOLUMEFRACTION);
                    assert(waterVolumeFraction.rows() == 1);
                }
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataScalar(PointId, WhatIter->GetComponentName(), waterVolumeFraction     .at(0,0));
            }
            break;
		default:
			std::cout << WhatIter->GetComponentEnum() << "\n";
			throw NuTo::MechanicsException("[NuTo::ElementBase::Visualize] unsupported datatype for visualization.");
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
void NuTo::ElementBase::SetNonlocalWeight(int rLocalIpNumber,	const ElementBase* rNonlocalElement, int rNonlocalIp, double rWeight)
{
	this->mElementData->SetNonlocalWeight(rLocalIpNumber,rNonlocalElement,rNonlocalIp,rWeight);
}

//! @brief returns a vector of weights for an ip
//! @param rIp local Ip
//! @param rNonlocalElement nonlocal element (must be in the range of the nonlocal element size stored at the element data level)
//! @return weights for each integration point of the nonlocal element
const std::vector<double>& NuTo::ElementBase::GetNonlocalWeights(int rIp, int rNonlocalElement)const
{
    return mElementData->GetNonlocalWeights(rIp,rNonlocalElement);
}

//! @brief returns a vector of the nonlocal elements
//! @retrun nonlocal elements
const std::vector<const NuTo::ElementBase*>& NuTo::ElementBase::GetNonlocalElements()const
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
//! @retrun nonlocal elements
int NuTo::ElementBase::GetNumNonlocalElements()const
{
	return this->mElementData->GetNumNonlocalElements();
}

//! @brief integrates the stress over the element
//! @param rStress integrated stress
void NuTo::ElementBase::GetIntegratedStress(FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rStress)
{
    // define local variables storing the element contribution would be wrong, since at the moment you pass it to the ptr_map, you will run into problems.
    // My theory (TT): In this case, the d'tor of the ElementOutputIpData would be called twice, trying to delete the same underlying data array twice.
    // --> use the ptr_maps as intended with pointers, or better the insert operators

    NuTo::Element::eOutput keyIP_DATA = Element::IP_DATA;
    boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;

    // first brackets: the ptr_map
    // second brackets: (KEY, arguments for c'tor of ElementOutputIpData)
    boost::assign::ptr_map_insert<ElementOutputIpData>(elementOutput)(keyIP_DATA, IpData::ENGINEERING_STRESS);

    this->Evaluate(elementOutput);

    // optional: extract the data from the ptr_map
    ElementOutputBase& elementOutputIpData = elementOutput.at(keyIP_DATA);

    const NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& ipStress = elementOutputIpData.GetFullMatrixDouble();
    const Eigen::VectorXd ipVolume = this->GetIntegrationPointVolume();

    rStress.Resize(ipStress.GetNumRows(), 1);
    for (int countIP = 0; countIP < ipStress.GetNumColumns(); countIP++)
    {
        rStress += (ipStress.col(countIP) * (ipVolume[countIP]));
    }
}

//! @brief integrates the strain over the element
//! @param rStrain integrated strain
void NuTo::ElementBase::GetIntegratedStrain(FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rStrain)
{
    // see [NuTo::ElementBase::GetIntegratedStress]
    NuTo::Element::eOutput keyIP_DATA(Element::IP_DATA);

    boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;
    boost::assign::ptr_map_insert<ElementOutputIpData>(elementOutput)(keyIP_DATA, IpData::ENGINEERING_STRAIN);

    this->Evaluate(elementOutput);
    ElementOutputBase& elementOutputIpData = elementOutput.at(keyIP_DATA);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rIpStrain(elementOutputIpData.GetFullMatrixDouble());
    Eigen::VectorXd ipVolume = this->GetIntegrationPointVolume();

    rStrain.Resize(rIpStrain.GetNumRows(), 1);
    for (int countIP = 0; countIP < rIpStrain.GetNumColumns(); countIP++)
    {
        rStrain += rIpStrain.col(countIP) * (ipVolume[countIP]);
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

//! @brief cast the base pointer to an Element1D, otherwise throws an exception
const NuTo::Element1D* NuTo::ElementBase::AsElement1D()const
{
    throw NuTo::MechanicsException("[NuTo::ElementBase::AsElement1D] Element is not of type Element1D.");
}

//! @brief cast the base pointer to an Element1D, otherwise throws an exception
NuTo::Element1D* NuTo::ElementBase::AsElement1D()
{
    throw NuTo::MechanicsException("[NuTo::ElementBase::AsElement1D] Element is not of type Element1D.");
}

//! @brief cast the base pointer to an Element2D, otherwise throws an exception
const NuTo::Element2D* NuTo::ElementBase::AsElement2D()const
{
    throw NuTo::MechanicsException("[NuTo::ElementBase::AsElement2D] Element is not of type Element2D.");
}

//! @brief cast the base pointer to an Element2D, otherwise throws an exception
NuTo::Element2D* NuTo::ElementBase::AsElement2D()
{
    throw NuTo::MechanicsException("[NuTo::ElementBase::AsElement2D] Element is not of type Element2D.");
}

//! @brief cast the base pointer to an Element3D, otherwise throws an exception
const NuTo::Element3D* NuTo::ElementBase::AsElement3D()const
{
    throw NuTo::MechanicsException("[NuTo::ElementBase::AsElement3D] Element is not of type Element3D.");
}

//! @brief cast the base pointer to an Element3D, otherwise throws an exception
NuTo::Element3D* NuTo::ElementBase::AsElement3D()
{
    throw NuTo::MechanicsException("[NuTo::ElementBase::AsElement3D] Element is not of type Element3D.");
}


//! @brief cast the base pointer to an BoundaryElement1D, otherwise throws an exception
const NuTo::BoundaryElement1D* NuTo::ElementBase::AsBoundaryElement1D()const
{
    throw NuTo::MechanicsException("[NuTo::ElementBase::AsBoundaryElement1D] BoundaryElement is not of type BoundaryElement1D.");
}

//! @brief cast the base pointer to an BoundaryElement1D, otherwise throws an exception
NuTo::BoundaryElement1D* NuTo::ElementBase::AsBoundaryElement1D()
{
    throw NuTo::MechanicsException("[NuTo::ElementBase::AsBoundaryElement1D] BoundaryElement is not of type BoundaryElement1D.");
}

//! @brief cast the base pointer to an BoundaryElement2D, otherwise throws an exception
const NuTo::BoundaryElement2D* NuTo::ElementBase::AsBoundaryElement2D()const
{
    throw NuTo::MechanicsException("[NuTo::ElementBase::AsBoundaryElement2D] BoundaryElement is not of type BoundaryElement2D.");
}

//! @brief cast the base pointer to an BoundaryElement2D, otherwise throws an exception
NuTo::BoundaryElement2D* NuTo::ElementBase::AsBoundaryElement2D()
{
    throw NuTo::MechanicsException("[NuTo::ElementBase::AsBoundaryElement2D] BoundaryElement is not of type BoundaryElement2D.");
}

//! @brief cast the base pointer to an BoundaryElement3D, otherwise throws an exception
const NuTo::BoundaryElement3D* NuTo::ElementBase::AsBoundaryElement3D()const
{
    throw NuTo::MechanicsException("[NuTo::ElementBase::AsBoundaryElement3D] BoundaryElement is not of type BoundaryElement3D.");
}

//! @brief cast the base pointer to an BoundaryElement3D, otherwise throws an exception
NuTo::BoundaryElement3D* NuTo::ElementBase::AsBoundaryElement3D()
{
    throw NuTo::MechanicsException("[NuTo::ElementBase::AsBoundaryElement3D] BoundaryElement is not of type BoundaryElement3D.");
}


//! @brief cast the base pointer to a Plane, otherwise throws an exception
const NuTo::Plane* NuTo::ElementBase::AsPlane()const
{
	throw NuTo::MechanicsException("[NuTo::ElementBase::AsElementPlane] Element is not of type Plane.");
}

//! @brief cast the base pointer to a Plane, otherwise throws an exception
NuTo::Plane* NuTo::ElementBase::AsPlane()
{
	throw NuTo::MechanicsException("[NuTo::ElementBase::AsElementPlane] Element is not of type Plane.");
}

//! @brief cast the base pointer to a Plane, otherwise throws an exception
const NuTo::Plane2D* NuTo::ElementBase::AsPlane2D()const
{
	throw NuTo::MechanicsException("[NuTo::ElementBase::AsPlane2D] Element is not of type Plane2D.");
}

//! @brief cast the base pointer to a Plane, otherwise throws an exception
NuTo::Plane2D* NuTo::ElementBase::AsPlane2D()
{
	throw NuTo::MechanicsException("[NuTo::ElementBase::AsPlane2D] Element is not of type Plane2D.");
}

//! @brief cast the base pointer to a Solid, otherwise throws an exception
const NuTo::Solid* NuTo::ElementBase::AsSolid()const
{
	throw NuTo::MechanicsException("[NuTo::ElementBase::AsElementSolid] Element is not of type ElementSolid.");
}

//! @brief cast the base pointer to a Solid, otherwise throws an exception
NuTo::Solid* NuTo::ElementBase::AsSolid()
{
	throw NuTo::MechanicsException("[NuTo::ElementBase::AsElementSolid] Element is not of type ElementSolid.");
}

//! @brief cast the base pointer to a Truss, otherwise throws an exception
const NuTo::Truss* NuTo::ElementBase::AsTruss()const
{
	throw NuTo::MechanicsException("[NuTo::ElementBase::AsElementTruss] Element is not of type ElementTruss.");
}

//! @brief cast the base pointer to a Truss, otherwise throws an exception
NuTo::Truss* NuTo::ElementBase::AsTruss()
{
	throw NuTo::MechanicsException("[NuTo::ElementBase::AsElementTruss] Element is not of type ElementTruss.");
}

//! @brief returns the Element Data Vector
//! this was necessary due to recursive problems for serialization (nonlocal data)
//! this method should only be called from the serialization routine of the structure
NuTo::ElementDataBase* NuTo::ElementBase::GetDataPtr()const
{
    return mElementData;
}

//! @brief sets the Element Data Vector
//! this was necessary due to recursive problems for serialization (nonlocal data)
//! this method should only be called from the serialization routine of the structure
void NuTo::ElementBase::SetDataPtr(NuTo::ElementDataBase* rElementData)
{
    if (mElementData!=0)
        delete mElementData;
    mElementData = rElementData;
}

void NuTo::ElementBase::ReorderNodes()
{
    const Eigen::MatrixX2i& reorderIndices = mInterpolationType->GetNodeRenumberingIndices();
    for (int i = 0; i < reorderIndices.rows(); ++i)
    {
        int i0 = reorderIndices(i,0);
        int i1 = reorderIndices(i,1);
        // swap nodes i0 and i1
        NodeBase* tmpNode0 = GetNode(i0);
        NodeBase* tmpNode1 = GetNode(i1);
        SetNode(i0, tmpNode1);
        SetNode(i1, tmpNode0);

        if (mStructure->GetVerboseLevel() > 5)
            mStructure->GetLogger() << "[NuTo::ElementBase::ReorderNodes] Swapped nodes " << i0 << " and " << i1<< ".\n";
    }
    mStructure->GetLogger() << "\n";
}
