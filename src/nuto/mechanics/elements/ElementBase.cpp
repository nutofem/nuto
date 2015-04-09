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
		IntegrationType::eIntegrationType rIntegrationType, IpData::eIpDataType rIpDataType)
        : mStructure(rStructure)
{
    //allocate element data
	ElementDataBase* ptrElementData;
	switch(rElementDataType)
	{
	case NuTo::ElementData::NOELEMENTDATA:
		throw MechanicsException("[NuTo::ElementWithDataBase::ElementWithDataBase] no elements without element data implemented.");
    break;
	case NuTo::ElementData::CONSTITUTIVELAWIP:
		ptrElementData = new NuTo::ElementDataConstitutiveIp(this,const_cast<StructureBase*>(rStructure)->GetPtrIntegrationType(rIntegrationType),rIpDataType);
    break;
	case NuTo::ElementData::CONSTITUTIVELAWIPNONLOCAL:
		ptrElementData = new NuTo::ElementDataConstitutiveIpNonlocal(this,const_cast<StructureBase*>(rStructure)->GetPtrIntegrationType(rIntegrationType),rIpDataType);
    break;
	case NuTo::ElementData::CONSTITUTIVELAWIPCRACK:
		ptrElementData = new NuTo::ElementDataConstitutiveIpCrack(this,const_cast<StructureBase*>(rStructure)->GetPtrIntegrationType(rIntegrationType),rIpDataType);
    break;

    default:
		throw MechanicsException("[NuTo::ElementWithDataBase::ElementWithDataBase] unsupported element data type.");
	}
	mElementData = ptrElementData;
}

NuTo::ElementBase::ElementBase(const StructureBase* rStructure, ElementData::eElementDataType rElementDataType,
		int rNumIp, IpData::eIpDataType rIpDataType)
        : mStructure(rStructure)
{
    //allocate element data
	ElementDataBase* ptrElementData;
	switch(rElementDataType)
	{
	case NuTo::ElementData::NOELEMENTDATA:
		throw MechanicsException("[NuTo::ElementWithDataBase::ElementWithDataBase] no elements without element data implemented.");
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

    default:
		throw MechanicsException("[NuTo::ElementWithDataBase::ElementWithDataBase] unsupported element data type.");
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
		message << "[NuTo::ElementWithDataBase::SetConstitutiveLaw] Constitutive Law " << mStructure->ConstitutiveLawGetId(rConstitutiveLaw)
				<<" does not match element type of element "<< mStructure->ElementGetId(this) <<"." <<std::endl;
	    throw MechanicsException(message.str());
	}
}

//! @brief returns true, if the constitutive law has been assigned
bool NuTo::ElementBase::HasConstitutiveLawAssigned(int rIp)const
{
	return mElementData->HasConstitutiveLawAssigned(rIp);
}

void NuTo::ElementBase::InterpolateCoordinatesFrom1D(double rLocalCoordinates, double rGlobalCoordinates[3]) const
{
    throw MechanicsException("[NuTo::ElementBase::InterpolateCoordinatesFrom1D] 1D geometry interpolation routine not implemented.");
}

void NuTo::ElementBase::InterpolateCoordinatesFrom2D(double rLocalCoordinates[2], double rGlobalCoordinates[3]) const
{
    throw MechanicsException("[NuTo::ElementBase::InterpolateCoordinatesFrom2D] 2D geometry interpolation routine not implemented.");
}

void NuTo::ElementBase::InterpolateCoordinatesFrom3D(double rLocalCoordinates[3], double rGlobalCoordinates[3]) const
{
    throw MechanicsException("[NuTo::ElementBase::InterpolateCoordinatesFrom3D] 3D geometry interpolation routine not implemented.");
}

void NuTo::ElementBase::InterpolateDisplacementsFrom1D(int rTimeDerivative, double rLocalCoordinates, double rGlobalDisplacements[3]) const
{
    throw MechanicsException("[NuTo::ElementBase::InterpolateDisplacementsFrom1D] 1D displacement interpolation routine not implemented.");
}

void NuTo::ElementBase::InterpolateDisplacementsFrom2D(int rTimeDerivative, double rLocalCoordinates[2], double rGlobalDisplacements[3]) const
{
    throw MechanicsException("[NuTo::ElementBase::InterpolateDisplacementsFrom2D] 2D displacement interpolation routine not implemented.");
}

void NuTo::ElementBase::InterpolateDisplacementsFrom3D(int rTimeDerivative, double rLocalCoordinates[3], double rGlobalDisplacements[3]) const
{
    throw MechanicsException("[NuTo::ElementBase::InterpolateDisplacementsFrom3D] 3D displacement interpolation routine not implemented.");
}

void NuTo::ElementBase::InterpolateTemperatureFrom1D(double rLocalCoordinates, double& rTemperature) const
{
    throw MechanicsException("[NuTo::ElementBase::InterpolateTemperatureFrom1D] 1D temperature interpolation routine not implemented.");
}

void NuTo::ElementBase::InterpolateTemperatureFrom2D(double rLocalCoordinates[2], double& rTemperature) const
{
    throw MechanicsException("[NuTo::ElementBase::InterpolateTemperatureFrom2D] 2D temperature interpolation routine not implemented.");
}

void NuTo::ElementBase::InterpolateTemperatureFrom3D(double rLocalCoordinates[3], double& rTemperature) const
{
    throw MechanicsException("[NuTo::ElementBase::InterpolateTemperatureFrom3D] 3D temperature interpolation routine not implemented.");
}

//! @brief ... interpolate three-dimensional global nonlocal eq strain from one-dimensional local point coordinates (element coordinates system)
//! @param rLocalCoordinates ... one-dimensional local point coordinates
//! @param rNonlocalEqStrain ... interpolated nonlocal eq strain
void NuTo::ElementBase::InterpolateNonlocalEqStrainFrom1D(double rLocalCoordinates, double& rNonlocalEqStrain) const
{
    throw MechanicsException("[NuTo::ElementBase::InterpolateNonlocalEqStrainFrom1D] 1D nonlocal eq strain interpolation routine not implemented.");
}


//! @brief ... interpolate three-dimensional global nonlocal eq strain from two-dimensional local point coordinates (element coordinates system)
//! @param rLocalCoordinates ... two-dimensional local point coordinates
//! @param rNonlocalEqStrain ... interpolated nonlocal eq strain
void NuTo::ElementBase::InterpolateNonlocalEqStrainFrom2D(double rLocalCoordinates[2], double& rNonlocalEqStrain) const
{
    throw MechanicsException("[NuTo::ElementBase::InterpolateNonlocalEqStrainFrom2D] 2D nonlocal eq strain interpolation routine not implemented.");
}

//! @brief ... interpolate three-dimensional global nonlocal eq strain from three-dimensional local point coordinates (element coordinates system)
//! @param rLocalCoordinates ... three-dimensional local point coordinates
//! @param rNonlocalEqStrain ... interpolated nonlocal eq strain
void NuTo::ElementBase::InterpolateNonlocalEqStrainFrom3D(double rLocalCoordinates[3], double& rNonlocalEqStrain) const
{
    throw MechanicsException("[NuTo::ElementBase::InterpolateNonlocalEqStrainFrom3D] 3D nonlocal eq strain interpolation routine not implemented.");
}

//! @brief ... interpolate three-dimensional global relative humidity from one-dimensional local point coordinates (element coordinates system)
//! @param rLocalCoordinates ... one-dimensional local point coordinates
//! @param rRelativeHumidity ... interpolated relative humidity
void NuTo::ElementBase::InterpolateRelativeHumidityFrom1D(double rLocalCoordinates, double& rRelativeHumidity) const
{
    throw MechanicsException("[NuTo::ElementBase::InterpolateRelativeHumidityFrom1D] 1D relative humidity interpolation routine not implemented.");
}

//! @brief ... interpolate three-dimensional global water volume fraction from one-dimensional local point coordinates (element coordinates system)
//! @param rLocalCoordinates ... one-dimensional local point coordinates
//! @param rWaterVolumeFraction ... interpolated water volume fraction
void NuTo::ElementBase::InterpolateWaterVolumeFractionFrom1D(double rLocalCoordinates, double& rWaterVolumeFraction) const
{
    throw MechanicsException("[NuTo::ElementBase::InterpolateWaterVolumeFractionFrom1D] 1D water volume fraction interpolation routine not implemented.");
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
		message << "[NuTo::ElementWithDataBase::SetIntegrationType] Integration Type does not match element type of element "<< mStructure->ElementGetId(this) <<"." <<std::endl;
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

//! @brief returns the natural coordinates of an given point
//! implemented with an exception for all elements, reimplementation required for those elements
//! @param rGlobCoords (input) ... pointer to the array of coordinates
//! @param rLocCoords (output) ... coordinates to be returned
//! @return True if coordinates are within the element, False otherwise
bool NuTo::ElementBase::GetLocalPointCoordinates(const double* rGlobCoords,  double* rLocCoords)const
{
	throw NuTo::MechanicsException("[NuTo::ElementBase::GetLocalPointCoordinates] not implemented for this element type.");
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

    if (NumVisualizationPoints == 0)
        // no visualisation
        return;

    // calculate global point coordinates and store point at the visualize object
    int dimension (VisualizationPointLocalCoordinates.size()/NumVisualizationPoints);
    assert(VisualizationPointLocalCoordinates.size() == NumVisualizationPoints * dimension);
    std::vector<unsigned int> PointIdVec;
    for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
    {
        double GlobalPointCoor[3];
        switch (dimension)
        {
        case 1:
            this->InterpolateCoordinatesFrom1D(VisualizationPointLocalCoordinates[PointCount], GlobalPointCoor);
            break;
        case 2:
            this->InterpolateCoordinatesFrom2D(&VisualizationPointLocalCoordinates[2*PointCount], GlobalPointCoor);
            break;
        case 3:
            this->InterpolateCoordinatesFrom3D(&VisualizationPointLocalCoordinates[3*PointCount], GlobalPointCoor);
            break;
        default:
            throw NuTo::MechanicsException("[NuTo::ElementBase::Visualize] invalid dimension of local coordinates");
        }
        unsigned int PointId = rVisualize.AddPoint(GlobalPointCoor);
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
                double GlobalDisplacements[3];
                switch (dimension)
                {
                case 1:
                    this->InterpolateDisplacementsFrom1D(0,VisualizationPointLocalCoordinates[PointCount], GlobalDisplacements);
                    break;
                case 2:
                    this->InterpolateDisplacementsFrom2D(0,&VisualizationPointLocalCoordinates[2*PointCount], GlobalDisplacements);
                    break;
                case 3:
                    this->InterpolateDisplacementsFrom3D(0,&VisualizationPointLocalCoordinates[3*PointCount], GlobalDisplacements);
                    break;
                default:
                    throw NuTo::MechanicsException("[NuTo::ElemenBase::Visualize] invalid dimension of local coordinates");
                }
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataVector(PointId, WhatIter->GetComponentName(), GlobalDisplacements);
            }
            break;
        case NuTo::VisualizeBase::VELOCITY:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                double GlobalVelocity[3];
                switch (dimension)
                {
                case 1:
                    this->InterpolateDisplacementsFrom1D(1,VisualizationPointLocalCoordinates[PointCount], GlobalVelocity);
                    break;
                case 2:
                    this->InterpolateDisplacementsFrom2D(1,&VisualizationPointLocalCoordinates[2*PointCount], GlobalVelocity);
                    break;
                case 3:
                    this->InterpolateDisplacementsFrom3D(1,&VisualizationPointLocalCoordinates[3*PointCount], GlobalVelocity);
                    break;
                default:
                    throw NuTo::MechanicsException("[NuTo::ElemenBase::Visualize] invalid dimension of local coordinates");
                }
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataVector(PointId, WhatIter->GetComponentName(), GlobalVelocity);
            }
            break;
        case NuTo::VisualizeBase::ACCELERATION:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                double GlobalAcceleration[3];
                switch (dimension)
                {
                case 1:
                    this->InterpolateDisplacementsFrom1D(2,VisualizationPointLocalCoordinates[PointCount], GlobalAcceleration);
                    break;
                case 2:
                    this->InterpolateDisplacementsFrom2D(2,&VisualizationPointLocalCoordinates[2*PointCount], GlobalAcceleration);
                    break;
                case 3:
                    this->InterpolateDisplacementsFrom3D(2,&VisualizationPointLocalCoordinates[3*PointCount], GlobalAcceleration);
                    break;
                default:
                    throw NuTo::MechanicsException("[NuTo::ElemenBase::Visualize] invalid dimension of local coordinates");
                }
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataVector(PointId, WhatIter->GetComponentName(), GlobalAcceleration);
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
                double nonlocalEqStrain = 0;
                if (GetSection()->GetIsNonlocalEqStrainDof())
                {
                    // calculate only if element has nonlocal eq strain dofs
                    switch (dimension)
                    {
                    case 1:
                        this->InterpolateNonlocalEqStrainFrom1D(VisualizationPointLocalCoordinates[PointCount], nonlocalEqStrain);
                        break;
                    case 2:
                        this->InterpolateNonlocalEqStrainFrom2D(&VisualizationPointLocalCoordinates[2*PointCount], nonlocalEqStrain);
                        break;
                    case 3:
                        this->InterpolateNonlocalEqStrainFrom3D(&VisualizationPointLocalCoordinates[3*PointCount], nonlocalEqStrain);
                        break;
                    default:
                        throw NuTo::MechanicsException("[NuTo::ElemenBase::Visualize] invalid dimension of local coordinates");
                    }
                }
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataScalar(PointId, WhatIter->GetComponentName(), nonlocalEqStrain);
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
                double temperature;
                switch (dimension)
                {
                case 1:
                    this->InterpolateTemperatureFrom1D(VisualizationPointLocalCoordinates[PointCount], temperature);
                    break;
                case 2:
                    this->InterpolateTemperatureFrom2D(&VisualizationPointLocalCoordinates[2*PointCount], temperature);
                    break;
                case 3:
                    this->InterpolateTemperatureFrom3D(&VisualizationPointLocalCoordinates[3*PointCount], temperature);
                    break;
                default:
                    throw NuTo::MechanicsException("[NuTo::ElemenBase::Visualize] invalid dimension of local coordinates");
                }

                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataVector(PointId, WhatIter->GetComponentName(), &temperature);
            }
            break;

        case NuTo::VisualizeBase::RELATIVE_HUMIDITY:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                double relativeHumidity = 0;
                if (GetSection()->GetIsRelativeHumidityDof())
                {
                    // calculate only if element has relative humidity dofs
                    switch (dimension)
                    {
                    case 1:
                        this->InterpolateRelativeHumidityFrom1D(VisualizationPointLocalCoordinates[PointCount], relativeHumidity);
                        break;
                    default:
                        throw NuTo::MechanicsException("[NuTo::ElemenBase::Visualize] invalid dimension of local coordinates");
                    }
                }
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataScalar(PointId, WhatIter->GetComponentName(), relativeHumidity);
            }
            break;
        case NuTo::VisualizeBase::WATER_VOLUME_FRACTION:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                double waterVolumeFraction = 0;
                if (GetSection()->GetIsWaterPhaseFractionDof())
                {
                    // calculate only if element has water volume fraction dofs
                    switch (dimension)
                    {
                    case 1:
                        this->InterpolateWaterVolumeFractionFrom1D(VisualizationPointLocalCoordinates[PointCount], waterVolumeFraction);
                        break;
                    default:
                        throw NuTo::MechanicsException("[NuTo::ElemenBase::Visualize] invalid dimension of local coordinates");
                    }
                }
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataScalar(PointId, WhatIter->GetComponentName(), waterVolumeFraction);
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
	// define variables storing the element contribution
	ElementOutputIpData elementOutputIpData(IpData::ENGINEERING_STRESS);
	NuTo::Element::eOutput keyIP_DATA(Element::IP_DATA);

	boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;
	elementOutput.insert(keyIP_DATA,&elementOutputIpData);

	this->Evaluate(elementOutput);
	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>&  rIpStress(elementOutputIpData.GetFullMatrixDouble());

    //this is certainly not the fastest approach, since the jacobian/derivates etc. is calculated twice, but it's not a critical routine
    std::vector<double> ipVolume;
    this->GetIntegrationPointVolume(ipVolume);

    rStress.Resize(rIpStress.GetNumRows(),1);
    for (int countIP=0; countIP<rIpStress.GetNumColumns(); countIP++)
    {
    	rStress+=(rIpStress.col(countIP)*(ipVolume[countIP]));
    }
}

//! @brief integrates the strain over the element
//! @param rStrain integrated strain
void NuTo::ElementBase::GetIntegratedStrain(FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rStrain)
{
	// define variables storing the element contribution
	ElementOutputIpData elementOutputIpData(IpData::ENGINEERING_STRAIN);
	NuTo::Element::eOutput keyIP_DATA(Element::IP_DATA);

	boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;
	elementOutput.insert(keyIP_DATA,&elementOutputIpData);

	this->Evaluate(elementOutput);
	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>&  rIpStrain(elementOutputIpData.GetFullMatrixDouble());

    //this is certainly not the fastest approach, since the jacobian/derivates etc. is calculated twice, but it's not a critical routine
    std::vector<double> ipVolume;
    this->GetIntegrationPointVolume(ipVolume);

    rStrain.Resize(rIpStrain.GetNumRows(),1);
    for (int countIP=0; countIP<rIpStrain.GetNumColumns(); countIP++)
    {
        rStrain+=rIpStrain.col(countIP)*(ipVolume[countIP]);
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

//! @brief cast the base pointer to a BoundaryGradientDamage1D, otherwise throws an exception
const NuTo::BoundaryGradientDamage1D* NuTo::ElementBase::AsBoundaryGradientDamage1D()const
{
	throw NuTo::MechanicsException("[NuTo::ElementBase::AsBoundaryGradientDamage1D] Element is not of type BoundaryGradientDamage1D.");
}

//! @brief cast the base pointer to a BoundaryGradientDamage1D, otherwise throws an exception
NuTo::BoundaryGradientDamage1D* NuTo::ElementBase::AsBoundaryGradientDamage1D()
{
	throw NuTo::MechanicsException("[NuTo::ElementBase::AsBoundaryGradientDamage1D] Element is not of type BoundaryGradientDamage1D.");
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

