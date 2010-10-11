
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION


#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constraints/ConstraintBase.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/ElementDataConstitutiveIp.h"
#include "nuto/mechanics/elements/ElementDataConstitutiveIpNonlocal.h"
#include "nuto/mechanics/groups/GroupBase.h"
#include "nuto/mechanics/loads/LoadBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/structures/StructureBase.h"

#include <eigen2/Eigen/QR>

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
    std::cout << "start serialize ElementBase " << mStructure->ElementGetId(this)  << " ptr " << mStructure <<std::endl;
#endif
    ar & BOOST_SERIALIZATION_NVP(mStructure)
    // the element data has to be saved on the main structure due to problems with a recursion on the stack (nonlocal data contains ptr to elements)
    // the idea is to first serialize all the elements in the table, and afterwards update the pointers of the element data in the element data routine
       & BOOST_SERIALIZATION_NVP(mElementData);
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

//! @brief sets the section of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need a section
//! @param rSection pointer to section
//! @return pointer to constitutive law
void NuTo::ElementBase::SetSection(const SectionBase* rSection)
{
    throw MechanicsException("[NuTo::ElementBase::SetSection] Section for this type of elements not required.");
}

//! @brief returns a pointer to the section of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need a section
//! @return pointer to section
const NuTo::SectionBase* NuTo::ElementBase::GetSection()const
{
    throw MechanicsException("[NuTo::ElementBase::GetSection] Section for this type of elements not required.");
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

void NuTo::ElementBase::InterpolateDisplacementsFrom1D(double rLocalCoordinates, double rGlobalDisplacements[3]) const
{
    throw MechanicsException("[NuTo::ElementBase::InterpolateDisplacementsFrom1D] 1D displacement interpolation routine not implemented.");
}

void NuTo::ElementBase::InterpolateDisplacementsFrom2D(double rLocalCoordinates[2], double rGlobalDisplacements[3]) const
{
    throw MechanicsException("[NuTo::ElementBase::InterpolateDisplacementsFrom2D] 2D displacement interpolation routine not implemented.");
}

void NuTo::ElementBase::InterpolateDisplacementsFrom3D(double rLocalCoordinates[3], double rGlobalDisplacements[3]) const
{
    throw MechanicsException("[NuTo::ElementBase::InterpolateDisplacementsFrom3D] 3D displacement interpolation routine not implemented.");
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

//! @brief returns the number of integration points
//! @return number of integration points
int NuTo::ElementBase::GetNumIntegrationPoints()const
{
    return this->mElementData->GetIntegrationType()->GetNumIntegrationPoints();
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point
//! @return weight
double NuTo::ElementBase::GetIntegrationPointWeight(int rIpNum)const
{
    return this->mElementData->GetIntegrationType()->GetIntegrationPointWeight(rIpNum);
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

#ifdef ENABLE_VISUALIZE
void NuTo::ElementBase::Visualize(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat) const
{
    // get visualization cells from integration type
    const IntegrationTypeBase* IntegrationType = this->mElementData->GetIntegrationType();
    unsigned int NumVisualizationPoints;
    std::vector<double> VisualizationPointLocalCoordinates;
    unsigned int NumVisualizationCells;
    std::vector<NuTo::CellBase::eCellTypes> VisualizationCellType;
    std::vector<unsigned int> VisualizationCellsIncidence;
    std::vector<unsigned int> VisualizationCellsIP;
    IntegrationType->GetVisualizationCells(
        NumVisualizationPoints,
        VisualizationPointLocalCoordinates,
        NumVisualizationCells,
        VisualizationCellType,
        VisualizationCellsIncidence,
        VisualizationCellsIP);

    // calculate global point coordinates and store point at the visualize object
    assert(VisualizationPointLocalCoordinates.size() == NumVisualizationPoints * IntegrationType->GetCoordinateDimension());
    std::vector<unsigned int> PointIdVec;
    for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
    {
        double GlobalPointCoor[3];
        switch (IntegrationType->GetCoordinateDimension())
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

    // store data
    boost::ptr_list<VisualizeComponentBase>::const_iterator WhatIter = rWhat.begin();
    while (WhatIter != rWhat.end())
    {
        switch (WhatIter->GetComponentEnum())
        {
        case NuTo::VisualizeBase::DAMAGE:
        {
        	FullMatrix<double> damage;
            this->GetIpData(NuTo::IpData::DAMAGE,damage);
            //this->GetCurrentIpData(IPDATA,damage);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataScalar(CellId, WhatIter->GetComponentName(), damage.mEigenMatrix.data()[theIp]);
            }
        }
        break;
        case NuTo::VisualizeBase::DISPLACEMENTS:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                double GlobalDisplacements[3];
                switch (IntegrationType->GetCoordinateDimension())
                {
                case 1:
                    this->InterpolateDisplacementsFrom1D(VisualizationPointLocalCoordinates[PointCount], GlobalDisplacements);
                    break;
                case 2:
                    this->InterpolateDisplacementsFrom2D(&VisualizationPointLocalCoordinates[2*PointCount], GlobalDisplacements);
                    break;
                case 3:
                    this->InterpolateDisplacementsFrom3D(&VisualizationPointLocalCoordinates[3*PointCount], GlobalDisplacements);
                    break;
                default:
                    throw NuTo::MechanicsException("[NuTo::ElemenBase::Visualize] invalid dimension of local coordinates");
                }
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataVector(PointId, WhatIter->GetComponentName(), GlobalDisplacements);
            }
            break;
        case NuTo::VisualizeBase::ENGINEERING_STRAIN:
        {
            FullMatrix<double> EngineeringStrains;
            this->GetIpData(NuTo::IpData::ENGINEERING_STRAIN,EngineeringStrains);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                const double* EngineeringStrainVector = &(EngineeringStrains.mEigenMatrix.data()[theIp*6]);
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
            FullMatrix<double> EngineeringStrains;
            this->GetIpData(NuTo::IpData::ENGINEERING_PLASTIC_STRAIN,EngineeringStrains);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                const double* EngineeringStrainVector = &(EngineeringStrains.mEigenMatrix.data()[theIp*6]);
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
            FullMatrix<double> EngineeringStress;
            this->GetIpData(NuTo::IpData::ENGINEERING_STRESS,EngineeringStress);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                const double* EngineeringStressVector = &(EngineeringStress.mEigenMatrix.data()[theIp*6]);
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
            FullMatrix<double> EngineeringStress;
            this->GetIpData(NuTo::IpData::ENGINEERING_STRESS,EngineeringStress);
            for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
                unsigned int theIp = VisualizationCellsIP[CellCount];
                const double* EngineeringStressVector = &(EngineeringStress.mEigenMatrix.data()[theIp*6]);
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
                rVisualize.SetCellDataVector(CellId, WhatIter->GetComponentName(), EngineeringStressTensor.part<Eigen::SelfAdjoint>().eigenvalues().data());
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
         default:
            throw NuTo::MechanicsException("[NuTo::ElementBase::Visualize] unsupported datatype for visualization.");
        }
        WhatIter++;
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

//! @brief returns a vector of the nonlocal elements
//! @param rConstitutive constitutive model for the nonlocale elements
//! @retrun nonlocal elements
int NuTo::ElementBase::GetNumNonlocalElements()const
{
	return this->mElementData->GetNumNonlocalElements();
}

//! @brief integrates the stress over the element
//! @param rStress integrated stress
void NuTo::ElementBase::GetIntegratedStress(FullMatrix<double>& rStress)const
{
    NuTo::FullMatrix<double> rIpStress;
    this->GetIpData(NuTo::IpData::ENGINEERING_STRESS, rIpStress);
    //this is certainly not the fastest approach, since the jacobian/derivates etc. is calculated twice, but it's not a critical routine
    std::vector<double> ipVolume;
    this->GetIntegrationPointVolume(ipVolume);

    rStress.Resize(rIpStress.GetNumRows(),1);
    for (int countIP=0; countIP<rIpStress.GetNumColumns(); countIP++)
    {
        rStress+=rIpStress.GetColumn(countIP)*(ipVolume[countIP]);
    }
}

//! @brief integrates the strain over the element
//! @param rStrain integrated strain
void NuTo::ElementBase::GetIntegratedStrain(FullMatrix<double>& rStrain)const
{
    NuTo::FullMatrix<double> rIpStrain;
    this->GetIpData(NuTo::IpData::ENGINEERING_STRAIN, rIpStrain);
    //this is certainly not the fastest approach, since the jacobian/derivates etc. is calculated twice, but it's not a critical routine
    std::vector<double> ipVolume;
    this->GetIntegrationPointVolume(ipVolume);

    rStrain.Resize(rIpStrain.GetNumRows(),1);
    for (int countIP=0; countIP<rIpStrain.GetNumColumns(); countIP++)
    {
        rStrain+=rIpStrain.GetColumn(countIP)*(ipVolume[countIP]);
    }
}

//! @brief cast the base pointer to a Plane, otherwise throws an exception
const NuTo::Plane* NuTo::ElementBase::AsPlane()const
{
	throw NuTo::MechanicsException("[NuTo::ElementBase::AsElementPlane] Element is not of type ElementPlane.");
}

//! @brief cast the base pointer to a Plane, otherwise throws an exception
NuTo::Plane* NuTo::ElementBase::AsPlane()
{
	throw NuTo::MechanicsException("[NuTo::ElementBase::AsElementPlane] Element is not of type ElementPlane.");
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
    mElementData = rElementData;
}

