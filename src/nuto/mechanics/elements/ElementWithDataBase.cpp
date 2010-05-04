#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constraints/ConstraintBase.h"
#include "nuto/mechanics/elements/ElementWithDataBase.h"
#include "nuto/mechanics/elements/ElementDataConstitutiveIp.h"
#include "nuto/mechanics/groups/GroupBase.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/loads/LoadBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/structures/StructureBase.h"

#ifdef ENABLE_VISUALIZE
   #include "nuto/visualize/CellBase.h"
   #include "nuto/visualize/VisualizeComponentBase.h"
#endif // ENABLE_VISUALIZE

//! @brief constructor
NuTo::ElementWithDataBase::ElementWithDataBase(const StructureBase* rStructure, ElementData::eElementDataType rElementDataType,
		IntegrationType::eIntegrationType rIntegrationType, IpData::eIpDataType rIpDataType) :
        NuTo::ElementBase(rStructure)
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
		throw MechanicsException("[NuTo::ElementWithDataBase::ElementWithDataBase] CONSTITUTIVELAWIPNONLOCAL to be implemented.");
    break;
	default:
		throw MechanicsException("[NuTo::ElementWithDataBase::ElementWithDataBase] unsupported element data type.");
	}

	mElementData = ptrElementData;


//	printf("Allocate %p",static_cast<void*>(ptrElementData));
}

// destructor
NuTo::ElementWithDataBase::~ElementWithDataBase()
{
	//std::cout << "NuTo::ElementWithDataBase::~ElementWithDataBase()" << std::endl;
    if(this->mElementData != 0)
    {
//    	printf("Delete %p",static_cast<void*>(mElementData));
        delete this->mElementData;
    }
}

//! @brief returns a pointer to the constitutive law for an integration point
//! @param integration point number (counting from zero)
//! @return pointer to constitutive law
const NuTo::ConstitutiveBase* NuTo::ElementWithDataBase::GetConstitutiveLaw(int rIp)const
{
    return mElementData->GetConstitutiveLaw(rIp);
}


//! @brief returns a pointer to the constitutive law for an integration point
//! @param integration point number (counting from zero)
//! @return pointer to constitutive law
NuTo::ConstitutiveBase* NuTo::ElementWithDataBase::GetConstitutiveLaw(int rIp)
{
    return mElementData->GetConstitutiveLaw(rIp);
}

//! @brief sets the constitutive law for an element
//! @param rConstitutiveLaw Pointer to constitutive law entry
void NuTo::ElementWithDataBase::SetConstitutiveLaw(ConstitutiveBase* rConstitutiveLaw)
{
    //check compatibility between element type and constitutive law
	if (rConstitutiveLaw->CheckElementCompatibility(this->GetEnumType()))
	{
		//printf("check element constitutive was positive.\n");
		printf("[NuTo::ElementWithDataBase::SetConstitutiveLaw] %p\n",mElementData);

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

//! @brief sets the integration type of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need an integration type
//! @param rIntegrationType pointer to integration type
void NuTo::ElementWithDataBase::SetIntegrationType(const NuTo::IntegrationTypeBase* rIntegrationType, NuTo::IpData::eIpDataType rIpDataType)
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
const NuTo::IntegrationTypeBase* NuTo::ElementWithDataBase::GetIntegrationType()const
{
    return mElementData->GetIntegrationType();
}

//! @brief returns the number of integration points
//! @return number of integration points
int NuTo::ElementWithDataBase::GetNumIntegrationPoints()const
{
    return this->mElementData->GetIntegrationType()->GetNumIntegrationPoints();
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point
//! @return weight
double NuTo::ElementWithDataBase::GetIntegrationPointWeight(int rIpNum)const
{
    return this->mElementData->GetIntegrationType()->GetIntegrationPointWeight(rIpNum);
}

//! @brief returns the static data of an integration point
//! @param rIp integration point
//! @return static data pointer
NuTo::ConstitutiveStaticDataBase* NuTo::ElementWithDataBase::GetStaticData(int rIp)const
{
	return this->mElementData->GetStaticData(rIp);
}


//! @brief returns the static data of an integration point
//! @param rIp integration point
//! @return static data pointer
NuTo::ConstitutiveStaticDataBase* NuTo::ElementWithDataBase::GetStaticData(int rIp)
{
	return this->mElementData->GetStaticData(rIp);
}

#ifdef ENABLE_VISUALIZE
void NuTo::ElementWithDataBase::Visualize(VisualizeUnstructuredGrid& rVisualize, const std::list<NuTo::VisualizeComponentBase*>& rWhat) const
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
            throw NuTo::MechanicsException("[NuTo::ElementWithDataBase::Visualize] invalid dimension of local coordinates");
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
            throw NuTo::MechanicsException("[NuTo::ElementWithDataBase::Visualize] unsupported visualization cell type");
        }
    }

    // store data
    std::list<VisualizeComponentBase*>::const_iterator WhatIter = rWhat.begin();
    while (WhatIter != rWhat.end())
    {
        switch ((*WhatIter)->GetComponentEnum())
        {
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
                    throw NuTo::MechanicsException("[NuTo::ElementWithDataBase::Visualize] invalid dimension of local coordinates");
                }
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataVector(PointId, (*WhatIter)->GetComponentName(), GlobalDisplacements);
            }
            break;
        case NuTo::VisualizeBase::ENGINEERING_STRAIN:
        {
            FullMatrix<double> EngineeringStrains;
            this->GetEngineeringStrain(EngineeringStrains);
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
                rVisualize.SetCellDataTensor(CellId, (*WhatIter)->GetComponentName(), EngineeringStrainTensor);
            }
        }
        break;
        case NuTo::VisualizeBase::ENGINEERING_STRESS:
        {
            FullMatrix<double> EngineeringStress;
            this->GetEngineeringStress(EngineeringStress);
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

                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataTensor(CellId, (*WhatIter)->GetComponentName(), EngineeringStressTensor);
            }
        }
        break;
        /*        case NuTo::VisualizeBase::NONLOCAL_WEIGHT:
        {
            // get nonlocal weights
        	ConstitutiveBase* constitutive(0);
        	std::vector<double> Weights;
        	Weights.resize(this->GetNumIntegrationPoints());

        	// get local number within nonlocal elements for current element
        	try
        	{
            	const std::vector<const NuTo::ElementBase*> nonlocalElements(mElementData->GetNonlocalElements(constitutive));
            	std::vector<double>::iterator it = std::find(nonlocalElements.begin(),nonlocalElements.end(),mStructure->ElementGetElementPtr((*WhatIter)->GetElementId()));
            	if (it!=nonlocalElements.end())
            	{
            		// get the weights
            		Weights = GetNonlocalWeights(constitutive, (*it));
            	} // otherwise they just remain zero
        	}
        	catch(...)
        	{

        	}
        	// get nonlocal weights
        	std::vector<double> Weights(GetNonlocalWeights(constitutive, nonlocalElement));

           for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
            {
        	   unsigned int theIp = VisualizationCellsIP[CellCount];
        	   // get nonlocal weight
                unsigned int CellId = CellIdVec[CellCount];
                rVisualize.SetCellDataScalar(CellId, (*WhatIter)->GetComponentName(), EngineeringStressTensor);
            }

        }
        break;
*/        default:
            throw NuTo::MechanicsException("[NuTo::ElementWithDataBase::Visualize] unsupported datatype for visualization.");
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
void NuTo::ElementWithDataBase::AddNonlocalIp(int rLocalIpNumber, const ConstitutiveBase* rConstitutive,
		const ElementWithDataBase* rNonlocalElement, int rNonlocalIp, double rWeight)
{
	this->mElementData->AddNonlocalIp(rLocalIpNumber,rConstitutive,rNonlocalElement,rNonlocalIp,rWeight);
}
