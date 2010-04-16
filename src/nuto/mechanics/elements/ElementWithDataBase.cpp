#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/elements/ElementWithDataBase.h"
#include "nuto/mechanics/elements/ElementDataConstitutive.h"
#include "nuto/mechanics/elements/ElementDataConstitutiveStaticData.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/visualize/CellBase.h"

//! @brief constructor
NuTo::ElementWithDataBase::ElementWithDataBase(const StructureBase* rStructure,
		ElementDataBase::eElementDataType rElementDataType, IntegrationTypeBase::eIntegrationType rIntegrationType):
        NuTo::ElementBase(rStructure)
{
    //allocate element data
	ElementDataBase* ptrElementData;
	switch(rElementDataType)
	{
	case NuTo::ElementDataBase::CONSTITUTIVELAWELEMENT_NOSTATICDATA:
		ptrElementData = new NuTo::ElementDataConstitutive(const_cast<StructureBase*>(rStructure)->GetPtrIntegrationType(rIntegrationType));
    break;
	case NuTo::ElementDataBase::CONSTITUTIVELAWELEMENT_STATICDATA:
		ptrElementData = new NuTo::ElementDataConstitutiveStaticData(const_cast<StructureBase*>(rStructure)->GetPtrIntegrationType(rIntegrationType));
    break;
	case NuTo::ElementDataBase::CONSTITUTIVELAWIP_NOSTATICDATA:
		throw MechanicsException("[NuTo::ElementWithDataBase::ElementWithDataBase] different constitutive laws for each IP not yet implemented.");
    break;
	case NuTo::ElementDataBase::CONSTITUTIVELAWIP_STATICDATA:
		throw MechanicsException("[NuTo::ElementWithDataBase::ElementWithDataBase] different constitutive laws for each IP not yet implemented.");
    break;
	default:
		throw MechanicsException("[NuTo::ElementWithDataBase::ElementWithDataBase] unsopported element data type.");
	}

	mElementData = ptrElementData;
}

// destructor
NuTo::ElementWithDataBase::~ElementWithDataBase()
{
    if(this->mElementData != 0)
    {
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
	printf("NuTo::ElementWithDataBase::SetConstitutiveLaw\n");
	if (rConstitutiveLaw->CheckElementCompatibility(this->GetEnumType()))
	{
		printf("check element constitutive was positive.\n");
		//set material
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
void NuTo::ElementWithDataBase::SetIntegrationType(const NuTo::IntegrationTypeBase* rIntegrationType)
{
    mElementData->SetIntegrationType(this, rIntegrationType);
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
void NuTo::ElementWithDataBase::Visualize(VisualizeUnstructuredGrid& rVisualize, const std::map<std::string,NuTo::VisualizeBase::eVisualizeWhat>& rWhat) const
{
    // get vizualization cells from integration type
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
    std::map<std::string,NuTo::VisualizeBase::eVisualizeWhat>::const_iterator WhatIter = rWhat.begin();
    while (WhatIter != rWhat.end())
    {
        switch (WhatIter->second)
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
                rVisualize.SetPointDataVector(PointId, WhatIter->first, GlobalDisplacements);
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
                rVisualize.SetCellDataTensor(CellId, WhatIter->first, EngineeringStrainTensor);
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
                rVisualize.SetCellDataTensor(CellId, WhatIter->first, EngineeringStressTensor);
            }
        }
        break;
        default:
            throw NuTo::MechanicsException("[NuTo::ElementWithDataBase::Visualize] unsupported datatype for visualization.");
        }
        WhatIter++;
    }
}
#endif // ENABLE_VISUALIZE
