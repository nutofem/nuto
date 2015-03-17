
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/MechanicsException.h"

#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/moistureTransport/RelativeHumidity.h"
#include "nuto/mechanics/constitutive/moistureTransport/WaterPhaseFraction.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/elements/ElementOutputBase.h"
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/sections/SectionTruss.h"

#include "nuto/mechanics/elements/BoundaryMoistureTransport1D.h"




//! @brief constructor
NuTo::BoundaryMoistureTransport1D::BoundaryMoistureTransport1D(const StructureBase* rStructure,
            Truss* rRealBoundaryElement,
            int rSurfaceEdge,
            BoundaryCondition::eType rBoundaryConditionType,
            ElementData::eElementDataType rElementDataType,
            IntegrationType::eIntegrationType rIntegrationType,
            IpData::eIpDataType rIpDataType
            ) :
            NuTo::ElementBase::ElementBase(rStructure, rElementDataType, rIntegrationType, rIpDataType)
{
    mSurfaceEdge = rSurfaceEdge;
    mRealBoundaryElement = rRealBoundaryElement;
    mBoundaryConditionType = rBoundaryConditionType;
}

NuTo::BoundaryMoistureTransport1D::BoundaryMoistureTransport1D(const StructureBase* rStructure,
            Truss* rRealBoundaryElement,
            NodeBase* rSurfaceNode,
            BoundaryCondition::eType rBoundaryConditionType,
            ElementData::eElementDataType rElementDataType,
            IntegrationType::eIntegrationType rIntegrationType,
            IpData::eIpDataType rIpDataType
            ) :
            NuTo::ElementBase::ElementBase(rStructure, rElementDataType, rIntegrationType, rIpDataType)
{
    mRealBoundaryElement = rRealBoundaryElement;
    mSurfaceEdge = CalculateSurfaceEdge(rSurfaceNode);
    mBoundaryConditionType = rBoundaryConditionType;
}


int NuTo::BoundaryMoistureTransport1D::CalculateSurfaceEdge(const NodeBase* rNode) const
{
    // check edge 0:
    std::vector<const NodeBase*> node(1);
    mRealBoundaryElement->GetSurfaceNodes(0, node);
    if (node[0] == rNode)
        return 0;

    // check edge 1:
    mRealBoundaryElement->GetSurfaceNodes(1, node);
    if (node[0] == rNode)
        return 1;

    throw MechanicsException("[NuTo::BoundaryMoistureTransport1D::CalculateSurfaceEdge] The element does not contain this boundary node");
}

//! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
//! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
void NuTo::BoundaryMoistureTransport1D::ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
{
    NodeBase* node = GetNode(0);
    if (node == rOldPtr)
    {
        node = rNewPtr;
    }
}


//! @brief Allocates static data for an integration point of an element
//! @param rConstitutiveLaw constitutive law, which is called to allocate the static data object
NuTo::ConstitutiveStaticDataBase* NuTo::BoundaryMoistureTransport1D::AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw)const
{
    return rConstitutiveLaw->AllocateStaticDataEngineeringStress_EngineeringStrain1D(this);
}

//! @brief calculates output data fo the elmement
//! @param eOutput ... coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
//!                    @param updateStaticData (with DummyOutput), IPData, globalrow/column dofs etc.
NuTo::Error::eError NuTo::BoundaryMoistureTransport1D::Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput)
{

    if (mStructure->GetHessianConstant(1)==false)
        throw MechanicsException("[NuTo::BoundaryMoistureTransport1D::Evaluate] only implemented for a constant Hessian for the first derivative (damping).");
    if (mStructure->GetHessianConstant(2)==false)
        throw MechanicsException("[NuTo::BoundaryMoistureTransport1D::Evaluate] only implemented for a constant Hessian for the second derivative (mass).");

    try
    {
        const SectionBase* sectionReal(mRealBoundaryElement->GetSection());
        if (sectionReal==0)
        {
            throw MechanicsException("[NuTo::BoundaryMoistureTransport1D::Evaluate] no section allocated for real boundary element.");
        }
        std::vector<double> localNodeCoordReal(mRealBoundaryElement->GetNumNodesGeometry());
        std::vector<double> nodalRelativeHumidity(mRealBoundaryElement->GetNumShapeFunctionsMoistureTransport());
        std::vector<double> nodalWaterPhaseFraction(mRealBoundaryElement->GetNumShapeFunctionsMoistureTransport());


        mRealBoundaryElement->CalculateLocalCoordinates(localNodeCoordReal);

        int numWaterPhaseFractionReal = mRealBoundaryElement->GetNumNodesField();
        int numWaterPhaseFractionDofsReal = sectionReal->GetIsWaterPhaseFractionDof() ? numWaterPhaseFractionReal : 0;

        int numRelativeHumidityReal = mRealBoundaryElement->GetNumNodesField();
        int numRelativeHumidityDofsReal = sectionReal->GetIsRelativeHumidityDof() ? numRelativeHumidityReal : 0;

        if (numWaterPhaseFractionDofsReal>0 || sectionReal->GetInputConstitutiveIsWaterPhaseFraction())
        {
            nodalWaterPhaseFraction.resize(numWaterPhaseFractionReal);
            mRealBoundaryElement->CalculateNodalWaterPhaseFraction(0,nodalWaterPhaseFraction);
        }
        if (numRelativeHumidityDofsReal>0 || sectionReal->GetInputConstitutiveIsRelativeHumidity())
        {
            nodalRelativeHumidity.resize(numRelativeHumidityReal);
            mRealBoundaryElement->CalculateNodalRelativeHumidity(0,nodalRelativeHumidity);
        }

        //allocate relative humidity
        RelativeHumidity relativeHumidity;

        //allocate water phase fraction
        WaterPhaseFraction waterPhaseFraction;

        ConstitutiveTangentLocal<1,1> tangentSurfaceMoistureTransportCoefficient;


        //allocate space for local shape functions
        std::vector<double> shapeFunctionsGeometry(mRealBoundaryElement->GetNumNodesGeometry());
        std::vector<double> shapeFunctionsMoistureTransport(mRealBoundaryElement->GetNumShapeFunctionsMoistureTransport());

        //allocate space for local derivative shape functions
        std::vector<double> derivativeShapeFunctionsGeometryNatural(mRealBoundaryElement->GetLocalDimension()*mRealBoundaryElement->GetNumNodesGeometry());
        std::vector<double> derivativeShapeFunctionsGeometryLocal(mRealBoundaryElement->GetLocalDimension()*mRealBoundaryElement->GetNumNodesGeometry());
        std::vector<double> derivativeShapeFunctionsMoistureTransportNatural(mRealBoundaryElement->GetLocalDimension()*mRealBoundaryElement->GetNumShapeFunctionsMoistureTransport());
        std::vector<double> derivativeShapeFunctionsMoistureTransportLocal(mRealBoundaryElement->GetLocalDimension()*mRealBoundaryElement->GetNumShapeFunctionsMoistureTransport());

        //define inputs and outputs
        std::map< NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase* > constitutiveInputList;
        std::map< NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase* > constitutiveOutputList;

        if (numWaterPhaseFractionDofsReal>0)
        {
            constitutiveInputList[NuTo::Constitutive::Input::WATER_PHASE_FRACTION] = &waterPhaseFraction;
        }
        if (numRelativeHumidityDofsReal>0)
        {
            constitutiveInputList[NuTo::Constitutive::Input::RELATIVE_HUMIDITY] = &relativeHumidity;
        }

        int numDofsReal = numRelativeHumidityDofsReal
                        + numWaterPhaseFractionDofsReal;


        //define outputs
        for (auto it = rElementOutput.begin(); it!=rElementOutput.end(); it++)
        {
            switch(it->first)
            {
            case Element::INTERNAL_GRADIENT:
            {
                break;
            }
            case Element::HESSIAN_0_TIME_DERIVATIVE:
            {
                it->second->GetFullMatrixDouble().Resize(numDofsReal, numDofsReal);
                it->second->SetSymmetry(true);
                it->second->SetConstant(true);
                if (numRelativeHumidityDofsReal || numWaterPhaseFractionDofsReal)
                {
                   if(numRelativeHumidityDofsReal != numWaterPhaseFractionDofsReal)
                   {
                        throw MechanicsException(std::string("[NuTo::BoundaryMoistureTransport1D::Evaluate] Number of relative humidity dofs must be equal to the number of water phase fraction dofs. Current number of dofs:\n")+
                                                 std::string("Relative humidity dofs: ")+std::to_string(numRelativeHumidityDofsReal)+
                                                 std::string("\nWater phase fraction dofs: ")+std::to_string(numWaterPhaseFractionDofsReal) + std::string("\n"));
                    }
                    constitutiveOutputList[NuTo::Constitutive::Output::BOUNDARY_SURFACE_MOISTURE_TRANSPORT_COEFFICIENT]  = &tangentSurfaceMoistureTransportCoefficient;
                }
            break;
            }

            case Element::HESSIAN_1_TIME_DERIVATIVE:
            {
                it->second->GetFullMatrixDouble().Resize(numDofsReal, numDofsReal);
                it->second->SetSymmetry(true);
                it->second->SetConstant(true);
                if (numRelativeHumidityDofsReal || numWaterPhaseFractionDofsReal)
                {
                    if(numRelativeHumidityDofsReal != numWaterPhaseFractionDofsReal)
                    {
                        throw MechanicsException(std::string("[NuTo::BoundaryMoistureTransport1D::Evaluate] Number of relative humidity dofs must be equal to the number of water phase fraction dofs. Current number of dofs:\n")+
                                                 std::string("Relative humidity dofs: ")+std::to_string(numRelativeHumidityDofsReal)+
                                                 std::string("\nWater phase fraction dofs: ")+std::to_string(numWaterPhaseFractionDofsReal) + std::string("\n"));
                    }
                }
                break;
            }
            case Element::HESSIAN_2_TIME_DERIVATIVE:
            {
                it->second->GetFullMatrixDouble().Resize(numDofsReal, numDofsReal);
                it->second->SetSymmetry(true);
                it->second->SetConstant(true);
                break;
            }
            case Element::UPDATE_STATIC_DATA:
                constitutiveOutputList[NuTo::Constitutive::Output::UPDATE_STATIC_DATA] = 0;
            break;
            case Element::GLOBAL_ROW_DOF:
                this->CalculateGlobalRowDofs(it->second->GetVectorInt(),
                        numRelativeHumidityDofsReal,
                        numWaterPhaseFractionDofsReal);
            break;
            case Element::GLOBAL_COLUMN_DOF:
                this->CalculateGlobalRowDofs(it->second->GetVectorInt(),
                        numRelativeHumidityDofsReal,
                        numWaterPhaseFractionDofsReal);
            break;
            default:
                throw MechanicsException("[NuTo::BoundaryMoistureTransport1D::Evaluate] element output not implemented.");
            }
        }
        // element has only one IP on the boundary

        //the ip is located on the boundary of the real element
        std::vector<double> naturalSurfaceCoordinates(1);

        mRealBoundaryElement->CalculateNaturalSurfaceCoordinates(mSurfaceEdge, 0, naturalSurfaceCoordinates);
        //allocate space for local ip coordinates
        double naturalIPCoordinate = naturalSurfaceCoordinates[0];

        mRealBoundaryElement->CalculateDerivativeShapeFunctionsGeometry(naturalIPCoordinate, derivativeShapeFunctionsGeometryNatural);

        //determinant of the Jacobian
        double detJReal = mRealBoundaryElement->DetJacobian(derivativeShapeFunctionsGeometryNatural,localNodeCoordReal);

        if (numRelativeHumidityDofsReal || numWaterPhaseFractionDofsReal)
        {
            if(numRelativeHumidityDofsReal != numWaterPhaseFractionDofsReal)
            {
                throw MechanicsException(std::string("[NuTo::Truss::Evaluate] Number of relative humidity dofs must be equal to the number of water phase fraction dofs. Current number of dofs:\n")+
                                         std::string("Relative humidity dofs: ")+std::to_string(numRelativeHumidityDofsReal)+
                                         std::string("\nWater phase fraction dofs: ")+std::to_string(numWaterPhaseFractionDofsReal) + std::string("\n"));
            }
            mRealBoundaryElement->CalculateShapeFunctionsMoistureTransport(naturalIPCoordinate,shapeFunctionsMoistureTransport);
            mRealBoundaryElement->CalculateRelativeHumidity(shapeFunctionsMoistureTransport,nodalRelativeHumidity,relativeHumidity);
            mRealBoundaryElement->CalculateWaterPhaseFraction(shapeFunctionsMoistureTransport,nodalWaterPhaseFraction,waterPhaseFraction);
            mRealBoundaryElement->CalculateDerivativeShapeFunctionsMoistureTransport(naturalIPCoordinate,derivativeShapeFunctionsMoistureTransportNatural);
            for (unsigned int count=0; count<derivativeShapeFunctionsMoistureTransportNatural.size(); count++)
            {
                derivativeShapeFunctionsMoistureTransportLocal[count] = derivativeShapeFunctionsMoistureTransportNatural[count]/detJReal;
            }

        }

        ConstitutiveBase* constitutivePtr = GetConstitutiveLaw(0);
        Error::eError error = constitutivePtr->Evaluate1D(this, 0, constitutiveInputList, constitutiveOutputList);

        if (error!=Error::SUCCESSFUL)
            return error;

        // calculate local IP coordinates
        double localIPcoordinate = 0;
        for (unsigned int iNode = 0; iNode < localNodeCoordReal.size(); ++iNode)
        {
            localIPcoordinate += localNodeCoordReal[iNode] * shapeFunctionsGeometry[iNode];
        }
        double area = sectionReal->GetArea() * sectionReal->AsSectionTruss()->GetAreaFactor(localIPcoordinate);

        for (auto it = rElementOutput.begin(); it!=rElementOutput.end(); it++)
        {
            switch(it->first)
            {
            case Element::INTERNAL_GRADIENT:
            {
                break;
            }
            case Element::HESSIAN_0_TIME_DERIVATIVE:
            {
                if (numWaterPhaseFractionDofsReal>0 || numRelativeHumidityDofsReal>0)
                {
                    mRealBoundaryElement->AddDetJNtCN(shapeFunctionsMoistureTransport,tangentSurfaceMoistureTransportCoefficient,area,0,0,it->second->GetFullMatrixDouble());
                    mRealBoundaryElement->AddDetJNtCN(shapeFunctionsMoistureTransport,tangentSurfaceMoistureTransportCoefficient,area,2,2,it->second->GetFullMatrixDouble());
                }
                //NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> test = it->second->GetFullMatrixDouble();
                //it->second->GetFullMatrixDouble().Info();         // Check Element Matrix if needed
                //int a=0;
                break;
            }
            case Element::HESSIAN_2_TIME_DERIVATIVE:
            {
                break;
            }
            case Element::GLOBAL_ROW_DOF:
                CalculateGlobalRowDofs(it->second->GetVectorInt(),numRelativeHumidityDofsReal, numWaterPhaseFractionDofsReal);
                break;
            case Element::GLOBAL_COLUMN_DOF:
                CalculateGlobalRowDofs(it->second->GetVectorInt(),numRelativeHumidityDofsReal, numWaterPhaseFractionDofsReal);
                break;
            default:

                std::cout << "I want to calculate other stuff as well! " << std::endl;
                break;
          }
        }

    }
    catch (NuTo::MechanicsException& e)
    {
        std::stringstream ss;
        ss << mStructure->ElementGetId(this);
        e.AddMessage("[NuTo::BoundaryMoistureTransport1D::Evaluate] Error evaluating element data of element"	+ ss.str() + ".");
        throw e;
    }

    return Error::SUCCESSFUL;
}

// build global row dofs
void NuTo::BoundaryMoistureTransport1D::CalculateGlobalRowDofs(
        std::vector<int>& rGlobalRowDofs,
        int rNumRelativeHumidityDofs,
        int rNumWaterPhaseFractionDofs) const
{
    rGlobalRowDofs.resize(rNumRelativeHumidityDofs+rNumWaterPhaseFractionDofs);



    // indices for the current dof type
    int iGlobalWaterPhaseFraction      = 0;
    int iGlobalRelativeHumidity        = 0;

    int numNodes(mRealBoundaryElement->GetNumNodes());

    for (int iNode = 0; iNode < numNodes; iNode++)
    {
        const NodeBase * node = mRealBoundaryElement->GetNode(iNode);

        int iStartWaterPhaseFraction = 0;

        if (node->GetNumWaterPhaseFraction()>0 && rNumWaterPhaseFractionDofs>0)
        {
            rGlobalRowDofs[iStartWaterPhaseFraction+iGlobalWaterPhaseFraction] = node->GetDofWaterPhaseFraction();
            iGlobalWaterPhaseFraction++;
        }

        int iStartRelativeHumidity = iStartWaterPhaseFraction + rNumWaterPhaseFractionDofs;
        if (node->GetNumRelativeHumidity()>0 && rNumRelativeHumidityDofs>0)
        {
            rGlobalRowDofs[iStartRelativeHumidity+iGlobalRelativeHumidity]  = node->GetDofRelativeHumidity();
            iGlobalRelativeHumidity++;
        }
    }

}



