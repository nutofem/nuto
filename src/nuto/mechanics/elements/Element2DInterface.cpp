#include "nuto/mechanics/elements/Element2DInterface.h"
#include "nuto/mechanics/constitutive/mechanics/InterfaceSlip.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixInt.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixDouble.h"
#include "nuto/mechanics/elements/ElementOutputVectorInt.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/elements/ElementOutputIpData.h"

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeBase.h"
#include "nuto/visualize/VisualizeUnstructuredGrid.h"
#include <boost/ptr_container/ptr_list.hpp>
#include <boost/assign/ptr_map_inserter.hpp>
#endif // ENABLE_VISUALIZE

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

NuTo::Element2DInterface::Element2DInterface(const NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase*>& rNodes, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType, InterpolationType* rInterpolationType) :
        NuTo::Element2D::Element2D(rStructure, rNodes, rElementDataType, rIpDataType, rInterpolationType)
{
}

NuTo::Error::eError NuTo::Element2DInterface::Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput)
{
    try
    {
        const SectionBase* section(GetSection());
        if (section == nullptr)
            throw MechanicsException("[NuTo::Element2DInterface::Evaluate] no section allocated for element.");

        const std::set<Node::eAttributes>& dofs = mInterpolationType->GetDofs();
        const std::set<Node::eAttributes>& activeDofs = mInterpolationType->GetActiveDofs();
        unsigned int numActiveDofs = mInterpolationType->GetNumActiveDofs();
        // extract all node values and store them
        std::map<Node::eAttributes, Eigen::MatrixXd> nodalValues;
        for (auto dof : dofs)
        {
            nodalValues[dof] = ExtractNodeValues(0, dof);
        }

        //****************************************//
        //    CONSTITUTIVE INPUT DECLARATION      //
        //****************************************//

        std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*> constitutiveInputList;

        InterfaceSlip interfaceSlip;

        //****************************************//
        //    CONSTITUTIVE OUTPUT DECLARATION     //
        //****************************************//

        std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*> constitutiveOutputList;

        ConstitutiveTangentLocal<2, 2> constitutiveMatrix;
        ConstitutiveTangentLocal<2, 1> interfaceStresses;

        //****************************************//
        //     FILL CONSTITUTIVE INPUT LIST       //
        //****************************************//

        for (auto dof : dofs)
        {
            if (mInterpolationType->IsConstitutiveInput(dof) == false)
                continue;
            switch (dof)
            {
            case Node::DISPLACEMENTS:
            {
                constitutiveInputList[NuTo::Constitutive::Input::INTERFACE_SLIP] = &(interfaceSlip);
            }
                break;
            default:
                throw MechanicsException("[NuTo::Element2DInterface::Evaluate] Constitutive input for " + Node::AttributeToString(dof) + " not implemented.");
            }
        }

        //****************************************//
        //     FILL CONSTITUTIVE OUTPUT LIST       //
        //****************************************//

        for (auto it = rElementOutput.begin(); it != rElementOutput.end(); it++)
        {
            switch (it->first)
            {
            case Element::INTERNAL_GRADIENT:
                it->second->GetFullVectorDouble().Resize(numActiveDofs);
                constitutiveOutputList[NuTo::Constitutive::Output::INTERFACE_STRESSES] = &(interfaceStresses);
                break;
            case Element::HESSIAN_0_TIME_DERIVATIVE:
            {
                it->second->GetFullMatrixDouble().Resize(numActiveDofs, numActiveDofs);
                it->second->GetFullMatrixDouble().setZero();
                it->second->SetSymmetry(true);
                it->second->SetConstant(true);
                for (auto dof : activeDofs)
                {
                    switch (dof)
                    {
                    case Node::DISPLACEMENTS:
                    {

                        constitutiveOutputList[NuTo::Constitutive::Output::INTERFACE_CONSTITUTIVE_MATRIX] = &(constitutiveMatrix);
                    }
                        break;
                    default:
                        throw MechanicsException("[NuTo::Element2DInterface::Evaluate] Constitutive output HESSIAN_0_TIME_DERIVATIVE for " + Node::AttributeToString(dof) + " not implemented.");

                    }
                }
            }
                break;
            case Element::HESSIAN_1_TIME_DERIVATIVE:
                for (auto dof : activeDofs)
                {
                    switch (dof)
                    {
                    case Node::DISPLACEMENTS:
                        break;
                    default:
                        throw MechanicsException("[NuTo::Element2DInterface::Evaluate] Constitutive output HESSIAN_1_TIME_DERIVATIVE for " + Node::AttributeToString(dof) + " not implemented.");
                    }
                }
                break;
            case Element::HESSIAN_2_TIME_DERIVATIVE:
                for (auto dof : activeDofs)
                {
                    switch (dof)
                    {
                    case Node::DISPLACEMENTS:
                        break;
                    default:
                        throw MechanicsException("[NuTo::Element2DInterface::Evaluate] Constitutive output HESSIAN_2_TIME_DERIVATIVE for " + Node::AttributeToString(dof) + " not implemented.");
                    }
                }
                break;
            case Element::LUMPED_HESSIAN_2_TIME_DERIVATIVE:
                break;
            case Element::UPDATE_STATIC_DATA:
                break;
            case Element::UPDATE_TMP_STATIC_DATA:
                break;
            case Element::IP_DATA:
                switch (it->second->GetIpDataType())
                {
                case NuTo::IpData::ENGINEERING_STRAIN:
                    break;
                case NuTo::IpData::ENGINEERING_STRESS:
                    break;
                default:
                    throw MechanicsException("[NuTo::Element2DInterface::Evaluate] this ip data type is not implemented.");
                }
                break;
            case Element::GLOBAL_ROW_DOF:
            {
                const Eigen::VectorXi& globalRowDofsEigen = CalculateGlobalRowDofs();

                std::vector<int> globalRowDofsStd(globalRowDofsEigen.data(), globalRowDofsEigen.data() + globalRowDofsEigen.rows());
                it->second->GetVectorInt() = globalRowDofsStd;
            }
                break;
            case Element::GLOBAL_COLUMN_DOF:
            {
                const Eigen::VectorXi& globalColumnDofsEigen = CalculateGlobalColumnDofs();
                std::vector<int> globalColumnDofsStd(globalColumnDofsEigen.data(), globalColumnDofsEigen.data() + globalColumnDofsEigen.rows());
                it->second->GetVectorInt() = globalColumnDofsStd;
            }
                break;
            default:
                throw MechanicsException("[NuTo::Element2DInterface::Evaluate] element output not implemented.");
            }
        }

        //****************************************//
        //     CALCULATE CONSTITUTIVE INPUTS      //
        //****************************************//

        double detJacobian;
        std::map<Node::eAttributes, Eigen::VectorXd> shapeFunctions;
        std::map<Node::eAttributes, Eigen::MatrixXd> BMatrix;

        // loop over the integration points
        for (int theIP = 0; theIP < GetNumIntegrationPoints(); theIP++)
        {

            detJacobian = CalculateDetJacobian(mInterpolationType->Get(Node::eAttributes::COORDINATES).GetDerivativeShapeFunctionsNatural(theIP), nodalValues.at(Node::eAttributes::COORDINATES));

            double gaussIntegrationFactor = mElementData->GetIntegrationPointWeight(theIP) * detJacobian;

            // calculate shape functions and their derivatives
            for (auto dof : dofs)
            {
                if (dof == Node::COORDINATES)
                    continue;

                shapeFunctions[dof] = mInterpolationType->Get(dof).GetShapeFunctions(theIP);
                BMatrix[dof] = mInterpolationType->Get(dof).GetDerivativeShapeFunctionsNatural(theIP);
            }

            // define constitutive inputs
            for (auto dof : dofs)
            {
                if (mInterpolationType->IsConstitutiveInput(dof) == false)
                    continue;
                switch (dof)
                {
                case Node::DISPLACEMENTS:
                {
                    interfaceSlip = CalculateInterfaceSlip(shapeFunctions.at(dof), nodalValues.at(dof));
                }
                    break;
                default:
                    throw MechanicsException("[NuTo::Element2DInterface::Evaluate] Constitutive input for " + Node::AttributeToString(dof) + " not implemented.");
                }
            }

            //****************************************//
            //      EVALUATE CONSTITUTIVE LAW         //
            //****************************************//

            ConstitutiveBase* constitutivePtr = GetConstitutiveLaw(theIP);

            try
            {
                Error::eError error = constitutivePtr->Evaluate2D(this, theIP, constitutiveInputList, constitutiveOutputList);
                if (error != Error::SUCCESSFUL)
                    return error;
            } catch (NuTo::MechanicsException &e)
            {
                e.AddMessage("[NuTo::Element2DInterfaceInterface::Evaluate] error evaluating the constitutive model.");
                throw e;
            }

            //****************************************//
            //      CALCULATE OUTPUT                  //
            //****************************************//

            for (auto it = rElementOutput.begin(); it != rElementOutput.end(); it++)
            {
                switch (it->first)
                {
                case Element::INTERNAL_GRADIENT:
                {
                    for (auto dof : activeDofs)
                    {
                        switch (dof)
                        {
                        case Node::DISPLACEMENTS:
                        {
                            AddInternalForceVector(interfaceStresses, BMatrix.at(Node::DISPLACEMENTS), it->second->GetFullVectorDouble(), gaussIntegrationFactor);
                        }
                            break;

                        default:
                            throw MechanicsException("[NuTo::Element2DInterfaceInterface::Evaluate] Element output INTERNAL_GRADIENT for " + Node::AttributeToString(dof) + " not implemented.");

                        }
                    }

                }
                    break;
                case Element::HESSIAN_0_TIME_DERIVATIVE:
                {

                    for (auto dof : activeDofs)
                    {
                        switch (dof)
                        {
                        case Node::DISPLACEMENTS:
                        {
                            AddElementStiffnessMatrix(constitutiveMatrix, BMatrix.at(Node::DISPLACEMENTS), it->second->GetFullMatrixDouble(), gaussIntegrationFactor);
                        }
                            break;

                        default:
                            throw MechanicsException("[NuTo::Element2DInterfaceInterface::Evaluate] Element output HESSIAN_0_TIME_DERIVATIVE for " + Node::AttributeToString(dof) + " not implemented.");

                        }
                    }
                }
                    break;
                case Element::HESSIAN_1_TIME_DERIVATIVE:
                case Element::HESSIAN_2_TIME_DERIVATIVE:
                case Element::LUMPED_HESSIAN_2_TIME_DERIVATIVE:
                case Element::UPDATE_STATIC_DATA:
                case Element::UPDATE_TMP_STATIC_DATA:
                case Element::IP_DATA:
                    break;
                case Element::GLOBAL_ROW_DOF:
                {
                    const Eigen::VectorXi& globalRowDofsEigen = CalculateGlobalRowDofs();
                    std::vector<int> globalRowDofsStd(globalRowDofsEigen.data(), globalRowDofsEigen.data() + globalRowDofsEigen.rows());
                    it->second->GetVectorInt() = globalRowDofsStd;
                }
                    break;
                case Element::GLOBAL_COLUMN_DOF:
                {
                    const Eigen::VectorXi& globalColumnDofsEigen = CalculateGlobalColumnDofs();
                    std::vector<int> globalColumnDofsStd(globalColumnDofsEigen.data(), globalColumnDofsEigen.data() + globalColumnDofsEigen.rows());
                    it->second->GetVectorInt() = globalColumnDofsStd;
                }
                    break;

                default:
                    throw MechanicsException("[NuTo::Element2DInterface::Evaluate] element output not implemented.");
                }
            }
        }            // loop over ip
    } catch (NuTo::MechanicsException& e)
    {
        std::stringstream ss;
        ss << mStructure->ElementGetId(this);
        e.AddMessage("[NuTo::Element2DInterface::Evaluate] Error evaluating element data of element " + ss.str() + ".");
        throw e;
    }

    return Error::SUCCESSFUL;

}

NuTo::Element::eElementType NuTo::Element2DInterface::GetEnumType() const
{
    return NuTo::Element::ELEMENT2DINTERFACE;
}

int NuTo::Element2DInterface::GetLocalDimension() const
{
    return 2;
}

const NuTo::InterfaceSlip NuTo::Element2DInterface::CalculateInterfaceSlip(const Eigen::MatrixXd& rShapeFunctions, const Eigen::MatrixXd& rNodeDisplacements)
{
    assert(rNodeDisplacements.rows() == 2 and "2d interface elements needs 2d coordinates");
    assert(rNodeDisplacements.cols() == 4 and "2d interface elements needs 4 nodes");
    assert(rShapeFunctions.rows() == 4 and "rShapeFunctions needs 4 interpolation functions");
    assert(rShapeFunctions.cols() == 1 and "rShapeFunctions needs 1d interpolation functions");

    InterfaceSlip interfaceSlip;

    Eigen::VectorXd interfaceSlipEigen(2, 1);

    const double N00 = rShapeFunctions(0, 0);
    const double N01 = rShapeFunctions(1, 0);

    interfaceSlipEigen(0, 0) = -N00 * rNodeDisplacements(0, 0) - N01 * rNodeDisplacements(0, 1) + N01 * rNodeDisplacements(0, 2) + N00 * rNodeDisplacements(0, 3);
    interfaceSlipEigen(1, 0) = -N00 * rNodeDisplacements(1, 0) - N01 * rNodeDisplacements(1, 1) + N01 * rNodeDisplacements(1, 2) + N00 * rNodeDisplacements(1, 3);

    Eigen::MatrixXd rotationMatrix = CalculateRotationMatrix();
    interfaceSlipEigen = rotationMatrix * interfaceSlipEigen;
    interfaceSlip.SetInterfaceSlip(interfaceSlipEigen);
    return interfaceSlip;
}

Eigen::MatrixXd NuTo::Element2DInterface::CalculateRotationMatrix()
{
    const Eigen::MatrixXd globalNodeCoordinates = ExtractNodeValues(0, Node::COORDINATES);

    unsigned int globalDimension = globalNodeCoordinates.rows();

    assert((globalDimension == 2 or globalDimension == 3) and "Need 2d or 3d coordinates");

    // the rotation matrix consists of three basis vectors that define the new coordinate system.
    Eigen::Vector3d basisVector00 = Eigen::Vector3d::Zero();
    Eigen::Vector3d basisVector01 = Eigen::Vector3d::Zero();
    Eigen::Vector3d basisVector02 = Eigen::Vector3d::Zero();

    // basisVector00 is aligned with the truss. Note: Trusses have to be straight!
    basisVector00.block(0, 0, globalDimension, 1) = globalNodeCoordinates.col(1) - globalNodeCoordinates.col(0);

    // check if basisVector00 is linear independent to unitVectorZ and calculate the cross product to determine basisVector01
    if (basisVector00.at(0, 0) > 1e-8 or basisVector00.at(1, 0) > 1e-8)
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

void NuTo::Element2DInterface::AddInternalForceVector(const ConstitutiveTangentLocal<2, 1>& rInterfaceStresses, const Eigen::MatrixXd& rBMatrix, NuTo::FullVector<double, Eigen::Dynamic>& rInternalForceVector, const double rGaussIntegrationFactor)
{
    const unsigned int globalDimension = GetStructure()->GetDimension();
    const unsigned int numberOfNodes = mInterpolationType->Get(Node::DISPLACEMENTS).GetNumNodes();
    const Eigen::MatrixXd transformationMatrix = CalculateTransformationMatrix(globalDimension, numberOfNodes);

    rInternalForceVector += transformationMatrix * rBMatrix.transpose() * rInterfaceStresses * rGaussIntegrationFactor;

}

void NuTo::Element2DInterface::AddElementStiffnessMatrix(const ConstitutiveTangentLocal<2, 2>& rConstitutiveMatrix, const Eigen::MatrixXd& rBMatrix, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rElementStiffnessMatrix, const double rGaussIntegrationFactor)
{
    const unsigned int globalDimension = GetStructure()->GetDimension();
    const unsigned int numberOfNodes = mInterpolationType->Get(Node::DISPLACEMENTS).GetNumNodes();
    const Eigen::MatrixXd transformationMatrix = CalculateTransformationMatrix(globalDimension, numberOfNodes);

    rElementStiffnessMatrix += transformationMatrix * rBMatrix.transpose() * rConstitutiveMatrix * rBMatrix * rGaussIntegrationFactor * transformationMatrix.transpose();

}

void NuTo::Element2DInterface::CheckElement()
{

    Eigen::MatrixXd elementCoords = ExtractNodeValues(0, Node::eAttributes::COORDINATES);
    assert(elementCoords.cols()==4 and "Only implemented for 4 nodes");
    if ( (elementCoords.col(0) - elementCoords.col(3)).norm() > 1.e-6 or (elementCoords.col(1) - elementCoords.col(2)).norm() > 1.e-6 )
        std::cout << "[NuTo::Element2DInterface::CheckElement] Node 0 and node 3 should coincide. Node 1 and node 2 should coincide. Make sure you know what you are doing." << std::endl;
}

double NuTo::Element2DInterface::CalculateDetJacobian(const Eigen::MatrixXd& rDerivativeShapeFunctions, const Eigen::MatrixXd& rNodeCoordinates) const
{
    return 0.5 * (rNodeCoordinates.col(0) - rNodeCoordinates.col(1)).norm();
}

#ifdef ENABLE_VISUALIZE

void NuTo::Element2DInterface::GetVisualizationCells(unsigned int& NumVisualizationPoints, std::vector<double>& VisualizationPointLocalCoordinates, unsigned int& NumVisualizationCells, std::vector<NuTo::CellBase::eCellTypes>& VisualizationCellType,
        std::vector<unsigned int>& VisualizationCellsIncidence, std::vector<unsigned int>& VisualizationCellsIP) const
{
    const IntegrationTypeBase* integrationType = this->mElementData->GetIntegrationType();
    if (integrationType == nullptr or integrationType->GetNumIntegrationPoints() != 2)
        throw MechanicsException("[NuTo::Element2DInterface::GetVisualizationCells] Integration type for this element is not IntegrationType1D2NGauss2Ip.");

    //  3   4   5
    //  x---x---x
    //  | 0 | 1 |   Visualization Cell
    //  x---x---x
    //  0   1   2

    NumVisualizationPoints = 4;
    // Point 0
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 1
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 2
    VisualizationPointLocalCoordinates.push_back(+1);
    VisualizationPointLocalCoordinates.push_back(+1);

    // Point 3
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(+1);

    NumVisualizationCells = 1;

    // cell 0
    VisualizationCellType.push_back(NuTo::CellBase::QUAD);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIP.push_back(0);

}

void NuTo::Element2DInterface::Visualize(VisualizeUnstructuredGrid& rVisualize, const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat)
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
        assert(dimension == 1 or dimension == 2 or dimension == 3);
        Eigen::Vector3d GlobalPointCoor = Eigen::Vector3d::Zero();

        GlobalPointCoor.head<2>() = ExtractNodeValues(0, Node::eAttributes::COORDINATES).col(PointCount);

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
            throw NuTo::MechanicsException("[NuTo::Element2DInterface::Visualize] unsupported visualization cell type");
        }
    }

    //determine the ipdata and determine the map
    boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;
    boost::ptr_list<VisualizeComponentBase>::const_iterator WhatIter = rWhat.begin();
    for (auto WhatIter = rWhat.begin(); WhatIter != rWhat.end(); WhatIter++)
    {
        switch (WhatIter->GetComponentEnum())
        {
        case NuTo::VisualizeBase::ENGINEERING_STRESS:
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
    for (auto itElementOutput = elementOutput.begin(); itElementOutput != elementOutput.end(); itElementOutput++)
    {
        switch (itElementOutput->second->GetIpDataType())
        {
        case NuTo::IpData::ENGINEERING_STRESS:
            break;
        default:
            throw MechanicsException("[NuTo::Element2DInterface::Visualize] other ipdatatypes not supported.");
        }
    }

    // store data
    for (auto WhatIter = rWhat.begin(); WhatIter != rWhat.end(); WhatIter++)
    {
        switch (WhatIter->GetComponentEnum())
        {
        case NuTo::VisualizeBase::DISPLACEMENTS:
            for (unsigned int PointCount = 0; PointCount < NumVisualizationPoints; PointCount++)
            {
                Eigen::Vector3d GlobalDisplacements = Eigen::Vector3d::Zero();
                GlobalDisplacements.head<2>() = ExtractNodeValues(0, Node::eAttributes::DISPLACEMENTS).col(PointCount);
                unsigned int PointId = PointIdVec[PointCount];
                rVisualize.SetPointDataVector(PointId, WhatIter->GetComponentName(), GlobalDisplacements.data());
            }
            break;
        case NuTo::VisualizeBase::ENGINEERING_STRESS:
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
        default:
            std::cout << WhatIter->GetComponentEnum() << "\n";
            throw NuTo::MechanicsException("[NuTo::Element2DInterface::Visualize] unsupported datatype for visualization.");
        }
    }
}
#endif // ENABLE_VISUALIZE

