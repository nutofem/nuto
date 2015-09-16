#include "nuto/mechanics/elements/Element2DInterface.h"
#include "nuto/mechanics/constitutive/mechanics/InterfaceSlip.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/elements/ElementOutputIpData.h"

NuTo::Element2DInterface::Element2DInterface(const NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase*>& rNodes, ElementData::eElementDataType rElementDataType,
        IpData::eIpDataType rIpDataType, InterpolationType* rInterpolationType) :
        NuTo::Element2D::Element2D(rStructure, rNodes, rElementDataType, rIpDataType, rInterpolationType)
{
}

NuTo::Error::eError NuTo::Element2DInterface::Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput)
{
    try
    {

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
                constitutiveInputList[NuTo::Constitutive::Input::INTERFACE_SLIP] = &(interfaceSlip);
                break;
            default:
                throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t Constitutive input for " + Node::AttributeToString(dof) + " not implemented.");
            }
        }

        //****************************************//
        //     FILL CONSTITUTIVE OUTPUT LIST      //
        //****************************************//

        for (auto it = rElementOutput.begin(); it != rElementOutput.end(); it++)
        {
            switch (it->first)
            {
            case Element::INTERNAL_GRADIENT:
            {
                it->second->GetFullVectorDouble().Resize(numActiveDofs);
                for (auto dof : activeDofs)
                {
                    switch (dof)
                    {
                    case Node::DISPLACEMENTS:
                        constitutiveOutputList[NuTo::Constitutive::Output::INTERFACE_STRESSES] = &(interfaceStresses);
                        break;
                    default:
                        throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t Constitutive output INTERNAL_GRADIENT for " + Node::AttributeToString(dof) + " not implemented.");
                    }
                }
            }
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
                        constitutiveOutputList[NuTo::Constitutive::Output::INTERFACE_CONSTITUTIVE_MATRIX] = &(constitutiveMatrix);
                        break;
                    default:
                        throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t Constitutive output HESSIAN_0_TIME_DERIVATIVE for " + Node::AttributeToString(dof) + " not implemented.");

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
                        throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t Constitutive output HESSIAN_1_TIME_DERIVATIVE for " + Node::AttributeToString(dof) + " not implemented.");
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
                        throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t Constitutive output HESSIAN_2_TIME_DERIVATIVE for " + Node::AttributeToString(dof) + " not implemented.");
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
                    throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t this ip data type is not implemented.");
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
                throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t element output not implemented.");
            }
        }

        //****************************************//
        //     CALCULATE CONSTITUTIVE INPUTS      //
        //****************************************//

        double detJacobian;
        std::map<Node::eAttributes, Eigen::VectorXd> shapeFunctions;

        // loop over the integration points
        for (int theIP = 0; theIP < GetNumIntegrationPoints(); theIP++)
        {

            detJacobian = CalculateDetJacobian(nodalValues.at(Node::eAttributes::COORDINATES));

            double gaussIntegrationFactor = mElementData->GetIntegrationPointWeight(theIP) * detJacobian;

            // calculate shape functions and their derivatives
            for (auto dof : dofs)
            {
                if (dof == Node::COORDINATES)
                    continue;

                shapeFunctions[dof] = mInterpolationType->Get(dof).GetShapeFunctions(theIP);
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
                    throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t Constitutive input for " + Node::AttributeToString(dof) + " not implemented.");
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
                e.AddMessage(std::string(__PRETTY_FUNCTION__) + ":\t error evaluating the constitutive model.");
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
                            AddInternalForceVector(interfaceStresses, shapeFunctions.at(Node::DISPLACEMENTS), it->second->GetFullVectorDouble(), gaussIntegrationFactor);
                        }
                            break;

                        default:
                            throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t Element output INTERNAL_GRADIENT for " + Node::AttributeToString(dof) + " not implemented.");

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
                            AddElementStiffnessMatrix(constitutiveMatrix, shapeFunctions.at(Node::DISPLACEMENTS), it->second->GetFullMatrixDouble(), gaussIntegrationFactor);
                        }
                            break;

                        default:
                            throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t Element output HESSIAN_0_TIME_DERIVATIVE for " + Node::AttributeToString(dof) + " not implemented.");

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
                    throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t element output not implemented.");
                }
            }
        }            // loop over ip
    } catch (NuTo::MechanicsException& e)
    {
        std::stringstream ss;
        ss << mStructure->ElementGetId(this);
        e.AddMessage(std::string(__PRETTY_FUNCTION__) + ":\t Error evaluating element data of element " + ss.str() + ".");
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

const NuTo::InterfaceSlip NuTo::Element2DInterface::CalculateInterfaceSlip(const Eigen::VectorXd& rShapeFunctions, const Eigen::MatrixXd& rNodeDisplacements)
{
    assert(rNodeDisplacements.rows() == 2 and "2d interface elements needs 2d coordinates");
    assert(rShapeFunctions.cols() == 1 and "rShapeFunctions needs 1d interpolation functions");

    const unsigned int globalDimension = GetStructure()->GetDimension();

    InterfaceSlip interfaceSlip;

    Eigen::VectorXd interfaceSlipEigen(globalDimension);
    interfaceSlipEigen.setZero();

    interfaceSlipEigen(0, 0) = rShapeFunctions.dot(rNodeDisplacements.row(0));
    interfaceSlipEigen(1, 0) = rShapeFunctions.dot(rNodeDisplacements.row(1));

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

void NuTo::Element2DInterface::AddInternalForceVector(const ConstitutiveTangentLocal<2, 1>& rInterfaceStresses, const Eigen::VectorXd& rShapefunctions,
        NuTo::FullVector<double, Eigen::Dynamic>& rInternalForceVector, const double rGaussIntegrationFactor)
{
    const unsigned int globalDimension = GetStructure()->GetDimension();
    const unsigned int numberOfNodes = mInterpolationType->Get(Node::DISPLACEMENTS).GetNumNodes();
    const Eigen::MatrixXd transformationMatrix = CalculateTransformationMatrix(globalDimension, numberOfNodes);

    Eigen::MatrixXd BMatrix(globalDimension, globalDimension * numberOfNodes);
    BMatrix.setZero();

    for (unsigned int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        const unsigned int index = globalDimension * iNode;
        BMatrix(0, index + 0) = rShapefunctions(iNode);
        BMatrix(1, index + 1) = rShapefunctions(iNode);
    }

    rInternalForceVector += transformationMatrix * BMatrix.transpose() * rInterfaceStresses * rGaussIntegrationFactor;

}

void NuTo::Element2DInterface::AddElementStiffnessMatrix(const ConstitutiveTangentLocal<2, 2>& rConstitutiveMatrix, const Eigen::VectorXd& rShapefunctions,
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rElementStiffnessMatrix, const double rGaussIntegrationFactor)
{
    const unsigned int globalDimension = GetStructure()->GetDimension();
    const unsigned int numberOfNodes = mInterpolationType->Get(Node::DISPLACEMENTS).GetNumNodes();
    const Eigen::MatrixXd transformationMatrix = CalculateTransformationMatrix(globalDimension, numberOfNodes);

    Eigen::MatrixXd BMatrix(globalDimension, globalDimension * numberOfNodes);
    BMatrix.setZero();

    for (unsigned int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        const unsigned int index = globalDimension * iNode;
        BMatrix(0, index + 0) = rShapefunctions(iNode);
        BMatrix(1, index + 1) = rShapefunctions(iNode);
    }

    rElementStiffnessMatrix += transformationMatrix * BMatrix.transpose() * rConstitutiveMatrix * BMatrix * rGaussIntegrationFactor * transformationMatrix.transpose();

}

void NuTo::Element2DInterface::CheckElement()
{
    Eigen::MatrixXd elementCoords = ExtractNodeValues(0, Node::eAttributes::COORDINATES);
    assert((elementCoords.cols() == 4 or elementCoords.cols() == 6) and "Only implemented for 4 or 6 nodes");
}

double NuTo::Element2DInterface::CalculateDetJacobian(const Eigen::MatrixXd& rNodeCoordinates) const
{
    // returns 0.5 * element length
    return 0.5 * (rNodeCoordinates.col(0) - rNodeCoordinates.col(0.5 * rNodeCoordinates.cols() - 1)).norm();
}

#ifdef ENABLE_VISUALIZE
void NuTo::Element2DInterface::GetVisualizationCells(unsigned int& NumVisualizationPoints, std::vector<double>& VisualizationPointLocalCoordinates, unsigned int& NumVisualizationCells,
        std::vector<NuTo::CellBase::eCellTypes>& VisualizationCellType, std::vector<unsigned int>& VisualizationCellsIncidence, std::vector<unsigned int>& VisualizationCellsIP) const
{

    const IntegrationTypeBase* integrationType = this->mElementData->GetIntegrationType();
    switch (integrationType->GetNumIntegrationPoints())
    {
    case 2:

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
        break;
    case 3:

        NumVisualizationPoints = 6;
        // Point 0
        VisualizationPointLocalCoordinates.push_back(-1);
        VisualizationPointLocalCoordinates.push_back(-1);

        // Point 1
        VisualizationPointLocalCoordinates.push_back(+0);
        VisualizationPointLocalCoordinates.push_back(-1);

        // Point 2
        VisualizationPointLocalCoordinates.push_back(+1);
        VisualizationPointLocalCoordinates.push_back(-1);

        // Point 3
        VisualizationPointLocalCoordinates.push_back(+1);
        VisualizationPointLocalCoordinates.push_back(+1);

        // Point 4
        VisualizationPointLocalCoordinates.push_back(+0);
        VisualizationPointLocalCoordinates.push_back(+1);

        // Point 5
        VisualizationPointLocalCoordinates.push_back(-1);
        VisualizationPointLocalCoordinates.push_back(+1);

        NumVisualizationCells = 2;

        // cell 0
        VisualizationCellType.push_back(NuTo::CellBase::QUAD);
        VisualizationCellsIncidence.push_back(0);
        VisualizationCellsIncidence.push_back(1);
        VisualizationCellsIncidence.push_back(4);
        VisualizationCellsIncidence.push_back(5);
        VisualizationCellsIP.push_back(0);

        // cell 0
        VisualizationCellType.push_back(NuTo::CellBase::QUAD);
        VisualizationCellsIncidence.push_back(1);
        VisualizationCellsIncidence.push_back(2);
        VisualizationCellsIncidence.push_back(3);
        VisualizationCellsIncidence.push_back(4);
        VisualizationCellsIP.push_back(1);

        break;
    default:
        throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t Integration type not valid for this element.");
    }


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
            throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t unsupported visualization cell type");
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
            throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t other ipdatatypes not supported.");
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
            throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t unsupported datatype for visualization.");
        }
    }

}
#endif // ENABLE_VISUALIZE
