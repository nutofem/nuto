/*
 * Element1D.cpp
 *
 *  Created on: 8 May 2015
 *      Author: ttitsche
 */

#include "nuto/mechanics/elements/Element1D.h"
#include "nuto/mechanics/constitutive/mechanics/Damage.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/thermal/HeatFlux1D.h"
#include "nuto/mechanics/constitutive/thermal/HeatFlux3D.h"
#include "nuto/mechanics/constitutive/thermal/Temperature.h"
#include "nuto/mechanics/constitutive/thermal/TemperatureGradient1D.h"
#include "nuto/mechanics/constitutive/thermal/TemperatureGradient3D.h"
#include "nuto/mechanics/constitutive/mechanics/LocalEqStrain.h"
#include "nuto/mechanics/constitutive/mechanics/NonlocalEqStrain.h"
#include "nuto/mechanics/constitutive/moistureTransport/RelativeHumidity.h"
#include "nuto/mechanics/constitutive/moistureTransport/WaterVolumeFraction.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentNonlocal.h"
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixInt.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixDouble.h"
#include "nuto/mechanics/elements/ElementOutputVectorInt.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionTruss.h"
#include "nuto/mechanics/structures/StructureBase.h"

#include "nuto/math/FullMatrix.h"

NuTo::Element1D::Element1D(const NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase*>& rNodes, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType, InterpolationType* rInterpolationType) :
        NuTo::ElementBase::ElementBase(rStructure, rElementDataType, rIpDataType, rInterpolationType),
        mNodes(rNodes),
        mSection(0)
{
}

NuTo::Error::eError NuTo::Element1D::Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput)
{
    if (mStructure->GetHessianConstant(1) == false)
        throw MechanicsException("[NuTo::Element1D::Evaluate] only implemented for a constant Hessian for the first derivative (damping).");
    if (mStructure->GetHessianConstant(2) == false)
        throw MechanicsException("[NuTo::Element1D::Evaluate] only implemented for a constant Hessian for the second derivative (mass).");

    try
    {
        const SectionBase* section(GetSection());
        if (section == 0)
            throw MechanicsException("[NuTo::Element1D::Evaluate] no section allocated for element.");

        const std::set<Node::eDof>& dofs = mInterpolationType->GetDofs();
        const std::set<Node::eDof>& activeDofs = mInterpolationType->GetActiveDofs();

        int numActiveDofs = mInterpolationType->GetNumActiveDofs();

        // extract all node values and store them
        std::map<Node::eDof, Eigen::MatrixXd> nodalValues;
        for (auto dof : dofs)
        {
            nodalValues[dof] = ExtractNodeValues(0, dof);
        }

        double total_mass(0.);

        /*****************************************\
         *    CONSTITUTIVE INPUT DECLARATION     *
         \*****************************************/

        std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*> constitutiveInputList;

        EngineeringStress1D engineeringStress1D;
        EngineeringStress3D engineeringStress3D;
        EngineeringStrain1D engineeringStrain1D;
        EngineeringStrain3D engineeringStrain3D;

        EngineeringStrain3D engineeringPlasticStrain3D;

        //allocate deformation gradient
        DeformationGradient1D deformationGradient;

        //allocate nonlocal eq strain (gauss point value, input of constitutive relation)
        NonlocalEqStrain nonlocalEqStrain;

        //allocate nonlocal parameter
        ConstitutiveTangentLocal<1, 1> nonlocalParameter;

        //allocate local eq strain output of constitutive relation
        LocalEqStrain localEqStrain;

        //allocate damage
        Damage damage;

        // Moisture Transport
        // ------------------

        RelativeHumidity relativeHumidity;
        RelativeHumidity relativeHumidityD1;
        RelativeHumidity relativeHumidityGradient;          //! ---> should only work for 1d
        WaterVolumeFraction waterVolumeFraction;
        WaterVolumeFraction waterVolumeFractionD1;
        WaterVolumeFraction waterVolumeFractionGradient;    //! ---> should only work for 1d

        /*****************************************\
         *    CONSTITUTIVE OUTPUT DECLARATION    *
         \*****************************************/

        std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*> constitutiveOutputList;

        //allocate tangents
        ConstitutiveTangentLocal<1, 1> tangentStressStrain;
        ConstitutiveTangentLocal<1, 1> tangentLocalEqStrainStrain;
        ConstitutiveTangentLocal<1, 1> tangentStressNonlocalEqStrain;

        // Moisture Transport
        // ------------------

        // Naming scheme for matrices: tangent_D_X_D_Y_HN_AB
        // D_X_D_Y:     defines the derivative of X with respect to Y
        // HN:          defines the matrix the tangent belongs to, for example H0 for stiffness or H1 for damping
        // AB:          defines the combination of (derivative) shapefunctions the tangent has to be multiplied with, for example BN means: derivative shapefunctions * tangent * shapefunctions

        // Internal Gradient
        ConstitutiveTangentLocal<1, 1> residualWaterPhaseN;
        ConstitutiveTangentLocal<1, 1> residualWaterPhaseB;
        ConstitutiveTangentLocal<1, 1> residualVaporPhaseN;
        ConstitutiveTangentLocal<1, 1> residualVaporPhaseB;
        // Hessian 0
        ConstitutiveTangentLocal<1, 1> tangent_D_Residual_RH_D_RH_H0_BB;
        ConstitutiveTangentLocal<1, 1> tangent_D_Residual_RH_D_RH_H0_NN;
        ConstitutiveTangentLocal<1, 1> tangent_D_Residual_RH_D_WV_H0_BN;
        ConstitutiveTangentLocal<1, 1> tangent_D_Residual_RH_D_WV_H0_NN;
        ConstitutiveTangentLocal<1, 1> tangent_D_Residual_WV_D_RH_H0_NN;
        ConstitutiveTangentLocal<1, 1> tangent_D_Residual_WV_D_WV_H0_BB;
        ConstitutiveTangentLocal<1, 1> tangent_D_Residual_WV_D_WV_H0_BN;
        ConstitutiveTangentLocal<1, 1> tangent_D_Residual_WV_D_WV_H0_NN;
        // Hessian 1
        ConstitutiveTangentLocal<1, 1> tangent_D_Residual_RH_D_RH_H1_NN;
        ConstitutiveTangentLocal<1, 1> tangent_D_Residual_RH_D_WV_H1_NN;
        ConstitutiveTangentLocal<1, 1> tangent_D_Residual_WV_D_WV_H1_NN;

        /*****************************************\
         *     FILL CONSTITUTIVE INPUT LIST      *
         \*****************************************/

        for (auto dof : dofs)
        {
            if (mInterpolationType->IsConstitutiveInput(dof) == false)
                continue;
            switch (dof)
            {
            case Node::DISPLACEMENTS:
            {
                constitutiveInputList[NuTo::Constitutive::Input::DEFORMATION_GRADIENT_1D] = &(deformationGradient);
            }
                break;
            case Node::NONLOCALEQSTRAIN:
            {
                constitutiveInputList[NuTo::Constitutive::Input::NONLOCAL_EQ_STRAIN] = &(nonlocalEqStrain);
            }
                break;
            case Node::RELATIVEHUMIDITY:
            {
                constitutiveInputList[NuTo::Constitutive::Input::RELATIVE_HUMIDITY] = &(relativeHumidity);
                constitutiveInputList[NuTo::Constitutive::Input::RELATIVE_HUMIDITY_D1] = &(relativeHumidityD1);
                constitutiveInputList[NuTo::Constitutive::Input::RELATIVE_HUMIDITY_GRADIENT] = &(relativeHumidityGradient);
            }
                break;
            case Node::WATERVOLUMEFRACTION:
            {
                constitutiveInputList[NuTo::Constitutive::Input::WATER_VOLUME_FRACTION] = &(waterVolumeFraction);
                constitutiveInputList[NuTo::Constitutive::Input::WATER_VOLUME_FRACTION_D1] = &(waterVolumeFractionD1);
                constitutiveInputList[NuTo::Constitutive::Input::WATER_VOLUME_FRACTION_GRADIENT] = &(waterVolumeFractionGradient);
            }
                break;
            default:
                throw MechanicsException("[NuTo::Element1D::Evaluate] Constitutive input for " + Node::DofToString(dof) + " not implemented.");
            }
        }

        /*****************************************\
         *    FILL CONSTITUTIVE OUTPUT LIST      *
         \*****************************************/

        for (auto it = rElementOutput.begin(); it != rElementOutput.end(); it++)
        {
            switch (it->first)
            {
            case Element::INTERNAL_GRADIENT:
                it->second->GetFullVectorDouble().Resize(numActiveDofs);
                //if the stiffness matrix is constant, the corresponding internal force is calculated via the Kd
                //on the global level
                if (mStructure->GetHessianConstant(0) == false)
                {
                    for (auto dof : activeDofs)
                    {
                        switch (dof)
                        {
                        case Node::DISPLACEMENTS:
                        {
                            constitutiveOutputList[NuTo::Constitutive::Output::ENGINEERING_STRESS_1D] = &(engineeringStress1D);
                        }
                            break;
                        case Node::NONLOCALEQSTRAIN:
                        {
                            constitutiveOutputList[NuTo::Constitutive::Output::LOCAL_EQ_STRAIN] = &localEqStrain;
                            constitutiveOutputList[NuTo::Constitutive::Output::NONLOCAL_PARAMETER_XI] = &nonlocalParameter;
                        }
                            break;
                        case Node::RELATIVEHUMIDITY:
                        {
                            if (activeDofs.find(Node::WATERVOLUMEFRACTION) != activeDofs.end())
                            {
                                constitutiveOutputList[NuTo::Constitutive::Output::RESIDUAL_VAPOR_PHASE_N] = &residualVaporPhaseN;
                                constitutiveOutputList[NuTo::Constitutive::Output::RESIDUAL_VAPOR_PHASE_B] = &residualVaporPhaseB;
                            }
                        }
                            break;

                        case Node::WATERVOLUMEFRACTION:
                        {
                            if (activeDofs.find(Node::RELATIVEHUMIDITY) != activeDofs.end())
                            {
                                constitutiveOutputList[NuTo::Constitutive::Output::RESIDUAL_WATER_PHASE_N] = &residualWaterPhaseN;
                                constitutiveOutputList[NuTo::Constitutive::Output::RESIDUAL_WATER_PHASE_B] = &residualWaterPhaseB;
                            }
                        }
                            break;
                        default:
                            throw MechanicsException("[NuTo::Element1D::Evaluate] Constitutive output INTERNAL_GRADIENT for " + Node::DofToString(dof) + " not implemented.");

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
                    {
                        constitutiveOutputList[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_1D] = &tangentStressStrain;
                        if (activeDofs.find(Node::NONLOCALEQSTRAIN) != activeDofs.end())
                        {
                            constitutiveOutputList[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_1D] = &tangentStressNonlocalEqStrain;
                        }
                    }
                        break;
                    case Node::NONLOCALEQSTRAIN:
                    {
                        constitutiveOutputList[NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_XI_D_STRAIN_1D] = &tangentLocalEqStrainStrain;
                        constitutiveOutputList[NuTo::Constitutive::Output::NONLOCAL_PARAMETER_XI] = &nonlocalParameter;
                    }
                        break;
                    case Node::RELATIVEHUMIDITY:
                    {
                        if (activeDofs.find(Node::WATERVOLUMEFRACTION) != activeDofs.end())
                        {
                            constitutiveOutputList[NuTo::Constitutive::Output::D_RESIDUAL_RH_D_RH_H0_BB] = &tangent_D_Residual_RH_D_RH_H0_BB;
                            constitutiveOutputList[NuTo::Constitutive::Output::D_RESIDUAL_RH_D_RH_H0_NN] = &tangent_D_Residual_RH_D_RH_H0_NN;
                            constitutiveOutputList[NuTo::Constitutive::Output::D_RESIDUAL_RH_D_WV_H0_BN] = &tangent_D_Residual_RH_D_WV_H0_BN;
                            constitutiveOutputList[NuTo::Constitutive::Output::D_RESIDUAL_RH_D_WV_H0_NN] = &tangent_D_Residual_RH_D_WV_H0_NN;
                        }
                    }
                        break;
                    case Node::WATERVOLUMEFRACTION:
                    {
                        if (activeDofs.find(Node::RELATIVEHUMIDITY) != activeDofs.end())
                        {
                            constitutiveOutputList[NuTo::Constitutive::Output::D_RESIDUAL_WV_D_RH_H0_NN] = &tangent_D_Residual_WV_D_RH_H0_NN;
                            constitutiveOutputList[NuTo::Constitutive::Output::D_RESIDUAL_WV_D_WV_H0_BB] = &tangent_D_Residual_WV_D_WV_H0_BB;
                            constitutiveOutputList[NuTo::Constitutive::Output::D_RESIDUAL_WV_D_WV_H0_BN] = &tangent_D_Residual_WV_D_WV_H0_BN;
                            constitutiveOutputList[NuTo::Constitutive::Output::D_RESIDUAL_WV_D_WV_H0_NN] = &tangent_D_Residual_WV_D_WV_H0_NN;
                        }
                    }
                        break;
                    default:
                        throw MechanicsException("[NuTo::Element1D::Evaluate] Constitutive output HESSIAN_0_TIME_DERIVATIVE for " + Node::DofToString(dof) + " not implemented.");
                    }
                }
            }
                break;
            case Element::HESSIAN_1_TIME_DERIVATIVE:
                it->second->GetFullMatrixDouble().Resize(numActiveDofs, numActiveDofs);
                it->second->SetSymmetry(true);
                it->second->SetConstant(true);
                for (auto dof : activeDofs)
                {
                    switch (dof)
                    {
                    case Node::RELATIVEHUMIDITY:
                    {
                        if (activeDofs.find(Node::WATERVOLUMEFRACTION) != activeDofs.end())
                        {
                            constitutiveOutputList[NuTo::Constitutive::Output::D_RESIDUAL_RH_D_RH_H1_NN] = &tangent_D_Residual_RH_D_RH_H1_NN;
                            constitutiveOutputList[NuTo::Constitutive::Output::D_RESIDUAL_RH_D_WV_H1_NN] = &tangent_D_Residual_RH_D_WV_H1_NN;
                        }
                    }
                        break;
                    case Node::WATERVOLUMEFRACTION:
                    {
                        if (activeDofs.find(Node::RELATIVEHUMIDITY) != activeDofs.end())
                        {
                            constitutiveOutputList[NuTo::Constitutive::Output::D_RESIDUAL_WV_D_WV_H1_NN] = &tangent_D_Residual_WV_D_WV_H1_NN;
                        }
                    }
                        break;
                    default:
                        throw MechanicsException("[NuTo::Element1D::Evaluate] Constitutive output HESSIAN_1_TIME_DERIVATIVE for " + Node::DofToString(dof) + " not implemented.");
                    }
                }
                break;
            case Element::HESSIAN_2_TIME_DERIVATIVE:
            {
                it->second->GetFullMatrixDouble().Resize(numActiveDofs, numActiveDofs);
                it->second->SetSymmetry(true);
                it->second->SetConstant(true);
                //there is only a constant mass part for the mechanics problem
            }
                break;
            case Element::LUMPED_HESSIAN_2_TIME_DERIVATIVE:
                it->second->GetFullVectorDouble().Resize(numActiveDofs);
                break;
            case Element::UPDATE_STATIC_DATA:
                constitutiveOutputList[NuTo::Constitutive::Output::UPDATE_STATIC_DATA] = 0;
                break;
            case Element::UPDATE_TMP_STATIC_DATA:
                constitutiveOutputList[NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA] = 0;
                break;
            case Element::IP_DATA:
                switch (it->second->GetIpDataType())
                {
                case NuTo::IpData::ENGINEERING_STRAIN:
                    it->second->GetFullMatrixDouble().Resize(6, GetNumIntegrationPoints());
                    //define outputs
                    constitutiveOutputList[NuTo::Constitutive::Output::ENGINEERING_STRAIN_3D] = &(engineeringStrain3D);
                    break;
                case NuTo::IpData::ENGINEERING_STRESS:
                    it->second->GetFullMatrixDouble().Resize(6, GetNumIntegrationPoints());
                    //define outputs
                    constitutiveOutputList[NuTo::Constitutive::Output::ENGINEERING_STRESS_3D] = &(engineeringStress3D);
                    break;
                case NuTo::IpData::ENGINEERING_PLASTIC_STRAIN:
                    it->second->GetFullMatrixDouble().Resize(6, GetNumIntegrationPoints());
                    //define outputs
                    constitutiveOutputList[NuTo::Constitutive::Output::ENGINEERING_PLASTIC_STRAIN_3D] = &(engineeringPlasticStrain3D);
                    break;
                case NuTo::IpData::DAMAGE:
                    it->second->GetFullMatrixDouble().Resize(1, GetNumIntegrationPoints());
                    //define outputs
                    constitutiveOutputList[NuTo::Constitutive::Output::DAMAGE] = &(damage);
                    break;
                case NuTo::IpData::LOCAL_EQ_STRAIN:
                    it->second->GetFullMatrixDouble().Resize(1, GetNumIntegrationPoints());
                    constitutiveOutputList[NuTo::Constitutive::Output::LOCAL_EQ_STRAIN] = &localEqStrain;
                    break;
                default:
                    throw MechanicsException("[NuTo::Plane::Evaluate] this ip data type is not implemented.");
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
                throw MechanicsException("[NuTo::Element1D::Evaluate] element output not implemented.");
            }
        }				//end for: constitutive output list

        /*****************************************\
         *     CALCULATE CONSTITUTIVE INPUTS     *
         \*****************************************/

        std::map<Node::eDof, Eigen::VectorXd> shapeFunctions;
        std::map<Node::eDof, Eigen::MatrixXd> derivativeShapeFunctions;

        // loop over the integration points
        for (int theIP = 0; theIP < GetNumIntegrationPoints(); theIP++)
        {

            // calculate Jacobian
            const Eigen::MatrixXd& derivativeShapeFunctionsGeometryNatural = mInterpolationType->Get(Node::COORDINATES).GetDerivativeShapeFunctionsNatural(theIP);

            double detJacobian = CalculateJacobian(derivativeShapeFunctionsGeometryNatural, nodalValues[Node::COORDINATES]);

            const Eigen::VectorXd shapeFunctionsAtIp = mInterpolationType->Get(Node::COORDINATES).GetShapeFunctions(theIP);
            const Eigen::VectorXd globalIPCoordinate = nodalValues[Node::COORDINATES] * shapeFunctionsAtIp;

            double factor = detJacobian * mSection->GetArea() * mSection->AsSectionTruss()->GetAreaFactor(globalIPCoordinate[0]) * (mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP));

            // calculate shape functions and their derivatives
            for (auto dof : dofs)
            {
                if (dof == Node::COORDINATES)
                    continue;
                const InterpolationBase& interpolationType = mInterpolationType->Get(dof);
                shapeFunctions[dof] = interpolationType.GetShapeFunctions(theIP);
                derivativeShapeFunctions[dof] = interpolationType.GetDerivativeShapeFunctionsNatural(theIP) / detJacobian;
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
                    deformationGradient = CalculateDeformationGradient(derivativeShapeFunctions.at(dof), nodalValues.at(dof));
                }
                    break;
                case Node::NONLOCALEQSTRAIN:
                {
                    nonlocalEqStrain(0, 0) = (nodalValues[Node::NONLOCALEQSTRAIN] * shapeFunctions[Node::NONLOCALEQSTRAIN])(0, 0);
                }
                    break;
                case Node::RELATIVEHUMIDITY:
                {
                    relativeHumidity(0, 0) = (nodalValues[Node::RELATIVEHUMIDITY] * shapeFunctions[Node::RELATIVEHUMIDITY])(0, 0);
                    relativeHumidityD1(0, 0) = (ExtractNodeValues(1, Node::RELATIVEHUMIDITY) * shapeFunctions[Node::RELATIVEHUMIDITY])(0, 0);
                    relativeHumidityGradient(0, 0) = (nodalValues[Node::RELATIVEHUMIDITY] * derivativeShapeFunctions[Node::RELATIVEHUMIDITY])(0, 0);
                }
                    break;
                case Node::WATERVOLUMEFRACTION:
                {
                    waterVolumeFraction(0, 0) = (nodalValues[Node::WATERVOLUMEFRACTION] * shapeFunctions[Node::WATERVOLUMEFRACTION])(0, 0);
                    waterVolumeFractionD1(0, 0) = (ExtractNodeValues(1, Node::WATERVOLUMEFRACTION) * shapeFunctions[Node::WATERVOLUMEFRACTION])(0, 0);
                    waterVolumeFractionGradient(0, 0) = (nodalValues[Node::WATERVOLUMEFRACTION] * derivativeShapeFunctions[Node::WATERVOLUMEFRACTION])(0, 0);
                }
                    break;
                default:
                    throw MechanicsException("[NuTo::Element1D::Evaluate] Constitutive input for " + Node::DofToString(dof) + " not implemented.");
                }
            }

            /*****************************************\
             *      EVALUATE CONSTITUTIVE LAW        *
             \*****************************************/

            ConstitutiveBase* constitutivePtr = GetConstitutiveLaw(theIP);
            try
            {
                Error::eError error = constitutivePtr->Evaluate1D(this, theIP, constitutiveInputList, constitutiveOutputList);
                if (error != Error::SUCCESSFUL)
                    return error;
            } catch (NuTo::MechanicsException &e)
            {
                e.AddMessage("[NuTo::Element1D::Evaluate] error evaluating the constitutive model.");
                throw e;
            }

            /*****************************************\
             *           CALCULATE OUTPUT            *
             \*****************************************/

            //calculate Kee detJacobian*(NtN+cBtB)
            FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> Kee;
            if (dofs.find(Node::NONLOCALEQSTRAIN) != dofs.end())
            {
                CalculateKee(shapeFunctions.at(Node::NONLOCALEQSTRAIN), derivativeShapeFunctions.at(Node::NONLOCALEQSTRAIN), nonlocalParameter, factor, Kee);
            }

            for (auto it = rElementOutput.begin(); it != rElementOutput.end(); it++)
            {
                switch (it->first)
                {
                case Element::INTERNAL_GRADIENT:
                {
                    //if the stiffness matrix is constant, the corresponding internal force is calculated via the Kd
                    //on the global level
                    if (mStructure->GetHessianConstant(0) == false)
                    {
                        for (auto dof : activeDofs)
                        {
                            int startIndex = mInterpolationType->Get(dof).GetLocalStartIndex();
                            switch (dof)
                            {
                            case Node::DISPLACEMENTS:
                            {
                                AddDetJBtSigma(derivativeShapeFunctions.at(dof), engineeringStress1D, factor, startIndex, it->second->GetFullVectorDouble());
                            }
                                break;
                            case Node::NONLOCALEQSTRAIN:
                            {
                                int startIndexNonlocalEqStrain = mInterpolationType->Get(Node::NONLOCALEQSTRAIN).GetLocalStartIndex();
                                AddDetJRnonlocalEqStrain(shapeFunctions.at(Node::NONLOCALEQSTRAIN), localEqStrain, Kee, nodalValues.at(Node::NONLOCALEQSTRAIN), factor / nonlocalParameter.GetValue(0), startIndexNonlocalEqStrain, it->second->GetFullVectorDouble());
                            }
                                break;
                            case Node::RELATIVEHUMIDITY:
                            {
                                if (activeDofs.find(Node::WATERVOLUMEFRACTION) != activeDofs.end())
                                {
                                    AddDetJBtX(derivativeShapeFunctions.at(Node::RELATIVEHUMIDITY), residualVaporPhaseB, factor, startIndex, it->second->GetFullVectorDouble());

                                    AddDetJNtX(shapeFunctions.at(Node::RELATIVEHUMIDITY), residualVaporPhaseN, factor, startIndex, it->second->GetFullVectorDouble());
                                }
                            }
                                break;
                            case Node::WATERVOLUMEFRACTION:
                            {
                                if (activeDofs.find(Node::RELATIVEHUMIDITY) != activeDofs.end())
                                {
                                    AddDetJBtX(derivativeShapeFunctions.at(Node::WATERVOLUMEFRACTION), residualWaterPhaseB, factor, startIndex, it->second->GetFullVectorDouble());

                                    AddDetJNtX(shapeFunctions.at(Node::WATERVOLUMEFRACTION), residualWaterPhaseN, factor, startIndex, it->second->GetFullVectorDouble());
                                }
                            }
                                break;
                            default:
                                throw MechanicsException("[NuTo::Element1D::Evaluate] Element output INTERNAL_GRADIENT for " + Node::DofToString(dof) + " not implemented.");

                            }
                        }
                    }
                }
                    break;
                case Element::HESSIAN_0_TIME_DERIVATIVE:
                {
                    //factor for the numerical integration
                    assert(mSection->GetArea() > 0);
                    if (tangentStressStrain.GetConstant() == false)
                        it->second->SetConstant(false);
                    if (tangentStressStrain.GetSymmetry() == false)
                        it->second->SetSymmetry(false);

                    for (auto dof : activeDofs)
                    {
                        int startIndex = mInterpolationType->Get(dof).GetLocalStartIndex();
                        switch (dof)
                        {
                        case Node::DISPLACEMENTS:
                        {
                            AddDetJBtCB(derivativeShapeFunctions.at(Node::DISPLACEMENTS), tangentStressStrain, factor, startIndex, startIndex, it->second->GetFullMatrixDouble());

                            if (activeDofs.find(Node::NONLOCALEQSTRAIN) != activeDofs.end())
                            {
                                int startIndexNonlocalEqStrain = mInterpolationType->Get(Node::NONLOCALEQSTRAIN).GetLocalStartIndex();
                                AddDetJBtdSigmadNonlocalEqStrainN(derivativeShapeFunctions.at(Node::DISPLACEMENTS), tangentStressNonlocalEqStrain, shapeFunctions.at(Node::NONLOCALEQSTRAIN), factor, startIndex, startIndexNonlocalEqStrain, it->second->GetFullMatrixDouble());
                            }
                        }
                            break;
                        case Node::NONLOCALEQSTRAIN:
                        {
                            int startIndexNonlocalEqStrain = mInterpolationType->Get(Node::NONLOCALEQSTRAIN).GetLocalStartIndex();
                            it->second->GetFullMatrixDouble().AddBlock(startIndexNonlocalEqStrain, startIndexNonlocalEqStrain, Kee);

                            if (activeDofs.find(Node::DISPLACEMENTS) != activeDofs.end())
                            {
                                int startIndexDisplacement = mInterpolationType->Get(Node::DISPLACEMENTS).GetLocalStartIndex();
                                AddDetJNtdLocalEqStraindEpsilonB(shapeFunctions.at(Node::NONLOCALEQSTRAIN), tangentLocalEqStrainStrain, derivativeShapeFunctions.at(Node::DISPLACEMENTS), factor, startIndexNonlocalEqStrain, startIndexDisplacement, it->second->GetFullMatrixDouble());
                            }

                        }
                            break;

                        case Node::RELATIVEHUMIDITY:
                        {
                            if (activeDofs.find(Node::WATERVOLUMEFRACTION) != activeDofs.end())
                            {
                                int indexRelHum = startIndex;
                                int indexWatVol = mInterpolationType->Get(Node::WATERVOLUMEFRACTION).GetLocalStartIndex();

                                auto& RelHumShapeFunction = shapeFunctions.at(Node::RELATIVEHUMIDITY);
                                auto& WatVolShapeFunction = shapeFunctions.at(Node::WATERVOLUMEFRACTION);

                                auto& RelHumDerivativeShapeFunction = derivativeShapeFunctions.at(Node::RELATIVEHUMIDITY);

                                // | - - |
                                // | - X |

                                AddDetJBtXB(RelHumDerivativeShapeFunction, RelHumDerivativeShapeFunction, tangent_D_Residual_RH_D_RH_H0_BB, factor, indexRelHum, indexRelHum, it->second->GetFullMatrixDouble());

                                AddDetJNtXN(RelHumShapeFunction, RelHumShapeFunction, tangent_D_Residual_RH_D_RH_H0_NN, factor, indexRelHum, indexRelHum, it->second->GetFullMatrixDouble());

                                // coupling terms

                                // | - - |
                                // | X - |

                                AddDetJBtXN(RelHumDerivativeShapeFunction, WatVolShapeFunction, tangent_D_Residual_RH_D_WV_H0_BN, factor, indexRelHum, indexWatVol, it->second->GetFullMatrixDouble());

                                AddDetJNtXN(RelHumShapeFunction, WatVolShapeFunction, tangent_D_Residual_RH_D_WV_H0_NN, factor, indexRelHum, indexWatVol, it->second->GetFullMatrixDouble());

                                //it->second->GetFullMatrixDouble().Info();         // Check Element Matrix if needed
                                //NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> test = it->second->GetFullMatrixDouble();
                                //test.Info(4,5,true);
                                //int a=0;
                                //a++;

                            }
                        }
                            break;
                        case Node::WATERVOLUMEFRACTION:
                        {
                            if (activeDofs.find(Node::RELATIVEHUMIDITY) != activeDofs.end())
                            {
                                int indexRelHum = mInterpolationType->Get(Node::RELATIVEHUMIDITY).GetLocalStartIndex();
                                int indexWatVol = startIndex;

                                auto& RelHumShapeFunction = shapeFunctions.at(Node::RELATIVEHUMIDITY);
                                auto& WatVolShapeFunction = shapeFunctions.at(Node::WATERVOLUMEFRACTION);

                                auto& WatVolDerivativeShapeFunction = derivativeShapeFunctions.at(Node::WATERVOLUMEFRACTION);

                                // | X - |
                                // | - - |

                                AddDetJBtXB(WatVolDerivativeShapeFunction, WatVolDerivativeShapeFunction, tangent_D_Residual_WV_D_WV_H0_BB, factor, indexWatVol, indexWatVol, it->second->GetFullMatrixDouble());

                                AddDetJBtXN(WatVolDerivativeShapeFunction, WatVolShapeFunction, tangent_D_Residual_WV_D_WV_H0_BN, factor, indexWatVol, indexWatVol, it->second->GetFullMatrixDouble());

                                AddDetJNtXN(WatVolShapeFunction, WatVolShapeFunction, tangent_D_Residual_WV_D_WV_H0_NN, factor, indexWatVol, indexWatVol, it->second->GetFullMatrixDouble());

                                // coupling terms

                                // | - X |
                                // | - - |

                                AddDetJNtXN(WatVolShapeFunction, RelHumShapeFunction, tangent_D_Residual_WV_D_RH_H0_NN, factor, indexWatVol, indexRelHum, it->second->GetFullMatrixDouble());

                                //it->second->GetFullMatrixDouble().Info();         // Check Element Matrix if needed
                                //NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> test = it->second->GetFullMatrixDouble();
                                //int a=0;
                            }
                        }
                            break;
                        default:
                            throw MechanicsException("[NuTo::Element1D::Evaluate] Element output HESSIAN_0_TIME_DERIVATIVE for " + Node::DofToString(dof) + " not implemented.");

                        }
                    }
                }
                    break;
                case Element::HESSIAN_1_TIME_DERIVATIVE:
                {
                    //factor for the numerical integration
                    assert(mSection->GetArea() > 0);
                    if (tangentStressStrain.GetConstant() == false)
                        it->second->SetConstant(false);
                    if (tangentStressStrain.GetSymmetry() == false)
                        it->second->SetSymmetry(false);

                    for (auto dof : activeDofs)
                    {
                        int startIndex = mInterpolationType->Get(dof).GetLocalStartIndex();
                        switch (dof)
                        {
                        case Node::RELATIVEHUMIDITY:
                        {
                            if (activeDofs.find(Node::WATERVOLUMEFRACTION) != activeDofs.end())
                            {
                                int indexRelHum = startIndex;
                                int indexWatVol = mInterpolationType->Get(Node::WATERVOLUMEFRACTION).GetLocalStartIndex();

                                auto RelHumShapeFunction = shapeFunctions.at(Node::RELATIVEHUMIDITY);
                                auto WatVolShapeFunction = shapeFunctions.at(Node::WATERVOLUMEFRACTION);

                                // | - - |
                                // | - X |

                                AddDetJNtXN(RelHumShapeFunction, RelHumShapeFunction, tangent_D_Residual_RH_D_RH_H1_NN, factor, indexRelHum, indexRelHum, it->second->GetFullMatrixDouble());

                                // coupling terms

                                // | - - |
                                // | X - |

                                AddDetJNtXN(RelHumShapeFunction, WatVolShapeFunction, tangent_D_Residual_RH_D_WV_H1_NN, factor, indexRelHum, indexWatVol, it->second->GetFullMatrixDouble());

                                //it->second->GetFullMatrixDouble().Info();         // Check Element Matrix if needed
                                //NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> test = it->second->GetFullMatrixDouble();
                                //test.Info(4,4,true);
                                //int a=0;
                                //a++;

                            }
                        }
                            break;
                        case Node::WATERVOLUMEFRACTION:
                        {
                            if (activeDofs.find(Node::RELATIVEHUMIDITY) != activeDofs.end())
                            {
                                int indexWatVol = startIndex;

                                auto WatVolShapeFunction = shapeFunctions.at(Node::WATERVOLUMEFRACTION);

                                // | X - |
                                // | - - |

                                AddDetJNtXN(WatVolShapeFunction, WatVolShapeFunction, tangent_D_Residual_WV_D_WV_H1_NN, factor, indexWatVol, indexWatVol, it->second->GetFullMatrixDouble());

                                //it->second->GetFullMatrixDouble().Info();         // Check Element Matrix if needed
                                //NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> test = it->second->GetFullMatrixDouble();
                                //int a=0;
                            }
                        }
                            break;
                        default:
                            throw MechanicsException("[NuTo::Element1D::Evaluate] Element output HESSIAN_0_TIME_DERIVATIVE for " + Node::DofToString(dof) + " not implemented.");

                        }
                    }
                }
                    break;

                case Element::HESSIAN_2_TIME_DERIVATIVE:
                {
                    for (auto dof : activeDofs)
                    {
                        switch (dof)
                        {
                        case Node::DISPLACEMENTS:
                        {

                            double density = constitutivePtr->GetParameterDouble(Constitutive::eConstitutiveParameter::DENSITY);
                            const Eigen::MatrixXd tmpMatrix = shapeFunctions.at(dof) * shapeFunctions.at(dof).transpose() * factor * density;

                            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& result(it->second->GetFullMatrixDouble());
                            result.block(0, 0, tmpMatrix.rows(), tmpMatrix.cols()) += tmpMatrix;

                        }
                            break;
                        case Node::RELATIVEHUMIDITY:
                        case Node::WATERVOLUMEFRACTION:
                        {
                        }
                            break;
                        case Node::TEMPERATURE:
                            break;
                        default:
                            throw MechanicsException("[NuTo::Element1D::Evaluate] Element output HESSIAN_2_TIME_DERIVATIVE for " + Node::DofToString(dof) + " not implemented.");

                        }
                    }
                }
                    break;
                case Element::LUMPED_HESSIAN_2_TIME_DERIVATIVE:
                    for (auto dof : activeDofs)
                    {
                        switch (dof)
                        {
                        case Node::DISPLACEMENTS:
                        {
                            //this->CalculateShapeFunctionsField(localIPCoord, shapeFunctionsField);
                            // calculate local mass matrix (the nonlocal terms are zero)
                            // don't forget to include determinant of the Jacobian and area
                            // detJ * area * density * HtH, :
                            double density = constitutivePtr->GetParameterDouble(Constitutive::eConstitutiveParameter::DENSITY);
                            double factor2 (density * factor);
                            FullVector<double,Eigen::Dynamic>& result(it->second->GetFullVectorDouble());
                            total_mass+=factor2;
                            //calculate for the translational dofs the diagonal entries
                            int size = shapeFunctions.at(dof).size();

                            for (int count = 0; count < size; count++)
                            {
                                result(count)+= shapeFunctions.at(dof)(count)*shapeFunctions.at(dof)(count)*factor2;
                            }

                            if (theIP+1==GetNumIntegrationPoints())
                            {
                                //calculate sum of diagonal entries (is identical for all directions, that's why only x direction is calculated
                                double sum_diagonal(0);
                                for (int count = 0; count < size; count++)
                                {
                                    sum_diagonal+= result(count);
                                }

                                //scale so that the sum of the diagonals represents the full mass
                                double scaleFactor = total_mass/sum_diagonal;
                                for (int count = 0; count < size; count++)
                                {
                                    result(count) *= scaleFactor;
                                }
                            }
                        }
                            break;
                        default:
                            throw MechanicsException("[NuTo::Element1D::Evaluate] Element output LUMPED_HESSIAN_2_TIME_DERIVATIVE for " + Node::DofToString(dof) + " not implemented.");
                        }
                    }

                    break;
                case Element::UPDATE_STATIC_DATA:
                case Element::UPDATE_TMP_STATIC_DATA:
                    break;
                case Element::IP_DATA:
                    switch (it->second->GetIpDataType())
                    {
                    case NuTo::IpData::ENGINEERING_STRAIN:
                        //error = constitutivePtr->GetEngineeringStrain(this, theIP, deformationGradient, engineeringStrain);
                        memcpy(&(it->second->GetFullMatrixDouble().data()[theIP * 6]), engineeringStrain3D.GetData(), 6 * sizeof(double));
                        break;
                    case NuTo::IpData::ENGINEERING_STRESS:
                        //error = constitutivePtr->GetEngineeringStressFromEngineeringStrain(this, theIP, deformationGradient, engineeringStress);
                        memcpy(&(it->second->GetFullMatrixDouble().data()[theIP * 6]), engineeringStress3D.GetData(), 6 * sizeof(double));
                        break;
                    case NuTo::IpData::ENGINEERING_PLASTIC_STRAIN:
                        //error = constitutivePtr->GetEngineeringPlasticStrain(this, theIP, deformationGradient, engineeringStrain);
                        memcpy(&(it->second->GetFullMatrixDouble().data()[theIP * 6]), engineeringPlasticStrain3D.GetData(), 6 * sizeof(double));
                        break;
                    case NuTo::IpData::DAMAGE:
                        //error = constitutivePtr->GetDamage(this, theIP, deformationGradient, rIpData.mEigenMatrix.data()[theIP]);
                        memcpy(&(it->second->GetFullMatrixDouble().data()[theIP]), damage.GetData(), sizeof(double));
                        break;
                    case NuTo::IpData::LOCAL_EQ_STRAIN:
                        memcpy(&(it->second->GetFullMatrixDouble().data()[theIP]), localEqStrain.data(), sizeof(double));
                        break;
                    default:
                        throw MechanicsException("[NuTo::Plane::GetIpData] Ip data not implemented.");
                    }
                    break;
                case Element::GLOBAL_ROW_DOF:
                case Element::GLOBAL_COLUMN_DOF:
                    break;
                default:
                    throw MechanicsException("[NuTo::Element1D::Evaluate] element output not implemented.");
                }
            }

        }

    } catch (NuTo::MechanicsException& e)
    {
        std::stringstream ss;
        ss << mStructure->ElementGetId(this);
        e.AddMessage("[NuTo::Element1D::Evaluate] Error evaluating element data of element " + ss.str() + ".");
        throw e;
    }
    return Error::SUCCESSFUL;
}

NuTo::Element::eElementType NuTo::Element1D::GetEnumType() const
{
    return NuTo::Element::ELEMENT1D;
}

int NuTo::Element1D::GetLocalDimension()const
{
    return 1;
}

NuTo::NodeBase* NuTo::Element1D::GetNode(int rLocalNodeNumber)
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->GetNumNodes());
    return mNodes[rLocalNodeNumber];
}

const NuTo::NodeBase* NuTo::Element1D::GetNode(int rLocalNodeNumber) const
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->GetNumNodes());
    return mNodes[rLocalNodeNumber];
}

NuTo::NodeBase* NuTo::Element1D::GetNode(int rLocalNodeNumber, Node::eDof rDofType)
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->Get(rDofType).GetNumNodes());
    int nodeIndex = mInterpolationType->Get(rDofType).GetNodeIndex(rLocalNodeNumber);
    assert(nodeIndex < (int )mNodes.size());
    return mNodes[nodeIndex];
}

const NuTo::NodeBase* NuTo::Element1D::GetNode(int rLocalNodeNumber, Node::eDof rDofType) const
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->Get(rDofType).GetNumNodes());
    int nodeIndex = mInterpolationType->Get(rDofType).GetNodeIndex(rLocalNodeNumber);
    assert(nodeIndex < (int )mNodes.size());
    return mNodes[nodeIndex];
}

void NuTo::Element1D::SetNode(int rLocalNodeNumber, NodeBase* rNode)
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->GetNumNodes());
    assert(rNode != nullptr);
    mNodes[rLocalNodeNumber] = rNode;
}

void NuTo::Element1D::ResizeNodes(int rNewNumNodes)
{
    if (rNewNumNodes == (int) mNodes.size())
        return;

    if (rNewNumNodes > (int) mNodes.size())
    {
        // just resize (enlarge)
        mNodes.resize(rNewNumNodes);
    } else
    {
        throw MechanicsException("[NuTo::Element1D::ResizeNodes] Resize that reduces the number of nodes is not implemented yet.");
    }
}

void NuTo::Element1D::ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
{
    for (unsigned int count = 0; count < mNodes.size(); count++)
    {
        if (this->mNodes[count] == rOldPtr)
        {
            this->mNodes[count] = rNewPtr;
            break;
        }
    }
}

void NuTo::Element1D::SetSection(const SectionBase* rSection)
{
    mSection = rSection;
}

const NuTo::SectionBase* NuTo::Element1D::GetSection() const
{
    return mSection;
}

NuTo::ConstitutiveStaticDataBase* NuTo::Element1D::AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw) const
{
    return rConstitutiveLaw->AllocateStaticDataEngineeringStress_EngineeringStrain1D(this);
}

const Eigen::VectorXd NuTo::Element1D::GetIntegrationPointVolume() const
{

    Eigen::MatrixXd localNodeCoord = this->ExtractNodeValues(0, Node::COORDINATES);

    const InterpolationBase& interpolationType = mInterpolationType->Get(Node::COORDINATES);

    Eigen::VectorXd volume(GetNumIntegrationPoints());

    for (int theIP = 0; theIP < GetNumIntegrationPoints(); theIP++)
    {
        Eigen::MatrixXd derivativeShapeFunctionsNatural = interpolationType.GetDerivativeShapeFunctionsNatural(theIP);

        double detJacobian = CalculateJacobian(derivativeShapeFunctionsNatural, localNodeCoord);

        //attention in 2D, this is just the area, but that is required for the nonlocal model
        volume[theIP] = detJacobian * mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP);
    }
    return volume;
}

const Eigen::MatrixXd NuTo::Element1D::ExtractNodeValues(int rTimeDerivative, Node::eDof rDofType) const
{
    const InterpolationBase& interpolationTypeDof = GetInterpolationType()->Get(rDofType);

    int numNodes = interpolationTypeDof.GetNumNodes();
    int numDofsPerNode = interpolationTypeDof.GetNumDofsPerNode();


    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nodalValues;
    nodalValues.resize(numDofsPerNode, numNodes);

    for (int iNode = 0; iNode < numNodes; ++iNode)
    {
        const NodeBase* node = mNodes[interpolationTypeDof.GetNodeIndex(iNode)];

        switch (rDofType)
        {
        case Node::COORDINATES:
             nodalValues.block(0, iNode, numDofsPerNode, 1) = node->GetCoordinates1D();
            break;

        case Node::DISPLACEMENTS:
            nodalValues.block(0, iNode, numDofsPerNode, 1) = node->GetDisplacements1D(rTimeDerivative);
            break;

        case Node::TEMPERATURE:
            nodalValues(0, iNode) = node->GetTemperature(rTimeDerivative);
            break;

        case Node::NONLOCALEQSTRAIN:
            nodalValues(0, iNode) = node->GetNonlocalEqStrain(rTimeDerivative);
            break;

        case Node::RELATIVEHUMIDITY:
            nodalValues(0, iNode) = node->GetRelativeHumidity(rTimeDerivative);
            break;

        case Node::WATERVOLUMEFRACTION:
            nodalValues(0, iNode) = node->GetWaterVolumeFraction(rTimeDerivative);
            break;

        default:
            throw MechanicsException("[NuTo::Element1D::ExtractNodeValues] Not implemented for " + Node::DofToString(rDofType) + ".");
        }
    }
    return nodalValues;
}

const Eigen::VectorXi NuTo::Element1D::CalculateGlobalRowDofs() const
{
    int numActiveDos = mInterpolationType->GetNumActiveDofs();
    Eigen::VectorXi globalRowDofs(numActiveDos);

    for (auto dof : mInterpolationType->GetActiveDofs())
    {
        const InterpolationBase& interpolationType = mInterpolationType->Get(dof);
        int index = interpolationType.GetLocalStartIndex(); //

        for (int iNodeDof = 0; iNodeDof < interpolationType.GetNumNodes(); ++iNodeDof)
        {
            const NodeBase* nodePtr = mNodes[interpolationType.GetNodeIndex(iNodeDof)];
            switch (dof)
            {
            case Node::DISPLACEMENTS:
            {
                globalRowDofs[index++] = nodePtr->GetDofDisplacement(0);
            }
                break;
            case Node::TEMPERATURE:
            {
                globalRowDofs[index++] = nodePtr->GetDofTemperature();
            }
                break;
            case Node::NONLOCALEQSTRAIN:
            {
                globalRowDofs[index++] = nodePtr->GetDofNonlocalEqStrain();
            }
                break;
            case Node::RELATIVEHUMIDITY:
            {
                globalRowDofs[index++] = nodePtr->GetDofRelativeHumidity();
            }
                break;
            case Node::WATERVOLUMEFRACTION:
            {
                globalRowDofs[index++] = nodePtr->GetDofWaterVolumeFraction();
            }
                break;
            default:
                throw MechanicsException("[NuTo::Element1D::CalculateGlobalRowDofs] Not implemented for " + Node::DofToString(dof) + ".");

            }
        }
    }
    return globalRowDofs;
}

const Eigen::VectorXi NuTo::Element1D::CalculateGlobalColumnDofs() const
{
    int numNonlocalElements = GetNumNonlocalElements();
    if (numNonlocalElements == 0)
        return CalculateGlobalRowDofs();
    else
    {
        throw NuTo::MechanicsException("[NuTo::Element1D::CalculateGlobalColumnDofs] not implemented for nonlocal integral element formulation.");
    }
}

const NuTo::DeformationGradient1D NuTo::Element1D::CalculateDeformationGradient(const Eigen::MatrixXd& rDerivativeShapeFunctionsLocal, const Eigen::MatrixXd& rNodalDisplacements) const
{
    DeformationGradient1D deformationGradient;

    assert(rDerivativeShapeFunctionsLocal.rows() == rNodalDisplacements.cols());
    assert(rDerivativeShapeFunctionsLocal.cols() == 1);
    assert(rNodalDisplacements.rows() == 1);

    deformationGradient.mDeformationGradient = 1 + (rNodalDisplacements * rDerivativeShapeFunctionsLocal).at(0, 0);
    return deformationGradient;
}

//! @brief adds to a matrix the product B^tCB, where B contains the derivatives of the shape functions and C is the constitutive tangent
//! eventually include also area/width of an element (that's the thermal solution)
//! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param ConstitutiveTangentBase constitutive tangent matrix
//! @param rFactor factor including determinant of Jacobian and IP weight
//! @param rRow row, where to start to add the submatrix
//! @param rCoefficientMatrix to be added to
void NuTo::Element1D::AddDetJBtCB(const Eigen::MatrixXd& rDerivativeShapeFunctionsGlobal, const ConstitutiveTangentLocal<1, 1>& rConstitutiveTangent, double rFactor, int rRow, int rCol, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCoefficientMatrix) const
{
    int numNodes = mInterpolationType->Get(Node::DISPLACEMENTS).GetNumNodes();
    assert(rDerivativeShapeFunctionsGlobal.rows() == numNodes);
    assert(rDerivativeShapeFunctionsGlobal.cols() == 1);

    Eigen::MatrixXd result = rFactor * rConstitutiveTangent.data()[0] * rDerivativeShapeFunctionsGlobal * rDerivativeShapeFunctionsGlobal.transpose();

    rCoefficientMatrix.block(rRow, rCol, numNodes, numNodes) += result;

}

//! @brief adds up the internal force vector
//! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param rEngineeringStress stress
//! @param factor factor including det Jacobian area and integration point weight
//! @param rRow start row (in case of a multifield problem)
//! @param rResult resforce vector
void NuTo::Element1D::AddDetJBtSigma(const Eigen::MatrixXd& rDerivativeShapeFunctionsLocal, const EngineeringStress1D& rEngineeringStress, double rFactor, int rRow, FullVector<double, Eigen::Dynamic>& rResult) const
{
//    assert(rResult.GetNumRows()==2*GetNumNodesField() && rResult.GetNumColumns()==1);
    assert(rResult.GetNumColumns() == 1);
    int numNodes = mInterpolationType->Get(Node::DISPLACEMENTS).GetNumNodes();
    assert(rDerivativeShapeFunctionsLocal.rows() == numNodes);
    assert(rDerivativeShapeFunctionsLocal.cols() == 1);

    Eigen::VectorXd result = rFactor * rEngineeringStress.GetData()[0] * rDerivativeShapeFunctionsLocal;

    rResult.block(rRow, 0, numNodes, 1) += result;

}

//! @brief add detJ B.T dSigma/dnonlocalEqStrain N
//! @param derivativeShapeFunctions of the ip for all shape functions
//! @param tangentStressNonlocalEqStrain derivative of the stress with respect to the nonlocal eq strain
//! @param rShapeFunctions of the ip for all shape functions
//! @param rFactor factor including detJ and area
//! @param rResult result
void NuTo::Element1D::AddDetJBtdSigmadNonlocalEqStrainN(const Eigen::MatrixXd& rDerivativeShapeFunctions, ConstitutiveTangentLocal<1, 1>& rTangentStressNonlocalEqStrain, Eigen::VectorXd rShapeFunctions, double rFactor, int rRow, int rCol,
        FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rResult) const
{
    double tmpfactor = rTangentStressNonlocalEqStrain(0) * rFactor;
    int numDofsDisplacement = mInterpolationType->Get(Node::DISPLACEMENTS).GetNumNodes();
    int numDofsNonlocalEqStrain = mInterpolationType->Get(Node::NONLOCALEQSTRAIN).GetNumNodes();

    assert(rDerivativeShapeFunctions.rows() == numDofsDisplacement);
    assert(rDerivativeShapeFunctions.cols() == 1);
    assert(rShapeFunctions.rows() == numDofsNonlocalEqStrain);
    assert(rShapeFunctions.cols() == 1);

    rResult.block(rRow, rCol, numDofsDisplacement, numDofsNonlocalEqStrain) += tmpfactor * rDerivativeShapeFunctions * rShapeFunctions.transpose();

}

//! @brief add detJ N_transpose dEqStrain/dEpsilon B
//! @param rShapeFunctions of the ip for the nonlocal eq strain dofs
//! @param rTangentLocalEqStrainStrain derivative of the local eq strains with respect to the strain
//! @param rderivativeShapeFunctions of the ip for the displacement dofs
//! @param rFactor factor including detJ and area
//! @param rResult result
void NuTo::Element1D::AddDetJNtdLocalEqStraindEpsilonB(Eigen::VectorXd rShapeFunctions, ConstitutiveTangentLocal<1, 1>& rTangentLocalEqStrainStrain, const Eigen::MatrixXd& rDerivativeShapeFunctions, double rFactor, int rRow, int rCol,
        FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rResult) const
{
    double tmpfactor = rTangentLocalEqStrainStrain(0) * rFactor;
    int numDofsDisplacement = mInterpolationType->Get(Node::DISPLACEMENTS).GetNumNodes();
    int numDofsNonlocalEqStrain = mInterpolationType->Get(Node::NONLOCALEQSTRAIN).GetNumNodes();

    assert(rDerivativeShapeFunctions.rows() == numDofsDisplacement);
    assert(rDerivativeShapeFunctions.cols() == 1);
    assert(rShapeFunctions.rows() == numDofsNonlocalEqStrain);
    assert(rShapeFunctions.cols() == 1);

    rResult.block(rRow, rCol, numDofsNonlocalEqStrain, numDofsDisplacement) -= tmpfactor * rShapeFunctions * rDerivativeShapeFunctions.transpose();
}

//! @brief adds up the constitutive Tangent times the Shape Functions
//! @param rShapeFunctions the shape functions
//! @param rConstitutiveTangent the result given by the constitutive law
//! @param factor factor including det Jacobian area and integration point weight
//! @param rRow start row in case of a multifield problem
//! @param rResult result vector
void NuTo::Element1D::AddDetJNtX(Eigen::VectorXd &rShapeFunctions, ConstitutiveTangentLocal<1, 1> &rConstitutiveTangent, double rFactor, int rRow, FullVector<double, Eigen::Dynamic> &rResult) const
{
    double tmpfactor = rConstitutiveTangent(0) * rFactor;
    int NumDofs = rShapeFunctions.rows();

    assert(rShapeFunctions.cols() == 1);

    Eigen::VectorXd result = tmpfactor * rShapeFunctions;

    rResult.block(rRow, 0, NumDofs, 1) += result;
}

//! @brief adds up the constitutive Tangent times the derivative shape functions
//! @param rDerivativeShapeFunctions the derivative shape functions
//! @param rConstitutiveTangent the result given by the constitutive law
//! @param factor factor including det Jacobian area and integration point weight
//! @param rRow start row in case of a multifield problem
//! @param rResult result vector
void NuTo::Element1D::AddDetJBtX(const Eigen::MatrixXd &rDerivativeShapeFunctions, ConstitutiveTangentLocal<1, 1> &rConstitutiveTangent, double rFactor, int rRow, FullVector<double, Eigen::Dynamic> &rResult) const
{
    double tmpfactor = rConstitutiveTangent(0) * rFactor;
    int NumDofs = rDerivativeShapeFunctions.rows();

    assert(rDerivativeShapeFunctions.cols() == 1);

    Eigen::VectorXd result = tmpfactor * rDerivativeShapeFunctions;

    rResult.block(rRow, 0, NumDofs, 1) += result;
}

//! @brief adds to a matrix the product N1^t X N2, where N1 and N2 contains the the shape functions and X is the constitutive tangent
//! @param rShapeFunction1 shape function 1
//! @param rShapeFunction2 shape function 2
//! @param ConstitutiveTangentBase constitutive tangent matrix
//! @param rFactor factor including area, determinant of Jacobian and IP weight
//! @param rRow row, where to start to add the submatrix
//! @param rCol col, where to start to add the submatrix
//! @param rResult result
void NuTo::Element1D::AddDetJNtXN(const Eigen::VectorXd &rShapeFunction1, const Eigen::VectorXd &rShapeFunction2, const ConstitutiveTangentLocal<1, 1> &rConstitutiveTangent, double rFactor, int rRow, int rCol, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> &rResult) const
{
    int NumRowDofs = rShapeFunction1.rows();
    int NumColDofs = rShapeFunction2.rows();

    // add assertions!!!

    rFactor *= rConstitutiveTangent(0, 0);

    Eigen::MatrixXd result = rFactor * rShapeFunction1 * rShapeFunction2.transpose();

    rResult.block(rRow, rCol, NumRowDofs, NumColDofs) += result;
}

//! @brief adds to a matrix the product B^t X N, where B is the derivative shape functions, N the shape function and X the constitutive tangent
//! @param rDerivativeShapeFunction derivative shape function
//! @param rShapeFunction shape function
//! @param ConstitutiveTangentBase constitutive tangent matrix
//! @param rFactor factor including area, determinant of Jacobian and IP weight
//! @param rRow row, where to start to add the submatrix
//! @param rCol col, where to start to add the submatrix
//! @param rResult result
void NuTo::Element1D::AddDetJBtXN(const Eigen::MatrixXd &rDerivativeShapeFunction, const Eigen::VectorXd &rShapeFunction, const ConstitutiveTangentLocal<1, 1> &rConstitutiveTangent, double rFactor, int rRow, int rCol, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> &rResult) const
{
    int NumRowDofs = rDerivativeShapeFunction.rows();
    int NumColDofs = rShapeFunction.rows();

    // add assertions!!!

    rFactor *= rConstitutiveTangent(0, 0);

    Eigen::MatrixXd result = rFactor * rDerivativeShapeFunction * rShapeFunction.transpose();

    rResult.block(rRow, rCol, NumRowDofs, NumColDofs) += result;
}

//! @brief adds to a matrix the product B1^t X B2, where B1 and B2 are the derivative shape functions and X is the constitutive tangent
//! @param rDerivativeShapeFunction1 derivative shape function 1
//! @param rDerivativeShapeFunction2 derivative shape function 2
//! @param ConstitutiveTangentBase constitutive tangent matrix
//! @param rFactor factor including area, determinant of Jacobian and IP weight
//! @param rRow row, where to start to add the submatrix
//! @param rCol col, where to start to add the submatrix
//! @param rResult result
void NuTo::Element1D::AddDetJBtXB(const Eigen::MatrixXd &rDerivativeShapeFunctions1, const Eigen::MatrixXd &rDerivativeShapeFunctions2, const ConstitutiveTangentLocal<1, 1> &rConstitutiveTangent, double rFactor, int rRow, int rCol, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> &rResult) const
{
    int NumRowDofs = rDerivativeShapeFunctions1.rows();
    int NumColDofs = rDerivativeShapeFunctions2.rows();

    // add assertions!!!

    rFactor *= rConstitutiveTangent(0, 0);

    Eigen::MatrixXd result = rFactor * rDerivativeShapeFunctions1 * rDerivativeShapeFunctions2.transpose();

    rResult.block(rRow, rCol, NumRowDofs, NumColDofs) += result;
}

//! @brief calculates the Kee matrix
//! @param rShapeFunctions of the ip for all shape functions
//! @param rDerivativeShapeFunctions of the ip for all shape functions
//! @param nonlocal gradient radius xi
//! @param rFactor multiplication factor (detJ area..)
//! @param Kee return matrix with detJ * (Nt N + cBtB)
void NuTo::Element1D::CalculateKee(Eigen::VectorXd rShapeFunctions, const Eigen::MatrixXd& rDerivativeShapeFunctions, ConstitutiveTangentLocal<1, 1>& rNonlocalParameter, double rFactor, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rKee) const
{
    assert(rShapeFunctions.cols() == 1);
    assert(rDerivativeShapeFunctions.cols() == 1);

    rKee.Resize(rShapeFunctions.rows(), rShapeFunctions.rows());

    rKee = rFactor * (rShapeFunctions * rShapeFunctions.transpose() / rNonlocalParameter(0, 0) + rDerivativeShapeFunctions * rDerivativeShapeFunctions.transpose());

}

//! @brief add Kee*nonlocalEqStrain-detJ*N.T*localEqStrain (detJ is already included in Kkk)
//! @param rShapeFunctions of the ip for all shape functions
//! @param rLocalEqStrain local eq. strain values
//! @param rKkk stiffness matrix Kkk
//! @param rNodeNonlocalEqStrain nodal nonlocal eq strain values
//! @param rFactor factor including detJ and area
//! @param rResult result
void NuTo::Element1D::AddDetJRnonlocalEqStrain(const Eigen::VectorXd& rShapeFunctions, LocalEqStrain& rLocalEqStrain, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rKee, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& rNodeNonlocalEqStrain, double rFactor,
        int startIndexNonlocalEqStrain, FullVector<double, Eigen::Dynamic>& rResult) const
{
    assert(rResult.GetNumRows() >= (int )(startIndexNonlocalEqStrain + rShapeFunctions.size()));
    assert(rShapeFunctions.size() == rNodeNonlocalEqStrain.size());

    // perform Kee * nodeNonlocalEqStrain
    rResult.segment(startIndexNonlocalEqStrain, rShapeFunctions.size()) += rKee * rNodeNonlocalEqStrain.transpose() - rLocalEqStrain.GetValue(0) * rFactor * rShapeFunctions;

}

double NuTo::Element1D::CalculateJacobian(const Eigen::MatrixXd& rDerivativeShapeFunctions, const Eigen::MatrixXd& rNodeCoordinates) const
{

    assert(rDerivativeShapeFunctions.rows() == rNodeCoordinates.cols());
    assert(rDerivativeShapeFunctions.cols() == 1);
    assert(rNodeCoordinates.rows() == 1);

    return (rNodeCoordinates * rDerivativeShapeFunctions).at(0, 0);
}

void NuTo::Element1D::CheckElement()
{

    int numIntegrationPoints = GetNumIntegrationPoints();
    // check number of integration points
    if (numIntegrationPoints < 1)
    {
        throw MechanicsException("[NuTo::Element1D::CheckElement] invalid integration type.");
    }

    int theIP = 0;
    const Eigen::MatrixXd& derivativeShapeFunctions = mInterpolationType->Get(Node::COORDINATES).GetDerivativeShapeFunctionsNatural(theIP);
    Eigen::MatrixXd nodeCoordinates = ExtractNodeValues(0, Node::COORDINATES);
    double detJacobian = CalculateJacobian(derivativeShapeFunctions, nodeCoordinates);
    if (detJacobian < 0)
    {
        this->ReorderNodes();
        // recalculate node coordinates after reordering
        nodeCoordinates = ExtractNodeValues(0, Node::COORDINATES);
    }

    double length = 0;
    for (int iIP = 0; iIP < numIntegrationPoints; ++iIP)
    {
        const Eigen::MatrixXd& derivativeShapeFunctions = mInterpolationType->Get(Node::COORDINATES).GetDerivativeShapeFunctionsNatural(iIP);
        double detJacobian = CalculateJacobian(derivativeShapeFunctions, nodeCoordinates);
        if (detJacobian <= 0)
            throw MechanicsException("[NuTo::Element1D::CheckElement] Determinant of the Jacobian <= zero, no inversion possible.");
        length += this->GetIntegrationPointWeight(iIP) * detJacobian;
    }

    // check element length
    if (length < 1e-14)
    {
        throw MechanicsException("[NuTo::Truss1D::CheckElement] element with zero length (check nodes).");
    }

}

#ifdef ENABLE_SERIALIZATION
template void NuTo::Element1D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Element1D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Element1D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Element1D::save(Archive & ar, const unsigned int version)const
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Element1D " << std::endl;
#endif
    ar & boost::serialization::make_nvp("Element1D_ElementBase",boost::serialization::base_object<ElementBase >(*this));
    ar & boost::serialization::make_nvp("mSection", const_cast<SectionBase*&>(mSection));

    const std::uintptr_t* mNodesAddress = reinterpret_cast<const std::uintptr_t*>(mNodes.data());
    int size = mNodes.size();
    ar & boost::serialization::make_nvp("mNodes_size", size);
    ar & boost::serialization::make_nvp("mNodes", boost::serialization::make_array(mNodesAddress, size));
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Element1D" << std::endl;
#endif
}

template void NuTo::Element1D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Element1D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Element1D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Element1D::load(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start deserialize Element1D " << std::endl;
#endif
    ar & boost::serialization::make_nvp("Element1D_ElementBase",boost::serialization::base_object<ElementBase >(*this));
    ar & boost::serialization::make_nvp("mSection", const_cast<SectionBase*&>(mSection));

    int size = 0;
    ar & boost::serialization::make_nvp("mNodes_size", size);
    std::uintptr_t* mNodesAddress = new std::uintptr_t[size];
    ar & boost::serialization::make_nvp("mNodes", boost::serialization::make_array(mNodesAddress, size));
    mNodes.assign(reinterpret_cast<NodeBase**>(&mNodesAddress[0]), reinterpret_cast<NodeBase**>(&mNodesAddress[size]));
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish deserialize Element1D" << std::endl;
#endif
}

void NuTo::Element1D::SetNodePtrAfterSerialization(const std::map<std::uintptr_t, std::uintptr_t>& mNodeMapCast)
{
    for(std::vector<NodeBase*>::iterator it = mNodes.begin(); it != mNodes.end(); it++)
    {
        std::map<std::uintptr_t, std::uintptr_t>::const_iterator itCast = mNodeMapCast.find(reinterpret_cast<std::uintptr_t>(*it));
        if (itCast!=mNodeMapCast.end())
        {
            *it = reinterpret_cast<NodeBase*>(itCast->second);
        }
        else
            throw MechanicsException("[NuTo::Element1D] The NodeBase-Pointer could not be updated.");
    }
}

BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Element1D)
#endif // ENABLE_SERIALIZATION

