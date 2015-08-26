/*
 * Element2D.cpp
 *
 *  Created on: 9 Apr 2015
 *      Author: ttitsche
 */

#include "nuto/mechanics/elements/Element2D.h"
#include "nuto/mechanics/constitutive/mechanics/Damage.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/thermal/HeatFlux2D.h"
#include "nuto/mechanics/constitutive/thermal/HeatFlux3D.h"
#include "nuto/mechanics/constitutive/thermal/Temperature.h"
#include "nuto/mechanics/constitutive/thermal/TemperatureGradient2D.h"
#include "nuto/mechanics/constitutive/thermal/TemperatureGradient3D.h"
#include "nuto/mechanics/constitutive/mechanics/LocalEqStrain.h"
#include "nuto/mechanics/constitutive/mechanics/NonlocalEqStrain.h"
#include "nuto/mechanics/constitutive/moistureTransport/RelativeHumidity.h"
#include "nuto/mechanics/constitutive/moistureTransport/RelativeHumidityGradient2D.h"
#include "nuto/mechanics/constitutive/moistureTransport/WaterVolumeFraction.h"
#include "nuto/mechanics/constitutive/moistureTransport/WaterVolumeFractionGradient2D.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentNonlocal.h"
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixInt.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixDouble.h"
#include "nuto/mechanics/elements/ElementOutputVectorInt.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/structures/StructureBase.h"

#include "nuto/math/FullMatrix.h"


NuTo::Element2D::Element2D(const NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase*>& rNodes, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType, InterpolationType* rInterpolationType) :
        NuTo::ElementBase::ElementBase(rStructure, rElementDataType, rIpDataType, rInterpolationType),
        mNodes(rNodes),
        mSection(0)
{
}

NuTo::Error::eError NuTo::Element2D::Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput)
{
    if (mStructure->GetHessianConstant(1) == false)
        throw MechanicsException("[NuTo::Element2D::Evaluate] only implemented for a constant Hessian for the first derivative (damping).");
    if (mStructure->GetHessianConstant(2) == false)
        throw MechanicsException("[NuTo::Element2D::Evaluate] only implemented for a constant Hessian for the second derivative (mass).");

    try
    {
        const SectionBase* section(GetSection());
        if (section == 0)
            throw MechanicsException("[NuTo::Element2D::Evaluate] no section allocated for element.");

        const std::set<Node::eAttributes>& dofs = mInterpolationType->GetDofs();
        const std::set<Node::eAttributes>& activeDofs = mInterpolationType->GetActiveDofs();

        int numActiveDofs = mInterpolationType->GetNumActiveDofs();

        // extract all node values and store them
        std::map<Node::eAttributes, Eigen::MatrixXd> nodalValues;
        for (auto dof : dofs)
        {
            ExtractNodeValues(nodalValues[dof], 0, dof);
        }





        /*****************************************\
         *    CONSTITUTIVE INPUT DECLARATION     *
        \*****************************************/

        std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*> constitutiveInputList;


        EngineeringStress2D engineeringStress2D;
        EngineeringStress3D engineeringStress3D;
        EngineeringStrain2D engineeringStrain2D;
        EngineeringStrain3D engineeringStrain3D;

        EngineeringStrain3D engineeringPlasticStrain3D;

        //allocate damage
        Damage damage;

        // Moisture Transport
        // ------------------

        RelativeHumidity                relativeHumidity;
        RelativeHumidity                relativeHumidityD1;
        RelativeHumidityGradient2D      relativeHumidityGradient;
        WaterVolumeFraction             waterVolumeFraction;
        WaterVolumeFraction             waterVolumeFractionD1;
        WaterVolumeFractionGradient2D   waterVolumeFractionGradient;

        //allocate nonlocal eq strain
        NonlocalEqStrain nonlocalEqStrain;

        //allocate local eq strain
        LocalEqStrain localEqStrain;




        /*****************************************\
         *    CONSTITUTIVE OUTPUT DECLARATION    *
        \*****************************************/

        std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*> constitutiveOutputList;

        //allocate transient nonlocal parameter
        ConstitutiveTangentLocal<1, 1> nonlocalParameter;
        ConstitutiveTangentLocal<3, 1> tangentStressNonlocalEqStrain;
        ConstitutiveTangentLocal<3, 1> tangentLocalEqStrainStrain;

        //allocate vector of tangent matrices
        int NumNonlocalElements = GetNumNonlocalElements();
        int NumActiveDofsNonlocal = numActiveDofs;
        int NumNonlocalIps = 0;

        if (NumNonlocalElements != 0)
        {
            // replace the local displacement dofs with the nonlocal displacement dofs
            NumActiveDofsNonlocal -= mInterpolationType->Get(Node::DISPLACEMENTS).GetNumDofs();
            for (const ElementBase* nonlocalElement : GetNonlocalElements())
            {
                NumNonlocalIps += nonlocalElement->GetNumIntegrationPoints();
                NumActiveDofsNonlocal += nonlocalElement->GetNumNodes(Node::DISPLACEMENTS) * 2;
            }
        }

        NuTo::ConstitutiveTangentNonlocal<3, 3> nonlocalTangentStressStrain;

        DeformationGradient2D deformationGradient;

        // Moisture Transport
        // ------------------

        // Naming scheme for matrices: tangent_D_X_D_Y_HN_AB
        // D_X_D_Y:     defines the derivative of X with respect to Y
        // HN:          defines the matrix the tangent belongs to, for example H0 for stiffness or H1 for damping
        // AB:          defines the combination of (derivative) shapefunctions the tangent has to be multiplied with, for example BN means: derivative shapefunctions * tangent * shapefunctions

        // Internal Gradient
        ConstitutiveTangentLocal<1,1> residualWaterPhaseN;
        ConstitutiveTangentLocal<2,1> residualWaterPhaseB;
        ConstitutiveTangentLocal<1,1> residualVaporPhaseN;
        ConstitutiveTangentLocal<2,1> residualVaporPhaseB;
        // Hessian 0
        ConstitutiveTangentLocal<1,1> tangent_D_Residual_RH_D_RH_H0_BB;
        ConstitutiveTangentLocal<1,1> tangent_D_Residual_RH_D_RH_H0_NN;
        ConstitutiveTangentLocal<2,1> tangent_D_Residual_RH_D_WV_H0_BN;
        ConstitutiveTangentLocal<1,1> tangent_D_Residual_RH_D_WV_H0_NN;
        ConstitutiveTangentLocal<1,1> tangent_D_Residual_WV_D_RH_H0_NN;
        ConstitutiveTangentLocal<1,1> tangent_D_Residual_WV_D_WV_H0_BB;
        ConstitutiveTangentLocal<2,1> tangent_D_Residual_WV_D_WV_H0_BN;
        ConstitutiveTangentLocal<1,1> tangent_D_Residual_WV_D_WV_H0_NN;
        // Hessian 1
        ConstitutiveTangentLocal<1,1> tangent_D_Residual_RH_D_RH_H1_NN;
        ConstitutiveTangentLocal<1,1> tangent_D_Residual_RH_D_WV_H1_NN;
        ConstitutiveTangentLocal<1,1> tangent_D_Residual_WV_D_WV_H1_NN;

        ConstitutiveTangentLocal<1,1> residualNormFactorDisplacements;
        ConstitutiveTangentLocal<1,1> residualNormFactorRelativeHumidity;
        ConstitutiveTangentLocal<1,1> residualNormFactorWaterVolumeFraction;

        //for the lumped mass calculation
        double total_mass = 0.;





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
                constitutiveInputList[NuTo::Constitutive::Input::DEFORMATION_GRADIENT_2D] = &(deformationGradient);
            }
                break;
            case Node::NONLOCALEQSTRAIN:
            {
                constitutiveInputList[NuTo::Constitutive::Input::NONLOCAL_EQ_STRAIN] = &(nonlocalEqStrain);
            }
                break;
            case Node::RELATIVEHUMIDITY:
            {
                constitutiveInputList[NuTo::Constitutive::Input::RELATIVE_HUMIDITY]                 = &(relativeHumidity);
                constitutiveInputList[NuTo::Constitutive::Input::RELATIVE_HUMIDITY_D1]              = &(relativeHumidityD1);
                constitutiveInputList[NuTo::Constitutive::Input::RELATIVE_HUMIDITY_GRADIENT]        = &(relativeHumidityGradient);
            }
                break;
            case Node::WATERVOLUMEFRACTION:
            {
                constitutiveInputList[NuTo::Constitutive::Input::WATER_VOLUME_FRACTION]             = &(waterVolumeFraction);
                constitutiveInputList[NuTo::Constitutive::Input::WATER_VOLUME_FRACTION_D1]          = &(waterVolumeFractionD1);
                constitutiveInputList[NuTo::Constitutive::Input::WATER_VOLUME_FRACTION_GRADIENT]    = &(waterVolumeFractionGradient);
            }
                break;
            default:
                throw MechanicsException("[NuTo::Element2D::Evaluate] Constitutive input for " + Node::AttributeToString(dof) + " not implemented.");
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
                            constitutiveOutputList[NuTo::Constitutive::Output::ENGINEERING_STRESS_2D] = &(engineeringStress2D);
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
                                constitutiveOutputList[NuTo::Constitutive::Output::RESIDUAL_VAPOR_PHASE_N]                                      = &residualVaporPhaseN;
                                constitutiveOutputList[NuTo::Constitutive::Output::RESIDUAL_VAPOR_PHASE_B]                                      = &residualVaporPhaseB;
                            }
                        }
                            break;

                        case Node::WATERVOLUMEFRACTION:
                        {
                            if (activeDofs.find(Node::RELATIVEHUMIDITY) != activeDofs.end())
                            {
                                constitutiveOutputList[NuTo::Constitutive::Output::RESIDUAL_WATER_PHASE_N]                                      = &residualWaterPhaseN;
                                constitutiveOutputList[NuTo::Constitutive::Output::RESIDUAL_WATER_PHASE_B]                                      = &residualWaterPhaseB;
                            }
                        }
                            break;
                        default:
                            throw MechanicsException("[NuTo::Element2D::Evaluate] Constitutive output INTERNAL_GRADIENT for " + Node::AttributeToString(dof) + " not implemented.");

                        }
                    }
                }
                break;
            case Element::HESSIAN_0_TIME_DERIVATIVE:
            {
                it->second->GetFullMatrixDouble().Resize(numActiveDofs, NumActiveDofsNonlocal);
                it->second->GetFullMatrixDouble().setZero();
                it->second->SetSymmetry(true);
                it->second->SetConstant(true);
                for (auto dof : activeDofs)
                {
                    switch (dof)
                    {
                    case Node::DISPLACEMENTS:
                    {
                        if (NumNonlocalElements == 0)
                        {
                            nonlocalTangentStressStrain.SetNumSubMatrices(1);
                            constitutiveOutputList[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_2D] = &(nonlocalTangentStressStrain.GetSubMatrix_3x3(0));
                        } else
                        {
                            nonlocalTangentStressStrain.SetNumSubMatrices(NumNonlocalIps);
                            constitutiveOutputList[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_2D] = &nonlocalTangentStressStrain;
                        }

                        if (activeDofs.find(Node::NONLOCALEQSTRAIN) != activeDofs.end())
                        {
                            constitutiveOutputList[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_2D] = &tangentStressNonlocalEqStrain;
                        }
                    }
                        break;
                    case Node::NONLOCALEQSTRAIN:
                    {
                        constitutiveOutputList[NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_XI_D_STRAIN_2D] = &tangentLocalEqStrainStrain;
                        constitutiveOutputList[NuTo::Constitutive::Output::NONLOCAL_PARAMETER_XI] = &nonlocalParameter;
                    }
                        break;
                    case Node::RELATIVEHUMIDITY:
                    {
                        if (activeDofs.find(Node::WATERVOLUMEFRACTION) != activeDofs.end())
                        {
                            constitutiveOutputList[NuTo::Constitutive::Output::D_RESIDUAL_RH_D_RH_H0_BB]                                            = &tangent_D_Residual_RH_D_RH_H0_BB;
                            constitutiveOutputList[NuTo::Constitutive::Output::D_RESIDUAL_RH_D_RH_H0_NN]                                            = &tangent_D_Residual_RH_D_RH_H0_NN;
                            constitutiveOutputList[NuTo::Constitutive::Output::D_RESIDUAL_RH_D_WV_H0_BN]                                            = &tangent_D_Residual_RH_D_WV_H0_BN;
                            constitutiveOutputList[NuTo::Constitutive::Output::D_RESIDUAL_RH_D_WV_H0_NN]                                            = &tangent_D_Residual_RH_D_WV_H0_NN;
                        }
                    }
                        break;
                    case Node::WATERVOLUMEFRACTION:
                    {
                        if (activeDofs.find(Node::RELATIVEHUMIDITY) != activeDofs.end())
                        {
                            constitutiveOutputList[NuTo::Constitutive::Output::D_RESIDUAL_WV_D_RH_H0_NN]                                            = &tangent_D_Residual_WV_D_RH_H0_NN;
                            constitutiveOutputList[NuTo::Constitutive::Output::D_RESIDUAL_WV_D_WV_H0_BB]                                            = &tangent_D_Residual_WV_D_WV_H0_BB;
                            constitutiveOutputList[NuTo::Constitutive::Output::D_RESIDUAL_WV_D_WV_H0_BN]                                            = &tangent_D_Residual_WV_D_WV_H0_BN;
                            constitutiveOutputList[NuTo::Constitutive::Output::D_RESIDUAL_WV_D_WV_H0_NN]                                            = &tangent_D_Residual_WV_D_WV_H0_NN;
                        }
                    }
                        break;
                    default:
                        throw MechanicsException("[NuTo::Element2D::Evaluate] Constitutive output HESSIAN_0_TIME_DERIVATIVE for " + Node::AttributeToString(dof) + " not implemented.");

                    }
                }
            }
                break;
            case Element::HESSIAN_1_TIME_DERIVATIVE:
                it->second->GetFullMatrixDouble().Resize(numActiveDofs, NumActiveDofsNonlocal);
                it->second->SetSymmetry(true);
                it->second->SetConstant(true);
                for (auto dof : activeDofs)
                {
                    switch (dof)
                    {
                    case Node::DISPLACEMENTS:
                    {
                        break;
                    }
                    case Node::RELATIVEHUMIDITY:
                    {
                        if (activeDofs.find(Node::WATERVOLUMEFRACTION) != activeDofs.end())
                        {
                            constitutiveOutputList[NuTo::Constitutive::Output::D_RESIDUAL_RH_D_RH_H1_NN]                                            = &tangent_D_Residual_RH_D_RH_H1_NN;
                            constitutiveOutputList[NuTo::Constitutive::Output::D_RESIDUAL_RH_D_WV_H1_NN]                                            = &tangent_D_Residual_RH_D_WV_H1_NN;
                        }
                    }
                        break;
                    case Node::WATERVOLUMEFRACTION:
                    {
                        if (activeDofs.find(Node::RELATIVEHUMIDITY) != activeDofs.end())
                        {
                            constitutiveOutputList[NuTo::Constitutive::Output::D_RESIDUAL_WV_D_WV_H1_NN]                                            = &tangent_D_Residual_WV_D_WV_H1_NN;
                        }
                    }
                        break;
                    default:
                    {
                        throw MechanicsException("[NuTo::Element2D::Evaluate] Constitutive output HESSIAN_1_TIME_DERIVATIVE for " + Node::AttributeToString(dof) + " not implemented.");
                    }
                    }
                }
                break;
            case Element::HESSIAN_2_TIME_DERIVATIVE:
            {
                it->second->GetFullMatrixDouble().Resize(numActiveDofs, NumActiveDofsNonlocal);
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
                default:
                    throw MechanicsException("[NuTo::Element2D::Evaluate] this ip data type is not implemented.");
                }
                break;
            case Element::RESIDUAL_NORM_FACTOR:
            {
                it->second->GetFullVectorDouble().Resize(numActiveDofs);
                for (auto dof : activeDofs)
                {
                    switch (dof)
                    {
                    case Node::DISPLACEMENTS:
                    {
                        constitutiveOutputList[NuTo::Constitutive::Output::RESIDUAL_NORM_FACTOR_DISPLACEMENTS]                                      = &residualNormFactorDisplacements;
                        break;
                    }
                    case Node::RELATIVEHUMIDITY:
                    {
                        if (activeDofs.find(Node::WATERVOLUMEFRACTION) != activeDofs.end())
                        {
                            constitutiveOutputList[NuTo::Constitutive::Output::RESIDUAL_NORM_FACTOR_RELATIVE_HUMIDITY]                              = &residualNormFactorRelativeHumidity;
                        }
                    }
                        break;
                    case Node::WATERVOLUMEFRACTION:
                    {
                        if (activeDofs.find(Node::RELATIVEHUMIDITY) != activeDofs.end())
                        {
                            constitutiveOutputList[NuTo::Constitutive::Output::RESIDUAL_NORM_FACTOR_WATER_VOLUME_FRACTION]                          = &residualNormFactorWaterVolumeFraction;
                        }
                    }
                        break;
                    default:
                    {
                        throw MechanicsException("[NuTo::Element2D::Evaluate] Constitutive output RESIDUAL_NORM_FACTOR for " + Node::AttributeToString(dof) + " not implemented.");
                    }
                    }
                }
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
                throw MechanicsException("[NuTo::Element2D::Evaluate] element output not implemented.");
            }
        }





        /*****************************************\
         *     CALCULATE CONSTITUTIVE INPUTS     *
        \*****************************************/


        Eigen::Matrix2d invJacobian;
        double detJacobian;

        std::map<Node::eAttributes, Eigen::VectorXd> shapeFunctions;
        std::map<Node::eAttributes, Eigen::MatrixXd> derivativeShapeFunctions;

        // loop over the integration points
        for (int theIP = 0; theIP < GetNumIntegrationPoints(); theIP++)
        {

            // calculate Jacobian
            const Eigen::MatrixXd& derivativeShapeFunctionsGeometryNatural = mInterpolationType->Get(Node::COORDINATES).GetDerivativeShapeFunctionsNatural(theIP);

            this->CalculateJacobian(derivativeShapeFunctionsGeometryNatural, nodalValues[Node::COORDINATES], invJacobian, detJacobian);

            double factor(mSection->GetThickness() * (detJacobian) * (mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP)));
            // calculate shape functions and their derivatives

            for (auto dof : dofs)
            {
                if (dof == Node::COORDINATES)
                    continue;
                const InterpolationBase& interpolationType = mInterpolationType->Get(dof);
                shapeFunctions[dof] = interpolationType.GetShapeFunctions(theIP);
                // this lazy product here is so much faster than any other implementation via a seperate method
                // possibly due to more efficient mallocs
                derivativeShapeFunctions[dof] = interpolationType.GetDerivativeShapeFunctionsNatural(theIP).lazyProduct(invJacobian);
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
//                    constitutiveInputList[NuTo::Constitutive::Input::DEFORMATION_GRADIENT_2D] = &(deformationGradient);
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
                    relativeHumidity(0,0)        = (nodalValues[Node::RELATIVEHUMIDITY] * shapeFunctions[Node::RELATIVEHUMIDITY])(0,0);
                    relativeHumidityD1(0,0)      = (ExtractNodeValues(1, Node::RELATIVEHUMIDITY) * shapeFunctions[Node::RELATIVEHUMIDITY])(0,0);

                    relativeHumidityGradient.setZero();
                    relativeHumidityGradient.AddBlock(0,0,derivativeShapeFunctions[Node::RELATIVEHUMIDITY].transpose() * nodalValues[Node::RELATIVEHUMIDITY].transpose());
                }
                    break;
                case Node::WATERVOLUMEFRACTION:
                {
                    waterVolumeFraction(0,0)         = (nodalValues[Node::WATERVOLUMEFRACTION] * shapeFunctions[Node::WATERVOLUMEFRACTION])(0,0);
                    waterVolumeFractionD1(0,0)       = (ExtractNodeValues(1, Node::WATERVOLUMEFRACTION) * shapeFunctions[Node::WATERVOLUMEFRACTION])(0);

                    waterVolumeFractionGradient.setZero();
                    waterVolumeFractionGradient.AddBlock(0,0, derivativeShapeFunctions[Node::WATERVOLUMEFRACTION].transpose() * nodalValues[Node::WATERVOLUMEFRACTION].transpose());

                }
                    break;
                default:
                    throw MechanicsException("[NuTo::Element2D::Evaluate] Constitutive input for " + Node::AttributeToString(dof) + " not implemented.");
                }
            }





            /*****************************************\
             *      EVALUATE CONSTITUTIVE LAW        *
            \*****************************************/

            ConstitutiveBase* constitutivePtr = GetConstitutiveLaw(theIP);
            try
            {
                Error::eError error = constitutivePtr->Evaluate2D(this, theIP, constitutiveInputList, constitutiveOutputList);
                if (error != Error::SUCCESSFUL)
                    return error;
            } catch (NuTo::MechanicsException &e)
            {
                e.AddMessage("[NuTo::Element2D::Evaluate] error evaluating the constitutive model.");
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
                        double factor(mSection->GetThickness() * detJacobian * (mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP)));
                        for (auto dof : activeDofs)
                        {
                            int startIndex = mInterpolationType->Get(dof).GetLocalStartIndex();
                            switch (dof)
                            {
                            case Node::DISPLACEMENTS:
                            {
                                AddDetJBtSigma(derivativeShapeFunctions.at(dof), engineeringStress2D, factor, startIndex, it->second->GetFullVectorDouble());
                            }
                                break;
                            case Node::NONLOCALEQSTRAIN:
                            {
                                int startIndexNonlocalEqStrain = mInterpolationType->Get(Node::NONLOCALEQSTRAIN).GetLocalStartIndex();
                                AddDetJRnonlocalEqStrain(shapeFunctions.at(Node::NONLOCALEQSTRAIN), localEqStrain, Kee, nodalValues.at(Node::NONLOCALEQSTRAIN), factor / nonlocalParameter.GetValue(0),
                                        startIndexNonlocalEqStrain, it->second->GetFullVectorDouble());
                            }
                                break;
                            case Node::RELATIVEHUMIDITY:
                            {
                                if(activeDofs.find(Node::WATERVOLUMEFRACTION) != activeDofs.end())
                                {
                                    AddDetJBtX( derivativeShapeFunctions.at(Node::RELATIVEHUMIDITY),
                                                residualVaporPhaseB,
                                                factor,
                                                startIndex,
                                                it->second->GetFullVectorDouble());

                                    AddDetJNtX( shapeFunctions.at(Node::RELATIVEHUMIDITY),
                                                residualVaporPhaseN,
                                                factor,
                                                startIndex,
                                                it->second->GetFullVectorDouble());
                                }
                            }
                                break;
                            case Node::WATERVOLUMEFRACTION:
                            {
                                if(activeDofs.find(Node::RELATIVEHUMIDITY) != activeDofs.end())
                                {
                                    AddDetJBtX( derivativeShapeFunctions.at(Node::WATERVOLUMEFRACTION),
                                                residualWaterPhaseB,
                                                factor,
                                                startIndex,
                                                it->second->GetFullVectorDouble());

                                    AddDetJNtX( shapeFunctions.at(Node::WATERVOLUMEFRACTION),
                                                residualWaterPhaseN,
                                                factor,
                                                startIndex,
                                                it->second->GetFullVectorDouble());
                                }
                            }
                                break;
                            default:
                                throw MechanicsException("[NuTo::Element2D::Evaluate] Element output INTERNAL_GRADIENT for " + Node::AttributeToString(dof) + " not implemented.");

                            }
                        }
                    }
                }
                    break;
                case Element::HESSIAN_0_TIME_DERIVATIVE:
                 {
                     //factor for the numerical integration
                     assert(mSection->GetThickness() > 0);
                     double factor(mSection->GetThickness() * detJacobian * (mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP)));



                     for (auto dof : activeDofs)
                    {
                        int startIndex = mInterpolationType->Get(dof).GetLocalStartIndex();
                        switch (dof)
                        {
                        case Node::DISPLACEMENTS:
                        {
                            if (NumNonlocalElements != 0)
                            {

                                if (NumNonlocalElements != 0)
                                {
                                    if (nonlocalTangentStressStrain.GetConstant() == false)
                                        it->second->SetConstant(false);

                                    if (nonlocalTangentStressStrain.GetSymmetry() == false)
                                        it->second->SetSymmetry(false);
                                } else
                                {

                                    if (nonlocalTangentStressStrain.GetSubMatrix_3x3(0).GetConstant() == false)
                                        it->second->SetConstant(false);

                                    if (nonlocalTangentStressStrain.GetSubMatrix_3x3(0).GetSymmetry() == false)
                                        it->second->SetSymmetry(false);
                                }

                                //Nonlocal Model where the material model is still local
                                if (nonlocalTangentStressStrain.GetLocalSolution())
                                {
                                    //same as local model, e.g. in the unloading range or

                                    // calculate element stiffness matrix
                                    // don't forget to include determinant of the Jacobian and area
                                    AddDetJBtCB(derivativeShapeFunctions.at(dof), derivativeShapeFunctions.at(dof), nonlocalTangentStressStrain.GetSubMatrix_3x3(0), factor,
                                            it->second->GetFullMatrixDouble(), startIndex, startIndex);
                                } else
                                {
                                    //nonlocal BMatrix of nonlocal other integration point
                                    //sum over nonlocal elements and their ips
                                    const std::vector<const ElementBase*>& nonlocalElements(GetNonlocalElements());
                                    int firstCol(0);
                                    int totalNonlocalIp(0);
                                    Eigen::MatrixXd nonlocalLocalNodeCoord;
                                    Eigen::MatrixXd nonlocalDerivativeShapeFunctionsNatural;
                                    Eigen::MatrixXd nonlocalDerivativeShapeFunctionsLocal;

                                    Eigen::Matrix2d nonlocalInvJacobian;
                                    double nonlocalDetJacobian;

                                    for (unsigned int theNonlocalElement = 0; theNonlocalElement < nonlocalElements.size(); theNonlocalElement++)
                                    {
                                        const Element2D* nonlocalElement = nonlocalElements[theNonlocalElement]->AsElement2D();
                                        const InterpolationType* interpolationType = nonlocalElement->GetInterpolationType();
                                        //calculate local coordinates
                                        nonlocalLocalNodeCoord = nonlocalElement->ExtractNodeValues(0, Node::COORDINATES);

                                        //get weights
                                        const std::vector<double>& weights(GetNonlocalWeights(theIP, theNonlocalElement));

                                        for (int theNonlocalIp = 0; theNonlocalIp < nonlocalElement->GetNumIntegrationPoints(); theNonlocalIp++, totalNonlocalIp++)
                                        {
                                            if (weights[theNonlocalIp] == 0.)
                                                continue;

                                            // calculate derivatives of shape functions in natural coordinate system
                                            nonlocalDerivativeShapeFunctionsNatural = interpolationType->Get(Node::COORDINATES).GetDerivativeShapeFunctionsNatural(theNonlocalIp);

                                            CalculateJacobian(nonlocalDerivativeShapeFunctionsNatural, nonlocalLocalNodeCoord, nonlocalInvJacobian, nonlocalDetJacobian);

                                            nonlocalDerivativeShapeFunctionsLocal = interpolationType->Get(Node::DISPLACEMENTS).GetDerivativeShapeFunctionsNatural(theNonlocalIp).lazyProduct(
                                                    nonlocalInvJacobian);

                                            // calculate element stiffness matrix
                                            // don't forget to include determinant of the Jacobian and area
                                            AddDetJBtCB(derivativeShapeFunctions.at(Node::DISPLACEMENTS), nonlocalDerivativeShapeFunctionsLocal,
                                                    nonlocalTangentStressStrain.GetSubMatrix_3x3(totalNonlocalIp), factor, it->second->GetFullMatrixDouble(), 0, firstCol);
                                        }
                                        firstCol += interpolationType->Get(Node::DISPLACEMENTS).GetNumDofs();
                                    }
                                    assert(totalNonlocalIp == NumNonlocalIps);
                                }
                            } else
                            {
                                // calculate element stiffness matrix
                                // don't forget to include determinant of the Jacobian and area
                                AddDetJBtCB(derivativeShapeFunctions.at(dof), derivativeShapeFunctions.at(dof), nonlocalTangentStressStrain.GetSubMatrix_3x3(0), factor,
                                        it->second->GetFullMatrixDouble(), startIndex, startIndex);
                            }

                            if (activeDofs.find(Node::NONLOCALEQSTRAIN) != activeDofs.end())
                            {
                                int startIndexNonlocalEqStrain = mInterpolationType->Get(Node::NONLOCALEQSTRAIN).GetLocalStartIndex();
                                AddDetJBtdSigmadNonlocalEqStrainN(derivativeShapeFunctions.at(Node::DISPLACEMENTS), tangentStressNonlocalEqStrain, shapeFunctions.at(Node::NONLOCALEQSTRAIN), factor,
                                        startIndex, startIndexNonlocalEqStrain, it->second->GetFullMatrixDouble());
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
                                AddDetJNtdLocalEqStraindEpsilonB(shapeFunctions.at(Node::NONLOCALEQSTRAIN), tangentLocalEqStrainStrain, derivativeShapeFunctions.at(Node::DISPLACEMENTS), factor,
                                        startIndexNonlocalEqStrain, startIndexDisplacement, it->second->GetFullMatrixDouble());
                            }

                        }
                            break;
                        case Node::RELATIVEHUMIDITY:
                        {
                                                        if(activeDofs.find(Node::WATERVOLUMEFRACTION) != activeDofs.end())
                            {
                                int indexRelHum = startIndex;
                                int indexWatVol = mInterpolationType->Get(Node::WATERVOLUMEFRACTION).GetLocalStartIndex();

                                auto& RelHumShapeFunction = shapeFunctions.at(Node::RELATIVEHUMIDITY);
                                auto& WatVolShapeFunction = shapeFunctions.at(Node::WATERVOLUMEFRACTION);

                                auto& RelHumDerivativeShapeFunction = derivativeShapeFunctions.at(Node::RELATIVEHUMIDITY);

                                // | - - |
                                // | - X |


                                AddDetJBtXB(RelHumDerivativeShapeFunction,
                                            RelHumDerivativeShapeFunction,
                                            tangent_D_Residual_RH_D_RH_H0_BB,
                                            factor,
                                            indexRelHum,
                                            indexRelHum,
                                            it->second->GetFullMatrixDouble());

                                AddDetJNtXN(RelHumShapeFunction,
                                            RelHumShapeFunction,
                                            tangent_D_Residual_RH_D_RH_H0_NN,
                                            factor,
                                            indexRelHum,
                                            indexRelHum,
                                            it->second->GetFullMatrixDouble());

                                // coupling terms

                                // | - - |
                                // | X - |

                                AddDetJBtXN(RelHumDerivativeShapeFunction,
                                            WatVolShapeFunction,
                                            tangent_D_Residual_RH_D_WV_H0_BN,
                                            factor,
                                            indexRelHum,
                                            indexWatVol,
                                            it->second->GetFullMatrixDouble());

                                AddDetJNtXN(RelHumShapeFunction,
                                            WatVolShapeFunction,
                                            tangent_D_Residual_RH_D_WV_H0_NN,
                                            factor,
                                            indexRelHum,
                                            indexWatVol,
                                            it->second->GetFullMatrixDouble());


                                // Maybe replace the tangent check for constness by a constitutive law function that tells if the build matrix is constant or not
                                if(it->second->GetConstant() &&
                                   (tangent_D_Residual_RH_D_RH_H0_BB.GetConstant() == false ||
                                    tangent_D_Residual_RH_D_RH_H0_NN.GetConstant() == false ||
                                    tangent_D_Residual_RH_D_WV_H0_BN.GetConstant() == false ||
                                    tangent_D_Residual_RH_D_WV_H0_NN.GetConstant() == false))
                                {
                                    it->second->SetConstant(false);
                                }
/*
                                //it->second->GetFullMatrixDouble().Info();         // Check Element Matrix if needed
                                //NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> test = it->second->GetFullMatrixDouble();
                                //test.Info(4,5,true);
                                //int a=0;
                                //a++;
*/
                            }
                        }
                            break;
                        case Node::WATERVOLUMEFRACTION:
                        {

                            if(activeDofs.find(Node::RELATIVEHUMIDITY) != activeDofs.end())
                            {
                                int indexRelHum = mInterpolationType->Get(Node::RELATIVEHUMIDITY).GetLocalStartIndex();
                                int indexWatVol = startIndex;

                                auto& RelHumShapeFunction = shapeFunctions.at(Node::RELATIVEHUMIDITY);
                                auto& WatVolShapeFunction = shapeFunctions.at(Node::WATERVOLUMEFRACTION);

                                auto& WatVolDerivativeShapeFunction = derivativeShapeFunctions.at(Node::WATERVOLUMEFRACTION);

                                // | X - |
                                // | - - |

                                AddDetJBtXB(WatVolDerivativeShapeFunction,
                                            WatVolDerivativeShapeFunction,
                                            tangent_D_Residual_WV_D_WV_H0_BB,
                                            factor,
                                            indexWatVol,
                                            indexWatVol,
                                            it->second->GetFullMatrixDouble());


                                AddDetJBtXN(WatVolDerivativeShapeFunction,
                                            WatVolShapeFunction,
                                            tangent_D_Residual_WV_D_WV_H0_BN,
                                            factor,
                                            indexWatVol,
                                            indexWatVol,
                                            it->second->GetFullMatrixDouble());

                                AddDetJNtXN(WatVolShapeFunction,
                                            WatVolShapeFunction,
                                            tangent_D_Residual_WV_D_WV_H0_NN,
                                            factor,
                                            indexWatVol,
                                            indexWatVol,
                                            it->second->GetFullMatrixDouble());


                                // coupling terms

                                // | - X |
                                // | - - |

                                AddDetJNtXN(WatVolShapeFunction,
                                            RelHumShapeFunction,
                                            tangent_D_Residual_WV_D_RH_H0_NN,
                                            factor,
                                            indexWatVol,
                                            indexRelHum,
                                            it->second->GetFullMatrixDouble());


                                // Maybe replace the tangent check for constness by a constitutive law function that tells if the build matrix is constant or not
                                if(it->second->GetConstant() &&
                                   (tangent_D_Residual_WV_D_RH_H0_NN.GetConstant() == false ||
                                    tangent_D_Residual_WV_D_WV_H0_BB.GetConstant() == false ||
                                    tangent_D_Residual_WV_D_WV_H0_BN.GetConstant() == false ||
                                    tangent_D_Residual_WV_D_WV_H0_NN.GetConstant() == false))
                                {
                                    it->second->SetConstant(false);
                                }

                                //it->second->GetFullMatrixDouble().Info();         // Check Element Matrix if needed
                                //NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> test = it->second->GetFullMatrixDouble();
                                //int a=0;
                            }
                        }
                            break;
                        default:
                            throw MechanicsException("[NuTo::Element2D::Evaluate] Element output HESSIAN_0_TIME_DERIVATIVE for " + Node::AttributeToString(dof) + " not implemented.");

                        }
                    }
                }
                    break;
                case Element::HESSIAN_1_TIME_DERIVATIVE:
                {
                //factor for the numerical integration
                //assert(mSection->GetArea() > 0);
                //if (tangentStressStrain.GetConstant() == false)
                //    it->second->SetConstant(false);
                //if (tangentStressStrain.GetSymmetry() == false)
                //    it->second->SetSymmetry(false);

                for (auto dof : activeDofs)
                {
                    int startIndex = mInterpolationType->Get(dof).GetLocalStartIndex();
                    switch (dof)
                    {
                    case Node::DISPLACEMENTS:
                    {
                        break;
                    }


                    case Node::RELATIVEHUMIDITY:
                    {

                        if(activeDofs.find(Node::WATERVOLUMEFRACTION) != activeDofs.end())
                        {
                            int indexRelHum = startIndex;
                            int indexWatVol = mInterpolationType->Get(Node::WATERVOLUMEFRACTION).GetLocalStartIndex();

                            auto RelHumShapeFunction = shapeFunctions.at(Node::RELATIVEHUMIDITY);
                            auto WatVolShapeFunction = shapeFunctions.at(Node::WATERVOLUMEFRACTION);

                            // | - - |
                            // | - X |

                            AddDetJNtXN(RelHumShapeFunction,
                                        RelHumShapeFunction,
                                        tangent_D_Residual_RH_D_RH_H1_NN,
                                        factor,
                                        indexRelHum,
                                        indexRelHum,
                                        it->second->GetFullMatrixDouble());

                            // coupling terms

                            // | - - |
                            // | X - |

                            AddDetJNtXN(RelHumShapeFunction,
                                        WatVolShapeFunction,
                                        tangent_D_Residual_RH_D_WV_H1_NN,
                                        factor,
                                        indexRelHum,
                                        indexWatVol,
                                        it->second->GetFullMatrixDouble());


                            // Maybe replace the tangent check for constness by a constitutive law function that tells if the build matrix is constant or not
                            if(it->second->GetConstant() &&
                               (tangent_D_Residual_RH_D_RH_H1_NN.GetConstant() == false ||
                                tangent_D_Residual_RH_D_WV_H1_NN.GetConstant() == false))
                            {
                                it->second->SetConstant(false);
                            }

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

                        if(activeDofs.find(Node::RELATIVEHUMIDITY) != activeDofs.end())
                        {
                            int indexWatVol = startIndex;

                            auto WatVolShapeFunction = shapeFunctions.at(Node::WATERVOLUMEFRACTION);


                            // | X - |
                            // | - - |

                            AddDetJNtXN(WatVolShapeFunction,
                                        WatVolShapeFunction,
                                        tangent_D_Residual_WV_D_WV_H1_NN,
                                        factor,
                                        indexWatVol,
                                        indexWatVol,
                                        it->second->GetFullMatrixDouble());

                            // Maybe replace the tangent check for constness by a constitutive law function that tells if the build matrix is constant or not
                            if(it->second->GetConstant() &&
                               tangent_D_Residual_WV_D_WV_H1_NN.GetConstant() == false)
                            {
                                it->second->SetConstant(false);
                            }

                            //it->second->GetFullMatrixDouble().Info();         // Check Element Matrix if needed
                            //NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> test = it->second->GetFullMatrixDouble();
                            //int a=0;
                        }
                    }
                        break;
                    default:
                        throw MechanicsException("[NuTo::Element2D::Evaluate] Element output HESSIAN_1_TIME_DERIVATIVE for " + Node::AttributeToString(dof) + " not implemented.");

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

                            double factor(mSection->GetThickness() * detJacobian
                                                                   * (mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP))
                                                                   * constitutivePtr->GetParameterDouble(Constitutive::eConstitutiveParameter::DENSITY));
                            const Eigen::MatrixXd tmpMatrix = shapeFunctions.at(dof) * shapeFunctions.at(dof).transpose() * factor;
                            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& result(it->second->GetFullMatrixDouble());
                            for (int count = 0; count < tmpMatrix.rows(); count++)
                            {
                                for (int count2 = 0; count2 < tmpMatrix.cols(); count2++)
                                {
                                    result(2 * count, 2 * count2) += tmpMatrix(count, count2);
                                    result(2 * count + 1, 2 * count2 + 1) += tmpMatrix(count, count2);
                                }
                            }
                        }
                            break;
                        case Node::RELATIVEHUMIDITY:
                        case Node::WATERVOLUMEFRACTION:
                        case Node::TEMPERATURES:
                            break;
                        default:
                            throw MechanicsException("[NuTo::Element2D::Evaluate] Element output HESSIAN_2_TIME_DERIVATIVE for " + Node::AttributeToString(dof) + " not implemented.");

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

                            // calculate local mass matrix (the nonlocal terms are zero)
                            // don't forget to include determinant of the Jacobian and area
                            // detJ * area * density * HtH, :
                            double factor(mSection->GetThickness() * detJacobian
                                                                   * (mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP))
                                                                   * constitutivePtr->GetParameterDouble(Constitutive::eConstitutiveParameter::DENSITY));
                            FullVector<double, Eigen::Dynamic>& result(it->second->GetFullVectorDouble());
                            total_mass += factor;
                            const Eigen::VectorXd& shapeFunction = shapeFunctions.at(dof);
                            //calculate for the translational dofs the diagonal entries
                            for (int count = 0; count < shapeFunction.rows(); count++)
                            {
                                result(2 * count) += shapeFunction[count] * shapeFunction[count] * factor;
                            }

                            if (theIP + 1 == GetNumIntegrationPoints())
                            {
                                //calculate sum of diagonal entries (is identical for all directions, that's why only x direction is calculated
                                double sum_diagonal(0);
                                for (int count = 0; count < shapeFunction.rows(); count++)
                                {
                                    sum_diagonal += result(2 * count);
                                }

                                //scale so that the sum of the diagonals represents the full mass
                                double scaleFactor = total_mass / sum_diagonal;
                                for (int count = 0; count < shapeFunction.rows(); count++)
                                {
                                    result(2 * count) *= scaleFactor;
                                    result(2 * count + 1) = result(2 * count);
                                }
                            }
                        }
                            break;
                        default:
                            throw MechanicsException("[NuTo::Element2D::Evaluate] Element output LUMPED_HESSIAN_2_TIME_DERIVATIVE for " + Node::AttributeToString(dof) + " not implemented.");
                        }
                    }

                    break;

                case Element::RESIDUAL_NORM_FACTOR:
                {
                    double factor = 1.0;
                    for (auto dof : activeDofs)
                    {
                        int startIndex = mInterpolationType->Get(dof).GetLocalStartIndex();
                        switch (dof)
                        {
                        case Node::DISPLACEMENTS:
                        {
                            AddDetJNtX(shapeFunctions.at(dof),
                                       residualNormFactorDisplacements,
                                       factor,
                                       startIndex,
                                       it->second->GetFullVectorDouble());
                        }
                            break;
                        case Node::RELATIVEHUMIDITY:
                        {
                            if(activeDofs.find(Node::WATERVOLUMEFRACTION) != activeDofs.end())
                            {
                                AddDetJNtX(shapeFunctions.at(dof),
                                           residualNormFactorRelativeHumidity,
                                           factor,
                                           startIndex,
                                           it->second->GetFullVectorDouble());
                            }
                        }
                            break;
                        case Node::WATERVOLUMEFRACTION:
                        {
                            if(activeDofs.find(Node::RELATIVEHUMIDITY) != activeDofs.end())
                            {
                                AddDetJNtX(shapeFunctions.at(dof),
                                           residualNormFactorWaterVolumeFraction,
                                           factor,
                                           startIndex,
                                           it->second->GetFullVectorDouble());
                            }
                        }
                            break;
                        default:
                            throw MechanicsException("[NuTo::Element2D::Evaluate] Element output INTERNAL_GRADIENT for " + Node::AttributeToString(dof) + " not implemented.");

                        }
                    }
                    auto test = it->second->GetFullVectorDouble();
                    break;
                }
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
                    default:
                        throw MechanicsException("[NuTo::Element2D::GetIpData] Ip data not implemented.");
                    }
                    break;
                case Element::GLOBAL_ROW_DOF:
                case Element::GLOBAL_COLUMN_DOF:
                    break;
                default:
                    throw MechanicsException("[NuTo::Element2D::Evaluate] element output not implemented.");
                }
            }

        }

    } catch (NuTo::MechanicsException& e)
    {
        std::stringstream ss;
        ss << mStructure->ElementGetId(this);
        e.AddMessage("[NuTo::Element2D::Evaluate] Error evaluating element data of element" + ss.str() + ".");
        throw e;
    }
    return Error::SUCCESSFUL;
}

NuTo::Element::eElementType NuTo::Element2D::GetEnumType() const
{
    return NuTo::Element::ELEMENT2D;
}


int NuTo::Element2D::GetLocalDimension()const
{
    return 2;
}

NuTo::NodeBase* NuTo::Element2D::GetNode(int rLocalNodeNumber)
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->GetNumNodes());
    return mNodes[rLocalNodeNumber];
}

const NuTo::NodeBase* NuTo::Element2D::GetNode(int rLocalNodeNumber) const
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->GetNumNodes());
    return mNodes[rLocalNodeNumber];
}

NuTo::NodeBase* NuTo::Element2D::GetNode(int rLocalNodeNumber, Node::eAttributes rDofType)
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber <  (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->Get(rDofType).GetNumNodes());
    int nodeIndex = mInterpolationType->Get(rDofType).GetNodeIndex(rLocalNodeNumber);
    assert(nodeIndex < (int )mNodes.size());
    return mNodes[nodeIndex];
}

const NuTo::NodeBase* NuTo::Element2D::GetNode(int rLocalNodeNumber, Node::eAttributes rDofType) const
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->Get(rDofType).GetNumNodes());
    int nodeIndex = mInterpolationType->Get(rDofType).GetNodeIndex(rLocalNodeNumber);
    assert(nodeIndex < (int )mNodes.size());
    return mNodes[nodeIndex];
}

void NuTo::Element2D::SetNode(int rLocalNodeNumber, NodeBase* rNode)
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->GetNumNodes());
    assert(rNode != nullptr);
    mNodes[rLocalNodeNumber] = rNode;
}

void NuTo::Element2D::ResizeNodes(int rNewNumNodes)
{
    if (rNewNumNodes == (int) mNodes.size())
        return;

    if (rNewNumNodes > (int) mNodes.size())
    {
        // just resize (enlarge)
        mNodes.resize(rNewNumNodes);
    } else
    {
        throw MechanicsException("[NuTo::Element2D::ResizeNodes] Resize that reduces the number of nodes is not implemented yet.");
    }
}

void NuTo::Element2D::ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
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

void NuTo::Element2D::SetSection(const SectionBase* rSection)
{
    mSection = rSection;
}

const NuTo::SectionBase* NuTo::Element2D::GetSection() const
{
    return mSection;
}

NuTo::ConstitutiveStaticDataBase* NuTo::Element2D::AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw) const
{
    return rConstitutiveLaw->AllocateStaticDataEngineeringStress_EngineeringStrain2D(this);
}

const Eigen::VectorXd NuTo::Element2D::GetIntegrationPointVolume() const
{

    Eigen::MatrixXd localNodeCoord = this->ExtractNodeValues(0, Node::COORDINATES);

    const InterpolationBase& interpolationType = mInterpolationType->Get(Node::COORDINATES);

    double detJac;
    Eigen::Matrix2d dummyJacobian;

    Eigen::VectorXd volume(GetNumIntegrationPoints());

    for (int theIP = 0; theIP < GetNumIntegrationPoints(); theIP++)
    {
        Eigen::MatrixXd derivativeShapeFunctionsNatural = interpolationType.GetDerivativeShapeFunctionsNatural(theIP);

        CalculateJacobian(derivativeShapeFunctionsNatural, localNodeCoord, dummyJacobian, detJac);

        //attention in 2D, this is just the area, but that is required for the nonlocal model
        volume[theIP] = detJac * mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP);
    }
    return volume;
}

const Eigen::MatrixXd NuTo::Element2D::ExtractNodeValues(int rTimeDerivative, Node::eAttributes rDofType) const
{
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nodalValues;
    ExtractNodeValues(nodalValues, rTimeDerivative, rDofType);
    return nodalValues;
}

void NuTo::Element2D::ExtractNodeValues(Eigen::MatrixXd& rNodeValues, int rTimeDerivative, Node::eAttributes rDofType) const
{
    const InterpolationBase& interpolationTypeDof = GetInterpolationType()->Get(rDofType);

    int numNodes = interpolationTypeDof.GetNumNodes();
    int numDofs = interpolationTypeDof.GetNumDofs();
    int numDofsPerNode = numDofs / numNodes;

    rNodeValues.resize(numDofsPerNode, numNodes);

    for (int iNode = 0; iNode < numNodes; ++iNode)
    {
        const NodeBase* node = mNodes[interpolationTypeDof.GetNodeIndex(iNode)];

        switch (rDofType)
        {
        case Node::COORDINATES:
            rNodeValues.block<2, 1>(0, iNode) = node->GetCoordinates2D();
            break;

        case Node::DISPLACEMENTS:
            rNodeValues.block<2, 1>(0, iNode) = node->GetDisplacements2D(rTimeDerivative);
            break;

        case Node::TEMPERATURES:
            rNodeValues(0, iNode) = node->GetTemperature(rTimeDerivative);
            break;

        case Node::NONLOCALEQPLASTICSTRAIN:
            rNodeValues.block<2, 1>(0, iNode) = node->GetNonlocalEqPlasticStrain(rTimeDerivative);
            break;

        case Node::NONLOCALEQSTRAIN:
            rNodeValues(0, iNode) = node->GetNonlocalEqStrain(rTimeDerivative);
            break;

        case Node::WATERVOLUMEFRACTION:
            rNodeValues(0,iNode) = node->GetWaterVolumeFraction(rTimeDerivative);
            break;

        case Node::RELATIVEHUMIDITY:
            rNodeValues(0,iNode) = node->GetRelativeHumidity(rTimeDerivative);
            break;

        default:
            throw MechanicsException("[NuTo::Element2D::ExtractNodeValues] Not implemented for " + Node::AttributeToString(rDofType) + ".");
        }
    }
}

const Eigen::VectorXi NuTo::Element2D::CalculateGlobalRowDofs() const
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
                globalRowDofs[index++] = nodePtr->GetDofDisplacement(1);
            }
                break;
            case Node::TEMPERATURES:
            {
                globalRowDofs[index++] = nodePtr->GetDofTemperature();
            }
                break;
            case Node::NONLOCALEQSTRAIN:
            {
                globalRowDofs[index++] = nodePtr->GetDofNonlocalEqStrain();
            }
                break;
            case Node::WATERVOLUMEFRACTION:
            {
                globalRowDofs[index++] = nodePtr->GetDofWaterVolumeFraction();
            }
                break;
            case Node::RELATIVEHUMIDITY:
            {
                globalRowDofs[index++] = nodePtr->GetDofRelativeHumidity();
            }
                break;
            default:
                throw MechanicsException("[NuTo::Element2D::CalculateGlobalRowDofs] Not implemented for " + Node::AttributeToString(dof) + ".");

            }
        }
    }
    return globalRowDofs;
}

const Eigen::VectorXi NuTo::Element2D::CalculateGlobalColumnDofs() const
{
    int numNonlocalElements = GetNumNonlocalElements();
    if (numNonlocalElements == 0)
        return CalculateGlobalRowDofs();
    else
    {

        int totalDofs = 0;
        for (const ElementBase* nonlocalElement : GetNonlocalElements())
            totalDofs += nonlocalElement->GetInterpolationType()->Get(Node::DISPLACEMENTS).GetNumDofs();

        Eigen::VectorXi globalColumnDofs(totalDofs);

        int shift = 0;
        for (const ElementBase* nonlocalElement : GetNonlocalElements())
        {
            for (int nodeCount = 0; nodeCount < nonlocalElement->GetNumNodes(Node::DISPLACEMENTS); nodeCount++)
            {
                const NodeBase *nodePtr = nonlocalElement->GetNode(nodeCount, Node::DISPLACEMENTS);
                if (nodePtr->GetNumDisplacements() > 0 && totalDofs > 0)
                {
                    globalColumnDofs[shift + 2 * nodeCount] = nodePtr->GetDofDisplacement(0);
                    globalColumnDofs[shift + 2 * nodeCount + 1] = nodePtr->GetDofDisplacement(1);
                }
            }
            shift += nonlocalElement->GetInterpolationType()->Get(Node::DISPLACEMENTS).GetNumDofs();
        }
        return globalColumnDofs;
    }
}

const NuTo::DeformationGradient2D NuTo::Element2D::CalculateDeformationGradient(const Eigen::MatrixXd& rDerivativeShapeFunctionsLocal, const Eigen::MatrixXd& rNodalDisplacements) const
{
    DeformationGradient2D deformationGradient;

    assert(rDerivativeShapeFunctionsLocal.rows() == rNodalDisplacements.cols());
    assert(rDerivativeShapeFunctionsLocal.cols() == rNodalDisplacements.rows());

    Eigen::Matrix2d deformationGradientEigen = rNodalDisplacements.lazyProduct(rDerivativeShapeFunctionsLocal);
    deformationGradientEigen.noalias() += Eigen::Matrix2d::Identity();
    Eigen::Map<Eigen::Matrix2d>(deformationGradient.mDeformationGradient) = deformationGradientEigen;

    return deformationGradient;
}

//! @brief adds to a matrix the product B^tCBnonlocal, where B contains the derivatives of the shape functions and C is the constitutive tangent and Bnonlocal is the nonlocal B matrix
//! eventually include also area/width of an element (that's the mechanics solution)
//! @param rLocalDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param rNonlocalDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param ConstitutiveTangentBase constitutive tangent matrix
//! @param rFactor factor including determinant of Jacobian and IP weight
//! @param rCoefficientMatrix to be added to
//! &param rFirstCol first column of the coefficient matrix to be modified (corresponding to the current nonlocal element)
void NuTo::Element2D::AddDetJBtCB(const Eigen::MatrixXd& rLocalDerivativeShapeFunctionsLocal, const Eigen::MatrixXd& rNonlocalDerivativeShapeFunctionsLocal,
        const ConstitutiveTangentLocal<3, 3>& rConstitutiveTangent, double rFactor, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCoefficientMatrix, int rFirstRow, int rFirstCol) const
{
//    assert(rCoefficientMatrix.GetNumRows()==2*GetNumNodesField() && rFirstCol + (int)rNonlocalDerivativeShapeFunctionsLocal.size()<=rCoefficientMatrix.GetNumColumns());
    int numNodes = mInterpolationType->Get(Node::DISPLACEMENTS).GetNumNodes();
    assert(rLocalDerivativeShapeFunctionsLocal.rows() == numNodes);
    assert(rLocalDerivativeShapeFunctionsLocal.cols() == 2);
    const double *C = rConstitutiveTangent.data();
    double x1, x2, y1, y2, x1x2, y2x1, x2y1, y2y1;
    for (int theNode1 = 0; theNode1 < numNodes; theNode1++)
    {
        int node1mul2 = 2 * theNode1;
        int node1mul2plus1 = node1mul2 + 1;

        x1 = rFactor * rLocalDerivativeShapeFunctionsLocal(theNode1, 0);
        y1 = rFactor * rLocalDerivativeShapeFunctionsLocal(theNode1, 1);
        // >> is division by two, but faster
        for (unsigned int theNode2 = 0; theNode2 < rNonlocalDerivativeShapeFunctionsLocal.size() >> 1; theNode2++)
        {
            int node2mul2 = 2 * theNode2;
            int node2mul2plus1 = node2mul2 + 1;

            x2 = rNonlocalDerivativeShapeFunctionsLocal(theNode2, 0);
            y2 = rNonlocalDerivativeShapeFunctionsLocal(theNode2, 1);

            x1x2 = x2 * x1;
            y2x1 = y2 * x1;
            x2y1 = x2 * y1;
            y2y1 = y2 * y1;

            rCoefficientMatrix(rFirstRow + node1mul2, rFirstCol + node2mul2) += x1x2 * C[0] + x2y1 * C[2] + y2x1 * C[6] + y2y1 * C[8];
            rCoefficientMatrix(rFirstRow + node1mul2, rFirstCol + node2mul2plus1) += x1x2 * C[6] + x2y1 * C[8] + y2x1 * C[3] + y2y1 * C[5];
            rCoefficientMatrix(rFirstRow + node1mul2plus1, rFirstCol + node2mul2) += x1x2 * C[2] + x2y1 * C[1] + y2x1 * C[8] + y2y1 * C[7];
            rCoefficientMatrix(rFirstRow + node1mul2plus1, rFirstCol + node2mul2plus1) += x1x2 * C[8] + x2y1 * C[7] + y2x1 * C[5] + y2y1 * C[4];
        }
    }
}

//! @brief adds to a matrix the product B^tCB, where B contains the derivatives of the shape functions and C is the constitutive tangent
//! eventually include also area/width of an element (that's the thermal solution)
//! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param ConstitutiveTangentBase constitutive tangent matrix
//! @param rFactor factor including determinant of Jacobian and IP weight
//! @param rRow row, where to start to add the submatrix
//! @param rCoefficientMatrix to be added to
void NuTo::Element2D::AddDetJBtCB(const Eigen::MatrixXd& rDerivativeShapeFunctionsGlobal, const ConstitutiveTangentLocal<2, 2>& rConstitutiveTangent, double rFactor, int rRow, int rCol,
        FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCoefficientMatrix) const
{
    const double *C = rConstitutiveTangent.data();
    int numNodes = mInterpolationType->Get(Node::DISPLACEMENTS).GetNumNodes();
    assert(rDerivativeShapeFunctionsGlobal.rows() == numNodes);
    assert(rDerivativeShapeFunctionsGlobal.cols() == 2);

    double x1, x2, y1, y2;
    for (int theNode1 = 0; theNode1 < numNodes; theNode1++)
    {
        int node1mul2 = theNode1 + theNode1;
        int node1mul2plus1 = node1mul2 + 1;

        x1 = rFactor * rDerivativeShapeFunctionsGlobal(theNode1, 0);
        y1 = rFactor * rDerivativeShapeFunctionsGlobal(theNode1, 1);
        node1mul2 += rRow;
        node1mul2plus1 += rRow;
        for (int theNode2 = 0; theNode2 < numNodes; theNode2++)
        {
            int node2mul2 = theNode2 + theNode2;
            int node2mul2plus1 = node2mul2 + 1;
            node2mul2 += rCol;
            node2mul2plus1 += rCol;

            x2 = rDerivativeShapeFunctionsGlobal(theNode2, 0);
            y2 = rDerivativeShapeFunctionsGlobal(theNode2, 1);

            assert(rCoefficientMatrix.GetNumRows() > node1mul2plus1 && rCoefficientMatrix.GetNumColumns() > node1mul2plus1);
            assert(rCoefficientMatrix.GetNumRows() > node2mul2plus1 && rCoefficientMatrix.GetNumColumns() > node2mul2plus1);

            rCoefficientMatrix(node1mul2, node2mul2) += x1 * C[0] * x2;
            rCoefficientMatrix(node1mul2, node2mul2plus1) += x1 * C[2] * y2;
            rCoefficientMatrix(node1mul2plus1, node2mul2) += y1 * C[1] * x2;
            rCoefficientMatrix(node1mul2plus1, node2mul2plus1) += y1 * C[3] * y2;
        }
    }
}

//! @brief adds up the internal force vector
//! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param rEngineeringStress stress
//! @param factor factor including det Jacobian area and integration point weight
//! @param rRow start row (in case of a multifield problem)
//! @param rResult resforce vector
void NuTo::Element2D::AddDetJBtSigma(const Eigen::MatrixXd& rDerivativeShapeFunctionsLocal, const EngineeringStress2D& rEngineeringStress, double rFactor, int rRow,
        FullVector<double, Eigen::Dynamic>& rResult) const
{
//    assert(rResult.GetNumRows()==2*GetNumNodesField() && rResult.GetNumColumns()==1);
    assert(rResult.GetNumColumns() == 1);
    const double *s = rEngineeringStress.GetData();
    int numNodes = mInterpolationType->Get(Node::DISPLACEMENTS).GetNumNodes();
    assert(rDerivativeShapeFunctionsLocal.rows() == numNodes);
    assert(rDerivativeShapeFunctionsLocal.cols() == 2);

    double x1, y1;
    for (int theNode1 = 0; theNode1 < numNodes; theNode1++)
    {
        int node1mul2 = 2 * theNode1;
        int node1mul2plus1 = node1mul2 + 1;

        x1 = rFactor * rDerivativeShapeFunctionsLocal(theNode1, 0);
        y1 = rFactor * rDerivativeShapeFunctionsLocal(theNode1, 1);

        rResult(rRow + node1mul2) += x1 * s[0] + y1 * s[2];
        rResult(rRow + node1mul2plus1) += y1 * s[1] + x1 * s[2];
    }
}

//! @brief calculates the Kee matrix
//! @param rShapeFunctions of the ip for all shape functions
//! @param rDerivativeShapeFunctions of the ip for all shape functions
//! @param nonlocal gradient radius xi
//! @param rFactor multiplication factor (detJ area..)
//! @param Kee return matrix with detJ * (Nt N + cBtB)
void NuTo::Element2D::CalculateKee(Eigen::VectorXd rShapeFunctions, const Eigen::MatrixXd& rDerivativeShapeFunctions, ConstitutiveTangentLocal<1, 1>& rNonlocalParameter, double rFactor,
        FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rKee) const
{
    //resize and set to zero
    rKee.Resize(rShapeFunctions.rows(), rShapeFunctions.rows());

    rKee = rFactor * (1/rNonlocalParameter(0, 0)*rShapeFunctions * rShapeFunctions.transpose() + rDerivativeShapeFunctions * rDerivativeShapeFunctions.transpose());
}
//! @brief add Kee*nonlocalEqStrain-detJ*N.T*localEqStrain (detJ is already included in Kee)
//! @param rShapeFunctions of the ip for all shape functions
//! @param rLocalEqStrain local eq. strain values
//! @param rKee stiffness matrix Kee
//! @param rNodeNonlocalEqStrain nodal nonlocal eq strain values
//! @param rFactor factor including detJ and area
//! @param rResult result
void NuTo::Element2D::AddDetJRnonlocalEqStrain(const Eigen::VectorXd& rShapeFunctions, LocalEqStrain& rLocalEqStrain, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rKee,
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& rNodeNonlocalEqStrain, double rFactor, int startIndexNonlocalEqStrain, FullVector<double, Eigen::Dynamic>& rResult) const
{
    assert(rResult.GetNumRows() >= (int )(startIndexNonlocalEqStrain + rShapeFunctions.size()));
    assert(rShapeFunctions.size() == rNodeNonlocalEqStrain.size());

    // perform Kee * nodeNonlocalEqStrain
    rResult.segment(startIndexNonlocalEqStrain, rShapeFunctions.size()) += rKee * rNodeNonlocalEqStrain.transpose() - rLocalEqStrain.GetValue(0) * rFactor * rShapeFunctions;

}

//! @brief add detJ B.T dSigma/dnonlocalEqStrain N
//! @param derivativeShapeFunctions of the ip for all shape functions
//! @param tangentStressNonlocalEqStrain derivative of the stress with respect to the nonlocal eq strain
//! @param rShapeFunctions of the ip for all shape functions
//! @param rFactor factor including detJ and area
//! @param rResult result
void NuTo::Element2D::AddDetJBtdSigmadNonlocalEqStrainN(const Eigen::MatrixXd& rDerivativeShapeFunctions, ConstitutiveTangentLocal<3, 1>& rTangentStressNonlocalEqStrain,
        Eigen::VectorXd rShapeFunctions, double rFactor, int rRow, int rCol, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rResult) const
{
    Eigen::VectorXd tmpfactor = rTangentStressNonlocalEqStrain * rFactor;

    int numNodesDisplacement = mInterpolationType->Get(Node::DISPLACEMENTS).GetNumNodes();
    int numNodesNonlocalEqStrain = mInterpolationType->Get(Node::NONLOCALEQSTRAIN).GetNumNodes();

    assert(rDerivativeShapeFunctions.rows() == numNodesDisplacement);
    assert(rDerivativeShapeFunctions.cols() == GetStructure()->GetDimension());
    assert(rShapeFunctions.rows() == numNodesNonlocalEqStrain);
    assert(rShapeFunctions.cols() == 1);

    for (int theNode = 0; theNode < numNodesDisplacement; theNode++)
    {
        double dphi_dx = rDerivativeShapeFunctions(theNode, 0);
        double dphi_dy = rDerivativeShapeFunctions(theNode, 1);

        rResult.block(rRow + 2 * theNode, rCol, 1, numNodesNonlocalEqStrain) += (dphi_dx * tmpfactor(0) + dphi_dy * tmpfactor(2)) * rShapeFunctions.transpose();
        rResult.block(rRow + 2 * theNode + 1, rCol, 1, numNodesNonlocalEqStrain) += (dphi_dx * tmpfactor(2) + dphi_dy * tmpfactor(1)) * rShapeFunctions.transpose();
    }
}

//! @brief add detJ N_transpose dEqStrain/dEpsilon B
//! @param rShapeFunctions of the ip for the nonlocal eq strain dofs
//! @param rTangentLocalEqStrainStrain derivative of the local eq strains with respect to the strain
//! @param rderivativeShapeFunctions of the ip for the displacement dofs
//! @param rFactor factor including detJ and area
//! @param rResult result
void NuTo::Element2D::AddDetJNtdLocalEqStraindEpsilonB(Eigen::VectorXd rShapeFunctions, ConstitutiveTangentLocal<3, 1>& rTangentLocalEqStrainStrain, const Eigen::MatrixXd& rDerivativeShapeFunctions,
        double rFactor, int rRow, int rCol, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rResult) const
{
    Eigen::VectorXd tmpfactor = rTangentLocalEqStrainStrain * rFactor;

    int numNodesDisplacement = mInterpolationType->Get(Node::DISPLACEMENTS).GetNumNodes();
    int numNodesNonlocalEqStrain = mInterpolationType->Get(Node::NONLOCALEQSTRAIN).GetNumNodes();

    assert(rDerivativeShapeFunctions.rows() == numNodesDisplacement);
    assert(rDerivativeShapeFunctions.cols() == GetStructure()->GetDimension());
    assert(rShapeFunctions.rows() == numNodesNonlocalEqStrain);
    assert(rShapeFunctions.cols() == 1);

    for (int theNode = 0; theNode < numNodesDisplacement; theNode++)
    {
        double dphi_dx = rDerivativeShapeFunctions(theNode, 0);
        double dphi_dy = rDerivativeShapeFunctions(theNode, 1);

        rResult.block(rRow, rCol + 2 * theNode, numNodesNonlocalEqStrain, 1) -= (dphi_dx * tmpfactor(0) + dphi_dy * tmpfactor(2)) * rShapeFunctions;
        rResult.block(rRow, rCol + 2 * theNode + 1, numNodesNonlocalEqStrain, 1) -= (dphi_dx * tmpfactor(2) + dphi_dy * tmpfactor(1)) * rShapeFunctions;
    }

}




//! @brief adds up the constitutive Tangent times the derivative shape functions
//! @param rDerivativeShapeFunctions the derivative shape functions
//! @param rConstitutiveTangent the result given by the constitutive law
//! @param factor factor including det Jacobian area and integration point weight
//! @param rRow start row in case of a multifield problem
//! @param rResult result vector
void NuTo::Element2D::AddDetJBtX(const Eigen::MatrixXd &rDerivativeShapeFunctions, ConstitutiveTangentLocal<2,1> &rConstitutiveTangent, double rFactor, int rRow, FullVector<double, Eigen::Dynamic> &rResult) const
{
    Eigen::VectorXd tmpfactor = rConstitutiveTangent * rFactor;
    int NumDofs      = rDerivativeShapeFunctions.rows();

    assert(rDerivativeShapeFunctions.cols() == GetStructure()->GetDimension());

    Eigen::VectorXd result = rDerivativeShapeFunctions * tmpfactor;
    rResult.block(rRow, 0, NumDofs, 1) += result;
}








//! @brief adds up the constitutive Tangent times the Shape Functions
//! @param rShapeFunctions the shape functions
//! @param rConstitutiveTangent the result given by the constitutive law
//! @param factor factor including det Jacobian area and integration point weight
//! @param rRow start row in case of a multifield problem
//! @param rResult result vector
void NuTo::Element2D::AddDetJNtX(Eigen::VectorXd &rShapeFunctions, ConstitutiveTangentLocal<1, 1> &rConstitutiveTangent, double rFactor, int rRow, FullVector<double, Eigen::Dynamic> &rResult) const
{
    Eigen::VectorXd tmpfactor = rConstitutiveTangent * rFactor;
    int NumDofs      = rShapeFunctions.rows();

    assert(rShapeFunctions.cols() == 1);

    Eigen::VectorXd result =  rShapeFunctions * tmpfactor;

    rResult.block(rRow, 0, NumDofs, 1) += result;
}


//! @brief adds to a matrix the product B1^t X B2, where B1 and B2 are the derivative shape functions and X is the constitutive tangent
//! @param rDerivativeShapeFunction1 derivative shape function 1
//! @param rDerivativeShapeFunction2 derivative shape function 2
//! @param ConstitutiveTangentBase constitutive tangent matrix
//! @param rFactor factor including area, determinant of Jacobian and IP weight
//! @param rRow row, where to start to add the submatrix
//! @param rCol col, where to start to add the submatrix
//! @param rResult result
void NuTo::Element2D::AddDetJBtXB(const Eigen::MatrixXd &rDerivativeShapeFunctions1,
                                  const Eigen::MatrixXd &rDerivativeShapeFunctions2,
                                  const ConstitutiveTangentLocal<1, 1> &rConstitutiveTangent,
                                  double rFactor,
                                  int rRow,
                                  int rCol,
                                  FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> &rResult) const
{
    int NumRowDofs  = rDerivativeShapeFunctions1.rows();
    int NumColDofs  = rDerivativeShapeFunctions2.rows();

    assert(rDerivativeShapeFunctions1.cols() == GetStructure()->GetDimension());
    assert(rDerivativeShapeFunctions2.cols() == GetStructure()->GetDimension());

    rFactor *= rConstitutiveTangent(0,0);

    Eigen::MatrixXd result =  rDerivativeShapeFunctions1 * rFactor * rDerivativeShapeFunctions2.transpose();

    rResult.block(rRow,rCol,NumRowDofs,NumColDofs)+=result;
}


//! @brief adds to a matrix the product B^t X N, where B is the derivative shape functions, N the shape function and X the constitutive tangent
//! @param rDerivativeShapeFunction derivative shape function
//! @param rShapeFunction shape function
//! @param ConstitutiveTangentBase constitutive tangent matrix
//! @param rFactor factor including area, determinant of Jacobian and IP weight
//! @param rRow row, where to start to add the submatrix
//! @param rCol col, where to start to add the submatrix
//! @param rResult result
void NuTo::Element2D::AddDetJBtXN(const Eigen::MatrixXd &rDerivativeShapeFunction,
                                  const Eigen::VectorXd &rShapeFunction,
                                  const ConstitutiveTangentLocal<2, 1> &rConstitutiveTangent,
                                  double rFactor,
                                  int rRow,
                                  int rCol,
                                  FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> &rResult) const
{
    int NumRowDofs  = rDerivativeShapeFunction.rows();
    int NumColDofs  = rShapeFunction.rows();

    Eigen::VectorXd tmpfactor = rConstitutiveTangent * rFactor;

    assert(rShapeFunction.cols() == 1);
    assert(rDerivativeShapeFunction.cols() == GetStructure()->GetDimension() );


    Eigen::MatrixXd result =  rDerivativeShapeFunction * tmpfactor * rShapeFunction.transpose();

    rResult.block(rRow,rCol,NumRowDofs,NumColDofs)+=result;
}

//! @brief adds to a matrix the product N1^t X N2, where N1 and N2 contains the the shape functions and X is the constitutive tangent
//! @param rShapeFunction1 shape function 1
//! @param rShapeFunction2 shape function 2
//! @param ConstitutiveTangentBase constitutive tangent matrix
//! @param rFactor factor including area, determinant of Jacobian and IP weight
//! @param rRow row, where to start to add the submatrix
//! @param rCol col, where to start to add the submatrix
void NuTo::Element2D::AddDetJNtXN(const Eigen::VectorXd &rShapeFunction1,
                                  const Eigen::VectorXd &rShapeFunction2,
                                  const ConstitutiveTangentLocal<1,1> &rConstitutiveTangent,
                                  double rFactor,
                                  int rRow,
                                  int rCol,
                                  FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> &rResult) const
{
    int NumRowDofs  = rShapeFunction1.rows();
    int NumColDofs  = rShapeFunction2.rows();

    assert(rShapeFunction1.cols() == 1);
    assert(rShapeFunction2.cols() == 1);

    rFactor *= rConstitutiveTangent(0,0);

    Eigen::MatrixXd result =  rShapeFunction1 * rFactor * rShapeFunction2.transpose();

    rResult.block(rRow,rCol,NumRowDofs,NumColDofs)+=result;
}


void NuTo::Element2D::CalculateJacobian(const Eigen::MatrixXd& rDerivativeShapeFunctions, const Eigen::MatrixXd& rNodeCoordinates, Eigen::Matrix2d& rInvJacobian, double& rDetJac) const
{
    /*       jacobian
     j0, j2,
     j1, j3 */

    assert(rDerivativeShapeFunctions.rows() == rNodeCoordinates.cols());
    assert(rDerivativeShapeFunctions.cols() == 2);
    assert(rNodeCoordinates.rows() == 2);

    Eigen::Matrix2d jacobian = rNodeCoordinates.lazyProduct(rDerivativeShapeFunctions);
    rDetJac = jacobian.determinant();

    if (rDetJac == 0)
        throw MechanicsException("[NuTo::Element2D::CalculateJacobian] Determinant of the Jacobian is zero, no inversion possible.");

    rInvJacobian = jacobian.inverse();
}

double NuTo::Element2D::CalculateArea() const
{
//    Eigen::MatrixXd nodeCoordinates = ExtractNodeValues(0, Node::COORDINATES);
//    int numIntegrationPoints = GetNumIntegrationPoints();
//    Eigen::Matrix2d invJacobian;
//    double detJacobian;
//    double area = 0;
//    for (int iIP = 0; iIP < numIntegrationPoints; ++iIP)
//    {
//        const Eigen::MatrixXd& derivativeShapeFunctions = mInterpolationType->Get(Node::COORDINATES).GetDerivativeShapeFunctionsNatural(iIP);
//        CalculateJacobian(derivativeShapeFunctions, nodeCoordinates, invJacobian, detJacobian);
//        if (detJacobian <= 0)
//            throw MechanicsException("[NuTo::Element2D::CheckElement] Determinant of the Jacobian <= zero, no inversion possible.");
//        area += this->GetIntegrationPointWeight(iIP) * detJacobian;
//    }
//
//    std::cout << area << std::endl;

    double area = 0;

    const InterpolationBase& interpolationTypeCoords = mInterpolationType->Get(Node::COORDINATES);

    int numNodesCoordinates = interpolationTypeCoords.GetNumNodes();

    Eigen::Vector2d coordinates1, coordinates2;
    coordinates2 = mNodes[interpolationTypeCoords.GetNodeIndex(numNodesCoordinates - 1)]->GetCoordinates2D();

    for (int theNode = 0; theNode < numNodesCoordinates; theNode++)
    {
        coordinates1 = coordinates2;
        coordinates2 = mNodes[interpolationTypeCoords.GetNodeIndex(theNode)]->GetCoordinates2D();
        area += coordinates1[0] * coordinates2[1] - coordinates1[1] * coordinates2[0];
    }

//    std::cout << .5 * area << std::endl;

    return 0.5 * area;
}

void NuTo::Element2D::CheckElement()
{

    int numIntegrationPoints = GetNumIntegrationPoints();
    // check number of integration points
    if (numIntegrationPoints < 1)
    {
        throw MechanicsException("[NuTo::Element2D::CheckElement] invalid integration type.");
    }

    int theIP = 0;
    const Eigen::MatrixXd& derivativeShapeFunctions = mInterpolationType->Get(Node::COORDINATES).GetDerivativeShapeFunctionsNatural(theIP);
    Eigen::MatrixXd nodeCoordinates = ExtractNodeValues(0, Node::COORDINATES);
    Eigen::Matrix2d invJacobian;
    double detJacobian;

    CalculateJacobian(derivativeShapeFunctions, nodeCoordinates, invJacobian, detJacobian);
    if (detJacobian < 0)
    {
        this->ReorderNodes();
        // recalculate node coordinates after reordering
        nodeCoordinates = ExtractNodeValues(0, Node::COORDINATES);
    }

    double volume = 0;
    for (int iIP = 0; iIP < numIntegrationPoints; ++iIP)
    {
        const Eigen::MatrixXd& derivativeShapeFunctions = mInterpolationType->Get(Node::COORDINATES).GetDerivativeShapeFunctionsNatural(iIP);
        CalculateJacobian(derivativeShapeFunctions, nodeCoordinates, invJacobian, detJacobian);
        if (detJacobian <= 0)
            throw MechanicsException("[NuTo::Element2D::CheckElement] Determinant of the Jacobian <= zero, no inversion possible.");
        volume += this->GetIntegrationPointWeight(iIP) * detJacobian;
    }

    // check element volume
    if (volume < 1e-14)
    {
        throw MechanicsException("[NuTo::Element2D::CheckElement] element with zero volume (check nodes).");
    }
}




