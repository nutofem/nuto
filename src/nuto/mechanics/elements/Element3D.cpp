/*
 * Element3D.h
 *
 *  Created on: 19 May 2015
 *      Author: ttitsche
 */

#include "nuto/mechanics/elements/Element3D.h"
#include "nuto/mechanics/constitutive/mechanics/Damage.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/thermal/HeatFlux3D.h"
#include "nuto/mechanics/constitutive/thermal/Temperature.h"
#include "nuto/mechanics/constitutive/thermal/TemperatureGradient3D.h"
#include "nuto/mechanics/constitutive/mechanics/LocalEqStrain.h"
#include "nuto/mechanics/constitutive/mechanics/NonlocalEqStrain.h"
#include "nuto/mechanics/constitutive/mechanics/LocalEqTotalInelasticStrain.h"
#include "nuto/mechanics/constitutive/moistureTransport/RelativeHumidity.h"
#include "nuto/mechanics/constitutive/moistureTransport/RelativeHumidityGradient3D.h"
#include "nuto/mechanics/constitutive/moistureTransport/WaterVolumeFraction.h"
#include "nuto/mechanics/constitutive/moistureTransport/WaterVolumeFractionGradient3D.h"
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

NuTo::Element3D::Element3D(const NuTo::StructureBase* rStructure, const std::vector<NuTo::NodeBase*>& rNodes, ElementData::eElementDataType rElementDataType, IpData::eIpDataType rIpDataType,
        InterpolationType* rInterpolationType) :
        NuTo::ElementBase::ElementBase(rStructure, rElementDataType, rIpDataType, rInterpolationType),
        mNodes(rNodes),
        mSection(0)
{
}

NuTo::Error::eError NuTo::Element3D::Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput)
{
    if (mStructure->GetHessianConstant(1) == false)
        throw MechanicsException("[NuTo::Element3D::Evaluate] only implemented for a constant Hessian for the first derivative (damping).");
    if (mStructure->GetHessianConstant(2) == false)
        throw MechanicsException("[NuTo::Element3D::Evaluate] only implemented for a constant Hessian for the second derivative (mass).");

    try
    {
        const SectionBase* section(GetSection());
        if (section == 0)
            throw MechanicsException("[NuTo::Element3D::Evaluate] no section allocated for element.");

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

        //define inputs and outputs
        std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*> constitutiveInputList;
        std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*> constitutiveOutputList;

        EngineeringStress3D engineeringStress3D;
        EngineeringStrain3D engineeringStrain3D;

        EngineeringStrain3D engineeringPlasticStrain3D;

        NuTo::ConstitutiveTangentLocal<6, 6> tangentStressStrain;

        DeformationGradient3D deformationGradient;

        //allocate damage
        Damage damage;


        // Moisture Transport
        // ------------------

        RelativeHumidity                relativeHumidity;
        RelativeHumidity                relativeHumidityD1;
        RelativeHumidityGradient3D      relativeHumidityGradient;
        WaterVolumeFraction             waterVolumeFraction;
        WaterVolumeFraction             waterVolumeFractionD1;
        WaterVolumeFractionGradient3D   waterVolumeFractionGradient;


        //allocate nonlocal eq strain
        NonlocalEqStrain nonlocalEqStrain;

        //allocate local eq strain
        LocalEqStrain localEqStrain;

        //allocate local eq total inelastic strain
        LocalEqTotalInelasticStrain localEqTotlalInelasticStrain;

        //allocate transient nonlocal parameter
        ConstitutiveTangentLocal<1, 1> nonlocalParameter;
        ConstitutiveTangentLocal<6, 1> tangentStressNonlocalEqStrain;
        ConstitutiveTangentLocal<6, 1> tangentLocalEqStrainStrain;

        NuTo::ConstitutiveTangentNonlocal<6, 6> nonlocalTangentStressStrain;

        // Moisture Transport
        // ------------------

        // Naming scheme for matrices: tangent_D_X_D_Y_HN_AB
        // D_X_D_Y:     defines the derivative of X with respect to Y
        // HN:          defines the matrix the tangent belongs to, for example H0 for stiffness or H1 for damping
        // AB:          defines the combination of (derivative) shapefunctions the tangent has to be multiplied with, for example BN means: derivative shapefunctions * tangent * shapefunctions

        // Internal Gradient
        ConstitutiveTangentLocal<1,1> residualWaterPhaseN;
        ConstitutiveTangentLocal<3,1> residualWaterPhaseB;
        ConstitutiveTangentLocal<1,1> residualVaporPhaseN;
        ConstitutiveTangentLocal<3,1> residualVaporPhaseB;
        // Hessian 0
        ConstitutiveTangentLocal<1,1> tangent_D_Residual_RH_D_RH_H0_BB;
        ConstitutiveTangentLocal<1,1> tangent_D_Residual_RH_D_RH_H0_NN;
        ConstitutiveTangentLocal<3,1> tangent_D_Residual_RH_D_WV_H0_BN;
        ConstitutiveTangentLocal<1,1> tangent_D_Residual_RH_D_WV_H0_NN;
        ConstitutiveTangentLocal<1,1> tangent_D_Residual_WV_D_RH_H0_NN;
        ConstitutiveTangentLocal<1,1> tangent_D_Residual_WV_D_WV_H0_BB;
        ConstitutiveTangentLocal<3,1> tangent_D_Residual_WV_D_WV_H0_BN;
        ConstitutiveTangentLocal<1,1> tangent_D_Residual_WV_D_WV_H0_NN;
        // Hessian 1
        ConstitutiveTangentLocal<1,1> tangent_D_Residual_RH_D_RH_H1_NN;
        ConstitutiveTangentLocal<1,1> tangent_D_Residual_RH_D_WV_H1_NN;
        ConstitutiveTangentLocal<1,1> tangent_D_Residual_WV_D_WV_H1_NN;

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
                constitutiveInputList[NuTo::Constitutive::Input::DEFORMATION_GRADIENT_3D] = &(deformationGradient);
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
                throw MechanicsException("[NuTo::Element3D::Evaluate] Constitutive input for " + Node::AttributeToString(dof) + " not implemented.");
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
                            constitutiveOutputList[NuTo::Constitutive::Output::ENGINEERING_STRESS_3D] = &(engineeringStress3D);
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
                            throw MechanicsException("[NuTo::Element3D::Evaluate] Constitutive output INTERNAL_GRADIENT for " + Node::AttributeToString(dof) + " not implemented.");

                        }
                    }
                }
                break;
			case Element::INTERNAL_GRADIENT_ELASTIC:
				it->second->GetFullVectorDouble().Resize(numActiveDofs);
				for (auto dof : activeDofs)
				{
					switch (dof)
					{
					case Node::DISPLACEMENTS:
					{
						constitutiveOutputList[NuTo::Constitutive::Output::ENGINEERING_STRESS_ELASTIC_3D] = &(engineeringStress3D);
                    }
					break;
					default:
						throw MechanicsException(
								"[NuTo::Element3D::Evaluate] Constitutive output INTERNAL_GRADIENT for " + Node::AttributeToString(dof)
						+ " not implemented.");

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
                        constitutiveOutputList[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_3D] = &tangentStressStrain;
                        if (activeDofs.find(Node::NONLOCALEQSTRAIN) != activeDofs.end())
                        {
                            constitutiveOutputList[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_3D] = &tangentStressNonlocalEqStrain;
                        }
                    }
                        break;
                    case Node::NONLOCALEQSTRAIN:
                    {
                        constitutiveOutputList[NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_XI_D_STRAIN_3D] = &tangentLocalEqStrainStrain;
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
                        throw MechanicsException("[NuTo::Element3D::Evaluate] Constitutive output HESSIAN_0_TIME_DERIVATIVE for " + Node::AttributeToString(dof) + " not implemented.");
                    }
                }
            }
                break;
			case Element::HESSIAN_0_TIME_DERIVATIVE_ELASTIC:
				{
					it->second->GetFullMatrixDouble().Resize(numActiveDofs,numActiveDofs);
	                it->second->GetFullMatrixDouble().setZero();
					it->second->SetSymmetry(true);
					it->second->SetConstant(true);
	                for (auto dof : activeDofs)
	                {
	                    switch (dof)
	                    {
	                    case Node::DISPLACEMENTS:
	                    {
	                        constitutiveOutputList[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_ELASTIC_3D] = &tangentStressStrain;
	                    }
	                        break;
	                    default:
	                        throw MechanicsException(
	                                "[NuTo::Element3D::Evaluate] Constitutive output HESSIAN_0_TIME_DERIVATIVE for " + Node::AttributeToString(dof)
	                                        + " not implemented.");

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
                        throw MechanicsException("[NuTo::Element3D::Evaluate] Constitutive output HESSIAN_1_TIME_DERIVATIVE for " + Node::AttributeToString(dof) + " not implemented.");
                    }
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
			case Element::FATIGUE_SAVE_STATIC_DATA:
				constitutiveOutputList[NuTo::Constitutive::Output::FATIGUE_SAVE_STATIC_DATA] = 0;
			break;
			case Element::FATIGUE_RESTORE_STATIC_DATA:
				constitutiveOutputList[NuTo::Constitutive::Output::FATIGUE_RESTORE_STATIC_DATA] = 0;
			break;
			case Element::FATIGUE_EXTRAPOLATE_STATIC_DATA:
				constitutiveOutputList[NuTo::Constitutive::Output::FATIGUE_EXTRAPOLATE_STATIC_DATA] = 0;
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
                case NuTo::IpData::TOTAL_INELASTIC_EQ_STRAIN:
                    it->second->GetFullMatrixDouble().Resize(1, GetNumIntegrationPoints());
                    //define outputs
                    constitutiveOutputList[NuTo::Constitutive::Output::LOCAL_EQ_TOTAL_INELASTIC_STRAIN] = &(localEqTotlalInelasticStrain);
                    break;
                default:
                    throw MechanicsException("[NuTo::Element3D::Evaluate] this ip data type is not implemented.");
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
                throw MechanicsException("[NuTo::Element3D::Evaluate] element output not implemented.");
            }
        }






        /*****************************************\
         *     CALCULATE CONSTITUTIVE INPUTS     *
        \*****************************************/

        Eigen::Matrix3d invJacobian;
        double detJacobian;

        std::map<Node::eAttributes, Eigen::VectorXd> shapeFunctions;
        std::map<Node::eAttributes, Eigen::MatrixXd> derivativeShapeFunctions;

        // loop over the integration points
        for (int theIP = 0; theIP < GetNumIntegrationPoints(); theIP++)
        {
            // calculate Jacobian
            const Eigen::MatrixXd& derivativeShapeFunctionsGeometryNatural = mInterpolationType->Get(Node::COORDINATES).GetDerivativeShapeFunctionsNatural(theIP);

            this->CalculateJacobian(derivativeShapeFunctionsGeometryNatural, nodalValues[Node::COORDINATES], invJacobian, detJacobian);

            double factor = detJacobian * (mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP));

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
                    throw MechanicsException("[NuTo::Element3D::Evaluate] Constitutive input for " + Node::AttributeToString(dof) + " not implemented.");
                }
            }



            /*****************************************\
             *      EVALUATE CONSTITUTIVE LAW        *
             \*****************************************/

            ConstitutiveBase* constitutivePtr = GetConstitutiveLaw(theIP);
            try
            {
                Error::eError error = constitutivePtr->Evaluate3D(this, theIP, constitutiveInputList, constitutiveOutputList);
                if (error != Error::SUCCESSFUL)
                    return error;
            } catch (NuTo::MechanicsException &e)
            {
                e.AddMessage("[NuTo::Element3D::Evaluate] error evaluating the constitutive model.");
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
                        double factor(fabs(detJacobian * (mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP))));
                        for (auto dof : activeDofs)
                        {
                            int startIndex = mInterpolationType->Get(dof).GetLocalStartIndex();
                            switch (dof)
                            {
                            case Node::DISPLACEMENTS:
                            {
                                AddDetJBtSigma(derivativeShapeFunctions.at(dof), engineeringStress3D, factor, startIndex, it->second->GetFullVectorDouble());
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
                                throw MechanicsException("[NuTo::Element3D::Evaluate] Element output INTERNAL_GRADIENT for " + Node::AttributeToString(dof) + " not implemented.");

                            }
                        }
                    }
                }
                    break;
				case Element::INTERNAL_GRADIENT_ELASTIC:
				{
					// Jacobian
                    double factor(fabs(detJacobian * (mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP))));
                    for (auto dof : activeDofs)
                    {
                        int startIndex = mInterpolationType->Get(dof).GetLocalStartIndex();
                        switch (dof)
                        {
                        case Node::DISPLACEMENTS:
                        {
                            AddDetJBtSigma(derivativeShapeFunctions.at(dof), engineeringStress3D, factor, startIndex, it->second->GetFullVectorDouble());
                        }
                            break;
                        default:
                            throw MechanicsException(
                                    "[NuTo::Element3D::Evaluate] Element output INTERNAL_GRADIENT_ELASTIC for " + Node::AttributeToString(dof)
                                            + " not implemented.");

                        }
                    }
				}
					break;
                case Element::HESSIAN_0_TIME_DERIVATIVE:
                {
                    //factor for the numerical integration
                    double factor(detJacobian * (mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP)));
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
                            AddDetJBtCB(derivativeShapeFunctions.at(dof), tangentStressStrain, factor, startIndex, startIndex, it->second->GetFullMatrixDouble());

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
                            throw MechanicsException("[NuTo::Element3D::Evaluate] Element output HESSIAN_0_TIME_DERIVATIVE for " + Node::AttributeToString(dof) + " not implemented.");

                        }
                    }
                }
                    break;
				case Element::HESSIAN_0_TIME_DERIVATIVE_ELASTIC:
                {
                    //factor for the numerical integration
                    double factor(detJacobian * (mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP)));
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
                            AddDetJBtCB(derivativeShapeFunctions.at(dof), tangentStressStrain, factor, startIndex, startIndex,
                                    it->second->GetFullMatrixDouble());
                        }
                            break;
                        default:
                            throw MechanicsException(
                                    "[NuTo::Element3D::Evaluate] Element output HESSIAN_0_TIME_DERIVATIVE_ELASTIC for " + Node::AttributeToString(dof)
                                            + " not implemented.");

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
                        throw MechanicsException("[NuTo::Element3D::Evaluate] Element output HESSIAN_1_TIME_DERIVATIVE for " + Node::AttributeToString(dof) + " not implemented.");

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

                            double factor(detJacobian * (mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP))
                                                      * constitutivePtr->GetParameterDouble(Constitutive::eConstitutiveParameter::DENSITY));
                            const Eigen::MatrixXd tmpMatrix = shapeFunctions.at(dof) * shapeFunctions.at(dof).transpose() * factor;
                            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& result(it->second->GetFullMatrixDouble());
                            for (int count = 0; count < tmpMatrix.rows(); count++)
                            {
                                for (int count2 = 0; count2 < tmpMatrix.cols(); count2++)
                                {
                                    result(3 * count, 3 * count2) += tmpMatrix(count, count2);
                                    result(3 * count + 1, 3 * count2 + 1) += tmpMatrix(count, count2);
                                    result(3 * count + 2, 3 * count2 + 2) += tmpMatrix(count, count2);
                                }
                            }
                        }
                            break;
                        case Node::RELATIVEHUMIDITY:
                        case Node::WATERVOLUMEFRACTION:
                            break;
                        default:
                            throw MechanicsException("[NuTo::Element3D::Evaluate] Element output HESSIAN_2_TIME_DERIVATIVE for " + Node::AttributeToString(dof) + " not implemented.");

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
                            double factor(detJacobian * (mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP))
                                                      * constitutivePtr->GetParameterDouble(Constitutive::eConstitutiveParameter::DENSITY));
                            FullVector<double, Eigen::Dynamic>& result(it->second->GetFullVectorDouble());
                            total_mass += factor;
                            const Eigen::VectorXd& shapeFunction = shapeFunctions.at(dof);
                            //calculate for the translational dofs the diagonal entries
                            for (int count = 0; count < shapeFunction.rows(); count++)
                            {
                                result(3 * count) += shapeFunction[count] * shapeFunction[count] * factor;
                            }

                            if (theIP + 1 == GetNumIntegrationPoints())
                            {
                                //calculate sum of diagonal entries (is identical for all directions, that's why only x direction is calculated
                                double sum_diagonal(0);
                                for (int count = 0; count < shapeFunction.rows(); count++)
                                {
                                    sum_diagonal += result(3 * count);
                                }

                                //scale so that the sum of the diagonals represents the full mass
                                double scaleFactor = total_mass / sum_diagonal;
                                for (int count = 0; count < shapeFunction.rows(); count++)
                                {
                                    result(3 * count) *= scaleFactor;
                                    result(3 * count + 1) = result(3 * count);
                                    result(3 * count + 2) = result(3 * count);
                                }
                            }
                        }
                            break;
                        default:
                            throw MechanicsException("[NuTo::Element3D::Evaluate] Element output LUMPED_HESSIAN_2_TIME_DERIVATIVE for " + Node::AttributeToString(dof) + " not implemented.");
                        }
                    }

                    break;
                case Element::UPDATE_STATIC_DATA:
                case Element::UPDATE_TMP_STATIC_DATA:
                    break;
				case Element::FATIGUE_SAVE_STATIC_DATA:
				break;
				case Element::FATIGUE_RESTORE_STATIC_DATA:
				break;
				case Element::FATIGUE_EXTRAPOLATE_STATIC_DATA:
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
                        //error = constitutivePtr->GetDamage(this, theIP, deformationGradient, rIpData.mEigenMatrix.data()[theIP]);
                        memcpy(&(it->second->GetFullMatrixDouble().data()[theIP]), localEqStrain.data(), sizeof(double));
                        break;
                    case NuTo::IpData::TOTAL_INELASTIC_EQ_STRAIN:
                        //error = constitutivePtr->GetDamage(this, theIP, deformationGradient, rIpData.mEigenMatrix.data()[theIP]);
                        memcpy(&(it->second->GetFullMatrixDouble().data()[theIP]), localEqTotlalInelasticStrain.GetData(), sizeof(double));
                        break;
                    default:
                        throw MechanicsException("[NuTo::Element3D::GetIpData] Ip data not implemented.");
                    }
                    break;
                case Element::GLOBAL_ROW_DOF:
                case Element::GLOBAL_COLUMN_DOF:
                    break;
                default:
                    throw MechanicsException("[NuTo::Element3D::Evaluate] element output not implemented.");
                }
            }

        }

    } catch (NuTo::MechanicsException& e)
    {
        std::stringstream ss;
        ss << mStructure->ElementGetId(this);
        e.AddMessage("[NuTo::Element3D::Evaluate] Error evaluating element data of element" + ss.str() + ".");
        throw e;
    }
    return Error::SUCCESSFUL;
}

NuTo::Element::eElementType NuTo::Element3D::GetEnumType() const
{
    return NuTo::Element::ELEMENT3D;
}

int NuTo::Element3D::GetLocalDimension()const
{
    return 3;
}

NuTo::NodeBase* NuTo::Element3D::GetNode(int rLocalNodeNumber)
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->GetNumNodes());
    return mNodes[rLocalNodeNumber];
}

const NuTo::NodeBase* NuTo::Element3D::GetNode(int rLocalNodeNumber) const
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->GetNumNodes());
    return mNodes[rLocalNodeNumber];
}

NuTo::NodeBase* NuTo::Element3D::GetNode(int rLocalNodeNumber, Node::eAttributes rDofType)
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->Get(rDofType).GetNumNodes());
    int nodeIndex = mInterpolationType->Get(rDofType).GetNodeIndex(rLocalNodeNumber);
    assert(nodeIndex < (int )mNodes.size());
    return mNodes[nodeIndex];
}

const NuTo::NodeBase* NuTo::Element3D::GetNode(int rLocalNodeNumber, Node::eAttributes rDofType) const
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->Get(rDofType).GetNumNodes());
    int nodeIndex = mInterpolationType->Get(rDofType).GetNodeIndex(rLocalNodeNumber);
    assert(nodeIndex < (int )mNodes.size());
    return mNodes[nodeIndex];
}

void NuTo::Element3D::SetNode(int rLocalNodeNumber, NodeBase* rNode)
{
    assert(rLocalNodeNumber >= 0);
    assert(rLocalNodeNumber < (int )mNodes.size());
    assert(rLocalNodeNumber < mInterpolationType->GetNumNodes());
    assert(rNode != nullptr);
    mNodes[rLocalNodeNumber] = rNode;
}

void NuTo::Element3D::ResizeNodes(int rNewNumNodes)
{
    if (rNewNumNodes == (int) mNodes.size())
        return;

    if (rNewNumNodes > (int) mNodes.size())
    {
        // just resize (enlarge)
        mNodes.resize(rNewNumNodes);
    } else
    {
        throw MechanicsException("[NuTo::Element3D::ResizeNodes] Resize that reduces the number of nodes is not implemented yet.");
    }
}

void NuTo::Element3D::ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
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

void NuTo::Element3D::SetSection(const SectionBase* rSection)
{
    mSection = rSection;
}

const NuTo::SectionBase* NuTo::Element3D::GetSection() const
{
    return mSection;
}

NuTo::ConstitutiveStaticDataBase* NuTo::Element3D::AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw) const
{
    return rConstitutiveLaw->AllocateStaticDataEngineeringStress_EngineeringStrain3D(this);
}

const Eigen::VectorXd NuTo::Element3D::GetIntegrationPointVolume() const
{

    Eigen::MatrixXd localNodeCoord = this->ExtractNodeValues(0, Node::COORDINATES);

    const InterpolationBase& interpolationType = mInterpolationType->Get(Node::COORDINATES);

    double detJac;
    Eigen::Matrix3d dummyJacobian;

    Eigen::VectorXd volume(GetNumIntegrationPoints());

    for (int theIP = 0; theIP < GetNumIntegrationPoints(); theIP++)
    {
        Eigen::MatrixXd derivativeShapeFunctionsNatural = interpolationType.GetDerivativeShapeFunctionsNatural(theIP);

        CalculateJacobian(derivativeShapeFunctionsNatural, localNodeCoord, dummyJacobian, detJac);

        volume[theIP] = detJac * mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP);
    }
    return volume;
}

const Eigen::MatrixXd NuTo::Element3D::ExtractNodeValues(int rTimeDerivative, Node::eAttributes rDofType) const
{
    const InterpolationBase& interpolationTypeDof = GetInterpolationType()->Get(rDofType);

    int numNodes = interpolationTypeDof.GetNumNodes();
    int numDofs = interpolationTypeDof.GetNumDofs();
    int numDofsPerNode = numDofs / numNodes;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nodalValues;
    nodalValues.resize(numDofsPerNode, numNodes);

    for (int iNode = 0; iNode < numNodes; ++iNode)
    {
        const NodeBase* node = mNodes[interpolationTypeDof.GetNodeIndex(iNode)];

        switch (rDofType)
        {
        case Node::COORDINATES:
            nodalValues.block<3, 1>(0, iNode) = node->GetCoordinates3D();
            break;

        case Node::DISPLACEMENTS:
            nodalValues.block<3, 1>(0, iNode) = node->GetDisplacements3D(rTimeDerivative);
            break;

        case Node::TEMPERATURES:
            nodalValues(0, iNode) = node->GetTemperature(rTimeDerivative);
            break;

        case Node::NONLOCALEQSTRAIN:
            nodalValues(0, iNode) = node->GetNonlocalEqStrain(rTimeDerivative);
            break;

        case Node::WATERVOLUMEFRACTION:
            nodalValues(0, iNode) = node->GetWaterVolumeFraction(rTimeDerivative);
            break;

        case Node::RELATIVEHUMIDITY:
            nodalValues(0, iNode) = node->GetRelativeHumidity(rTimeDerivative);
            break;

        default:
            throw MechanicsException("[NuTo::Element3D::ExtractNodeValues] Not implemented for " + Node::AttributeToString(rDofType) + ".");
        }
    }
    return nodalValues;
}

void NuTo::Element3D::ExtractNodeValues(Eigen::MatrixXd& rNodeValues, int rTimeDerivative, Node::eAttributes rDofType) const
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
            rNodeValues.block<3, 1>(0, iNode) = node->GetCoordinates3D();
            break;

        case Node::DISPLACEMENTS:
            rNodeValues.block<3, 1>(0, iNode) = node->GetDisplacements3D(rTimeDerivative);
            break;

        case Node::TEMPERATURES:
            rNodeValues(0, iNode) = node->GetTemperature(rTimeDerivative);
            break;

        case Node::NONLOCALEQSTRAIN:
            rNodeValues(0, iNode) = node->GetNonlocalEqStrain(rTimeDerivative);
            break;

        case Node::WATERVOLUMEFRACTION:
            rNodeValues(0, iNode) = node->GetWaterVolumeFraction(rTimeDerivative);
            break;

        case Node::RELATIVEHUMIDITY:
            rNodeValues(0, iNode) = node->GetRelativeHumidity(rTimeDerivative);
            break;

        default:
            throw MechanicsException("[NuTo::Element3D::ExtractNodeValues] Not implemented for " + Node::AttributeToString(rDofType) + ".");
        }
    }
}

const Eigen::VectorXi NuTo::Element3D::CalculateGlobalRowDofs() const
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
                globalRowDofs[index++] = nodePtr->GetDofDisplacement(2);
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
                throw MechanicsException("[NuTo::Element3D::CalculateGlobalRowDofs] Not implemented for " + Node::AttributeToString(dof) + ".");

            }
        }
    }
    return globalRowDofs;
}

const Eigen::VectorXi NuTo::Element3D::CalculateGlobalColumnDofs() const
{
    int numNonlocalElements = GetNumNonlocalElements();
    if (numNonlocalElements == 0)
        return CalculateGlobalRowDofs();
    else
    {
        throw NuTo::MechanicsException("[NuTo::Element3D::CalculateGlobalColumnDofs] not implemented for nonlocal integral element formulation.");
    }
}

const NuTo::DeformationGradient3D NuTo::Element3D::CalculateDeformationGradient(const Eigen::MatrixXd& rDerivativeShapeFunctionsLocal, const Eigen::MatrixXd& rNodalDisplacements) const
{
    DeformationGradient3D deformationGradient;

    assert(rDerivativeShapeFunctionsLocal.rows() == rNodalDisplacements.cols());
    assert(rDerivativeShapeFunctionsLocal.cols() == rNodalDisplacements.rows());

    Eigen::Matrix3d deformationGradientEigen = rNodalDisplacements.lazyProduct(rDerivativeShapeFunctionsLocal);
    deformationGradientEigen.noalias() += Eigen::Matrix3d::Identity();
    Eigen::Map<Eigen::Matrix3d>(deformationGradient.mDeformationGradient) = deformationGradientEigen;

    return deformationGradient;
}

//! @brief adds to a matrix the product B^tCB, where B contains the derivatives of the shape functions and C is the constitutive tangent
//! eventually include also area/width of an element (that's the thermal solution)
//! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param ConstitutiveTangentBase constitutive tangent matrix
//! @param rFactor factor including determinant of Jacobian and IP weight
//! @param rRow row, where to start to add the submatrix
//! @param rCoefficientMatrix to be added to
void NuTo::Element3D::AddDetJBtCB(const Eigen::MatrixXd& rDerivativeShapeFunctionsGlobal, const ConstitutiveTangentLocal<6, 6>& rConstitutiveTangent, double rFactor, int rRow, int rCol,
        FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCoefficientMatrix) const
{

	Eigen::Matrix<double, 6 , Eigen::Dynamic> Bmat;
	int numDofs = 3*GetNumNodes();
	Bmat.setZero(6,numDofs);
	int theColumn(0);
    for (int theNode1 = 0; theNode1 < GetNumNodes(); theNode1++, theColumn+=3)
    {
    	Bmat(0,theColumn)  = rDerivativeShapeFunctionsGlobal(theNode1, 0);
    	Bmat(1,theColumn+1)= rDerivativeShapeFunctionsGlobal(theNode1, 1);
    	Bmat(2,theColumn+2)= rDerivativeShapeFunctionsGlobal(theNode1, 2);
    	Bmat(3,theColumn+1)= rDerivativeShapeFunctionsGlobal(theNode1, 2);
    	Bmat(3,theColumn+2)= rDerivativeShapeFunctionsGlobal(theNode1, 1);
    	Bmat(4,theColumn)  = rDerivativeShapeFunctionsGlobal(theNode1, 2);
    	Bmat(4,theColumn+2)= rDerivativeShapeFunctionsGlobal(theNode1, 0);
    	Bmat(5,theColumn)  = rDerivativeShapeFunctionsGlobal(theNode1, 1);
    	Bmat(5,theColumn+1)= rDerivativeShapeFunctionsGlobal(theNode1, 0);
    }

    const Eigen::Matrix<double, 6 , 6>& Cmod(rConstitutiveTangent*rFactor);
    rCoefficientMatrix.block(rRow, rCol, numDofs, numDofs).noalias() += Bmat.transpose()*Cmod*Bmat;

}

//! @brief adds up the internal force vector
//! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param rEngineeringStress stress
//! @param factor factor including det Jacobian area and integration point weight
//! @param rRow start row (in case of a multifield problem)
//! @param rResult resforce vector
void NuTo::Element3D::AddDetJBtSigma(const Eigen::MatrixXd& rDerivativeShapeFunctionsLocal, const EngineeringStress3D& rEngineeringStress, double rFactor, int rRow,
        FullVector<double, Eigen::Dynamic>& rResult) const
{
    const double *s = rEngineeringStress.GetData();

    double x1, y1, z1;
    for (int theNode1 = 0; theNode1 < rDerivativeShapeFunctionsLocal.rows(); theNode1++)
    {
        int node1mul3 = 3 * theNode1;
        int node1mul3plus1 = node1mul3 + 1;
        int node1mul3plus2 = node1mul3plus1 + 1;

        assert(rResult.GetNumRows() > node1mul3plus2);
        x1 = rFactor * rDerivativeShapeFunctionsLocal(theNode1, 0);
        y1 = rFactor * rDerivativeShapeFunctionsLocal(theNode1, 1);
        z1 = rFactor * rDerivativeShapeFunctionsLocal(theNode1, 2);

        rResult(rRow + node1mul3) += x1 * s[0] + y1 * s[3] + z1 * s[5];
        rResult(rRow + node1mul3plus1) += y1 * s[1] + x1 * s[3] + z1 * s[4];
        rResult(rRow + node1mul3plus2) += z1 * s[2] + y1 * s[4] + x1 * s[5];
    }
}

//! @brief calculates the Kee matrix
//! @param rShapeFunctions of the ip for all shape functions
//! @param rDerivativeShapeFunctions of the ip for all shape functions
//! @param nonlocal gradient radius xi
//! @param rFactor multiplication factor (detJ area..)
//! @param Kee return matrix with detJ * (Nt N + cBtB)
void NuTo::Element3D::CalculateKee(Eigen::VectorXd rShapeFunctions, const Eigen::MatrixXd& rDerivativeShapeFunctions, ConstitutiveTangentLocal<1, 1>& rNonlocalParameter, double rFactor,
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
void NuTo::Element3D::AddDetJRnonlocalEqStrain(const Eigen::VectorXd& rShapeFunctions, LocalEqStrain& rLocalEqStrain, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rKee,
        Eigen::MatrixXd& rNodeNonlocalEqStrain, double rFactor, int startIndexNonlocalEqStrain, FullVector<double, Eigen::Dynamic>& rResult) const
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
void NuTo::Element3D::AddDetJBtdSigmadNonlocalEqStrainN(const Eigen::MatrixXd& rDerivativeShapeFunctions, ConstitutiveTangentLocal<6, 1>& rTangentStressNonlocalEqStrain,
        Eigen::VectorXd rShapeFunctions, double rFactor, int rRow, int rCol, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rResult) const
{


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
        double dphi_dz = rDerivativeShapeFunctions(theNode, 2);

        rResult.block(rRow + 3 * theNode    , rCol, 1, numNodesNonlocalEqStrain) += rFactor*(dphi_dx * rTangentStressNonlocalEqStrain(0) + dphi_dy * rTangentStressNonlocalEqStrain(3) + dphi_dz * rTangentStressNonlocalEqStrain(5)) * rShapeFunctions.transpose();
        rResult.block(rRow + 3 * theNode + 1, rCol, 1, numNodesNonlocalEqStrain) += rFactor*(dphi_dy * rTangentStressNonlocalEqStrain(1) + dphi_dx * rTangentStressNonlocalEqStrain(3) + dphi_dz * rTangentStressNonlocalEqStrain(4)) * rShapeFunctions.transpose();
        rResult.block(rRow + 3 * theNode + 2, rCol, 1, numNodesNonlocalEqStrain) += rFactor*(dphi_dz * rTangentStressNonlocalEqStrain(2) + dphi_dy * rTangentStressNonlocalEqStrain(4) + dphi_dx * rTangentStressNonlocalEqStrain(5)) * rShapeFunctions.transpose();
    }
}
//! @brief add detJ N_transpose dEqStrain/dEpsilon B
//! @param rShapeFunctions of the ip for the nonlocal eq strain dofs
//! @param rTangentLocalEqStrainStrain derivative of the local eq strains with respect to the strain
//! @param rderivativeShapeFunctions of the ip for the displacement dofs
//! @param rFactor factor including detJ and area
//! @param rResult result
void NuTo::Element3D::AddDetJNtdLocalEqStraindEpsilonB(Eigen::VectorXd rShapeFunctions, ConstitutiveTangentLocal<6, 1>& rTangentLocalEqStrainStrain, const Eigen::MatrixXd& rDerivativeShapeFunctions,
        double rFactor, int rRow, int rCol, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rResult) const
{

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
        double dphi_dz = rDerivativeShapeFunctions(theNode, 2);

        rResult.block(rRow, rCol+ 3 * theNode    ,numNodesNonlocalEqStrain,1) -= rFactor*(dphi_dx * rTangentLocalEqStrainStrain(0) + dphi_dy * rTangentLocalEqStrainStrain(3) + dphi_dz * rTangentLocalEqStrainStrain(5)) * rShapeFunctions;
        rResult.block(rRow, rCol+ 3 * theNode + 1,numNodesNonlocalEqStrain,1) -= rFactor*(dphi_dy * rTangentLocalEqStrainStrain(1) + dphi_dx * rTangentLocalEqStrainStrain(3) + dphi_dz * rTangentLocalEqStrainStrain(4)) * rShapeFunctions;
        rResult.block(rRow, rCol+ 3 * theNode + 2,numNodesNonlocalEqStrain,1) -= rFactor*(dphi_dz * rTangentLocalEqStrainStrain(2) + dphi_dy * rTangentLocalEqStrainStrain(4) + dphi_dx * rTangentLocalEqStrainStrain(5)) * rShapeFunctions;

    }
}

//! @brief adds up the constitutive Tangent times the derivative shape functions
//! @param rDerivativeShapeFunctions the derivative shape functions
//! @param rConstitutiveTangent the result given by the constitutive law
//! @param factor factor including det Jacobian area and integration point weight
//! @param rRow start row in case of a multifield problem
//! @param rResult result vector
void NuTo::Element3D::AddDetJBtX(const Eigen::MatrixXd &rDerivativeShapeFunctions,
                                 ConstitutiveTangentLocal<3, 1> &rConstitutiveTangent,
                                 double rFactor,
                                 int rRow,
                                 FullVector<double, Eigen::Dynamic> &rResult) const
{
    Eigen::VectorXd tmpfactor = rConstitutiveTangent * rFactor;
    int NumDofs      = rDerivativeShapeFunctions.rows();

    assert(rDerivativeShapeFunctions.cols() == GetStructure()->GetDimension());

    Eigen::VectorXd result = rDerivativeShapeFunctions * tmpfactor;
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
void NuTo::Element3D::AddDetJNtXN(const Eigen::VectorXd &rShapeFunction1,
                                  const Eigen::VectorXd &rShapeFunction2,
                                  const ConstitutiveTangentLocal<1, 1> &rConstitutiveTangent,
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

//! @brief adds to a matrix the product B^t X N, where B is the derivative shape functions, N the shape function and X the constitutive tangent
//! @param rDerivativeShapeFunction derivative shape function
//! @param rShapeFunction shape function
//! @param ConstitutiveTangentBase constitutive tangent matrix
//! @param rFactor factor including area, determinant of Jacobian and IP weight
//! @param rRow row, where to start to add the submatrix
//! @param rCol col, where to start to add the submatrix
//! @param rResult result
void NuTo::Element3D::AddDetJBtXN(const Eigen::MatrixXd &rDerivativeShapeFunction,
                                  const Eigen::VectorXd &rShapeFunction,
                                  const ConstitutiveTangentLocal<3, 1> &rConstitutiveTangent,
                                  double rFactor,
                                  int rRow,
                                  int rCol,
                                  FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> &rResult) const
{
    int NumRowDofs  = rDerivativeShapeFunction.rows();
    int NumColDofs  = rShapeFunction.rows();

    Eigen::VectorXd tmpfactor = rConstitutiveTangent * rFactor;

    assert(rShapeFunction.cols() == 1);
    assert(rDerivativeShapeFunction.cols() == GetStructure()->GetDimension());


    Eigen::MatrixXd result =  rDerivativeShapeFunction * tmpfactor * rShapeFunction.transpose();

    rResult.block(rRow,rCol,NumRowDofs,NumColDofs)+=result;
}

//! @brief adds to a matrix the product B1^t X B2, where B1 and B2 are the derivative shape functions and X is the constitutive tangent
//! @param rDerivativeShapeFunction1 derivative shape function 1
//! @param rDerivativeShapeFunction2 derivative shape function 2
//! @param ConstitutiveTangentBase constitutive tangent matrix
//! @param rFactor factor including area, determinant of Jacobian and IP weight
//! @param rRow row, where to start to add the submatrix
//! @param rCol col, where to start to add the submatrix
//! @param rResult result
void NuTo::Element3D::AddDetJBtXB(const Eigen::MatrixXd &rDerivativeShapeFunctions1,
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


//! @brief adds up the constitutive Tangent times the Shape Functions
//! @param rShapeFunctions the shape functions
//! @param rConstitutiveTangent the result given by the constitutive law
//! @param factor factor including det Jacobian area and integration point weight
//! @param rRow start row in case of a multifield problem
//! @param rResult result vector
void NuTo::Element3D::AddDetJNtX(Eigen::VectorXd &rShapeFunctions,
                                 ConstitutiveTangentLocal<1, 1> &rConstitutiveTangent,
                                 double rFactor,
                                 int rRow,
                                 FullVector<double, Eigen::Dynamic> &rResult) const
{
    Eigen::VectorXd tmpfactor = rConstitutiveTangent * rFactor;
    int NumDofs      = rShapeFunctions.rows();

    assert(rShapeFunctions.cols() == 1);

    Eigen::VectorXd result =  rShapeFunctions * tmpfactor;

    rResult.block(rRow, 0, NumDofs, 1) += result;
}

void NuTo::Element3D::CalculateJacobian(const Eigen::MatrixXd& rDerivativeShapeFunctions, const Eigen::MatrixXd& rNodeCoordinates, Eigen::Matrix3d& rInvJacobian, double& rDetJac) const
{
    /*       jacobian
     j0, j2,
     j1, j3 */

    assert(rDerivativeShapeFunctions.rows() == rNodeCoordinates.cols());
    assert(rDerivativeShapeFunctions.cols() == 3);
    assert(rNodeCoordinates.rows() == 3);

    Eigen::Matrix3d jacobian = rNodeCoordinates.lazyProduct(rDerivativeShapeFunctions);
    rDetJac = jacobian.determinant();

    if (rDetJac == 0)
        throw MechanicsException("[NuTo::Element3D::CalculateJacobian] Determinant of the Jacobian is zero, no inversion possible.");

    rInvJacobian = jacobian.inverse();
}

void NuTo::Element3D::CheckElement()
{

    int numIntegrationPoints = GetNumIntegrationPoints();
    // check number of integration points
    if (numIntegrationPoints < 1)
    {
        throw MechanicsException("[NuTo::Element3D::CheckElement] invalid integration type.");
    }

    int theIP = 0;
    const Eigen::MatrixXd& derivativeShapeFunctions = mInterpolationType->Get(Node::COORDINATES).GetDerivativeShapeFunctionsNatural(theIP);
    Eigen::MatrixXd nodeCoordinates = ExtractNodeValues(0, Node::COORDINATES);
    Eigen::Matrix3d invJacobian;
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
            throw MechanicsException("[NuTo::Element3D::CheckElement] Determinant of the Jacobian <= zero, no inversion possible.");
        volume += this->GetIntegrationPointWeight(iIP) * detJacobian;
    }

    // check element volume
    if (volume < 1e-14)
    {
        throw MechanicsException("[NuTo::Element3D::CheckElement] element with zero volume (check nodes).");
    }
}

