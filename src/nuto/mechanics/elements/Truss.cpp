// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/math/FullMatrix.h"

#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"
#include "nuto/mechanics/constitutive/mechanics/Damage.h"
#include "nuto/mechanics/constitutive/mechanics/Displacements1D.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress1D.h"
#include "nuto/mechanics/constitutive/mechanics/LocalEqPlasticStrain.h"
#include "nuto/mechanics/constitutive/mechanics/LocalEqStrain.h"
#include "nuto/mechanics/constitutive/mechanics/NonlocalEqPlasticStrain.h"
#include "nuto/mechanics/constitutive/mechanics/NonlocalEqStrain.h"
#include "nuto/mechanics/constitutive/thermal/HeatFlux1D.h"
#include "nuto/mechanics/constitutive/thermal/HeatFlux3D.h"
#include "nuto/mechanics/constitutive/thermal/Temperature.h"
#include "nuto/mechanics/constitutive/thermal/TemperatureGradient1D.h"
#include "nuto/mechanics/constitutive/thermal/TemperatureGradient3D.h"
#include "nuto/mechanics/constitutive/moistureTransport/RelativeHumidity.h"
#include "nuto/mechanics/constitutive/moistureTransport/WaterPhaseFraction.h"
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixInt.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixDouble.h"
#include "nuto/mechanics/elements/ElementOutputVectorInt.h"
#include "nuto/mechanics/elements/Truss.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionTruss.h"
#include "nuto/mechanics/structures/StructureBase.h"

//! @brief constructor
NuTo::Truss::Truss(const StructureBase* rStructure, ElementData::eElementDataType rElementDataType,
		IntegrationType::eIntegrationType rIntegrationType, IpData::eIpDataType rIpDataType) :
        NuTo::ElementBase::ElementBase(rStructure, rElementDataType, rIntegrationType, rIpDataType)
{
    mSection = 0;
}

//! @brief calculates output data fo the elmement
//! @param eOutput ... coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
//!                    @param updateStaticData (with DummyOutput), IPData, globalrow/column dofs etc.
NuTo::Error::eError NuTo::Truss::Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput)
{
    if (mStructure->GetHessianConstant(1)==false)
    	throw MechanicsException("[NuTo::Truss::Evaluate] only implemented for a constant Hessian for the first derivative (damping).");
    if (mStructure->GetHessianConstant(2)==false)
    	throw MechanicsException("[NuTo::Truss::Evaluate] only implemented for a constant Hessian for the second derivative (mass).");

    try
    {
		// get section information determining which input on the constitutive level should be used
		const SectionBase* section(GetSection());
		if (section==0)
			throw MechanicsException("[NuTo::Truss::Evaluate] no section allocated for element.");

		//calculate coordinates
		int numCoordinates(GetNumNodesGeometry());
		std::vector<double> localNodeCoord(numCoordinates);
		CalculateLocalCoordinates(localNodeCoord);

		//calculate local displacements, velocities and accelerations
		//the difference between disp and dispdof is a problem where the displacements are fixed, but enter the constitutive equation
		//for example in a two stage problem, first solve mechanics, then thermal and so on
		int numDisp(GetNumNodesField());
		int numDispDofs = section->GetIsDisplacementDof() ? numDisp : 0;

		int numTemp(GetNumNodesField());
		int numTempDofs(section->GetIsTemperatureDof() ? numTemp : 0);

		int numNonlocalEqPlasticStrain(2*GetNumNodesField());
		int numNonlocalEqPlasticStrainDofs(section->GetIsNonlocalEqPlasticStrainDof() ? numNonlocalEqPlasticStrain : 0);

		int numNonlocalTotalStrain(GetNumShapeFunctionsNonlocalTotalStrain());
		int numNonlocalTotalStrainDofs(section->GetIsNonlocalTotalStrainDof() ? numNonlocalTotalStrain : 0);

		int numNonlocalEqStrain = GetNumShapeFunctionsNonlocalEqStrain();
		int numNonlocalEqStrainDofs = section->GetIsNonlocalEqStrainDof() ? numNonlocalEqStrain : 0;

        int numWaterPhaseFraction = GetNumNodesField();
        int numWaterPhaseFractionDofs = section->GetIsWaterPhaseFractionDof() ? numWaterPhaseFraction : 0;

        int numRelativeHumidity = GetNumNodesField();
        int numRelativeHumidityDofs = section->GetIsRelativeHumidityDof() ? numRelativeHumidity : 0;

		// node values of the Dofs
        std::vector<double> localNodeDisp,nodeTemp,nodeNonlocalEqPlasticStrain,nodeNonlocalTotalStrain, nodeNonlocalEqStrain, nodeWaterPhaseFraction, nodeRelativeHumidity;

		//calculate local displacements, velocities and accelerations
		if (numDispDofs>0)
		{
			localNodeDisp.resize(numDisp);
			CalculateLocalDisplacements(0,localNodeDisp);
		}
		//calculate temperatures, temperature rate and temperature a
		if (numTempDofs>0 || section->GetInputConstitutiveIsTemperature())
		{
			nodeTemp.resize(numTemp);
			CalculateNodalTemperatures(0,nodeTemp);
		}
		if (numNonlocalEqPlasticStrainDofs>0 || section->GetInputConstitutiveIsNonlocalEqPlasticStrain())
		{
			nodeNonlocalEqPlasticStrain.resize(numNonlocalEqPlasticStrain);
			CalculateNodalNonlocalEqPlasticStrain(0,nodeNonlocalEqPlasticStrain);
		}
		if (numNonlocalTotalStrainDofs>0 || section->GetInputConstitutiveIsNonlocalTotalStrain())
		{
			nodeNonlocalTotalStrain.resize(numNonlocalTotalStrain);
			CalculateNodalNonlocalTotalStrain(0,nodeNonlocalTotalStrain);
		}
		if (numNonlocalEqStrainDofs>0 || section->GetInputConstitutiveIsNonlocalEqStrain())
		{
		    nodeNonlocalEqStrain.resize(numNonlocalEqStrain);
		    CalculateNodalNonlocalEqStrain(0, nodeNonlocalEqStrain);
		}
        if (numWaterPhaseFractionDofs>0 || section->GetInputConstitutiveIsWaterPhaseFraction())
        {
            nodeWaterPhaseFraction.resize(numWaterPhaseFraction);
            CalculateNodalWaterPhaseFraction(0,nodeWaterPhaseFraction);

        }
        if (numRelativeHumidityDofs>0 || section->GetInputConstitutiveIsRelativeHumidity())
        {
            nodeRelativeHumidity.resize(numRelativeHumidity);
            CalculateNodalRelativeHumidity(0,nodeRelativeHumidity);
        }

		//allocate space for local ip coordinates
		double localIPCoord;


		//allocate space for derivatives of shape functions
		std::vector<double> shapeFunctionsField(GetNumNodesField());
		std::vector<double> shapeFunctionsNonlocalTotalStrain(GetNumShapeFunctionsNonlocalTotalStrain());
        std::vector<double> shapeFunctionsNonlocalEqStrain(GetNumShapeFunctionsNonlocalEqStrain());
        std::vector<double> shapeFunctionsMoistureTransport(GetNumShapeFunctionsMoistureTransport());



		//allocate space for local shape functions
        std::vector<double> derivativeShapeFunctionsGeometryNatural(GetLocalDimension()*GetNumNodesGeometry());

        std::vector<double> derivativeShapeFunctionsFieldNatural(GetLocalDimension()*GetNumNodesField());
        std::vector<double> derivativeShapeFunctionsFieldLocal(GetLocalDimension()*GetNumNodesField());

        std::vector<double> derivativeShapeFunctionsNonlocalTotalStrainNatural(GetLocalDimension()*GetNumShapeFunctionsNonlocalTotalStrain());
        std::vector<double> derivativeShapeFunctionsNonlocalTotalStrainLocal(GetLocalDimension()*GetNumShapeFunctionsNonlocalTotalStrain());

        std::vector<double> derivativeShapeFunctionsNonlocalEqStrainNatural(GetLocalDimension()*GetNumShapeFunctionsNonlocalEqStrain());
        std::vector<double> derivativeShapeFunctionsNonlocalEqStrainLocal(GetLocalDimension()*GetNumShapeFunctionsNonlocalEqStrain());

        std::vector<double> derivativeShapeFunctionsMoistureTransportNatural(GetLocalDimension()*GetNumShapeFunctionsMoistureTransport());
        std::vector<double> derivativeShapeFunctionsMoistureTransportLocal(GetLocalDimension()*GetNumShapeFunctionsMoistureTransport());


		//allocate deformation gradient
		DeformationGradient1D deformationGradient;

		EngineeringStrain3D engineeringStrain3D;

		//allocate global engineering plastic strain
		EngineeringStrain3D engineeringPlasticStrain3D;

		//allocate  damage (output of constitutive relation)
		Damage damage;

		//allocate nonlocal eq plastic strain (gauss point value, input of constitutive relation)
		NonlocalEqPlasticStrain nonlocalEqPlasticStrain;

        //allocate nonlocal eq strain (gauss point value, input of constitutive relation)
        NonlocalEqStrain nonlocalEqStrain;

		//allocate local eq plastic strain output of constitutive relation
		LocalEqPlasticStrain localEqPlasticStrain;

        //allocate local eq strain output of constitutive relation
        LocalEqStrain localEqStrain;

		EngineeringStrain1D nonlocalTotalStrain;

		EngineeringStrain1D localTotalStrain;

		//allocate temperature
		Temperature temperature;

		//allocate temperature gradient
		TemperatureGradient1D temperatureGradient1D;
		TemperatureGradient3D temperatureGradient3D;

		//allocate global engineering stress
		EngineeringStress1D engineeringStress1D;
		EngineeringStress3D engineeringStress3D;

		//allocate global heat flux
		HeatFlux1D heatFlux1D;
		HeatFlux3D heatFlux3D;

        //allocate relative humidity
        RelativeHumidity relativeHumidity;

        //allocate water phase fraction
        WaterPhaseFraction waterPhaseFraction;


		//allocate tangents
		ConstitutiveTangentLocal<1,1> tangentStressStrain;
		ConstitutiveTangentLocal<1,1> tangentStressTemperature;
		ConstitutiveTangentLocal<1,2> tangentStressNonlocalEqPlasticStrain;
        ConstitutiveTangentLocal<1,1> tangentStressNonlocalEqStrain;
		ConstitutiveTangentLocal<1,1> tangentStressNonlocalTotalStrain;
		ConstitutiveTangentLocal<1,1> tangentHeatFluxTemperatureGradient;
		ConstitutiveTangentLocal<1,1> tangentHeatFluxTemperatureRate;
		ConstitutiveTangentLocal<2,1> tangentLocalEqPlasticStrainStrain;
        ConstitutiveTangentLocal<1,1> tangentLocalEqStrainStrain;

        // allocate Moisture Transport Outputs
        ConstitutiveTangentLocal<1,1> tangentVaporPhaseDiffusionCoefficient;
        ConstitutiveTangentLocal<1,1> tangentWaterPhaseDiffusionCoefficient;
        ConstitutiveTangentLocal<1,1> tangentPhaseMassExchangeRate;
        ConstitutiveTangentLocal<1,1> tangentPhaseMassExchangeRateTimesEquilibriumSorptionCurve;
        ConstitutiveTangentLocal<1,1> tangentWaterPhaseDensity;
        ConstitutiveTangentLocal<1,1> tangentVaporPhaseSaturationDensityTimesVaporPhaseVolumeFraction;
        ConstitutiveTangentLocal<1,1> tangentVaporPhaseSaturationDensityTimesRelativeHumidity;

		//define inputs and outputs
		std::map< NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase* > constitutiveInputList;
		std::map< NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase* > constitutiveOutputList;

		if (numDispDofs>0)
		{
			constitutiveInputList[NuTo::Constitutive::Input::DEFORMATION_GRADIENT_1D] = &deformationGradient;
		}

		if (numTempDofs>0)
		{
			constitutiveInputList[NuTo::Constitutive::Input::TEMPERATURE_GRADIENT_1D] = &temperatureGradient1D;
		}

		if (section->GetInputConstitutiveIsTemperature())
		{
			constitutiveInputList[NuTo::Constitutive::Input::TEMPERATURE] = &temperature;
		}

		if (numNonlocalEqPlasticStrainDofs>0)
		{
			constitutiveInputList[NuTo::Constitutive::Input::NONLOCAL_EQ_PLASTIC_STRAIN] = &nonlocalEqPlasticStrain;
		}

		if (numNonlocalTotalStrainDofs>0)
		{
			constitutiveInputList[NuTo::Constitutive::Input::NONLOCAL_TOTAL_STRAIN_1D] = &nonlocalTotalStrain;
		}
        if (numNonlocalEqStrainDofs>0)
        {
            constitutiveInputList[NuTo::Constitutive::Input::NONLOCAL_EQ_STRAIN] = &nonlocalEqStrain;
        }
        if (numWaterPhaseFractionDofs>0)
        {
            constitutiveInputList[NuTo::Constitutive::Input::WATER_PHASE_FRACTION] = &waterPhaseFraction;
        }
        if (numRelativeHumidityDofs>0)
        {
            constitutiveInputList[NuTo::Constitutive::Input::RELATIVE_HUMIDITY] = &relativeHumidity;
        }
        // sum of all dofs
        int numDofs = numDispDofs
                + numTempDofs
                + numNonlocalEqPlasticStrainDofs
                + numNonlocalTotalStrainDofs
                + numNonlocalEqStrainDofs
                + numRelativeHumidityDofs
                + numWaterPhaseFractionDofs;

		//define outputs
		for (auto it = rElementOutput.begin(); it!=rElementOutput.end(); it++)
		{
			switch(it->first)
			{
			case Element::INTERNAL_GRADIENT:
				it->second->GetFullVectorDouble().Resize(numDofs);
				//if the stiffness matrix is constant, the corresponding internal force is calculated via the Kd
				//on the global level
				if (mStructure->GetHessianConstant(0)==false)
				{
					if (numDispDofs>0)
					{
						constitutiveOutputList[NuTo::Constitutive::Output::ENGINEERING_STRESS_1D] = &engineeringStress1D;
					}
					if (numTempDofs>0)
					{
						constitutiveOutputList[NuTo::Constitutive::Output::HEAT_FLUX_1D] = &heatFlux1D;
					}
					if (numNonlocalEqPlasticStrainDofs>0)
					{
						constitutiveOutputList[NuTo::Constitutive::Output::LOCAL_EQ_PLASTIC_STRAIN] = &localEqPlasticStrain;
					}
					if (numNonlocalTotalStrainDofs>0)
					{
						constitutiveOutputList[NuTo::Constitutive::Output::ENGINEERING_STRAIN_1D] = &localTotalStrain;
					}
                    if (numNonlocalEqStrainDofs>0)
                    {
                        constitutiveOutputList[NuTo::Constitutive::Output::LOCAL_EQ_STRAIN] = &localEqStrain;
                    }
				}
			break;
			case Element::HESSIAN_0_TIME_DERIVATIVE:
				{
					it->second->GetFullMatrixDouble().Resize(numDofs, numDofs);
					it->second->SetSymmetry(true);
					it->second->SetConstant(true);
					if (numDispDofs>0)
					{
						constitutiveOutputList[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_1D] = &tangentStressStrain;
						//mixed terms
						if (numTempDofs>0)
							constitutiveOutputList[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_TEMPERATURE_1D] = &tangentStressTemperature;
						if (numNonlocalEqPlasticStrainDofs>0)
							constitutiveOutputList[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_PLASTIC_STRAIN_1D] = &tangentStressNonlocalEqPlasticStrain;
						if (numNonlocalTotalStrainDofs>0)
							constitutiveOutputList[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_TOTAL_STRAIN_1D] = &tangentStressNonlocalTotalStrain;
                        if (numNonlocalEqStrainDofs>0)
                            constitutiveOutputList[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_1D] = &tangentStressNonlocalEqStrain;
                    }
					if (numTempDofs>0)
					{
						constitutiveOutputList[NuTo::Constitutive::Output::D_HEAT_FLUX_D_TEMPERATURE_RATE_1D] = &tangentHeatFluxTemperatureGradient;
						//mixed terms
						//if (numDisp)
						//    constitutiveOutputList.insert(std::pair<NuTo::Constitutive::eOutput, ConstitutiveOutputBase*>(NuTo::Constitutive::eOutput::D_HEAT_FLUX_D_ENGINEERING_STRAIN_3D, &tangentHeatFluxEngineeringStrain[timeDerivative]));
					}
					if (numNonlocalEqPlasticStrainDofs>0 && numDispDofs>0)
					{
                        constitutiveOutputList[NuTo::Constitutive::Output::D_LOCAL_EQ_PLASTIC_STRAIN_D_STRAIN_1D] = &tangentLocalEqPlasticStrainStrain;
					}

                    if (numNonlocalEqStrainDofs>0 && numDispDofs>0)
                    {
                        constitutiveOutputList[NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_D_STRAIN_1D] = &tangentLocalEqStrainStrain;
                    }
                    if (numRelativeHumidityDofs || numWaterPhaseFractionDofs)
                    {
                        if(numRelativeHumidityDofs != numWaterPhaseFractionDofs)
                        {
                            throw MechanicsException(std::string("[NuTo::Truss::Evaluate] Number of relative humidity dofs must be equal to the number of water phase fraction dofs. Current number of dofs:\n")+
                                                     std::string("Relative humidity dofs: ")+std::to_string(numRelativeHumidityDofs)+
                                                     std::string("\nWater phase fraction dofs: ")+std::to_string(numWaterPhaseFractionDofs) + std::string("\n"));
                        }
                        constitutiveOutputList[NuTo::Constitutive::Output::VAPOR_PHASE_DIFFUSION_COEFFICIENT]                           = &tangentVaporPhaseDiffusionCoefficient;
                        constitutiveOutputList[NuTo::Constitutive::Output::WATER_PHASE_DIFFUSION_COEFFICIENT]                           = &tangentWaterPhaseDiffusionCoefficient;
                        constitutiveOutputList[NuTo::Constitutive::Output::PHASE_MASS_EXCHANGE_RATE]                                    = &tangentPhaseMassExchangeRate;
                        constitutiveOutputList[NuTo::Constitutive::Output::PHASE_MASS_EXCHANGE_RATE_TIMES_EQUILIBRIUM_SORPTION_CURVE]   = &tangentPhaseMassExchangeRateTimesEquilibriumSorptionCurve;
                    }
				}
			break;
			case Element::HESSIAN_1_TIME_DERIVATIVE:
			{
                it->second->GetFullMatrixDouble().Resize(numDofs, numDofs);
                it->second->SetSymmetry(true);
                it->second->SetConstant(true);
				if (numDispDofs>0)
				{
					// Rayleigh damping should be introduced on the global level
				}
				if (numTempDofs>0)
				{
					constitutiveOutputList[NuTo::Constitutive::Output::D_HEAT_FLUX_D_TEMPERATURE_RATE_1D] = &tangentHeatFluxTemperatureRate;
				}
                if (numRelativeHumidityDofs || numWaterPhaseFractionDofs)
                {
                    if(numRelativeHumidityDofs != numWaterPhaseFractionDofs)
                    {
                        throw MechanicsException(std::string("[NuTo::Truss::Evaluate] Number of relative humidity dofs must be equal to the number of water phase fraction dofs. Current number of dofs:\n")+
                                                 std::string("Relative humidity dofs: ")+std::to_string(numRelativeHumidityDofs)+
                                                 std::string("\nWater phase fraction dofs: ")+std::to_string(numWaterPhaseFractionDofs) + std::string("\n"));
                    }
                    constitutiveOutputList[NuTo::Constitutive::Output::WATER_PHASE_DENSITY]                                                 = &tangentWaterPhaseDensity;
                    constitutiveOutputList[NuTo::Constitutive::Output::VAPOR_PHASE_SATURATION_DENSITY_TIMES_VAPOR_PHASE_VOLUME_FRACTION]    = &tangentVaporPhaseSaturationDensityTimesVaporPhaseVolumeFraction;
                    constitutiveOutputList[NuTo::Constitutive::Output::VAPOR_PHASE_SATURATION_DENSITY_TIMES_RELATIVE_HUMIDITY]              = &tangentVaporPhaseSaturationDensityTimesRelativeHumidity;
                }

			}
			break;
			case Element::HESSIAN_2_TIME_DERIVATIVE:
			{
				it->second->GetFullMatrixDouble().Resize(numDispDofs+numTempDofs+numNonlocalEqPlasticStrainDofs+numNonlocalTotalStrainDofs,numDispDofs+numTempDofs+numNonlocalEqPlasticStrainDofs+numNonlocalTotalStrainDofs);
				it->second->SetSymmetry(true);
				it->second->SetConstant(true);
				//there is only a constant mass part for the mechanics problem
			}
			break;
			case Element::UPDATE_STATIC_DATA:
				constitutiveOutputList[NuTo::Constitutive::Output::UPDATE_STATIC_DATA] = 0;
			break;
			case Element::UPDATE_TMP_STATIC_DATA:
				constitutiveOutputList[NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA] = 0;
			break;
			case Element::IP_DATA:
				switch(it->second->GetIpDataType())
				{
				case NuTo::IpData::ENGINEERING_STRAIN:
					it->second->GetFullMatrixDouble().Resize(6,GetNumIntegrationPoints());
					 //define outputs
					constitutiveOutputList[NuTo::Constitutive::Output::ENGINEERING_STRAIN_3D] = &engineeringStrain3D;
				break;
				case NuTo::IpData::ENGINEERING_STRESS:
					it->second->GetFullMatrixDouble().Resize(6,GetNumIntegrationPoints());
					 //define outputs
					 constitutiveOutputList[NuTo::Constitutive::Output::ENGINEERING_STRESS_3D] = &engineeringStress3D;
				break;
				case NuTo::IpData::ENGINEERING_PLASTIC_STRAIN:
					it->second->GetFullMatrixDouble().Resize(6,GetNumIntegrationPoints());
					 //define outputs
					 constitutiveOutputList[NuTo::Constitutive::Output::ENGINEERING_PLASTIC_STRAIN_3D] = &engineeringPlasticStrain3D;
					 break;
				break;
				case NuTo::IpData::DAMAGE:
					it->second->GetFullMatrixDouble().Resize(1,GetNumIntegrationPoints());
					//define outputs
					constitutiveOutputList[NuTo::Constitutive::Output::DAMAGE] = &damage;
				break;
				default:
					throw MechanicsException("[NuTo::Truss::Evaluate] this ip data type is not implemented.");
				}
			break;
			case Element::GLOBAL_ROW_DOF:
				this->CalculateGlobalRowDofs(it->second->GetVectorInt(),
				        numDispDofs,
				        numTempDofs,
				        numNonlocalEqPlasticStrainDofs,
				        numNonlocalTotalStrainDofs,
                        numNonlocalEqStrainDofs,
                        numRelativeHumidityDofs,
                        numWaterPhaseFractionDofs);
			break;
			case Element::GLOBAL_COLUMN_DOF:
				this->CalculateGlobalColumnDofs(it->second->GetVectorInt(),
				        numDispDofs,
				        numTempDofs,
				        numNonlocalEqPlasticStrainDofs,
				        numNonlocalTotalStrainDofs,
                        numNonlocalEqStrainDofs,
                        numRelativeHumidityDofs,
                        numWaterPhaseFractionDofs);
			break;
			default:
				throw MechanicsException("[NuTo::Truss::Evaluate] element output not implemented.");
			}
		}

		// loop over the integration points
		for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
		{
			GetLocalIntegrationPointCoordinates(theIP, localIPCoord);

			//derivative in natural coordinate system
			CalculateDerivativeShapeFunctionsGeometry(localIPCoord, derivativeShapeFunctionsGeometryNatural);

			//determinant of the Jacobian
			double detJ(this->DetJacobian(derivativeShapeFunctionsGeometryNatural,localNodeCoord));

            if (GetNumNodesGeometry()==GetNumNodesField())
            {
            	//isoparametric element
				//derivative in local coordinate system
				for (unsigned int count=0; count<derivativeShapeFunctionsGeometryNatural.size(); count++)
				{
					derivativeShapeFunctionsFieldLocal[count] = derivativeShapeFunctionsGeometryNatural[count]/detJ;
				}
            }
            else
            {
    			//derivative in natural coordinate system
    			CalculateDerivativeShapeFunctionsField(localIPCoord, derivativeShapeFunctionsFieldNatural);
				//derivative in local coordinate system
				for (unsigned int count=0; count<derivativeShapeFunctionsFieldNatural.size(); count++)
				{
					derivativeShapeFunctionsFieldLocal[count] = derivativeShapeFunctionsFieldNatural[count]/detJ;
				}
            }

			if (numDispDofs)
			{
				// determine deformation gradient from the local Displacements and the derivative of the shape functions
				CalculateDeformationGradient(derivativeShapeFunctionsFieldLocal, localNodeDisp, deformationGradient);
			}

			if (numTempDofs)
			{
				// determine temperature gradient from the local temperatures and the derivative of the shape functions
				// this is not included in the AddIpStiffness to avoid reallocation of the deformation gradient for each IP
				throw MechanicsException("[Nuto::Truss::Evaluate] temperature not yet implemented.");
	//TODO      CalculateTemperatureGradient(derivativeShapeFunctionsGlobal, localNodeCoord, nodeTemp, temperatureGradient);
			}

			if (section->GetInputConstitutiveIsTemperature())
			{
				// determine temperature
				throw MechanicsException("[Nuto::Truss::Evaluate] temperature not yet implemented.");
			}

			if (numNonlocalEqPlasticStrainDofs)
			{
				CalculateShapeFunctionsNonlocalTotalStrain(localIPCoord, shapeFunctionsNonlocalTotalStrain);
				CalculateNonlocalEqPlasticStrain(shapeFunctionsNonlocalTotalStrain, nodeNonlocalEqPlasticStrain, nonlocalEqPlasticStrain);
			}

			if (numNonlocalTotalStrainDofs)
			{
				CalculateShapeFunctionsNonlocalTotalStrain(localIPCoord, shapeFunctionsNonlocalTotalStrain);

				//derivative in natural coordinate system
				CalculateDerivativeShapeFunctionsNonlocalTotalStrain(localIPCoord, derivativeShapeFunctionsNonlocalTotalStrainNatural);

				//derivative in local coordinate system
				for (unsigned int count=0; count<derivativeShapeFunctionsNonlocalTotalStrainLocal.size(); count++)
				{
					derivativeShapeFunctionsNonlocalTotalStrainLocal[count] = derivativeShapeFunctionsNonlocalTotalStrainNatural[count]/detJ;
				}

				CalculateNonlocalTotalStrain(shapeFunctionsNonlocalTotalStrain, nodeNonlocalTotalStrain, nonlocalTotalStrain);
			}

            if (numNonlocalEqStrainDofs)
            {
                CalculateShapeFunctionsNonlocalEqStrain(localIPCoord, shapeFunctionsNonlocalEqStrain);
                CalculateNonlocalEqStrain(shapeFunctionsNonlocalEqStrain, nodeNonlocalEqStrain, nonlocalEqStrain);
                CalculateDerivativeShapeFunctionsNonlocalEqStrain(localIPCoord, derivativeShapeFunctionsNonlocalEqStrainNatural);
                for (unsigned int count=0; count<derivativeShapeFunctionsNonlocalEqStrainNatural.size(); count++)
                {
                    derivativeShapeFunctionsNonlocalEqStrainLocal[count] = derivativeShapeFunctionsNonlocalEqStrainNatural[count]/detJ;
                }

            }
            if (numRelativeHumidityDofs || numWaterPhaseFractionDofs)
            {
                if(numRelativeHumidityDofs != numWaterPhaseFractionDofs)
                {
                    throw MechanicsException(std::string("[NuTo::Truss::Evaluate] Number of relative humidity dofs must be equal to the number of water phase fraction dofs. Current number of dofs:\n")+
                                             std::string("Relative humidity dofs: ")+std::to_string(numRelativeHumidityDofs)+
                                             std::string("\nWater phase fraction dofs: ")+std::to_string(numWaterPhaseFractionDofs) + std::string("\n"));
                }
                CalculateShapeFunctionsMoistureTransport(localIPCoord,shapeFunctionsMoistureTransport);
                CalculateRelativeHumidity(shapeFunctionsMoistureTransport,nodeRelativeHumidity,relativeHumidity);
                CalculateWaterPhaseFraction(shapeFunctionsMoistureTransport,nodeWaterPhaseFraction,waterPhaseFraction);
                CalculateDerivativeShapeFunctionsMoistureTransport(localIPCoord,derivativeShapeFunctionsMoistureTransportNatural);
                for (unsigned int count=0; count<derivativeShapeFunctionsMoistureTransportNatural.size(); count++)
                {
                    derivativeShapeFunctionsMoistureTransportLocal[count] = derivativeShapeFunctionsMoistureTransportNatural[count]/detJ;
                }

            }

			ConstitutiveBase* constitutivePtr = GetConstitutiveLaw(theIP);
			Error::eError error = constitutivePtr->Evaluate1D(this, theIP,
					constitutiveInputList, constitutiveOutputList);
			if (error!=Error::SUCCESSFUL)
				return error;



			double areaFactor = 1.;

			if (numNonlocalEqStrainDofs)
			{
			    // calculate global IP coordinates
	            double coordsIP[3];
	            GetGlobalIntegrationPointCoordinates(theIP, coordsIP);

	            areaFactor = mSection->AsSectionTruss()->GetAreaFactor(coordsIP[0]);
//	            std::cout << coordsIP[0] << ": " << areaFactor << std::endl;
			}

			double factor = detJ * mSection->GetArea() * areaFactor *
					       (mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP));


			FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> Kkk;
			if (numNonlocalEqPlasticStrainDofs > 0 || numNonlocalTotalStrainDofs > 0)
			{
				//calculate Kkk detJ*(cBtB+NtN)
				//the nonlocal radius is in a gradient formulation is different from the nonlocal radius in an integral formulation
				CalculateKkk(shapeFunctionsNonlocalTotalStrain,derivativeShapeFunctionsNonlocalTotalStrainLocal,constitutivePtr->GetNonlocalRadius(),factor,Kkk);
			} else if (numNonlocalEqStrainDofs > 0)
            {
                //calculate Kkk detJ*(cBtB+NtN)
                //the nonlocal radius is in a gradient formulation is different from the nonlocal radius in an integral formulation
                CalculateKkk(shapeFunctionsNonlocalEqStrain,derivativeShapeFunctionsNonlocalEqStrainLocal,constitutivePtr->GetNonlocalRadius(),factor,Kkk);
            }

			//calculate output
			for (auto it = rElementOutput.begin(); it!=rElementOutput.end(); it++)
			{
				switch(it->first)
				{
				case Element::INTERNAL_GRADIENT:
				{
					//if the stiffness matrix is constant, the corresponding internal force is calculated via the Kd
					//on the global level
					if (mStructure->GetHessianConstant(0)==false)
					{
						if (numDispDofs>0)
						{
							AddDetJBtSigma(derivativeShapeFunctionsFieldLocal,engineeringStress1D, factor, 0, it->second->GetFullVectorDouble());
						}
						if (numTempDofs>0)
						{
							throw MechanicsException("[Nuto::Truss::Evaluate] temperature not yet implemented.");
		//TODO				    AddDetJBtCB(derivativeShapeFunctions,tangentStressStrain, factor, 0, it->second->GetFullMatrixDouble());
						}
						if (numNonlocalEqPlasticStrainDofs>0)
						{
							//add Kkk*nonlocalEqPlasticStrain+detJ*F
							AddDetJRnonlocalPlasticStrain(shapeFunctionsField,localEqPlasticStrain, Kkk, nodeNonlocalEqPlasticStrain, factor, numDispDofs+numTempDofs, it->second->GetFullVectorDouble());
						}
						if (numNonlocalTotalStrainDofs>0)
						{
							//add Kkk*nonlocalTotalStrain+detJ*F
							AddDetJRnonlocalTotalStrain(shapeFunctionsNonlocalTotalStrain,localTotalStrain, Kkk, nodeNonlocalTotalStrain, factor, numDispDofs+numTempDofs+numNonlocalEqPlasticStrainDofs, it->second->GetFullVectorDouble());
						}
                        if (numNonlocalEqStrainDofs>0)
                        {
                            //add Kkk*nonlocalEqStrain-detJ*F
                            AddDetJRnonlocalEqStrain(shapeFunctionsNonlocalEqStrain,localEqStrain, Kkk, nodeNonlocalEqStrain, factor, numDispDofs+numTempDofs+numNonlocalEqPlasticStrainDofs+numNonlocalTotalStrainDofs, it->second->GetFullVectorDouble());
                        }
					}
				}
				break;
				case Element::HESSIAN_0_TIME_DERIVATIVE:
					{
						if (numDispDofs>0)
						{
							//derivative of F(sigma) with respect to all unknowns
							AddDetJBtCB(derivativeShapeFunctionsFieldLocal, tangentStressStrain, factor, 0,0, it->second->GetFullMatrixDouble());
							if (tangentStressStrain.GetSymmetry()==false)
								it->second->SetSymmetry(false);
							if (tangentStressStrain.GetConstant()==false)
								it->second->SetConstant(false);
							if (numTempDofs>0)
								throw MechanicsException("[NuTo::Truss::Evaluate] mixed terms not yet implemented.");
							if (numNonlocalEqPlasticStrainDofs>0)
							{
								AddDetJBtdSigmadNonlocalEqPlasticStrainN(derivativeShapeFunctionsFieldLocal, tangentStressNonlocalEqPlasticStrain,shapeFunctionsField, factor, 0,numDispDofs+numTempDofs, it->second->GetFullMatrixDouble());
							}
							if (numNonlocalTotalStrainDofs>0)
							{
								AddDetJBtdSigmadNonlocalTotalStrainN(derivativeShapeFunctionsNonlocalTotalStrainLocal, tangentStressNonlocalTotalStrain,shapeFunctionsNonlocalTotalStrain, factor, 0,numDispDofs+numTempDofs+numNonlocalEqPlasticStrainDofs, it->second->GetFullMatrixDouble());
							}
							if (numNonlocalEqStrainDofs > 0)
							{
							    AddDetJBtdSigmadNonlocalEqStrainN(derivativeShapeFunctionsFieldLocal, tangentStressNonlocalEqStrain, shapeFunctionsNonlocalEqStrain, factor, 0, numDispDofs+numTempDofs+numNonlocalEqPlasticStrainDofs+numNonlocalTotalStrainDofs, it->second->GetFullMatrixDouble());
							}

						}
						if (numTempDofs>0)
						{
							//derivative of F(temperatureGradient) with respect to all unknowns
							throw MechanicsException("[Nuto::Truss::Evaluate] temperature not yet implemented.");
							//TODO AddDetJBtCB(derivativeShapeFunctionsGlobal, tangentHeatFluxTemperatureGradient[timeDerivative], factor, numDisp,numDisp, it->second->GetFullMatrixDouble());
							if (tangentHeatFluxTemperatureGradient.GetSymmetry()==false)
								it->second->SetSymmetry(false);
							if (tangentHeatFluxTemperatureGradient.GetConstant()==false)
								it->second->SetConstant(false);
							if (numDispDofs>0)
								throw MechanicsException("[NuTo::Truss::Evaluate] mixed terms not yet implemented.");
						}
						if (numNonlocalEqPlasticStrainDofs>0)
						{
							//derivative of F(nonlocalEqPlasticStrain) with respect to all unknowns
							if (numDispDofs>0)
							{
								AddDetJNtdLocalEqPlasticStraindEpsilonB(shapeFunctionsField,tangentLocalEqPlasticStrainStrain,derivativeShapeFunctionsFieldLocal, factor, numDispDofs+numTempDofs, 0, it->second->GetFullMatrixDouble());
							}
							it->second->GetFullMatrixDouble().AddBlock(numDispDofs+numTempDofs,numDispDofs+numTempDofs,Kkk);
							it->second->GetFullMatrixDouble().AddBlock(numDispDofs+numTempDofs+this->GetNumNodes(),numDispDofs+numTempDofs+this->GetNumNodes(),Kkk);
						}
						if (numNonlocalTotalStrainDofs>0)
						{
							//derivative of F(nonlocaltotalStrain) with respect to all unknowns
							if (numDispDofs>0)
							{
								AddDetJNtB(shapeFunctionsNonlocalTotalStrain,derivativeShapeFunctionsNonlocalTotalStrainLocal, factor, numDispDofs+numTempDofs+numNonlocalEqPlasticStrainDofs, 0, it->second->GetFullMatrixDouble());
							}
							it->second->GetFullMatrixDouble().AddBlock(numDispDofs+numTempDofs+numNonlocalEqPlasticStrainDofs,numDispDofs+numTempDofs+numNonlocalEqPlasticStrainDofs,Kkk);
						}
                        if (numNonlocalEqStrainDofs>0)
                        {
                            int prevDofs = numDispDofs + numTempDofs + numNonlocalEqPlasticStrainDofs + numNonlocalTotalStrainDofs;

                            //derivative of F(nonlocalEqStrain) with respect to all unknowns
                            if (numDispDofs > 0)
                            {
                                AddDetJNtdLocalEqStraindEpsilonB(shapeFunctionsNonlocalEqStrain, tangentLocalEqStrainStrain, derivativeShapeFunctionsFieldLocal, factor, prevDofs, 0, it->second->GetFullMatrixDouble());
                            }
                            // add Kkk
                            it->second->GetFullMatrixDouble().AddBlock(prevDofs,prevDofs,Kkk);

                        }
                        if (numWaterPhaseFractionDofs>0 || numRelativeHumidityDofs>0)
                        {

                            // Assembly of the "Stiffness"-Matrix (see: B.Johannesson - A Numerical Approach for Non-Linear Moisture Flow in Porous Materials with Account to Sorption Hysteresis)
                            // | K_w + R_w      -R_w_eq         |
                            // |                                |
                            // | -R_w           K_phi + R_w_eq  |



                            // | X - |      X = K_w + R_w
                            // | - - |

                            AddDetJBtCB(derivativeShapeFunctionsMoistureTransportLocal,tangentWaterPhaseDiffusionCoefficient,factor,0,0,it->second->GetFullMatrixDouble());
                            AddDetJNtCN(shapeFunctionsMoistureTransport,tangentPhaseMassExchangeRate,factor,0,0,it->second->GetFullMatrixDouble());

                            // | - - |      X = -R_w
                            // | X - |

                            AddDetJNtCN(shapeFunctionsMoistureTransport,tangentPhaseMassExchangeRate,-1.0*factor,2,0,it->second->GetFullMatrixDouble());

                            // | - X |      X = -R_w_eq
                            // | - - |
                            AddDetJNtCN(shapeFunctionsMoistureTransport,tangentPhaseMassExchangeRateTimesEquilibriumSorptionCurve,-1.0*factor,0,2,it->second->GetFullMatrixDouble());

                            // | - - |      X = K_phi + R_w_eq
                            // | - X |
                            AddDetJBtCB(derivativeShapeFunctionsMoistureTransportLocal,tangentVaporPhaseDiffusionCoefficient,factor,2,2,it->second->GetFullMatrixDouble());
                            AddDetJNtCN(shapeFunctionsMoistureTransport,tangentPhaseMassExchangeRateTimesEquilibriumSorptionCurve,factor,2,2,it->second->GetFullMatrixDouble());

                            //it->second->GetFullMatrixDouble().Info();         // Check Element Matrix if needed
                            //NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> test = it->second->GetFullMatrixDouble();
                            //int a=0;
                        }
                        BlowLocalMatrixToGlobal(it->second->GetFullMatrixDouble());
					}
				break;
				case Element::HESSIAN_1_TIME_DERIVATIVE:
				{
					if (numDispDofs>0)
					{
						//no damping term, do Rayleigh damping on the global level
					}
					if (numTempDofs>0)
					{
						throw MechanicsException("[Nuto::Truss::Evaluate] temperature not yet implemented.");
						//TODO AddDetJBtCB(derivativeShapeFunctionsGlobal, tangentHeatFluxTemperatureGradient[timeDerivative], factor, numDisp,numDisp, it->second->GetFullMatrixDouble());
	/*					if (tangentHeatFluxTemperatureGradient[timeDerivative].GetSymmetry()==false)
							it->second->SetSymmetry(false);
						if (tangentHeatFluxTemperatureGradient[timeDerivative].GetConstant()==false)
							it->second->SetConstant(false);

						if (numDisp>0)
							throw MechanicsException("[NuTo::Truss::Evaluate] mixed terms not yet implemented.");
							*/
					}
					if (numNonlocalEqPlasticStrainDofs>0)
					{
						//no damping term, do Rayleigh damping on the global level
					}
                    if (numWaterPhaseFractionDofs>0 || numRelativeHumidityDofs>0)
                    {

                        // Assembly of the "Damping"-Matrix (see: B.Johannesson - A Numerical Approach for Non-Linear Moisture Flow in Porous Materials with Account to Sorption Hysteresis)
                        // | C_w     0     |
                        // |               |
                        // | -M_phi  C_phi |


                        // | X - |      X = C_w
                        // | - - |
                        AddDetJNtCN(shapeFunctionsMoistureTransport,tangentWaterPhaseDensity,factor,0,0,it->second->GetFullMatrixDouble());

                        // | - - |      X = -M_phi
                        // | X - |
                        AddDetJNtCN(shapeFunctionsMoistureTransport,tangentVaporPhaseSaturationDensityTimesRelativeHumidity,-1.0*factor,2,0,it->second->GetFullMatrixDouble());


                        // | - - |      X = C_phi
                        // | - X |
                        AddDetJNtCN(shapeFunctionsMoistureTransport,tangentVaporPhaseSaturationDensityTimesVaporPhaseVolumeFraction,factor,2,2,it->second->GetFullMatrixDouble());

                        //it->second->GetFullMatrixDouble().Info();         // Check Element Matrix if needed
                        //NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> test = it->second->GetFullMatrixDouble();
                        //int a=0;
                    }
					BlowLocalMatrixToGlobal(it->second->GetFullMatrixDouble());
				}
				break;
				case Element::HESSIAN_2_TIME_DERIVATIVE:
				{
					if (numDispDofs>0)
					{
						double density = constitutivePtr->GetDensity();
						this->CalculateShapeFunctionsField(localIPCoord, shapeFunctionsField);

						// calculate local mass matrix
						// don't forget to include determinant of the Jacobian and area
						// detJ * area * density * HtH, :
						double factor2 (density * factor);
						this->AddDetJHtH(shapeFunctionsField, factor2, it->second->GetFullMatrixDouble());

						//if (numTemp>0) no mixed terms
					}
					if (numTempDofs>0)
					{
						//no termperature terms
					}
					if (numNonlocalEqPlasticStrainDofs>0)
					{
						//no termperature terms
					}
					BlowLocalMatrixToGlobal(it->second->GetFullMatrixDouble());
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
						memcpy(&(it->second->GetFullMatrixDouble().data()[theIP*6]),engineeringStrain3D.GetData(),6*sizeof(double));
					break;
					case NuTo::IpData::ENGINEERING_STRESS:
						//error = constitutivePtr->GetEngineeringStressFromEngineeringStrain(this, theIP, deformationGradient, engineeringStress);
						memcpy(&(it->second->GetFullMatrixDouble().data()[theIP*6]),engineeringStress3D.GetData(),6*sizeof(double));
					break;
					case NuTo::IpData::ENGINEERING_PLASTIC_STRAIN:
						//error = constitutivePtr->GetEngineeringPlasticStrain(this, theIP, deformationGradient, engineeringStrain);
						memcpy(&(it->second->GetFullMatrixDouble().data()[theIP*6]),engineeringPlasticStrain3D.GetData(),6*sizeof(double));
					break;
					case NuTo::IpData::DAMAGE:
						//error = constitutivePtr->GetDamage(this, theIP, deformationGradient, rIpData.mEigenMatrix.data()[theIP]);
						memcpy(&(it->second->GetFullMatrixDouble().data()[theIP]),damage.GetData(),sizeof(double));
					break;
					default:
						throw MechanicsException("[NuTo::Truss::Evaluate] Ip data not implemented.");
					}
				break;
				case Element::GLOBAL_ROW_DOF:
				case Element::GLOBAL_COLUMN_DOF:
				break;
				default:
					throw MechanicsException("[NuTo::Truss::Evaluate] element output not implemented.");
				}
			}
		}
    }
    catch (NuTo::MechanicsException& e)
    {
        std::stringstream ss;
        ss << mStructure->ElementGetId(this);
    	e.AddMessage("[NuTo::Truss::Evaluate] Error evaluating element data of element"	+ ss.str() + ".");
        throw e;
    }

    return Error::SUCCESSFUL;
}

//! @brief add detJ transpose N dOmega/depsilon B
//! @param rShapeFunctions of the ip for all shape functions
//! @param rTangentLocalEqPlasticStrainStrain derivative of the local eq plastic strains with respect to the strain
//! @param rderivativeShapeFunctions of the ip for all shape functions
//! @param rFactor factor including detJ and area
//! @param rResult result
void NuTo::Truss::AddDetJNtdLocalEqPlasticStraindEpsilonB(const std::vector<double>& rShapeFunctions, ConstitutiveTangentLocal<2,1>& rTangentLocalEqPlasticStrainStrain,
		const std::vector<double>& rDerivativeShapeFunctions, double rFactor, int rRow, int rCol, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult)
{
    assert(rShapeFunctions.size()==rDerivativeShapeFunctions.size());
    for (int countNonlocalEqPlasticStrain=0; countNonlocalEqPlasticStrain<2; countNonlocalEqPlasticStrain++ )
    {
		double tmpfactor = rTangentLocalEqPlasticStrainStrain(countNonlocalEqPlasticStrain,0)*rFactor;
		//these for loops are not optimal, we might use a map to an eigenvector and perform the multiplication directly
		for (unsigned int count=0; count<rShapeFunctions.size(); count++)
		{
			for (unsigned int count2=0; count2<rDerivativeShapeFunctions.size(); count2++)
			{
				rResult(rRow+count+rShapeFunctions.size()*countNonlocalEqPlasticStrain,rCol+count2)-=
						tmpfactor*rShapeFunctions[count]*(rDerivativeShapeFunctions[count2]);
			}
		}
    }
}

//! @brief add detJ transpose N dOmega/depsilon B
//! @param rShapeFunctions of the ip for all shape functions
//! @param rTangentLocalEqStrain derivative of the local eq plastic strains with respect to the strain
//! @param rderivativeShapeFunctions of the ip for all shape functions
//! @param rFactor factor including detJ and area
//! @param rResult result
void NuTo::Truss::AddDetJNtdLocalEqStraindEpsilonB(const std::vector<double>& rShapeFunctions, ConstitutiveTangentLocal<1,1>& rTangentLocalEqStrainStrain,
        const std::vector<double>& rDerivativeShapeFunctions, double rFactor, int rRow, int rCol, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult)
{

    double tmpfactor = rTangentLocalEqStrainStrain(0)*rFactor;
    //these for loops are not optimal, we might use a map to an eigenvector and perform the multiplication directly
    // ttitsche:: I did the benchmark: the conversion / mapping to Eigen::Matrices is too slow.
    for (unsigned int iN=0; iN<rShapeFunctions.size(); iN++) // N-Matrix
    {
        for (unsigned int iB=0; iB<rDerivativeShapeFunctions.size(); iB++) // B-Matrix
        {
            rResult(rRow+iN,rCol+iB) -= tmpfactor*rShapeFunctions[iN]*rDerivativeShapeFunctions[iB];
        }
    }


}

//! @brief add detJ transpose N dOmega/depsilon B
//! @param rShapeFunctions of the ip for all shape functions
//! @param rTangentLocalEqPlasticStrainStrain derivative of the local eq plastic strains with respect to the strain
//! @param rderivativeShapeFunctions of the ip for all shape functions
//! @param rFactor factor including detJ and area
//! @param rResult result
void NuTo::Truss::AddDetJNtB(const std::vector<double>& rShapeFunctions,
		const std::vector<double>& rDerivativeShapeFunctions, double rFactor, int rRow, int rCol,
		FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult)
{
    //assert(rShapeFunctions.size()==rDerivativeShapeFunctions.size());
	//these for loops are not optimal, we might use a map to an eigenvector and perform the multiplication directly
	for (unsigned int count=0; count<rShapeFunctions.size(); count++)
	{
		for (unsigned int count2=0; count2<rDerivativeShapeFunctions.size(); count2++)
		{
			rResult(rRow+count,rCol+count2)-=
					rFactor*rShapeFunctions[count]*(rDerivativeShapeFunctions[count2]);
		}
    }
}

//! @brief adds to a matrix the product N^t C N, where N contains the the shape functions and C is the constitutive tangent
//! eventually include also area/width of an element
//! @param rShapeFunctions shape functions
//! @param ConstitutiveTangentBase constitutive tangent matrix
//! @param rFactor factor including area, determinant of Jacobian and IP weight
//! @param rRow row, where to start to add the submatrix
//! @param rCol col, where to start to add the submatrix
//! @param rCoefficientMatrix to be added to
void NuTo::Truss::AddDetJNtCN(const std::vector<double>& rShapeFunctions,
                              const ConstitutiveTangentLocal<1,1>& rConstitutiveTangent, double rFactor,
                              int rRow, int rCol,
                              FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoefficientMatrix)const
{
    rFactor *=rConstitutiveTangent(0,0);
    for (int node1=0; node1<GetNumNodes(); node1++)
    {
        for (int node2=0; node2<GetNumNodes(); node2++)
        {
            rCoefficientMatrix(rRow+node1,rCol+node2)+=rFactor*rShapeFunctions[node1]*rShapeFunctions[node2];
        }
    }
}


//! @brief add Kkk*kappa+detJ*F (detJ is already included in Kkk)
//! @param derivativeShapeFunctions of the ip for all shape functions
//! @param tangentStressNonlocalEqPlasticStrain derivative of the stress with respect to the nonlocal eq plastic strain
//! @param rShapeFunctions of the ip for all shape functions
//! @param rFactor factor including detJ and area
//! @param rResult result
void NuTo::Truss::AddDetJBtdSigmadNonlocalEqPlasticStrainN(const std::vector<double>& rDerivativeShapeFunctions, ConstitutiveTangentLocal<1,2>& rTangentStressNonlocalEqPlasticStrain,
		const std::vector<double>& rShapeFunctions, double rFactor, int rRow, int rCol, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult)
{
    assert(rShapeFunctions.size()==rDerivativeShapeFunctions.size());

    for (int countNonlocalEqPlasticStrain=0; countNonlocalEqPlasticStrain<2; countNonlocalEqPlasticStrain++ )
    {
    	double tmpfactor = rTangentStressNonlocalEqPlasticStrain(countNonlocalEqPlasticStrain)*rFactor;
		//these tow for loops are not optimal, we might use a map to an eigenvector and perform the multiplication directly
		for (unsigned int count=0; count<rDerivativeShapeFunctions.size(); count++)
		{
			for (unsigned int count2=0; count2<rShapeFunctions.size(); count2++)
			{
				rResult(rRow+count,rCol+count2+rShapeFunctions.size()*countNonlocalEqPlasticStrain)+=tmpfactor*rDerivativeShapeFunctions[count]*rShapeFunctions[count2];
			}
		}
    }
}

//! @brief add Kkk*kappa+detJ*F (detJ is already included in Kkk)
//! @param derivativeShapeFunctions of the ip for all shape functions
//! @param tangentStressNonlocalEqPlasticStrain derivative of the stress with respect to the nonlocal eq plastic strain
//! @param rShapeFunctions of the ip for all shape functions
//! @param rFactor factor including detJ and area
//! @param rResult result
void NuTo::Truss::AddDetJBtdSigmadNonlocalTotalStrainN(const std::vector<double>& rDerivativeShapeFunctions, ConstitutiveTangentLocal<1,1>& rTangentStressNonlocalTotalStrain,
		const std::vector<double>& rShapeFunctions, double rFactor, int rRow, int rCol, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult)
{
	double tmpfactor = rTangentStressNonlocalTotalStrain(0,0)*rFactor;
	//these tow for loops are not optimal, we might use a map to an eigenvector and perform the multiplication directly
	for (unsigned int count=0; count<rDerivativeShapeFunctions.size(); count++)
	{
		for (unsigned int count2=0; count2<rShapeFunctions.size(); count2++)
		{
			rResult(rRow+count,rCol+count2)+=tmpfactor*rDerivativeShapeFunctions[count]*rShapeFunctions[count2];
		}
    }
}

//! @brief add detJ B.T dSigma/dnonlocalEqStrain N
//! @param derivativeShapeFunctions of the ip for all shape functions
//! @param tangentStressNonlocalEqStrain derivative of the stress with respect to the nonlocal eq strain
//! @param rShapeFunctions of the ip for all shape functions
//! @param rFactor factor including detJ and area
//! @param rResult result
void NuTo::Truss::AddDetJBtdSigmadNonlocalEqStrainN(const std::vector<double>& rDerivativeShapeFunctions, ConstitutiveTangentLocal<1,1>& rTangentStressNonlocalEqStrain,
        const std::vector<double>& rShapeFunctions, double rFactor, int rRow, int rCol, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult)
{
    double tmpfactor = rTangentStressNonlocalEqStrain(0)*rFactor;
    //these tow for loops are not optimal, we might use a map to an eigenvector and perform the multiplication directly
    for (unsigned int iB=0; iB<rDerivativeShapeFunctions.size(); iB++) // B-Matrix
    {
        for (unsigned int iN=0; iN<rShapeFunctions.size(); iN++)  // N-Matrix
        {
            rResult(rRow+iB,rCol+iN)+=tmpfactor*rDerivativeShapeFunctions[iB]*rShapeFunctions[iN];
        }
    }
}

//! @brief add Kkk*omega+detJ*F (detJ is already included in Kkk)
//! @param rShapeFunctions of the ip for all shape functions
//! @param rLocalEqPlasticStrain local eq. plastic strain values
//! @param rKkk stiffness matrix Kkk
//! @param rNodeNonlocalEqPlasticStrain nodal nonlocal eq plastic strain values
//! @param rFactor factor including detJ and area
//! @param rResult result
void NuTo::Truss::AddDetJRnonlocalPlasticStrain(const std::vector<double>& rShapeFunctions,const LocalEqPlasticStrain& rLocalEqPlasticStrain,
		const FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rKkk,
		const std::vector<double>& rNodeNonlocalEqPlasticStrain, double rFactor, int rRow, FullVector<double,Eigen::Dynamic>& rResult) const
{
	assert(rResult.GetNumRows()>=(int)(rRow+2*rShapeFunctions.size()));
	assert(2*rShapeFunctions.size()==rNodeNonlocalEqPlasticStrain.size());

	//first component
	FullVector<double,Eigen::Dynamic> tmpResult = rKkk*Eigen::Map<Eigen::VectorXd>(const_cast<double*>(&(rNodeNonlocalEqPlasticStrain[0])),rShapeFunctions.size());
	double localEqPlasticStrainMulFactor = rLocalEqPlasticStrain(0)*rFactor;
	for (unsigned int count=0; count<rShapeFunctions.size(); count++)
	{
		tmpResult(count)-=rShapeFunctions[count]*localEqPlasticStrainMulFactor;
	}
	rResult.segment(rRow,rShapeFunctions.size()) +=tmpResult;

	//second component
	tmpResult = rKkk*Eigen::Map<Eigen::VectorXd>(const_cast<double*>(&(rNodeNonlocalEqPlasticStrain[rShapeFunctions.size()])),rShapeFunctions.size());
	localEqPlasticStrainMulFactor = rLocalEqPlasticStrain(1)*rFactor;
	for (unsigned int count=0; count<rShapeFunctions.size(); count++)
	{
		tmpResult(count)-=rShapeFunctions[count]*localEqPlasticStrainMulFactor;
		//std::cout << "tmpResult \n" << tmpResult << std::endl;
	}
	rResult.segment(rRow+rShapeFunctions.size(),rShapeFunctions.size()) +=tmpResult;
}

//! @brief add Kkk*omega+detJ*F (detJ is already included in Kkk)
//! @param rShapeFunctions of the ip for all shape functions
//! @param rLocalTotalStrain local total strain values
//! @param rKkk stiffness matrix Kkk
//! @param rNodeNonlocalTotalStrain nodal nonlocal total strain values
//! @param rFactor factor including detJ and area
//! @param rResult result
void NuTo::Truss::AddDetJRnonlocalTotalStrain(const std::vector<double>& rShapeFunctions,const EngineeringStrain1D& rLocalTotalStrain,
		const FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rKkk,
		const std::vector<double>& rNodeNonlocalTotalStrain, double rFactor, int rRow, FullVector<double,Eigen::Dynamic>& rResult) const
{
	assert(rResult.GetNumRows()>=(int)(rRow+rShapeFunctions.size()));
	assert(rShapeFunctions.size()==rNodeNonlocalTotalStrain.size());

	//first component
	FullVector<double,Eigen::Dynamic> tmpResult = rKkk*Eigen::Map<Eigen::VectorXd>(const_cast<double*>(&(rNodeNonlocalTotalStrain[0])),rShapeFunctions.size());
	double localTotalStrainMulFactor = rLocalTotalStrain(0)*rFactor;
	for (unsigned int count=0; count<rShapeFunctions.size(); count++)
	{
		tmpResult(count)-=rShapeFunctions[count]*localTotalStrainMulFactor;
	}
	rResult.segment(rRow,rShapeFunctions.size()) +=tmpResult;

}

//! @brief add Kkk*nonlocalEqStrain-detJ*N.T*localEqStrain (detJ is already included in Kkk)
//! @param rShapeFunctions of the ip for all shape functions
//! @param rLocalEqStrain local eq. strain values
//! @param rKkk stiffness matrix Kkk
//! @param rNodeNonlocalEqStrain nodal nonlocal eq strain values
//! @param rFactor factor including detJ and area
//! @param rResult result
void NuTo::Truss::AddDetJRnonlocalEqStrain(const std::vector<double>& rShapeFunctions,const LocalEqStrain& rLocalEqStrain, const FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rKkk,
        const std::vector<double>& rNodeNonlocalEqStrain, double rFactor, int rRow, FullVector<double,Eigen::Dynamic>& rResult) const
{
    assert(rResult.GetNumRows()>=(int)(rRow+rShapeFunctions.size()));
    assert(rShapeFunctions.size()==rNodeNonlocalEqStrain.size());

    // perform Kkk * nodeNonlocalEqStrain
    FullVector<double,Eigen::Dynamic> tmpResult = rKkk*Eigen::Map<Eigen::VectorXd>(const_cast<double*>(&(rNodeNonlocalEqStrain[0])),rShapeFunctions.size());

    double localEqStrainMulFactor = rLocalEqStrain(0)*rFactor;
    for (unsigned int count=0; count<rShapeFunctions.size(); count++)
    {
        tmpResult(count) -= rShapeFunctions[count]*localEqStrainMulFactor;
    }
    rResult.segment(rRow,rShapeFunctions.size()) +=tmpResult;
}


//! @brief calculates the Kkk matrix
//! @param shapeFunctions of the ip for all shape functions
//! @param derivativeShapeFunctions of the ip for all shape functions
//! @param c nonlocal gradient radius
//! @param factor multiplication factor (detJ area..)
//! @param Kkk return matrix with detJ * NtT+cBtB
void NuTo::Truss::CalculateKkk(const std::vector<double>& shapeFunctions,const std::vector<double>& derivativeShapeFunctions,double nonlocalGradientRadius,double factor,
		FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& Kkk)
{
	//resize and set to zero
	Kkk.Resize(shapeFunctions.size(),shapeFunctions.size());
    //add NtN
	for (unsigned int count=0; count<shapeFunctions.size();count++)
	{
		for (unsigned int count2=0; count2<shapeFunctions.size();count2++)
		{
			Kkk(count,count2)=factor*(
					      shapeFunctions[count]*shapeFunctions[count2]+
					      nonlocalGradientRadius*derivativeShapeFunctions[count]*derivativeShapeFunctions[count2]);
		}
	}
}

//! @brief stores the temperatures of the nodes
//! @param time derivative (0 temperature, 1 temperature rate, 2 second time derivative of temperature)
//! @param temperature vector with already correct size allocated
//! this can be checked with an assertation
void NuTo::Truss::CalculateNodalTemperatures(int rTimeDerivative, std::vector<double>& rTemperatures)const
{
    assert(rTimeDerivative>=0);
    assert(rTimeDerivative<3);
	assert((int)rTemperatures.size()==GetNumNodesField());
    for (int count=0; count<GetNumNodesField(); count++)
    {
        if (GetNode(count)->GetNumTemperatures()!=1)
            throw MechanicsException("[NuTo::Truss::CalculateNodalTemperatures] Temperature is required as input to the constitutive model, but the node does not have this data.");
        rTemperatures[count] = GetNode(count)->GetTemperature(rTimeDerivative);
    }
}

//! @brief stores the nonlocal eq plastic strain of the nodes
//! @param time derivative (0 damage, 1 damage rate, 2 second time derivative of damage)
//! @param nonlocal eq plastic strain vector with already correct size allocated (2*nodes)
//! this can be checked with an assertation
void NuTo::Truss::CalculateNodalNonlocalEqPlasticStrain(int rTimeDerivative, std::vector<double>& rNodalNonlocalEqPlasticStrain)const
{
    assert(rTimeDerivative>=0);
    assert(rTimeDerivative<3);
	assert((int)rNodalNonlocalEqPlasticStrain.size()==2*GetNumNodesField());
	double nonlocalEqPlasticStrain[2];
    for (int count=0; count<GetNumNodesField(); count++)
    {
        if (GetNode(count)->GetNumNonlocalEqPlasticStrain()!=2)
            throw MechanicsException("[NuTo::Truss::CalculateNodalNonlocalEqPlasticStrain] Damage is required as input to the constitutive model, but the node does not have this data.");
        GetNode(count)->GetNonlocalEqPlasticStrain(nonlocalEqPlasticStrain);
        rNodalNonlocalEqPlasticStrain[count] = nonlocalEqPlasticStrain[0];
        rNodalNonlocalEqPlasticStrain[count+GetNumNodesField()] = nonlocalEqPlasticStrain[1];
    }
}

//! @brief stores the nonlocal total strain of the nodes
//! @param time derivative (0 damage, 1 damage rate, 2 second time derivative of damage)
//! @param nonlocal total strain vector with already correct size allocated (1*nodes)
//! this can be checked with an assertation
void NuTo::Truss::CalculateNodalNonlocalTotalStrain(int rTimeDerivative, std::vector<double>& rNodalNonlocalTotalStrain)const
{
    assert(rTimeDerivative>=0);
    assert(rTimeDerivative<3);
	assert((int)rNodalNonlocalTotalStrain.size()==GetNumShapeFunctionsNonlocalTotalStrain());
	double nonlocalTotalStrain[1];
    for (int count=0; count<GetNumShapeFunctionsNonlocalTotalStrain(); count++)
    {
        if (GetNodeNonlocalTotalStrain(count)->GetNumNonlocalTotalStrain()!=1)
            throw MechanicsException("[NuTo::Truss::CalculateNodalNonlocalTotalStrain] nonlocal strain is required as input to the constitutive model, but the node does not have this data.");
        GetNodeNonlocalTotalStrain(count)->GetNonlocalTotalStrain1D(nonlocalTotalStrain);
        rNodalNonlocalTotalStrain[count] = nonlocalTotalStrain[0];
    }
}

//! @brief stores the nonlocal eq strain of the nodes
//! @param time derivative (0 damage, 1 damage rate, 2 second time derivative of damage)
//! @param nonlocal eq strain vector with already correct size allocated (1*nodes)
//! this can be checked with an assertation
void NuTo::Truss::CalculateNodalNonlocalEqStrain(int rTimeDerivative, std::vector<double>& rNodalNonlocalEqStrain)const
{
    assert(rTimeDerivative>=0);
    assert(rTimeDerivative<3);
    assert((int)rNodalNonlocalEqStrain.size()==GetNumShapeFunctionsNonlocalEqStrain());
    for (int iNode=0; iNode < GetNumShapeFunctionsNonlocalEqStrain(); iNode++)
    {
        if (GetNodeNonlocalEqStrain(iNode)->GetNumNonlocalEqStrain()!=1)
            throw MechanicsException("[NuTo::Truss::CalculateNodalNonlocalEqStrain] Nonlocal eq strain is required as input to the constitutive model, but the node does not have this data.");
        rNodalNonlocalEqStrain[iNode] = GetNodeNonlocalEqStrain(iNode)->GetNonlocalEqStrain(rTimeDerivative);
    }
}

//! @brief stores the relative humidity of the nodes
//! @param time derivative
//! @param relative humidity vector with already correct size allocated (1*nodes)
//! this can be checked with an assertation
void NuTo::Truss::CalculateNodalRelativeHumidity(int rTimeDerivative, std::vector<double>& rNodalRelativeHumidity)const
{
    assert(rTimeDerivative>=0);
    assert(rTimeDerivative<3);
    assert((int)rNodalRelativeHumidity.size()==GetNumShapeFunctionsMoistureTransport());
    for (int iNode=0; iNode < GetNumShapeFunctionsMoistureTransport(); iNode++)
    {
        if (GetNode(iNode)->GetNumRelativeHumidity()!=1)
            throw MechanicsException("[NuTo::Truss::CalculateNodalRelativeHumidity] Relative humidity is required as input to the constitutive model, but the node does not have this data.");
        rNodalRelativeHumidity[iNode] = GetNode(iNode)->GetRelativeHumidity(rTimeDerivative);
    }
}

//! @brief stores the water phase fraction of the nodes
//! @param time derivative
//! @param water phase fraction vector with already correct size allocated (1*nodes)
//! this can be checked with an assertation
void NuTo::Truss::CalculateNodalWaterPhaseFraction(int rTimeDerivative, std::vector<double>& rNodalWaterPhaseFraction)const
{
    assert(rTimeDerivative>=0);
    assert(rTimeDerivative<3);
    assert((int)rNodalWaterPhaseFraction.size()==GetNumShapeFunctionsMoistureTransport());
    for (int iNode=0; iNode < GetNumShapeFunctionsMoistureTransport(); iNode++)
    {
        if (GetNode(iNode)->GetNumWaterPhaseFraction()!=1)
            throw MechanicsException("[NuTo::Truss::CalculateNodalWaterPhaseFraction] Water phase fraction is required as input to the constitutive model, but the node does not have this data.");
        rNodalWaterPhaseFraction[iNode] = GetNode(iNode)->GetWaterPhaseFraction(rTimeDerivative);
    }
}

//! @brief calculates deformation gradient1D
//! @param rRerivativeShapeFunctions derivatives of the shape functions
//! @param rLocalDisp local displacements
//! @param rConstitutiveInput (return value)
void NuTo::Truss::CalculateDeformationGradient(const std::vector<double>& rDerivativeShapeFunctions,
        const std::vector<double>& rLocalDisp,
        DeformationGradient1D& rDeformationGradient)const
{
    assert((int)rLocalDisp.size()==GetNumNodesField() && (int)rDerivativeShapeFunctions.size()==GetNumNodesField());
    rDeformationGradient.mDeformationGradient = 0;

    //normally, the inverse Jacobian should be calculated, but for a truss element, it is sufficient to use the inverse of the Jacobian determinant
    for (int count=0; count<GetNumNodes(); count++)
    {
        rDeformationGradient.mDeformationGradient+=rLocalDisp[count]*rDerivativeShapeFunctions[count];
    }
    rDeformationGradient.mDeformationGradient+=1.;
}

//! @brief returns the number of shape functions for the interpolation of the nonlocal total strains
//! this is required for the calculation of the derivatives of the shape functions
//! whose size is GetLocalDimension*GetNumShapeFunctions
//! @return local dimension
int NuTo::Truss::GetNumShapeFunctionsNonlocalTotalStrain()const
{
    return GetNumNodesGeometry();
}

//! @brief returns the number of shape functions for the interpolation of the nonlocal eq strains
//! this is required for the calculation of the derivatives of the shape functions
//! whose size is GetLocalDimension*GetNumShapeFunctions
//! @return local dimension
int NuTo::Truss::GetNumShapeFunctionsNonlocalEqStrain()const
{
    return GetNumNodesGeometry();
}

//! @brief returns the number of shape functions for the interpolation of the moisture transport
//! this is required for the calculation of the derivatives of the shape functions
//! whose size is GetLocalDimension*GetNumShapeFunctions
//! @return local dimension
int NuTo::Truss::GetNumShapeFunctionsMoistureTransport() const
{
    return GetNumNodesGeometry();
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point
//! @param rCoordinates local coordinates (return value)
void NuTo::Truss::GetLocalIntegrationPointCoordinates(int rIpNum, double& rCoordinates)const
{
    this->mElementData->GetIntegrationType()->GetLocalIntegrationPointCoordinates1D(rIpNum, rCoordinates);
    return;
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point
//! @param rCoordinates local coordinates (return value)
void  NuTo::Truss::GetGlobalIntegrationPointCoordinates(int rIpNum, double rCoordinates[3])const
{
    double naturalCoordinates;
    double nodeCoordinates[3];
    std::vector<double> shapeFunctions(GetNumNodes());
    GetLocalIntegrationPointCoordinates(rIpNum, naturalCoordinates);
    CalculateShapeFunctionsGeometry(naturalCoordinates, shapeFunctions);
    rCoordinates[0] = 0.;
    rCoordinates[1] = 0.;
    rCoordinates[2] = 0.;

    nodeCoordinates[0] = 0;
    nodeCoordinates[1] = 0;
    nodeCoordinates[2] = 0;
    for (int theNode=0; theNode<GetNumNodes(); theNode++)
    {
    	const NodeBase *nodePtr(GetNode(theNode));
    	switch (nodePtr->GetNumCoordinates())
    	{
    	case 1:
    		nodePtr->GetCoordinates1D(nodeCoordinates);
    	break;
    	case 2:
    		nodePtr->GetCoordinates2D(nodeCoordinates);
    	break;
    	case 3:
    		nodePtr->GetCoordinates3D(nodeCoordinates);
    	break;
    	default:
    		throw MechanicsException("[NuTo::Truss::GetGlobalIntegrationPointCoordinates] Node has to have 1, 2 or 3 coordinates.");
    	}

    	for (int theCoordinate=0; theCoordinate<nodePtr->GetNumCoordinates(); theCoordinate++)
    	{
    		rCoordinates[theCoordinate]+=shapeFunctions[theNode]*nodeCoordinates[theCoordinate];
    	}
    }
    return;
}

//! @brief calculates the shape functions, uses the geometry shape functions unless implemented differently
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes, size should already be correct, but can be checked with an assert
void NuTo::Truss::CalculateShapeFunctionsField(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const
{
    CalculateShapeFunctionsGeometry(rLocalCoordinates, rShapeFunctions);
}

//! @brief calculates the shape functions, uses the geometry shape functions unless implemented differently
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Truss::CalculateShapeFunctionsNonlocalTotalStrain(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const
{
    CalculateShapeFunctionsGeometry(rLocalCoordinates, rShapeFunctions);
}

//! @brief calculates the shape functions, uses the geometry shape functions unless implemented differently
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Truss::CalculateShapeFunctionsNonlocalEqStrain(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const
{
    CalculateShapeFunctionsGeometry(rLocalCoordinates, rShapeFunctions);
}

//! @brief calculates the shape functions, uses the geometry shape functions unless implemented differently
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Truss::CalculateShapeFunctionsMoistureTransport(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const
{
    CalculateShapeFunctionsGeometry(rLocalCoordinates,rShapeFunctions);
}

//! @brief calculates the derivative of the shape functions, uses the derivatives of the geometry shape function unless implemented differently
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
//! first all the directions for a single node, and then for the next node
void NuTo::Truss::CalculateDerivativeShapeFunctionsField(const double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const
{
    CalculateDerivativeShapeFunctionsGeometry(rLocalCoordinates, rDerivativeShapeFunctions);
}

//! @brief calculates the derivative of the shape functions, uses the derivatives of the geometry shape function unless implemented differently
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
//! first all the directions for a single node, and then for the next node
void NuTo::Truss::CalculateDerivativeShapeFunctionsNonlocalTotalStrain(const double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const
{
    CalculateDerivativeShapeFunctionsGeometry(rLocalCoordinates, rDerivativeShapeFunctions);
}

//! @brief calculates the derivative of the shape functions, uses the derivatives of the geometry shape function unless implemented differently
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
//! first all the directions for a single node, and then for the next node
void NuTo::Truss::CalculateDerivativeShapeFunctionsNonlocalEqStrain(const double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const
{
    CalculateDerivativeShapeFunctionsGeometry(rLocalCoordinates, rDerivativeShapeFunctions);
}

//! @brief calculates the derivative of the shape functions, uses the derivatives of the geometry shape function unless implemented differently
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes, size should already be correct, but can be checked with an assert
//! first all the directions for a single node, and then for the next node
void NuTo::Truss::CalculateDerivativeShapeFunctionsMoistureTransport(const double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const
{
    CalculateDerivativeShapeFunctionsGeometry(rLocalCoordinates,rDerivativeShapeFunctions);
}


//! @brief returns the nonlocal eq plastic strain interpolated from the nodal values
//! @param shapeFunctionsGlobal shape functions
//! @param rNodeDamage nonlocal eq plastic strain values of the nodes
//! @param rNonlocalEqentPlasticStrain return value
void NuTo::Truss::CalculateNonlocalEqPlasticStrain(const std::vector<double>& shapeFunctions,
		const std::vector<double>& rNodeEquivalentPlasticStrain, NonlocalEqPlasticStrain& rNonlocalEqentPlasticStrain)const
{
    assert(2*shapeFunctions.size()==rNodeEquivalentPlasticStrain.size());
    rNonlocalEqentPlasticStrain(0) = 0.;
    rNonlocalEqentPlasticStrain(1) = 0.;
	for (unsigned int count=0; count<shapeFunctions.size(); count++)
	{
		rNonlocalEqentPlasticStrain(0)+=shapeFunctions[count]*rNodeEquivalentPlasticStrain[count];
		rNonlocalEqentPlasticStrain(1)+=shapeFunctions[count]*rNodeEquivalentPlasticStrain[count+shapeFunctions.size()];
		//std::cout << "count " << count << " damage " << rNodeDamage[count] << " shape " << shapeFunctions[count] << std::endl;;
	}
	return;
}

//! @brief returns the nonlocal total strain interpolated from the nodal values
//! @param shapeFunctionsGlobal shape functions
//! @param rNodeNonlocalTotalStrain nonlocal total strain values of the nodes
//! @param rNonlocalTotalStrain return value
void NuTo::Truss::CalculateNonlocalTotalStrain(const std::vector<double>& shapeFunctions,
		const std::vector<double>& rNodeNonlocalTotalStrain, EngineeringStrain1D& rNonlocalTotalStrain)const
{
    assert(shapeFunctions.size()==rNodeNonlocalTotalStrain.size());
    rNonlocalTotalStrain(0) = 0.;
	for (unsigned int count=0; count<shapeFunctions.size(); count++)
	{
		rNonlocalTotalStrain(0)+=shapeFunctions[count]*rNodeNonlocalTotalStrain[count];
	}
	return;
}

//! @brief returns the nonlocal eq strain interpolated from the nodal values
//! @param shapeFunctions shape functions
//! @param rNodeNonlocalEqStrain nonlocal eq strain values of the nodes
//! @param rNonlocalEqStrain return value
void NuTo::Truss::CalculateNonlocalEqStrain(const std::vector<double>& shapeFunctions,
        const std::vector<double>& rNodeNonlocalEqStrain, NonlocalEqStrain& rNonlocalEqStrain)const
{
    assert(shapeFunctions.size()==rNodeNonlocalEqStrain.size());
    rNonlocalEqStrain(0) = 0.;

    for (unsigned int count=0; count<shapeFunctions.size(); count++)
    {
        rNonlocalEqStrain(0)+=shapeFunctions[count]*rNodeNonlocalEqStrain[count];
    }
}

//! @brief returns the relative humidity interpolated from the nodal values
//! @param shapeFunctions shape functions
//! @param rNodeRelativeHumidity relative humidity values of the nodes
//! @param rRelativeHumidity return value
void NuTo::Truss::CalculateRelativeHumidity(const std::vector<double>& rShapeFunctions,
                                            const std::vector<double>& rNodeRelativeHumidity,
                                            RelativeHumidity& rRelativeHumidity)const
{
    assert(rShapeFunctions.size()==rNodeRelativeHumidity.size());
    rRelativeHumidity(0) = 0.0;
    for (unsigned int count=0; count<rShapeFunctions.size(); count++)
    {
        rRelativeHumidity(0)+=rShapeFunctions[count]*rNodeRelativeHumidity[count];
    }
}

//! @brief returns the water phase fraction interpolated from the nodal values
//! @param shapeFunctions shape functions
//! @param rNodeWaterPhaseFraction water phase fraction values of the nodes
//! @param rWaterPhaseFraction return value
void NuTo::Truss::CalculateWaterPhaseFraction(const std::vector<double>& rShapeFunctions,
                                              const std::vector<double>& rNodeWaterPhaseFraction,
                                              WaterPhaseFraction& rWaterPhaseFraction)const
{
    assert(rShapeFunctions.size()==rNodeWaterPhaseFraction.size());
    rWaterPhaseFraction(0)=0.0;
    for (unsigned int count = 0; count < rShapeFunctions.size(); count++)
    {
        rWaterPhaseFraction(0)+=rShapeFunctions[count]*rNodeWaterPhaseFraction[count];
    }
}


//! @brief calculates the integration point data with the current displacements applied
//! @param rIpDataType data type to be stored for each integration point
//! @param rIpData return value with dimension (dim of data type) x (numIp)
/*NuTo::Error::eError NuTo::Truss::GetIpData(NuTo::IpData::eIpStaticDataType rIpDataType, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rIpData)const
{
	//calculate local coordinates
    std::vector<double> localNodeCoord(GetNumLocalDofs());
    CalculateLocalCoordinates(localNodeCoord);

    //calculate local displacements
    std::vector<double> localNodeDisp(GetNumLocalDofs());
    CalculateLocalDisplacements(localNodeDisp);

    //allocate space for local ip coordinates
    double localIPCoord;

    //allocate space for local shape functions
    std::vector<double> derivativeShapeFunctions(GetLocalDimension()*GetNumShapeFunctions());

    //allocate deformation gradient
    DeformationGradient1D deformationGradient;

    //allocate global engineering strain
    EngineeringStrain3D engineeringStrain;

    //allocate global engineering stress
	EngineeringStress3D engineeringStress;

    //material pointer
    const ConstitutiveEngineeringStressStrain *constitutivePtr;

    //allocate and initialize result matrix
    switch (rIpDataType)
    {
    case NuTo::IpData::ENGINEERING_STRAIN:
    case NuTo::IpData::ENGINEERING_STRESS:
    case NuTo::IpData::ENGINEERING_PLASTIC_STRAIN:
       	rIpData.Resize(6,GetNumIntegrationPoints());
    break;
    case NuTo::IpData::DAMAGE:
       	rIpData.Resize(1,GetNumIntegrationPoints());
    break;
    default:
    	throw MechanicsException("[NuTo::Plane::GetIpData] Ip data not implemented.");
    }

    //store the data
    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, localIPCoord);

        CalculateDerivativeShapeFunctions(localIPCoord, derivativeShapeFunctions);

        // determine deformation gradient from the local Displacements and the derivative of the shape functions
        CalculateDeformationGradient(derivativeShapeFunctions, localNodeCoord, localNodeDisp, deformationGradient);

        //call material law to calculate engineering strain
        constitutivePtr = GetConstitutiveLaw(theIP)->AsConstitutiveEngineeringStressStrain();

        Error::eError error;
        switch (rIpDataType)
        {
        case NuTo::IpData::ENGINEERING_STRAIN:
            error = constitutivePtr->GetEngineeringStrain(this, theIP, deformationGradient, engineeringStrain);
            memcpy(&(rIpData.mEigenMatrix.data()[theIP*6]),engineeringStrain.GetData(),6*sizeof(double));
        break;
        case NuTo::IpData::ENGINEERING_STRESS:
        	error = constitutivePtr->GetEngineeringStressFromEngineeringStrain(this, theIP, deformationGradient, engineeringStress);
            memcpy(&(rIpData.mEigenMatrix.data()[theIP*6]),engineeringStress.GetData(),6*sizeof(double));
        break;
        case NuTo::IpData::ENGINEERING_PLASTIC_STRAIN:
        	error = constitutivePtr->GetEngineeringPlasticStrain(this, theIP, deformationGradient, engineeringStrain);
            memcpy(&(rIpData.mEigenMatrix.data()[theIP*6]),engineeringStrain.GetData(),6*sizeof(double));
        break;
        case NuTo::IpData::DAMAGE:
        	error = constitutivePtr->GetDamage(this, theIP, deformationGradient, rIpData.mEigenMatrix.data()[theIP]);
        break;
        default:
        	throw MechanicsException("[NuTo::Plane::GetIpData] Ip data not implemented.");
        }
        if (error!=Error::SUCCESSFUL)
        	return error;
    }

    return Error::SUCCESSFUL;
}
*/

//! @brief sets the section of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need a section
//! @param rSection pointer to section
//! @return pointer to constitutive law
void NuTo::Truss::SetSection(const SectionBase* rSection)
{
    mSection = rSection;
}

//! @brief returns a pointer to the section of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need a section
//! @return pointer to section
const NuTo::SectionBase* NuTo::Truss::GetSection()const
{
    return mSection;
}

//! @brief Allocates static data for an integration point of an element
//! @param rConstitutiveLaw constitutive law, which is called to allocate the static data object
//! actually, both - the element type and the constitutive law are required to determine the static data object actually required
NuTo::ConstitutiveStaticDataBase* NuTo::Truss::AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw)const
{
    return rConstitutiveLaw->AllocateStaticDataEngineeringStress_EngineeringStrain1D(this);
}

//! @brief returns determinant of the Jacobian
//! @param derivativeShapeFunctions derivatives of the shape functions
//! @param localCoord local coordinates
//! @return determinant of the Jacobian
double NuTo::Truss::DetJacobian(const std::vector<double>& derivativeShapeFunctions,const std::vector<double>& localCoord)const
{
    assert((int)localCoord.size()==GetNumNodes() && (int)derivativeShapeFunctions.size()==GetNumNodes());
    double detJ(0);
    for (int count=0; count<GetNumNodes(); count++)
        detJ+=derivativeShapeFunctions[count]*localCoord[count];
    return detJ;
}

//! @brief adds to a matrix the product B^tCB, where B contains the derivatives of the shape functions and C is the constitutive tangent
//! eventually include also area/width of an element
//! @param rDerivativeShapeFunctions derivatives of the shape functions
//! @param ConstitutiveTangentBase constitutive tangent matrix
//! @param rFactor factor including area, determinant of Jacobian and IP weight
//! @param rRow row, where to start to add the submatrix
//! @param rCol col, where to start to add the submatrix
//! @param rCoefficientMatrix to be added to
void NuTo::Truss::AddDetJBtCB(const std::vector<double>& rDerivativeShapeFunctions,
                              const ConstitutiveTangentLocal<1,1>& rConstitutiveTangent, double rFactor,
                              int rRow, int rCol,
                              FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoefficientMatrix)const
{
    rFactor *=rConstitutiveTangent(0,0);
    for (int node1=0; node1<GetNumNodes(); node1++)
    {
        for (int node2=0; node2<GetNumNodes(); node2++)
        {
            rCoefficientMatrix(rRow+node1,rCol+node2)+=rFactor*rDerivativeShapeFunctions[node1]*rDerivativeShapeFunctions[node2];
        }
    }
}

//! @brief adds to a matrix the product factor * H^tH, where H contains the shape functions
//! @param rShapeFunctions ... shape functions
//! @param rFactor factor including area, determinant of Jacobian, IP weight and, eventually, the density
//! @param rCoefficientMatrix to be added to
void NuTo::Truss::AddDetJHtH(const std::vector<double>& rShapeFunctions, double rFactor, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoefficientMatrix)const
{
    for (int node1=0; node1<GetNumNodes(); node1++)
    {
        for (int node2=0; node2<GetNumNodes(); node2++)
        {
            rCoefficientMatrix(node1,node2)+=rFactor*rShapeFunctions[node1]*rShapeFunctions[node2];
        }
    }
}


//! @brief adds up the internal force vector
//! @param derivativeShapeFunctions derivatives of the shape functions
//! @param rEngineeringStress stress
//! @param factor factor including det Jacobian area and integration point weight
//! @param rRow start row in case of a multifield problem
//! @param rResult resforce vector
void NuTo::Truss::AddDetJBtSigma(const std::vector<double>& rDerivativeShapeFunctions,
                                 const EngineeringStress1D& rEngineeringStress, double rFactor,
                                 int rRow,
                                 FullVector<double,Eigen::Dynamic>& rResult)const
{
    rFactor*=rEngineeringStress.GetData()[0];
    for (int node1=0; node1<GetNumNodes(); node1++)
    {
        rResult(rRow+node1)+=rFactor*rDerivativeShapeFunctions[node1];
    }
}

//! @brief adds up the internal force vector
//! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param rHeatFlux stress
//! @param factor factor including det Jacobian area and integration point weight
//! @param rRow start row (in case of a multifield problem)
//! @param rResult resforce vector
void NuTo::Truss::AddDetJBtHeatFlux(const std::vector<double>& rDerivativeShapeFunctions,
                                 const HeatFlux1D& rHeatFlux,
                                 double rFactor,
                                 int rRow,
                                 FullVector<double,Eigen::Dynamic>& rResult)const
{
    rFactor*=rHeatFlux.GetData()[0];
    for (int node1=0; node1<GetNumNodes(); node1++)
    {
        rResult(rRow+node1)+=rFactor*rDerivativeShapeFunctions[node1];
    }
}

//! @brief calculates the volume of an integration point (weight * detJac)
//! @param rVolume  vector for storage of the ip volumes (area in 2D, length in 1D)
void NuTo::Truss::GetIntegrationPointVolume(std::vector<double>& rVolume)const
{
    //calculate local coordinates
    std::vector<double> localNodeCoord(GetNumNodesGeometry());
    CalculateLocalCoordinates(localNodeCoord);

    //allocate space for local ip coordinates
    double localIPCoord;

    //allocate space for local shape functions
    std::vector<double> derivativeShapeFunctions(GetNumNodesGeometry());

	rVolume.resize(GetNumIntegrationPoints());

    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, localIPCoord);

        CalculateDerivativeShapeFunctionsGeometry(localIPCoord, derivativeShapeFunctions);

		//attention in 1D, this is just the length, but that is required for the nonlocal model
		rVolume[theIP] = DetJacobian(derivativeShapeFunctions,localNodeCoord)
                       *(mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP));
    }
}

//! @brief returns a pointer to the i-th node of the element
//! @param local node number
//! @return pointer to the node
NuTo::NodeBase* NuTo::Truss::GetNodeNonlocalTotalStrain(int rLocalNodeNumber)
{
    return GetNode(rLocalNodeNumber);
}

//! @brief returns a pointer to the i-th node of the element
//! @param local node number
//! @return pointer to the node
const  NuTo::NodeBase* NuTo::Truss::GetNodeNonlocalTotalStrain(int rLocalNodeNumber)const
{
    return GetNode(rLocalNodeNumber);
}

//! @brief returns a pointer to the i-th node of the element
//! @param local node number
//! @return pointer to the node
NuTo::NodeBase* NuTo::Truss::GetNodeNonlocalEqStrain(int rLocalNodeNumber)
{
    return GetNode(rLocalNodeNumber);
}

//! @brief returns a pointer to the i-th node of the element
//! @param local node number
//! @return pointer to the node
const  NuTo::NodeBase* NuTo::Truss::GetNodeNonlocalEqStrain(int rLocalNodeNumber)const
{
    return GetNode(rLocalNodeNumber);
}


//! @brief cast the base pointer to a Truss, otherwise throws an exception
const NuTo::Truss* NuTo::Truss::AsTruss()const
{
	return this;
}

//! @brief cast the base pointer to a Truss, otherwise throws an exception
NuTo::Truss* NuTo::Truss::AsTruss()
{
	return this;
}

//! @brief ... extract global dofs from nodes (mapping of local column ordering of the element matrices to the global dof ordering)
//! @param rGlobalColumnDofs ... vector of global column dofs
//! @param rNumXxxxDofs ... number of Xxxx dofs
void NuTo::Truss::CalculateGlobalColumnDofs(std::vector<int>& rGlobalColumnDofs,int rNumDispDofs, int rNumTempDofs, int rNumNonlocalEqPlasticStrainDofs,int rNumNonlocalTotalStrainDofs, int rNumNonlocalEqStrainDofs, int rNumRelativeHumidityDofs, int rNumWaterPhaseFractionDofs) const
{
    CalculateGlobalRowDofs(
            rGlobalColumnDofs,
            rNumDispDofs,
            rNumTempDofs,
            rNumNonlocalEqPlasticStrainDofs,
            rNumNonlocalTotalStrainDofs,
            rNumNonlocalEqStrainDofs,
            rNumRelativeHumidityDofs,
            rNumWaterPhaseFractionDofs);
}


#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::Truss::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Truss::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Truss::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Truss::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Truss::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Truss::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Truss::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Truss" << std::endl;
#endif
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementBase)
           & BOOST_SERIALIZATION_NVP(mSection);
    }
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Truss" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Truss)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::Truss)
#endif // ENABLE_SERIALIZATION
