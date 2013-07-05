// $Id: Truss.cpp 627 2013-05-22 07:43:22Z unger3 $

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
#include "nuto/mechanics/constitutive/mechanics/NonlocalEqPlasticStrain.h"
#include "nuto/mechanics/constitutive/thermal/HeatFlux1D.h"
#include "nuto/mechanics/constitutive/thermal/HeatFlux3D.h"
#include "nuto/mechanics/constitutive/thermal/Temperature.h"
#include "nuto/mechanics/constitutive/thermal/TemperatureGradient1D.h"
#include "nuto/mechanics/constitutive/thermal/TemperatureGradient3D.h"
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixInt.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixDouble.h"
#include "nuto/mechanics/elements/ElementOutputVectorInt.h"
#include "nuto/mechanics/elements/BoundaryGradientDamage1D.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/structures/StructureBase.h"

//! @brief constructor
NuTo::BoundaryGradientDamage1D::BoundaryGradientDamage1D(const StructureBase* rStructure,
    		std::vector<NuTo::NodeBase* >& rNodes,
    		Truss* rRealBoundaryElement,
    		int rEdgeRealBoundaryElement,
    		ElementData::eElementDataType rElementDataType,
    		IntegrationType::eIntegrationType rIntegrationType,
    		IpData::eIpDataType rIpDataType
    		)
{

}

//! @brief calculates output data fo the elmement
//! @param eOutput ... coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
//!                    @param updateStaticData (with DummyOutput), IPData, globalrow/column dofs etc.
NuTo::Error::eError NuTo::BoundaryGradientDamage1D::Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput)
{
/*

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
		int numCoordinates(GetNumShapeFunctions());
		std::vector<double> localNodeCoord(numCoordinates);
		CalculateLocalCoordinates(localNodeCoord);

		//calculate local displacements, velocities and accelerations
		//the difference between disp and dispdof is a problem where the displacements are fixed, but enter the constitutive equation
		//for example in a two stage problem, first solve mechanics, then thermal and so on
		int numDisp(GetNumShapeFunctions());
		int numDispDofs = section->GetIsDisplacementDof() ? numDisp : 0;

		int numTemp(GetNumShapeFunctions());
		int numTempDofs(section->GetIsTemperatureDof() ? numTemp : 0);

		int numNonlocalEqPlasticStrain(2*GetNumShapeFunctions());
		int numNonlocalEqPlasticStrainDofs(section->GetIsNonlocalEqPlasticStrainDof() ? numNonlocalEqPlasticStrain : 0);

		int numNonlocalTotalStrain(GetNumShapeFunctions());
		int numNonlocalTotalStrainDofs(section->GetIsNonlocalTotalStrainDof() ? numNonlocalTotalStrain : 0);

		std::vector<double> localNodeDisp,nodeTemp,nodeNonlocalEqPlasticStrain,nodeNonlocalTotalStrain;

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

		//allocate space for local ip coordinates
		double localIPCoord;

		//allocate space for local shape functions
		std::vector<double> derivativeShapeFunctionsNatural(GetLocalDimension()*GetNumShapeFunctions());  //allocate space for derivatives of shape functions
		std::vector<double> derivativeShapeFunctionsLocal(GetLocalDimension()*GetNumShapeFunctions());    //allocate space for derivatives of shape functions
		std::vector<double> shapeFunctions(GetNumShapeFunctions());                                       //allocate space for derivatives of shape functions

		//allocate deformation gradient
		DeformationGradient1D deformationGradient;

		EngineeringStrain3D engineeringStrain3D;

		//allocate global engineering plastic strain
		EngineeringStrain3D engineeringPlasticStrain3D;

		//allocate  damage (output of constitutive relation)
		Damage damage;

		//allocate nonlocal eq plastic strain (nodal dof value, input of constitutive relation)
		NonlocalEqPlasticStrain nonLocalEqPlasticStrain;

		//allocate local eq plastic strain output of constitutive relation
		LocalEqPlasticStrain localEqPlasticStrain;

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

		//allocate tangents
		ConstitutiveTangentLocal<1,1> tangentStressStrain;
		ConstitutiveTangentLocal<1,1> tangentStressTemperature;
		ConstitutiveTangentLocal<1,2> tangentStressNonlocalEqPlasticStrain;
		ConstitutiveTangentLocal<1,1> tangentStressNonlocalTotalStrain;
		ConstitutiveTangentLocal<1,1> tangentHeatFluxTemperatureGradient;
		ConstitutiveTangentLocal<1,1> tangentHeatFluxTemperatureRate;
		ConstitutiveTangentLocal<2,1> tangentLocalEqPlasticStrainStrain;

		//define inputs and outputs
		std::map< NuTo::Constitutive::eInput, const ConstitutiveInputBase* > constitutiveInputList;
		std::map< NuTo::Constitutive::eOutput, ConstitutiveOutputBase* > constitutiveOutputList;

		if (numDispDofs>0)
		{
			constitutiveInputList[NuTo::Constitutive::eInput::DEFORMATION_GRADIENT_1D] = &deformationGradient;
		}

		if (numTempDofs>0)
		{
			constitutiveInputList[NuTo::Constitutive::eInput::TEMPERATURE_GRADIENT_1D] = &temperatureGradient1D;
		}

		if (section->GetInputConstitutiveIsTemperature())
		{
			constitutiveInputList[NuTo::Constitutive::eInput::TEMPERATURE] = &temperature;
		}

		if (numNonlocalEqPlasticStrainDofs>0)
		{
			constitutiveInputList[NuTo::Constitutive::eInput::NONLOCAL_EQ_PLASTIC_STRAIN] = &nonLocalEqPlasticStrain;
		}

		if (numNonlocalTotalStrainDofs>0)
		{
			constitutiveInputList[NuTo::Constitutive::eInput::NONLOCAL_TOTAL_STRAIN_1D] = &nonlocalTotalStrain;
		}

		//define outputs
		for (auto it = rElementOutput.begin(); it!=rElementOutput.end(); it++)
		{
			switch(it->first)
			{
			case Element::INTERNAL_GRADIENT:
				it->second->GetFullVectorDouble().Resize(numDispDofs+numTempDofs+numNonlocalEqPlasticStrainDofs+numNonlocalTotalStrainDofs);
				//if the stiffness matrix is constant, the corresponding internal force is calculated via the Kd
				//on the global level
				if (mStructure->GetHessianConstant(0)==false)
				{
					if (numDispDofs>0)
					{
						constitutiveOutputList[NuTo::Constitutive::eOutput::ENGINEERING_STRESS_1D] = &engineeringStress1D;
					}
					if (numTempDofs>0)
					{
						constitutiveOutputList[NuTo::Constitutive::eOutput::HEAT_FLUX_1D] = &heatFlux1D;
					}
					if (numNonlocalEqPlasticStrainDofs>0)
					{
						constitutiveOutputList[NuTo::Constitutive::eOutput::LOCAL_EQ_PLASTIC_STRAIN] = &localEqPlasticStrain;
					}
					if (numNonlocalTotalStrainDofs>0)
					{
						constitutiveOutputList[NuTo::Constitutive::eOutput::ENGINEERING_STRAIN_1D] = &localTotalStrain;
					}
				}
			break;
			case Element::HESSIAN_0_TIME_DERIVATIVE:
				{
					it->second->GetFullMatrixDouble().Resize(numDispDofs+numTempDofs+numNonlocalEqPlasticStrainDofs+numNonlocalTotalStrainDofs,numDispDofs+numTempDofs+numNonlocalEqPlasticStrainDofs+numNonlocalTotalStrainDofs);
					it->second->SetSymmetry(true);
					it->second->SetConstant(true);
					if (numDispDofs>0)
					{
						constitutiveOutputList[NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_1D] = &tangentStressStrain;
						//mixed terms
						if (numTempDofs>0)
							constitutiveOutputList[NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_TEMPERATURE_1D] = &tangentStressTemperature;
						if (numNonlocalEqPlasticStrainDofs>0)
							constitutiveOutputList[NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_PLASTIC_STRAIN_1D] = &tangentStressNonlocalEqPlasticStrain;
						if (numNonlocalTotalStrainDofs>0)
							constitutiveOutputList[NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_NONLOCAL_TOTAL_STRAIN_1D] = &tangentStressNonlocalTotalStrain;
					}
					if (numTempDofs>0)
					{
						constitutiveOutputList[NuTo::Constitutive::eOutput::D_HEAT_FLUX_D_TEMPERATURE_RATE_1D] = &tangentHeatFluxTemperatureGradient;
						//mixed terms
						//if (numDisp)
						//    constitutiveOutputList.insert(std::pair<NuTo::Constitutive::eOutput, ConstitutiveOutputBase*>(NuTo::Constitutive::eOutput::D_HEAT_FLUX_D_ENGINEERING_STRAIN_3D, &tangentHeatFluxEngineeringStrain[timeDerivative]));
					}
					if (numNonlocalEqPlasticStrainDofs>0)
					{
						//mixed terms
						if (numDispDofs)
							constitutiveOutputList[NuTo::Constitutive::eOutput::D_LOCAL_EQ_PLASTIC_STRAIN_D_STRAIN_1D] = &tangentLocalEqPlasticStrainStrain;
					}
				}
			break;
			case Element::HESSIAN_1_TIME_DERIVATIVE:
			{
				it->second->GetFullMatrixDouble().Resize(numDispDofs+numTempDofs+numNonlocalEqPlasticStrainDofs+numNonlocalTotalStrainDofs,numDispDofs+numTempDofs+numNonlocalEqPlasticStrainDofs+numNonlocalTotalStrainDofs);
				it->second->SetSymmetry(true);
				it->second->SetConstant(true);
				if (numDispDofs>0)
				{
					// Rayleigh damping should be introduced on the global level
				}
				if (numTempDofs>0)
				{
					constitutiveOutputList[NuTo::Constitutive::eOutput::D_HEAT_FLUX_D_TEMPERATURE_RATE_1D] = &tangentHeatFluxTemperatureRate;
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
				constitutiveOutputList[NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA] = 0;
			break;
			case Element::UPDATE_TMP_STATIC_DATA:
				constitutiveOutputList[NuTo::Constitutive::eOutput::UPDATE_TMP_STATIC_DATA] = 0;
			break;
			case Element::IP_DATA:
				switch(it->second->GetIpDataType())
				{
				case NuTo::IpData::ENGINEERING_STRAIN:
					it->second->GetFullMatrixDouble().Resize(6,GetNumIntegrationPoints());
					 //define outputs
					constitutiveOutputList[NuTo::Constitutive::eOutput::ENGINEERING_STRAIN_3D] = &engineeringStrain3D;
				break;
				case NuTo::IpData::ENGINEERING_STRESS:
					it->second->GetFullMatrixDouble().Resize(6,GetNumIntegrationPoints());
					 //define outputs
					 constitutiveOutputList[NuTo::Constitutive::eOutput::ENGINEERING_STRESS_3D] = &engineeringStress3D;
				break;
				case NuTo::IpData::ENGINEERING_PLASTIC_STRAIN:
					it->second->GetFullMatrixDouble().Resize(6,GetNumIntegrationPoints());
					 //define outputs
					 constitutiveOutputList[NuTo::Constitutive::eOutput::ENGINEERING_PLASTIC_STRAIN_3D] = &engineeringPlasticStrain3D;
					 break;
				break;
				case NuTo::IpData::DAMAGE:
					it->second->GetFullMatrixDouble().Resize(1,GetNumIntegrationPoints());
					//define outputs
					constitutiveOutputList[NuTo::Constitutive::eOutput::DAMAGE] = &damage;
				break;
				default:
					throw MechanicsException("[NuTo::Truss::Evaluate] this ip data type is not implemented.");
				}
			break;
			case Element::GLOBAL_ROW_DOF:
				this->CalculateGlobalRowDofs(it->second->GetVectorInt(),numDispDofs,numTempDofs,numNonlocalEqPlasticStrainDofs,numNonlocalTotalStrainDofs);
			break;
			case Element::GLOBAL_COLUMN_DOF:
				this->CalculateGlobalColumnDofs(it->second->GetVectorInt(),numDispDofs,numTempDofs,numNonlocalEqPlasticStrainDofs,numNonlocalTotalStrainDofs);
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
			CalculateDerivativeShapeFunctions(localIPCoord, derivativeShapeFunctionsNatural);

			//determinant of the Jacobian
			double detJ(this->DetJacobian(derivativeShapeFunctionsNatural,localNodeCoord));

			//derivative in local coordinate system
			for (unsigned int count=0; count<derivativeShapeFunctionsNatural.size(); count++)
			{
				derivativeShapeFunctionsLocal[count] = derivativeShapeFunctionsNatural[count]/detJ;
			}

			if (numDispDofs)
			{
				// determine deformation gradient from the local Displacements and the derivative of the shape functions
				CalculateDeformationGradient(derivativeShapeFunctionsLocal, localNodeCoord, localNodeDisp, deformationGradient);
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
				CalculateShapeFunctions(localIPCoord, shapeFunctions);
				CalculateNonlocalEqPlasticStrain(shapeFunctions, nodeNonlocalEqPlasticStrain, nonLocalEqPlasticStrain);
			}

			if (numNonlocalTotalStrainDofs)
			{
				CalculateShapeFunctions(localIPCoord, shapeFunctions);
				CalculateNonlocalTotalStrain(shapeFunctions, nodeNonlocalTotalStrain, nonlocalTotalStrain);
			}

			ConstitutiveBase* constitutivePtr = GetConstitutiveLaw(theIP);
			Error::eError error = constitutivePtr->Evaluate1D(this, theIP,
					constitutiveInputList, constitutiveOutputList);
			if (error!=Error::SUCCESSFUL)
				return error;

			double factor (detJ*mSection->GetArea()*
					       (mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP)));

			FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> Kkk;
			if (numNonlocalEqPlasticStrainDofs>0 || numNonlocalTotalStrainDofs>0)
			{
				//calculate Kkk detJ*(cBtB+NtN)
				//the nonlocal radius is in a gradient formulation is different from the nonlocal radius in an integral formulation
				CalculateKkk(shapeFunctions,derivativeShapeFunctionsLocal,constitutivePtr->GetNonlocalRadius(),factor,Kkk);
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
							AddDetJBtSigma(derivativeShapeFunctionsLocal,engineeringStress1D, factor, 0, it->second->GetFullVectorDouble());
						}
						if (numTempDofs>0)
						{
							throw MechanicsException("[Nuto::Truss::Evaluate] temperature not yet implemented.");
		//TODO				    AddDetJBtCB(derivativeShapeFunctions,tangentStressStrain, factor, 0, it->second->GetFullMatrixDouble());
						}
						if (numNonlocalEqPlasticStrainDofs>0)
						{
							//add Kkk*nonlocalEqPlasticStrain+detJ*F
							AddDetJRnonlocalplasticStrain(shapeFunctions,localEqPlasticStrain, Kkk, nodeNonlocalEqPlasticStrain, factor, numDispDofs+numTempDofs, it->second->GetFullVectorDouble());
						}
						if (numNonlocalTotalStrainDofs>0)
						{
							//add Kkk*nonlocalTotalStrain+detJ*F
							AddDetJRnonlocalTotalStrain(shapeFunctions,localTotalStrain, Kkk, nodeNonlocalTotalStrain, factor, numDispDofs+numTempDofs+numNonlocalEqPlasticStrainDofs, it->second->GetFullVectorDouble());
						}
					}
				}
				break;
				case Element::HESSIAN_0_TIME_DERIVATIVE:
					{
						if (numDispDofs>0)
						{
							//derivative of F(sigma) with respect to all unknowns
							AddDetJBtCB(derivativeShapeFunctionsLocal, tangentStressStrain, factor, 0,0, it->second->GetFullMatrixDouble());
							if (tangentStressStrain.GetSymmetry()==false)
								it->second->SetSymmetry(false);
							if (tangentStressStrain.GetConstant()==false)
								it->second->SetConstant(false);
							if (numTempDofs>0)
								throw MechanicsException("[NuTo::Truss::Evaluate] mixed terms not yet implemented.");
							if (numNonlocalEqPlasticStrainDofs>0)
							{
								AddDetJBtdSigmadNonlocalEqPlasticStrainN(derivativeShapeFunctionsLocal, tangentStressNonlocalEqPlasticStrain,shapeFunctions, factor, 0,numDispDofs+numTempDofs, it->second->GetFullMatrixDouble());
							}
							if (numNonlocalTotalStrainDofs>0)
							{
								AddDetJBtdSigmadNonlocalTotalStrainN(derivativeShapeFunctionsLocal, tangentStressNonlocalTotalStrain,shapeFunctions, factor, 0,numDispDofs+numTempDofs+numNonlocalEqPlasticStrainDofs, it->second->GetFullMatrixDouble());
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
								AddDetJNtdLocalEqPlasticStraindEpsilonB(shapeFunctions,tangentLocalEqPlasticStrainStrain,derivativeShapeFunctionsLocal, factor, numDispDofs+numTempDofs, 0, it->second->GetFullMatrixDouble());
							}
							it->second->GetFullMatrixDouble().AddBlock(numDispDofs+numTempDofs,numDispDofs+numTempDofs,Kkk);
							it->second->GetFullMatrixDouble().AddBlock(numDispDofs+numTempDofs+this->GetNumNodes(),numDispDofs+numTempDofs+this->GetNumNodes(),Kkk);
						}
						if (numNonlocalTotalStrainDofs>0)
						{
							//derivative of F(nonlocaltotalStrain) with respect to all unknowns
							if (numDispDofs>0)
							{
								AddDetJNtB(shapeFunctions,derivativeShapeFunctionsLocal, factor, numDispDofs+numTempDofs+numNonlocalEqPlasticStrainDofs, 0, it->second->GetFullMatrixDouble());
							}
							it->second->GetFullMatrixDouble().AddBlock(numDispDofs+numTempDofs+numNonlocalEqPlasticStrainDofs,numDispDofs+numTempDofs+numNonlocalEqPlasticStrainDofs,Kkk);
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
					}
					if (numNonlocalEqPlasticStrainDofs>0)
					{
						//no damping term, do Rayleigh damping on the global level
					}
					BlowLocalMatrixToGlobal(it->second->GetFullMatrixDouble());
				}
				break;
				case Element::HESSIAN_2_TIME_DERIVATIVE:
				{
					if (numDispDofs>0)
					{
						double density = constitutivePtr->GetDensity();
						this->CalculateShapeFunctions(localIPCoord, shapeFunctions);

						// calculate local mass matrix
						// don't forget to include determinant of the Jacobian and area
						// detJ * area * density * HtH, :
						double factor2 (density * factor);
						this->AddDetJHtH(shapeFunctions, factor2, it->second->GetFullMatrixDouble());

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
    catch (NuTo::MechanicsException e)
    {
        std::stringstream ss;
        ss << mStructure->ElementGetId(this);
    	e.AddMessage("[NuTo::Truss::Evaluate] Error evaluating element data of element"	+ ss.str() + ".");
        throw e;
    }
*/
    return Error::SUCCESSFUL;
}

//! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
//! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
void NuTo::BoundaryGradientDamage1D::ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
{
    for (int count=0; count<(int)mNodes.size(); count++)
    {
        if (this->mNodes[count]==rOldPtr)
        {
            this->mNodes[count]=rNewPtr;
            break;
        }
    }
}


//! @brief cast the base pointer to a BoundaryGradientDamage1D, otherwise throws an exception
const NuTo::BoundaryGradientDamage1D* NuTo::BoundaryGradientDamage1D::AsBoundaryGradientDamage1D()const
{
	return this;
}

//! @brief cast the base pointer, otherwise throws an exception
NuTo::BoundaryGradientDamage1D* NuTo::BoundaryGradientDamage1D::AsBoundaryGradientDamage1D()
{
	return this;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::BoundaryGradientDamage1D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::BoundaryGradientDamage1D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::BoundaryGradientDamage1D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::BoundaryGradientDamage1D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::BoundaryGradientDamage1D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::BoundaryGradientDamage1D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::BoundaryGradientDamage1D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize BoundaryGradientDamage1D" << std::endl;
#endif
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementBase);
           //& BOOST_SERIALIZATION_NVP(mSection);
    }
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize BoundaryGradientDamage1D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::BoundaryGradientDamage1D)
#endif // ENABLE_SERIALIZATION
