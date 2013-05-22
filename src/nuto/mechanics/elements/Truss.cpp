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
#include "nuto/mechanics/elements/Truss.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
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
//! @param rTangentLocalEqPlasticStrainStrain derivative of the local eq plastic strains with respect to the strain
//! @param rderivativeShapeFunctions of the ip for all shape functions
//! @param rFactor factor including detJ and area
//! @param rResult result
void NuTo::Truss::AddDetJNtB(const std::vector<double>& rShapeFunctions,
		const std::vector<double>& rDerivativeShapeFunctions, double rFactor, int rRow, int rCol,
		FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult)
{
    assert(rShapeFunctions.size()==rDerivativeShapeFunctions.size());
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
    assert(rShapeFunctions.size()==rDerivativeShapeFunctions.size());

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

//! @brief add Kkk*omega+detJ*F (detJ is already included in Koo)
//! @param rShapeFunctions of the ip for all shape functions
//! @param rLocalEqPlasticStrain local eq. plastic strain values
//! @param rKkk stiffness matrix Kkk
//! @param rNodeNonlocalEqPlasticStrain nodal nonlocal eq plastic strain values
//! @param rFactor factor including detJ and area
//! @param rResult result
void NuTo::Truss::AddDetJRnonlocalplasticStrain(const std::vector<double>& rShapeFunctions,const LocalEqPlasticStrain& rLocalEqPlasticStrain,
		const FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rKkk,
		const std::vector<double>& rNodeNonlocalEqPlasticStrain, double rFactor, int rRow, FullVector<double,Eigen::Dynamic>& rResult)
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

//! @brief add Kkk*omega+detJ*F (detJ is already included in Koo)
//! @param rShapeFunctions of the ip for all shape functions
//! @param rLocalTotalStrain local total strain values
//! @param rKkk stiffness matrix Kkk
//! @param rNodeNonlocalTotalStrain nodal nonlocal total strain values
//! @param rFactor factor including detJ and area
//! @param rResult result
void NuTo::Truss::AddDetJRnonlocalTotalStrain(const std::vector<double>& rShapeFunctions,const EngineeringStrain1D& rLocalTotalStrain,
		const FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rKkk,
		const std::vector<double>& rNodeNonlocalTotalStrain, double rFactor, int rRow, FullVector<double,Eigen::Dynamic>& rResult)
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
	assert((int)rTemperatures.size()==GetNumShapeFunctions());
    for (int count=0; count<GetNumShapeFunctions(); count++)
    {
        if (GetNode(count)->GetNumTemperatures()!=1)
            throw MechanicsException("[NuTo::Truss::CalculateTemperatures] Temperature is required as input to the constitutive model, but the node does not have this data.");
        GetNode(count)->GetTemperature(rTimeDerivative,&(rTemperatures[count]));
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
	assert((int)rNodalNonlocalEqPlasticStrain.size()==2*GetNumShapeFunctions());
	double nonlocalEqPlasticStrain[2];
    for (int count=0; count<GetNumShapeFunctions(); count++)
    {
        if (GetNode(count)->GetNumNonlocalEqPlasticStrain()!=2)
            throw MechanicsException("[NuTo::Truss::CalculateNodalDamage] Damage is required as input to the constitutive model, but the node does not have this data.");
        GetNode(count)->GetNonlocalEqPlasticStrain(nonlocalEqPlasticStrain);
        rNodalNonlocalEqPlasticStrain[count] = nonlocalEqPlasticStrain[0];
        rNodalNonlocalEqPlasticStrain[count+GetNumShapeFunctions()] = nonlocalEqPlasticStrain[1];
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
	assert((int)rNodalNonlocalTotalStrain.size()==GetNumShapeFunctions());
	double nonlocalTotalStrain[1];
    for (int count=0; count<GetNumShapeFunctions(); count++)
    {
        if (GetNode(count)->GetNumNonlocalTotalStrain()!=1)
            throw MechanicsException("[NuTo::Truss::CalculateNodalNonlocalStrain] nonlocal strain is required as input to the constitutive model, but the node does not have this data.");
        GetNode(count)->GetNonlocalTotalStrain1D(nonlocalTotalStrain);
        rNodalNonlocalTotalStrain[count] = nonlocalTotalStrain[0];
    }
}

//! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the stiffness matrix
/*NuTo::Error::eError NuTo::Truss::CalculateCoefficientMatrix_0(NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoefficientMatrix,
        std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn, bool& rSymmetry)const
{
    // get section information determining which input on the constitutive level should be used
    const SectionBase* section(GetSection());
    if (section==0)
        throw MechanicsException("[NuTo::Truss::CalculateCoefficientMatrix_0] no section allocated for element.");

    //calculate local coordinates
    std::vector<double> localNodeCoord(GetNumShapeFunctions());
    CalculateLocalCoordinates(localNodeCoord);

    //calculate local displacements
    int numDisp(GetNumShapeFunctions());
    std::vector<double> localNodeDisp(numDisp);
    if (section->GetInputConstitutiveIsDeformationGradient())
    {
        CalculateLocalDisplacements(localNodeDisp);
    }
    int numDisp(0);
    if (section->GetIsDisplacementDof())
        numDisp = numDisp;

    //calculate temperatures
    int numTemp(GetNumShapeFunctions());
    std::vector<double> nodeTemp(numTemp);
    if (section->GetInputConstitutiveIsTemperatureGradient() || section->GetInputConstitutiveIsTemperature())
    {
        throw MechanicsException("[NuTo::Truss::CalculateCoefficientMatrix_0] temperature not fully implemented.");
    	//CalculateTemperatures(nodeTemp);
    }
    int numTemp(0);
    if (section->GetIsTemperatureDof())
        numTemp=numTemp;

    //allocate space for local ip coordinates
    double localIPCoord;

    //allocate space for local shape functions
    std::vector<double> derivativeShapeFunctions(GetLocalDimension()*GetNumShapeFunctions());

    //allocate deformation gradient
    DeformationGradient1D deformationGradient;

    //allocate deformation gradient
    TemperatureGradient1D temperatureGradient;

    //allocate deformation gradient
    ConstitutiveTangentLocal1x1 tangentStressStrain;

    //allocate deformation gradient
    ConstitutiveTangentLocal1x1 tangentFluxGradTemp;

    //material pointer
    const ConstitutiveEngineeringStressStrain *constitutivePtr;

    //allocate and initialize result matrix
    rCoefficientMatrix.Resize(GetNumLocalDofs(),GetNumLocalDofs());
    bool areAllIpsSymmetric=(true);
    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, localIPCoord);

        CalculateDerivativeShapeFunctions(localIPCoord, derivativeShapeFunctions);

        // determine deformation gradient from the local Displacements and the derivative of the shape functions
        // this is not included in the AddIpStiffness to avoid reallocation of the deformation gradient for each IP
        CalculateDeformationGradient(derivativeShapeFunctions, localNodeCoord, localNodeDisp, deformationGradient);

        //call material law to calculate tangent
        constitutivePtr = GetConstitutiveLaw(theIP)->AsConstitutiveEngineeringStressStrain();
        if (constitutivePtr==0)
            throw MechanicsException("[NuTo::Solid::GetEngineeringStress] Constitutive law can not deal with engineering stresses and strains");
        Error::eError error = constitutivePtr->GetTangent_EngineeringStress_EngineeringStrain(this, theIP,
                deformationGradient, &tangentStressStrain);
        if (error!=Error::SUCCESSFUL)
        	return error;

        areAllIpsSymmetric &= tangentStressStrain.GetSymmetry();

        // calculate local stiffness matrix
        // don't forget to include determinant of the Jacobian and area
        // theoretically, the factor  is
        // detJ * area*BtCB, but B includes 1/detJ, which finally gives:
        double factor (mSection->GetArea()
                       /DetJacobian(derivativeShapeFunctions,localNodeCoord)
                       *(mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP)));

        AddDetJBtCB(derivativeShapeFunctions,tangentStressStrain, factor, rCoefficientMatrix);
    }

    // eventually blow local matrix to full matrix - only relevant for
    // truss in 2D and 3D
    BlowLocalMatrixToGlobal(rCoefficientMatrix);

    // symmetry flag
    rSymmetry = areAllIpsSymmetric;

    // calculate list of global dofs related to the entries in the element stiffness matrix
    this->CalculateGlobalRowDofs(rGlobalDofsRow);
	this->CalculateGlobalColumnDofs(rGlobalDofsColumn);

    return Error::SUCCESSFUL;
}
*/

//! @brief calculates the gradient of the internal potential
//! for a mechanical problem, this corresponds to the internal force vector
/*NuTo::Error::eError NuTo::Truss::CalculateGradientInternalPotential(NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult,
        std::vector<int>& rGlobalDofs)const
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

    //allocate global engineering stress
    EngineeringStress1D engineeringStress;

    //material pointer
    const ConstitutiveEngineeringStressStrain *constitutivePtr;

    //allocate and initialize result matrix
    rResult.Resize(GetNumLocalDofs(),1);
    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, localIPCoord);

        CalculateDerivativeShapeFunctions(localIPCoord, derivativeShapeFunctions);

        // determine deformation gradient from the local Displacements and the derivative of the shape functions
        // this is not included in the AddIpStiffness to avoid reallocation of the deformation gradient for each IP
        CalculateDeformationGradient(derivativeShapeFunctions, localNodeCoord, localNodeDisp, deformationGradient);

        //call material law to calculate stress
        constitutivePtr = dynamic_cast<const ConstitutiveEngineeringStressStrain*>(GetConstitutiveLaw(theIP));
        if (constitutivePtr==0)
            throw MechanicsException("[NuTo::Truss::CalculateGradientInternalPotential] Constitutive law can not deal with engineering stresses and strains");
        Error::eError error = constitutivePtr->GetEngineeringStressFromEngineeringStrain(this, theIP,
                deformationGradient, engineeringStress);
        if (error!=Error::SUCCESSFUL)
        	return error;

        // calculate local stiffness matrix
        // theoretically, the factor  is
        // detJ * area*BtSigma, but B includes 1/detJ, which finally gives:
        double factor (mSection->GetArea()
                       *(mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP)));

        AddDetJBtSigma(derivativeShapeFunctions,engineeringStress, factor, rResult);
    }

    // eventually blow local matrix to full matrix - only relevant for
    // truss in 2D and 3D
    BlowLocalVectorToGlobal(rResult);

    // calculate list of global dofs related to the entries in the element stiffness matrix
    this->CalculateGlobalRowDofs(rGlobalDofs);

    return Error::SUCCESSFUL;
}
*/

//! @brief Update the static data of an element
/*NuTo::Error::eError NuTo::Truss::UpdateStaticData(NuTo::Element::eUpdateType rUpdateType)
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

    //material pointer
    const ConstitutiveEngineeringStressStrain *constitutivePtr;

    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, localIPCoord);

        CalculateDerivativeShapeFunctions(localIPCoord, derivativeShapeFunctions);

        // determine deformation gradient from the local Displacements and the derivative of the shape functions
        // this is not included in the AddIpStiffness to avoid reallocation of the deformation gradient for each IP
        CalculateDeformationGradient(derivativeShapeFunctions, localNodeCoord, localNodeDisp, deformationGradient);

        //call material law to update static data
        constitutivePtr = GetConstitutiveLaw(theIP)->AsConstitutiveEngineeringStressStrain();
        if (constitutivePtr==0)
            throw MechanicsException("[NuTo::Truss::UpdateStaticData] Constitutive law can not deal with engineering stresses and strains");
        Error::eError error;
        switch(rUpdateType)
        {
        case NuTo::Element::STATICDATA:
        	error = constitutivePtr->UpdateStaticData_EngineeringStress_EngineeringStrain(this, theIP, deformationGradient);
        break;
        case NuTo::Element::TMPSTATICDATA:
        	error = constitutivePtr->UpdateTmpStaticData_EngineeringStress_EngineeringStrain(this, theIP, deformationGradient);
        break;
        default:
        	throw MechanicsException("[NuTo::Truss::UpdateStaticData] update static data flag not known (neither static not tmpstatic data");
        }
        if (error!=Error::SUCCESSFUL)
        	return error;
    }

    return Error::SUCCESSFUL;
}
*/

//! @brief calculates deformation gradient1D
//! @param rRerivativeShapeFunctions derivatives of the shape functions
//! @param rLocalDisp local displacements
//! @param rConstitutiveInput (return value)
void NuTo::Truss::CalculateDeformationGradient(const std::vector<double>& rDerivativeShapeFunctions,
        const std::vector<double>& rLocalCoord,
        const std::vector<double>& rLocalDisp,
        DeformationGradient1D& rDeformationGradient)const
{
    assert((int)rLocalDisp.size()==GetNumNodes() && (int)rDerivativeShapeFunctions.size()==GetNumNodes());
    rDeformationGradient.mDeformationGradient = 0;

    //normally, the inverse Jacobian should be calculated, but for a truss element, it is sufficient to use the inverse of the Jacobian determinant
    double factor(1./DetJacobian(rDerivativeShapeFunctions, rLocalCoord));
    for (int count=0; count<GetNumNodes(); count++)
    {
        rDeformationGradient.mDeformationGradient+=rLocalDisp[count]*rDerivativeShapeFunctions[count];
    }
    rDeformationGradient.mDeformationGradient*=factor;
    rDeformationGradient.mDeformationGradient+=1.;
}

//! @brief calculates the coefficient matrix for the 1-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the damping matrix
/*NuTo::Error::eError NuTo::Truss::CalculateCoefficientMatrix_1(NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult,
        std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn, bool& rSymmetry)const
{
    throw MechanicsException("[NuTo::Truss::CalculateCoefficientMatrix_1] to be implemented.");
}
*/

//! @brief calculates the coefficient matrix for the 2-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the Mass matrix
/*NuTo::Error::eError NuTo::Truss::CalculateCoefficientMatrix_2(NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoefficientMatrix,
        std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn, bool& rSymmetry)const
{
    //calculate local coordinates
    std::vector<double> localNodeCoord(this->GetNumLocalDofs());
    this->CalculateLocalCoordinates(localNodeCoord);

    //allocate space for local shape functions
    std::vector<double> derivativeShapeFunctions(this->GetLocalDimension()*this->GetNumShapeFunctions());
    std::vector<double> shapeFunctions(this->GetNumShapeFunctions());

    //allocate and initialize result matrix
    rCoefficientMatrix.Resize(this->GetNumLocalDofs(),this->GetNumLocalDofs());
    for (int theIP=0; theIP<this->GetNumIntegrationPoints(); theIP++)
    {
    	double localIPCoord;
        this->GetLocalIntegrationPointCoordinates(theIP, localIPCoord);

        this->CalculateShapeFunctions(localIPCoord, shapeFunctions);
        this->CalculateDerivativeShapeFunctions(localIPCoord, derivativeShapeFunctions);

        //call material law to calculate tangent
        const ConstitutiveBase* constitutivePtr = this->GetConstitutiveLaw(theIP);
        if (constitutivePtr==0)
        {
            throw MechanicsException("[NuTo::Truss::CalculateCoefficientMatrix_2] Constitutive law can not found at integration point.");
        }
        double density = constitutivePtr->GetDensity();

        // calculate local mass matrix
        // don't forget to include determinant of the Jacobian and area
        // detJ * area * density * HtH, :
        double factor (density * this->mSection->GetArea() * this->DetJacobian(derivativeShapeFunctions,localNodeCoord)
                       *(this->mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP)));
        this->AddDetJHtH(shapeFunctions, factor, rCoefficientMatrix);
    }

    // eventually blow local matrix to full matrix - only relevant for
    // truss in 2D and 3D
    this->BlowLocalMatrixToGlobal(rCoefficientMatrix);

    // calculate list of global dofs related to the entries in the element stiffness matrix
    this->CalculateGlobalRowDofs(rGlobalDofsRow);
    this->CalculateGlobalColumnDofs(rGlobalDofsColumn);
    rSymmetry = true;

    return Error::SUCCESSFUL;
}
*/

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
    CalculateShapeFunctions(naturalCoordinates, shapeFunctions);
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

//! @brief returns the nonlocal eq plastic strain interpolated from the nodal values
//! @param shapeFunctionsGlobal shape functions
//! @param rNodeDamage nonlocal eq plastic strain values of the nodes
//! @param rNonlocalEqentPlasticStrain return value
void NuTo::Truss::CalculateNonlocalEqPlasticStrain(const std::vector<double>& shapeFunctions,
		const std::vector<double>& rNodeEquivalentPlasticStrain, NonlocalEqPlasticStrain& rNonlocalEqentPlasticStrain)
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
		const std::vector<double>& rNodeNonlocalTotalStrain, EngineeringStrain1D& rNonlocalTotalStrain)
{
    assert(shapeFunctions.size()==rNodeNonlocalTotalStrain.size());
    rNonlocalTotalStrain(0) = 0.;
	for (unsigned int count=0; count<shapeFunctions.size(); count++)
	{
		rNonlocalTotalStrain(0)+=shapeFunctions[count]*rNodeNonlocalTotalStrain[count];
	}
	return;
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
    std::vector<double> localNodeCoord(GetNumLocalDofs());
    CalculateLocalCoordinates(localNodeCoord);

    //allocate space for local ip coordinates
    double localIPCoord;

    //allocate space for local shape functions
    std::vector<double> derivativeShapeFunctions(GetLocalDimension()*GetNumShapeFunctions());

	rVolume.resize(GetNumIntegrationPoints());

     for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, localIPCoord);

        CalculateDerivativeShapeFunctions(localIPCoord, derivativeShapeFunctions);

		//attention in 1D, this is just the length, but that is required for the nonlocal model
		rVolume[theIP] = DetJacobian(derivativeShapeFunctions,localNodeCoord)
                       *(mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP));
    }
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
