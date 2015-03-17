// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

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
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentNonlocal.h"
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixInt.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixDouble.h"
#include "nuto/mechanics/elements/ElementOutputVectorInt.h"
#include "nuto/mechanics/elements/Plane.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/structures/StructureBase.h"

#include "nuto/math/FullMatrix.h"

//! @brief constructor
NuTo::Plane::Plane(const StructureBase* rStructure, ElementData::eElementDataType rElementDataType,
        IntegrationType::eIntegrationType rIntegrationType, IpData::eIpDataType rIpDataType) :
        NuTo::ElementBase::ElementBase(rStructure, rElementDataType, rIntegrationType, rIpDataType)
{
    mSection = 0;
}


//! @brief calculates output data fo the elmement
//! @param eOutput ... coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
//!                    @param updateStaticData (with DummyOutput), IPData, globalrow/column dofs etc.
NuTo::Error::eError NuTo::Plane::Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput)
{
	if (mStructure->GetHessianConstant(1)==false)
    	throw MechanicsException("[NuTo::Plane::Evaluate] only implemented for a constant Hessian for the first derivative (damping).");
    if (mStructure->GetHessianConstant(2)==false)
    	throw MechanicsException("[NuTo::Plane::Evaluate] only implemented for a constant Hessian for the second derivative (mass).");

    try
    {
    	// get section information determining which input on the constitutive level should be used
		const SectionBase* section(GetSection());
		if (section==0)
			throw MechanicsException("[NuTo::Plane::Evaluate] no section allocated for element.");

		//calculate coordinates
		int numLocalCoordinates(2*GetNumNodesGeometry());
		std::vector<double> localNodeCoord(numLocalCoordinates);
		CalculateLocalCoordinates(localNodeCoord);

		//calculate local displacements
		int numLocalDisp(2*GetNumNodesField());
		std::vector<double> localNodeDisp(numLocalDisp);
		if (section->GetInputConstitutiveIsDeformationGradient())
		{
			CalculateLocalDisplacements(localNodeDisp);
		}
		int numLocalDispDofs(0);
		if (section->GetIsDisplacementDof())
			numLocalDispDofs = numLocalDisp;

		//calculate temperatures
		int numTemp(GetNumNodesField());
		std::vector<double> nodeTemp(numTemp);
		if (section->GetInputConstitutiveIsTemperatureGradient() || section->GetInputConstitutiveIsTemperature())
		{
			CalculateTemperatures(nodeTemp);
		}
		int numTempDofs(0);
		if (section->GetIsTemperatureDof())
			numTempDofs=numTemp;

		//allocate space for ip coordinates in natural coordinate system (-1,1)
		double naturalIPCoord[2];
		double nonlocalNaturalIPCoord[2];

		//allocate space for derivatives of shape functions in natural coordinate system
		std::vector<double> derivativeShapeFunctionsGeometryNatural(2*GetNumNodesGeometry());
		std::vector<double> derivativeShapeFunctionsFieldNatural(2*GetNumNodesField());
		//allocate space for derivatives of shape functions in local coordinate system
		std::vector<double> derivativeShapeFunctionsFieldLocal(2*GetNumNodesField());
		std::vector<double> shapeFunctionsField(GetNumNodesField());    //allocate space for shape functions

		std::vector<double> nonlocalDerivativeShapeFunctionsGeometryNatural;
		std::vector<double> nonlocalDerivativeShapeFunctionsFieldNatural;
		std::vector<double> nonlocalDerivativeShapeFunctionsFieldLocal;

		//allocate deformation gradient
		DeformationGradient2D deformationGradient;

		EngineeringStrain2D engineeringStrain2D;
		EngineeringStrain3D engineeringStrain3D;

		//allocate global engineering plastic strain
		EngineeringStrain3D engineeringPlasticStrain3D;

		//allocate damage
		Damage damage;

		//allocate temperature
	    Temperature temperature;

		//allocate temperature gradient
		TemperatureGradient2D temperatureGradient;

		//allocate global engineering stress
		EngineeringStress2D engineeringStress2D;
		EngineeringStress3D engineeringStress3D;

		//allocate global heat flux
		HeatFlux2D heatFlux2D;
		HeatFlux3D heatFlux3D;

		//allocate vector of tangent matrices
		int NumNonlocalElements(GetNumNonlocalElements());
		// in case of a local formulation, just use a nonlocal matrix with 1 entry

		// calculate number of DOFS involved in all nonlocal elements
		int NumNonlocalIps;
		int numNonlocalDispDofs(0);
		if (NumNonlocalElements==0)
		{
			NumNonlocalIps = 1;
			numNonlocalDispDofs = numLocalDispDofs;
		}
		else
		{
			//allocate and initialize result matrix
			const std::vector<const ElementBase*>& nonlocalElements(GetNonlocalElements());
			NumNonlocalIps=0;
			numNonlocalDispDofs=0;
			for (int theNonlocalElement=0; theNonlocalElement<NumNonlocalElements; theNonlocalElement++)
			{
				numNonlocalDispDofs += nonlocalElements[theNonlocalElement]->AsPlane()->GetNumNodesField()*2;
				NumNonlocalIps += nonlocalElements[theNonlocalElement]->GetNumIntegrationPoints();
			}
		}

		NuTo::ConstitutiveTangentNonlocal<3,3> nonlocalTangentStressStrain;
		//ConstitutiveTangentNonlocal<3,1> nonlocalTangentStressTemperature(NumNonlocalIps)[3];
		ConstitutiveTangentLocal<2,2> tangentHeatFluxTemperatureGradient;
		ConstitutiveTangentLocal<2,1> tangentHeatFluxTemperatureRate;

		//InvJacobian and determinant of Jacobian
		double invJacobian[4], detJac;

		//for the lumped mass calculation
		double total_mass(0.);

		//InvJacobian and determinant of Jacobian for nonlocal model
		 double nonlocalInvJacobian[4], nonlocalDetJac;

		//define inputs and outputs
		std::map< NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase* > constitutiveInputList;
		std::map< NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase* > constitutiveOutputList;

		if (numLocalDispDofs>0)
		{
			constitutiveInputList[NuTo::Constitutive::Input::DEFORMATION_GRADIENT_2D] = &deformationGradient;
		}

		if (numTempDofs>0)
		{
			constitutiveInputList[NuTo::Constitutive::Input::TEMPERATURE_GRADIENT_2D] = &temperatureGradient;
		}

		if (section->GetInputConstitutiveIsTemperature())
		{
			constitutiveInputList[NuTo::Constitutive::Input::TEMPERATURE] = &temperature;
		}

		//define outputs
		for (auto it = rElementOutput.begin(); it!=rElementOutput.end(); it++)
		{
			switch(it->first)
			{
			case Element::INTERNAL_GRADIENT:
				it->second->GetFullVectorDouble().Resize(numLocalDispDofs+numTempDofs);
	    	    //if the stiffness matrix is constant, the corresponding internal force is calculated via the Kd
	    		//on the global level
	    		if (mStructure->GetHessianConstant(0)==false)
	    	    {
					if (numLocalDispDofs>0)
					{
						 constitutiveOutputList[NuTo::Constitutive::Output::ENGINEERING_STRESS_2D] = &engineeringStress2D;
					}
					if (numTempDofs>0)
					{
						constitutiveOutputList[NuTo::Constitutive::Output::HEAT_FLUX_2D] = &heatFlux2D;
					}
	    	    }
			break;
			case Element::HESSIAN_0_TIME_DERIVATIVE:
				{
					it->second->GetFullMatrixDouble().Resize(numLocalDispDofs+numTempDofs,numNonlocalDispDofs+numTempDofs);
					it->second->GetFullMatrixDouble().setZero();
					it->second->SetSymmetry(true);
					it->second->SetConstant(true);

					if (numLocalDispDofs>0)
					{
						if (NumNonlocalElements==0)
						{
							nonlocalTangentStressStrain.SetNumSubMatrices(1);
							constitutiveOutputList[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_2D] = &(nonlocalTangentStressStrain.GetSubMatrix_3x3(0));
						}
						else
						{
							nonlocalTangentStressStrain.SetNumSubMatrices(NumNonlocalIps);
							constitutiveOutputList[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_2D] = &nonlocalTangentStressStrain;
						}

						//mixed term
						if (numTempDofs)
							//constitutiveOutputList[NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_TEMPERATURE_2D] = &tangentStressTemperature[timeDerivative];
							throw MechanicsException("[NuTo::Plane::Evaluate] Mixed problem not implemented.");
					}
					if (numTempDofs>0)
					{
						constitutiveOutputList[NuTo::Constitutive::Output::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT_2D] = &tangentHeatFluxTemperatureGradient;
						//mixed term
						//if (numDispDofs)
						//    constitutiveOutputList.insert(std::pair<NuTo::Constitutive::eOutput, ConstitutiveOutputBase*>(NuTo::Constitutive::eOutput::D_HEAT_FLUX_D_ENGINEERING_STRAIN_3D, &tangentHeatFluxEngineeringStrain[timeDerivative]));
					}
				}
			break;
			case Element::HESSIAN_1_TIME_DERIVATIVE:
                it->second->GetFullMatrixDouble().Resize(numLocalDispDofs+numTempDofs,numNonlocalDispDofs+numTempDofs);
                it->second->SetSymmetry(true);
                it->second->SetConstant(true);
                if (numLocalDispDofs>0)
                {
                        // Rayleigh damping should be introduced on the global level
                }
                if (numTempDofs>0)
                {
                        constitutiveOutputList[NuTo::Constitutive::Output::D_HEAT_FLUX_D_TEMPERATURE_RATE_2D] = &tangentHeatFluxTemperatureRate;
                }
			break;
			case Element::HESSIAN_2_TIME_DERIVATIVE:
	        {
				it->second->GetFullMatrixDouble().Resize(numLocalDispDofs+numTempDofs,numNonlocalDispDofs+numTempDofs);
				it->second->SetSymmetry(true);
				it->second->SetConstant(true);
				//there is only a constant mass part for the mechanics problem
	        }
			break;
			case Element::LUMPED_HESSIAN_2_TIME_DERIVATIVE:
				it->second->GetFullVectorDouble().Resize(numLocalDispDofs);
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
					throw MechanicsException("[NuTo::Plane::Evaluate] this ip data type is not implemented.");
				}
			break;
			case Element::GLOBAL_ROW_DOF:
				this->CalculateGlobalRowDofs(it->second->GetVectorInt(),numLocalDispDofs,numTempDofs);
			break;
			case Element::GLOBAL_COLUMN_DOF:
				this->CalculateGlobalColumnDofs(it->second->GetVectorInt(),numNonlocalDispDofs,numTempDofs);
			break;
			default:
				throw MechanicsException("[NuTo::Plane::Evaluate] element output not implemented.");
			}
		}

		// loop over the integration points
		for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
		{
			GetLocalIntegrationPointCoordinates(theIP, naturalIPCoord);

			CalculateDerivativeShapeFunctionsGeometryNatural(naturalIPCoord, derivativeShapeFunctionsGeometryNatural);
			CalculateJacobian(derivativeShapeFunctionsGeometryNatural,localNodeCoord, invJacobian, detJac);

//Das muss bei Masseberechnung z.B. nicht berechnet werden
			if (GetNumNodesGeometry()==GetNumNodesField())
			{
				//isoparametric elements
				CalculateDerivativeShapeFunctionsLocal(derivativeShapeFunctionsGeometryNatural,invJacobian,
														derivativeShapeFunctionsFieldLocal);
			}
			else
			{
				// sub or superparametric elements
				CalculateDerivativeShapeFunctionsFieldNatural(naturalIPCoord, derivativeShapeFunctionsFieldNatural);
				CalculateDerivativeShapeFunctionsLocal(derivativeShapeFunctionsFieldNatural,invJacobian,
														derivativeShapeFunctionsFieldLocal);
			}

			if (section->GetInputConstitutiveIsDeformationGradient())
			{
				// determine deformation gradient from the local Displacements and the derivative of the shape functions
				// this is not included in the AddIpStiffness to avoid reallocation of the deformation gradient for each IP
				CalculateDeformationGradient(derivativeShapeFunctionsFieldLocal, localNodeDisp, deformationGradient);
			}
			if (section->GetInputConstitutiveIsTemperatureGradient())
			{
				// determine deformation gradient from the local Displacements and the derivative of the shape functions
				// this is not included in the AddIpStiffness to avoid reallocation of the deformation gradient for each IP
				CalculateTemperatureGradient(derivativeShapeFunctionsFieldLocal, nodeTemp, temperatureGradient);
			}

	        ConstitutiveBase* constitutivePtr = GetConstitutiveLaw(theIP);
			try
			{
				Error::eError error = constitutivePtr->Evaluate2D(this, theIP,
						constitutiveInputList, constitutiveOutputList);
				if (error!=Error::SUCCESSFUL)
					return error;
			}
			catch (NuTo::MechanicsException &e)
			{
				e.AddMessage("[NuTo::Plane::Evaluate] error evaluating the constitutive model.");
				throw e;
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
						// Jacobian
						double factor(mSection->GetThickness()*fabs(detJac*(mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP))));
						if (numLocalDispDofs>0)
						{
							AddDetJBtSigma(derivativeShapeFunctionsFieldLocal,engineeringStress2D, factor, 0, it->second->GetFullVectorDouble());
						}
						if (numTempDofs>0)
						{
							AddDetJBtHeatFlux(derivativeShapeFunctionsFieldLocal,heatFlux2D, factor, numLocalDispDofs, it->second->GetFullVectorDouble());
						}
		    	    }
				}
				break;
				case Element::HESSIAN_0_TIME_DERIVATIVE:
					{
						//factor for the numerical integration
						assert(mSection->GetThickness()>0);
						double factor(mSection->GetThickness()*detJac*(mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP)));
						if (nonlocalTangentStressStrain.GetConstant()==false)
							it->second->SetConstant(false);
						if (nonlocalTangentStressStrain.GetSymmetry()==false)
							it->second->SetSymmetry(false);

						if (numLocalDispDofs>0)
						{
							if (NumNonlocalElements!=0)
							{
								//Nonlocal Model where the material model is still local
								if (nonlocalTangentStressStrain.GetLocalSolution())
								{
									//same as local model, e.g. in the unloading range or

									// calculate element stiffness matrix
									// don't forget to include determinant of the Jacobian and area
									AddDetJBtCB(derivativeShapeFunctionsFieldLocal, derivativeShapeFunctionsFieldLocal, nonlocalTangentStressStrain.GetSubMatrix_3x3(0), factor, it->second->GetFullMatrixDouble(), 0, 0);
								}
								else
								{
									//nonlocal BMatrix of nonlocal other integration point
									//sum over nonlocal elements and their ips
									const std::vector<const ElementBase*>& nonlocalElements(GetNonlocalElements());
									int firstCol(0);
									int totalNonlocalIp(0);
									for (unsigned int theNonlocalElement = 0; theNonlocalElement<nonlocalElements.size(); theNonlocalElement++)
									{
										const Plane* nonlocalElement(nonlocalElements[theNonlocalElement]->AsPlane());
										//calculate local coordinates
										std::vector<double> nonlocalLocalNodeCoord(nonlocalElement->GetNumNodesField()*2);
										nonlocalElement->CalculateLocalCoordinates(nonlocalLocalNodeCoord);

										//get weights
										const std::vector<double>& weights(GetNonlocalWeights(theIP,theNonlocalElement));

										for (int theNonlocalIp = 0; theNonlocalIp < nonlocalElement->GetNumIntegrationPoints(); theNonlocalIp++, totalNonlocalIp++)
										{
											if (weights[theNonlocalIp]==0.)
												continue;

											// get IP coordinates
											nonlocalElement->GetLocalIntegrationPointCoordinates(theNonlocalIp, nonlocalNaturalIPCoord);

											// calculate derivatives of shape functions in natural coordinate system
											nonlocalDerivativeShapeFunctionsGeometryNatural.resize(2*nonlocalElement->GetNumNodesGeometry());
											nonlocalElement->CalculateDerivativeShapeFunctionsGeometryNatural(nonlocalNaturalIPCoord, nonlocalDerivativeShapeFunctionsGeometryNatural);

											// calculate Jacobian in order to transform natural to local coordinate system
											nonlocalElement->CalculateJacobian(nonlocalDerivativeShapeFunctionsGeometryNatural,nonlocalLocalNodeCoord, nonlocalInvJacobian, nonlocalDetJac);

											// calculate derivates of shape functions in local coordinate system
											nonlocalDerivativeShapeFunctionsFieldNatural.resize(2*nonlocalElement->GetNumNodesField());
											nonlocalDerivativeShapeFunctionsFieldLocal.resize(nonlocalDerivativeShapeFunctionsFieldNatural.size());

											if (nonlocalElement->GetNumNodesGeometry() == nonlocalElement->GetNumNodesField())
											{
												//isoparametric element
												nonlocalElement->CalculateDerivativeShapeFunctionsLocal(nonlocalDerivativeShapeFunctionsGeometryNatural,nonlocalInvJacobian,
																						nonlocalDerivativeShapeFunctionsFieldLocal);
											}
											else
											{

												nonlocalElement->CalculateDerivativeShapeFunctionsFieldNatural(nonlocalNaturalIPCoord, nonlocalDerivativeShapeFunctionsFieldNatural);
												nonlocalElement->CalculateDerivativeShapeFunctionsLocal(nonlocalDerivativeShapeFunctionsFieldNatural,nonlocalInvJacobian,
																						nonlocalDerivativeShapeFunctionsFieldLocal);
											}

											// calculate element stiffness matrix
											// don't forget to include determinant of the Jacobian and area
											AddDetJBtCB(derivativeShapeFunctionsFieldLocal, nonlocalDerivativeShapeFunctionsFieldLocal, nonlocalTangentStressStrain.GetSubMatrix_3x3(totalNonlocalIp), factor, it->second->GetFullMatrixDouble(), 0, firstCol);
										}
										firstCol+=2*nonlocalElement->GetNumNodesField();
									}
									assert(totalNonlocalIp==NumNonlocalIps);
								}
							}
							else
							{
								// calculate element stiffness matrix
								// don't forget to include determinant of the Jacobian and area
								AddDetJBtCB(derivativeShapeFunctionsFieldLocal,derivativeShapeFunctionsFieldLocal, nonlocalTangentStressStrain.GetSubMatrix_3x3(0), factor, it->second->GetFullMatrixDouble(),0,0);
							}
							if (numTempDofs>0)
								 throw MechanicsException("[NuTo::Plane::Evaluate] mixed terms not yet implemented.");
						} // end disp dof
						if (numTempDofs>0)
						{
							AddDetJBtCB(derivativeShapeFunctionsFieldLocal, tangentHeatFluxTemperatureGradient, factor, numLocalDispDofs,numNonlocalDispDofs, it->second->GetFullMatrixDouble());
							if (tangentHeatFluxTemperatureGradient.GetSymmetry()==false)
								it->second->SetSymmetry(false);
							if (tangentHeatFluxTemperatureGradient.GetConstant()==false)
								it->second->SetConstant(false);
							if (numLocalDispDofs>0)
								throw MechanicsException("[NuTo::Plane::Evaluate] mixed terms not yet implemented.");
						}
					}
				break;
				case Element::HESSIAN_1_TIME_DERIVATIVE:
                {
					if (numLocalDispDofs>0)
					{
									//no damping term, do Rayleigh damping on the global level
					}
					if (numTempDofs>0)
					{
					throw MechanicsException("[Nuto::Plane::Evaluate] temperature not yet implemented.");
					}
                }
                break;
				case Element::HESSIAN_2_TIME_DERIVATIVE:
				{
                    if (numLocalDispDofs>0)
                    {
						this->CalculateShapeFunctionsField(naturalIPCoord, shapeFunctionsField);

						// calculate local mass matrix (the nonlocal terms are zero)
						// don't forget to include determinant of the Jacobian and area
						// detJ * area * density * HtH, :
				        double factor(mSection->GetThickness()*detJac*(mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP))*constitutivePtr->GetDensity());
				        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> tmpMatrix;
				        tmpMatrix = (factor*Eigen::Matrix<double,Eigen::Dynamic,1>::Map(&(shapeFunctionsField[0]),GetNumNodesField()))*Eigen::Matrix<double,1,Eigen::Dynamic>::Map(&(shapeFunctionsField[0]),GetNumNodesField());
				        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& result(it->second->GetFullMatrixDouble());
				        for (int count=0; count<GetNumNodesField(); count++)
				        {
				            for (int count2=0; count2<GetNumNodesField(); count2++)
				            {
				            	result(2*count,2*count2) += tmpMatrix(count,count2);
				            	result(2*count+1,2*count2+1) += tmpMatrix(count,count2);
				            }
				        }

						//if (numTemp>0) no mixed terms
                    }
					if (numTempDofs>0)
					{
							//no termperature terms
					}
				}
				break;
				case Element::LUMPED_HESSIAN_2_TIME_DERIVATIVE:
                    if (numLocalDispDofs>0)
                    {
						this->CalculateShapeFunctionsField(naturalIPCoord, shapeFunctionsField);
						// calculate local mass matrix (the nonlocal terms are zero)
						// don't forget to include determinant of the Jacobian and area
						// detJ * area * density * HtH, :
				        double factor(mSection->GetThickness()*detJac*(mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP))*constitutivePtr->GetDensity());
						FullVector<double,Eigen::Dynamic>& result(it->second->GetFullVectorDouble());
						total_mass+=factor;
            			//calculate for the translational dofs the diagonal entries
						for (int count=0; count<GetNumNodesField(); count++)
						{
							result(2*count)+=shapeFunctionsField[count]*shapeFunctionsField[count]*factor;
						}

						if (theIP+1==GetNumIntegrationPoints())
						{
							//calculate sum of diagonal entries (is identical for all directions, that's why only x direction is calculated
							double sum_diagonal(0);
							for (int count=0; count<GetNumNodesField(); count++)
							{
								sum_diagonal+= result(2*count);
							}

							//scale so that the sum of the diagonals represents the full mass
							double scaleFactor = total_mass/sum_diagonal;
							for (int count=0; count<GetNumNodesField(); count++)
							{
								result(2*count) *= scaleFactor;
								result(2*count+1) = result(2*count);
							}
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
						throw MechanicsException("[NuTo::Plane::GetIpData] Ip data not implemented.");
					}
				break;
				case Element::GLOBAL_ROW_DOF:
				case Element::GLOBAL_COLUMN_DOF:
				break;
				default:
					throw MechanicsException("[NuTo::Plane::Evaluate] element output not implemented.");
				}
			}
		}

		//resize the matrix in case of a local matrix
		/*if (numLocalDispDofs>0)
		{
			auto it = rElementOutput.find(Element::HESSIAN_0_TIME_DERIVATIVE);
		}
		*/
    }
    catch (NuTo::MechanicsException e)
    {
        std::stringstream ss;
        ss << mStructure->ElementGetId(this);
    	e.AddMessage("[NuTo::Plane::Evaluate] Error evaluating element data of element"	+ ss.str() + ".");
        throw e;
    }

    return Error::SUCCESSFUL;
}

//! @brief adds to a matrix the product B^tCBnonlocal, where B contains the derivatives of the shape functions and C is the constitutive tangent and Bnonlocal is the nonlocal B matrix
//! eventually include also area/width of an element (that's the mechanics solution)
//! @param rLocalDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param rNonlocalDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param ConstitutiveTangentBase constitutive tangent matrix
//! @param rFactor factor including determinant of Jacobian and IP weight
//! @param rCoefficientMatrix to be added to
//! &param rFirstCol first column of the coefficient matrix to be modified (corresponding to the current nonlocal element)
void NuTo::Plane::AddDetJBtCB(const std::vector<double>& rLocalDerivativeShapeFunctionsLocal,const std::vector<double>& rNonlocalDerivativeShapeFunctionsLocal,
                              const ConstitutiveTangentLocal<3,3>& rConstitutiveTangent, double rFactor,
                              FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoefficientMatrix, int rFirstRow, int rFirstCol)const
{
    assert(rCoefficientMatrix.GetNumRows()==2*GetNumNodesField() && rFirstCol + (int)rNonlocalDerivativeShapeFunctionsLocal.size()<=rCoefficientMatrix.GetNumColumns());
    assert((int)rLocalDerivativeShapeFunctionsLocal.size()==2*GetNumNodesField());
    const double *C = rConstitutiveTangent.data();
    double x1,x2,y1,y2,x1x2,y2x1,x2y1,y2y1;

    for (int theNode1=0; theNode1<GetNumNodes(); theNode1++)
    {
        int node1mul2 = 2*theNode1;
        int node1mul2plus1 = node1mul2+1;

        x1 = rFactor * rLocalDerivativeShapeFunctionsLocal[node1mul2];
        y1 = rFactor * rLocalDerivativeShapeFunctionsLocal[node1mul2plus1];
        // << is division by two, but faster
        for (unsigned int theNode2=0; theNode2<rNonlocalDerivativeShapeFunctionsLocal.size()>>1; theNode2++)
        {
            int node2mul2 = 2*theNode2;
            int node2mul2plus1 = node2mul2+1;

            x2 = rNonlocalDerivativeShapeFunctionsLocal[node2mul2];
            y2 = rNonlocalDerivativeShapeFunctionsLocal[node2mul2plus1];

            x1x2 = x2*x1;
            y2x1 = y2*x1;
            x2y1 = x2*y1;
            y2y1 = y2*y1;

            rCoefficientMatrix(rFirstRow+node1mul2,rFirstCol+node2mul2)           +=x1x2*C[0] +x2y1*C[2] + y2x1*C[6] +y2y1*C[8];
            rCoefficientMatrix(rFirstRow+node1mul2,rFirstCol+node2mul2plus1)     +=x1x2*C[6] +x2y1*C[8] + y2x1*C[3] +y2y1*C[5];
            rCoefficientMatrix(rFirstRow+node1mul2plus1,rFirstCol+node2mul2)     +=x1x2*C[2] +x2y1*C[1] + y2x1*C[8] +y2y1*C[7];
            rCoefficientMatrix(rFirstRow+node1mul2plus1,rFirstCol+node2mul2plus1)+=x1x2*C[8] +x2y1*C[7] + y2x1*C[5] +y2y1*C[4];
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
void NuTo::Plane::AddDetJBtCB(const std::vector<double>& rDerivativeShapeFunctionsGlobal,
                              const ConstitutiveTangentLocal<2,2>& rConstitutiveTangent, double rFactor,
                              int rRow, int rCol,
                              FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoefficientMatrix)const
{
    const double *C = rConstitutiveTangent.data();
    double x1,x2,y1,y2;
    for (int theNode1=0; theNode1<GetNumNodes(); theNode1++)
    {
        int node1mul2 = theNode1+theNode1;
        int node1mul2plus1 = node1mul2+1;

        assert((int)rDerivativeShapeFunctionsGlobal.size()>node1mul2plus1);
        x1 = rFactor * rDerivativeShapeFunctionsGlobal[node1mul2];
        y1 = rFactor * rDerivativeShapeFunctionsGlobal[node1mul2plus1];
        node1mul2 +=rRow;
        node1mul2plus1 +=rRow;
        for (int theNode2=0; theNode2<GetNumNodes(); theNode2++)
        {
            int node2mul2 = theNode2 + theNode2;
            int node2mul2plus1 = node2mul2+1;
            node2mul2 +=rCol;
            node2mul2plus1 +=rCol;

            assert((int)rDerivativeShapeFunctionsGlobal.size()>node2mul2plus1);
            x2 = rDerivativeShapeFunctionsGlobal[node2mul2];
            y2 = rDerivativeShapeFunctionsGlobal[node2mul2plus1];

            assert(rCoefficientMatrix.GetNumRows()>node1mul2plus1 && rCoefficientMatrix.GetNumColumns()>node1mul2plus1);
            assert(rCoefficientMatrix.GetNumRows()>node2mul2plus1 && rCoefficientMatrix.GetNumColumns()>node2mul2plus1);

            rCoefficientMatrix(node1mul2,node2mul2)          +=x1*C[0]*x2;
            rCoefficientMatrix(node1mul2,node2mul2plus1)     +=x1*C[2]*y2;
            rCoefficientMatrix(node1mul2plus1,node2mul2)     +=y1*C[1]*x2;
            rCoefficientMatrix(node1mul2plus1,node2mul2plus1)+=y1*C[3]*y2;
        }
    }
}


//! @brief adds up the internal force vector
//! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param rEngineeringStress stress
//! @param factor factor including det Jacobian area and integration point weight
//! @param rRow start row (in case of a multifield problem)
//! @param rResult resforce vector
void NuTo::Plane::AddDetJBtSigma(const std::vector<double>& rDerivativeShapeFunctionsLocal,
                                 const EngineeringStress2D& rEngineeringStress,
                                 double rFactor,
                                 int rRow,
                                 FullVector<double,Eigen::Dynamic>& rResult)const
{
    assert(rResult.GetNumRows()==2*GetNumNodesField() && rResult.GetNumColumns()==1);
    const double *s = rEngineeringStress.GetData();
    double x1,y1;
    for (int theNode1=0; theNode1<GetNumNodes(); theNode1++)
    {
        int node1mul2 = 2*theNode1;
        int node1mul2plus1 = node1mul2+1;

        x1 = rFactor * rDerivativeShapeFunctionsLocal[node1mul2];
        y1 = rFactor * rDerivativeShapeFunctionsLocal[node1mul2plus1];

        rResult(rRow + node1mul2)     +=x1*s[0]+y1*s[2];
        rResult(rRow + node1mul2plus1)+=y1*s[1]+x1*s[2];
    }
}

//! @brief adds up the internal force vector
//! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param rHeatFlux stress
//! @param factor factor including det Jacobian area and integration point weight
//! @param rRow start row (in case of a multifield problem)
//! @param rResult resforce vector
void NuTo::Plane::AddDetJBtHeatFlux(const std::vector<double>& rDerivativeShapeFunctionsGlobal,
                                 const HeatFlux2D& rHeatFlux,
                                 double rFactor,
                                 int rRow,
                                 FullVector<double,Eigen::Dynamic>& rResult)const
{
    const double *s = rHeatFlux.GetData();
    double x1,y1;
    for (int theNode1=0; theNode1<GetNumNodes(); theNode1++)
    {
        int node1mul2 = 2*theNode1;
        int node1mul2plus1 = node1mul2+1;

        assert(rResult.GetNumRows()>node1mul2plus1);
        x1 = rFactor * rDerivativeShapeFunctionsGlobal[node1mul2];
        y1 = rFactor * rDerivativeShapeFunctionsGlobal[node1mul2plus1];

        rResult(rRow + node1mul2)     +=x1*s[0];
        rResult(rRow + node1mul2plus1)+=y1*s[1];
    }
}


//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Plane::CalculateShapeFunctionsField(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)const
{
    CalculateShapeFunctionsGeometry(rNaturalCoordinates,rShapeFunctions);
}

//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for the nonlocal eq strain nodes
void NuTo::Plane::CalculateShapeFunctionsNonlocalEqStrain(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)const
{
    CalculateShapeFunctionsGeometry(rNaturalCoordinates,rShapeFunctions);
}

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Plane::CalculateDerivativeShapeFunctionsFieldNatural(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const
{
    CalculateDerivativeShapeFunctionsGeometryNatural(rNaturalCoordinates,rDerivativeShapeFunctions);
}

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::Plane::CalculateDerivativeShapeFunctionsNonlocalEqStrainNatural(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)const
{
    CalculateDerivativeShapeFunctionsGeometryNatural(rNaturalCoordinates,rDerivativeShapeFunctions);
}

//! @brief Calculates the the inverse of the Jacobian and its determinant
//! @param rDerivativeShapeFunctions Derivatives of the shape functions (dN1dx, dN1dy, dN1dz, dN2dx, ..
//! @param rNodeCoordinates Node coordinates (X1,Y1,Z1,X2,Y2,Z2,...
//! @param rInvJacobian inverse Jacobian matrix (return value)
//! @param rDetJac determinant of the Jacobian (return value)
void NuTo::Plane::CalculateJacobian(const std::vector<double>& rDerivativeShapeFunctions,
                                    const std::vector<double>& rNodeCoordinates,
                                    double rInvJacobian[4],
                                    double& rDetJac)const
{
    /*       jacobian
           j0, j2,
           j1, j3 */

    assert((int)rDerivativeShapeFunctions.size()==2*GetNumNodesGeometry() && (int)rNodeCoordinates.size()==2*GetNumNodesGeometry());
    double  j0(0.),j1(0.),j2(0.),j3(0.),x,y;

    int theDerivative(0);
    for (int count = 0; count < GetNumNodesGeometry(); count++)
    {
        x = rNodeCoordinates[theDerivative];
        y = rNodeCoordinates[theDerivative+1];

        j0 += rDerivativeShapeFunctions[theDerivative] * x;
        j1 += rDerivativeShapeFunctions[theDerivative] * y;
        theDerivative++;

        j2 += rDerivativeShapeFunctions[theDerivative] * x;
        j3 += rDerivativeShapeFunctions[theDerivative] * y;
        theDerivative++;
    }

    rDetJac = (j0*j3-j1*j2);

    if (rDetJac==0)
        throw MechanicsException("[NuTo::Plane::CalculateJacobian] Determinant of the Jacobian is zero, no inversion possible.");

    if (rInvJacobian!=0)
    {
        double invDeterminant(1./rDetJac);
        rInvJacobian[0]=j3*invDeterminant;
        rInvJacobian[1]=-j2*invDeterminant;
        rInvJacobian[2]=-j1*invDeterminant;
        rInvJacobian[3]=j0*invDeterminant;
    }
}

//! @brief calculates the derivative of the shape functions with respect to global coordinates
//! @param std::vector<double>& rDerivativeShapeFunctions derivatives of the shape functions
//! @param rJacInv inverse of the Jacobian
//! @param rDerivativeShapeFunctionsGlobal derivaties of the shape functions with respect to global coordinates
//! size should already be correct, but can be checked with an assert
//! first all the directions for a single node, and then for the next node
void NuTo::Plane::CalculateDerivativeShapeFunctionsLocal(const std::vector<double>& rDerivativeShapeFunctionsNatural, const double rJacInv[4], std::vector<double>& rDerivativeShapeFunctionsLocal)const
{
    int numShapeFunctions = rDerivativeShapeFunctionsLocal.size()/2;
    assert(rDerivativeShapeFunctionsLocal.size()==(size_t)2*numShapeFunctions);
    assert(rDerivativeShapeFunctionsNatural.size()==(size_t)2*numShapeFunctions);
    for (int count=0; count<numShapeFunctions; count++)
    {
        int mul2count = 2*count;
        int mul2countplus1 = mul2count+1;
        rDerivativeShapeFunctionsLocal[mul2count] =
            rDerivativeShapeFunctionsNatural[mul2count]     *rJacInv[0]+
            rDerivativeShapeFunctionsNatural[mul2countplus1]*rJacInv[2];

        rDerivativeShapeFunctionsLocal[mul2countplus1] =
            rDerivativeShapeFunctionsNatural[mul2count]     *rJacInv[1]+
            rDerivativeShapeFunctionsNatural[mul2countplus1]*rJacInv[3];
    }
}

//! @brief returns the number of nodes in this element(field interpolation)
//! @return number of nodes
int NuTo::Plane::GetNumNodesField()const
{
    return GetNumNodesGeometry();
}

//! @brief returns a pointer to the i-th node of the element
//! @param local node number
//! @return pointer to the node
NuTo::NodeBase* NuTo::Plane::GetNodeField(int rLocalNodeNumber)
{
    return GetNode(rLocalNodeNumber);
}

//! @brief returns a pointer to the i-th node of the element
//! @param local node number
//! @return pointer to the node
const NuTo::NodeBase* NuTo::Plane::GetNodeField(int rLocalNodeNumber)const
{
    return GetNode(rLocalNodeNumber);
}

//! @brief returns the number of nodes in this element(nonlocal eq strain interpolation)
//! @return number of nodes
int NuTo::Plane::GetNumNodesNonlocalEqStrain()const
{
    return GetNumNodesGeometry();
}

//! @brief returns a pointer to the i-th node of the element
//! @param local node number
//! @return pointer to the node
NuTo::NodeBase* NuTo::Plane::GetNodeNonlocalEqStrain(int rLocalNodeNumber)
{
    return GetNode(rLocalNodeNumber);
}

//! @brief returns a pointer to the i-th node of the element
//! @param local node number
//! @return pointer to the node
const NuTo::NodeBase* NuTo::Plane::GetNodeNonlocalEqStrain(int rLocalNodeNumber)const
{
    return GetNode(rLocalNodeNumber);
}

//! @brief calculates the deformation gradient in 2D
//! @param rRerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param rLocalDisp local displacements
//! @param rDeformationGradient (return value)
void NuTo::Plane::CalculateDeformationGradient(const std::vector<double>& rDerivativeShapeFunctionsLocal,
        const std::vector<double>& rLocalDisp,
        DeformationGradient2D& rDeformationGradient)const
{
    assert((int)rLocalDisp.size()==2*GetNumNodesField() && (int)rDerivativeShapeFunctionsLocal.size()==2*GetNumNodesField());

    rDeformationGradient.mDeformationGradient[0] = 1.;
    rDeformationGradient.mDeformationGradient[1] = 0.;
    rDeformationGradient.mDeformationGradient[2] = 0.;
    rDeformationGradient.mDeformationGradient[3] = 1.;

    int theDisp(0);
    double dNdX,dNdY;
    for (int count=0; count<GetNumNodes(); count++)
    {
        dNdX = rDerivativeShapeFunctionsLocal[theDisp];
        dNdY = rDerivativeShapeFunctionsLocal[theDisp+1];

        rDeformationGradient.mDeformationGradient[0]+=rLocalDisp[theDisp]* dNdX;
        rDeformationGradient.mDeformationGradient[1]+=rLocalDisp[theDisp]* dNdY;
        theDisp++;
        rDeformationGradient.mDeformationGradient[2]+=rLocalDisp[theDisp]* dNdX;
        rDeformationGradient.mDeformationGradient[3]+=rLocalDisp[theDisp]* dNdY;
        theDisp++;
    }
}

//! @brief calculates the temperature gradient in 3D
//! @param rRerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param rTemp nodal temperatures
//! @param rTemperatureGradient (return value)
void NuTo::Plane::CalculateTemperatureGradient(const std::vector<double>& rDerivativeShapeFunctionsLocal,
        const std::vector<double>& rTemp,
        TemperatureGradient2D& rTemperatureGradient)const
{
    assert((int)rTemp.size()==GetNumNodes() && (int)rDerivativeShapeFunctionsLocal.size()==2*GetNumNodes());

    rTemperatureGradient.mTemperatureGradient[0] = 0.;
    rTemperatureGradient.mTemperatureGradient[1] = 0.;

    int theTemp(0);
    double dNdX,dNdY;
    for (int theNode=0; theNode<GetNumNodes(); theNode++)
    {
        dNdX = rDerivativeShapeFunctionsLocal[theTemp];
        theTemp++;
        dNdY = rDerivativeShapeFunctionsLocal[theTemp];
        theTemp++;

        rTemperatureGradient.mTemperatureGradient[0]+=rTemp[theNode]* dNdX;
        rTemperatureGradient.mTemperatureGradient[1]+=rTemp[theNode]* dNdY;
    }
}


//! @brief stores the temperatures of the nodes
//! @param temperature vector with already correct size allocated
//! this can be checked with an assertation
void NuTo::Plane::CalculateTemperatures(std::vector<double>& rTemperatures)const
{
    assert((int)rTemperatures.size()==GetNumNodesField());
    for (int count=0; count<GetNumNodesField(); count++)
    {
        if (GetNode(count)->GetNumTemperatures()!=1)
            throw MechanicsException("[NuTo::Solid::CalculateTemperatures] Temperature is required as input to the constitutive model, but the node does not have this data.");
        rTemperatures[count] = GetNode(count)->GetTemperature();
    }
}

//! @brief calculates the coefficient matrix for the 1-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the damping matrix
/*NuTo::Error::eError NuTo::Plane::CalculateCoefficientMatrix_1(NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult,
        std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn, bool& rSymmetry)const
{
    throw MechanicsException("[NuTo::Plane::CalculateCoefficientMatrix_1] to be implemented.");
}
*/

//! @brief calculates the coefficient matrix for the 2-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the Mass matrix
/*NuTo::Error::eError NuTo::Plane::CalculateCoefficientMatrix_2(NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult,
        std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn, bool& rSymmetry)const
{
   //calculate local coordinates
    std::vector<double> localNodeCoord(GetNumLocalDofs());
    CalculateLocalCoordinates(localNodeCoord);

    //allocate space for ip coordinates in natural coordinate system (-1,1)
    double naturalIPCoord[2];

    //allocate space for derivatives of shape functions in natural coordinate system
    std::vector<double> derivativeShapeFunctionsNatural(2*GetNumShapeFunctions());
    //allocate space for shape functions in natural coordinate system
    std::vector<double> shapeFunctions(GetNumShapeFunctions());

    //InvJacobian and determinant of Jacobian
    double invJacobian[4], detJac;

    //material pointer
    const ConstitutiveEngineeringStressStrain *constitutivePtr;

    //allocate and initialize result matrix
    rResult.Resize(GetNumLocalDofs(),GetNumLocalDofs());
    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, naturalIPCoord);

        CalculateShapeFunctions(naturalIPCoord, shapeFunctions);

        CalculateDerivativeShapeFunctionsNatural(naturalIPCoord, derivativeShapeFunctionsNatural);

        CalculateJacobian(derivativeShapeFunctionsNatural,localNodeCoord, invJacobian, detJac);

        //call constitutive law and calculate stress
        constitutivePtr = GetConstitutiveLaw(theIP)->AsConstitutiveEngineeringStressStrain();
        if (constitutivePtr==0)
            throw MechanicsException("[NuTo::Plane::CalculateGradientInternalPotential] Constitutive law can not deal with engineering stresses and strains");

        //add to local mass vector
        assert(mSection->GetThickness()>0);

        //this global weight function is 1 for all standard problems, it's not equal to one for the multiscale method
        double factor(mSection->GetThickness()*detJac*(mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP))*constitutivePtr->GetDensity());
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> tmpMatrix;
        tmpMatrix = (factor*Eigen::Matrix<double,Eigen::Dynamic,1>::Map(&(shapeFunctions[0]),GetNumShapeFunctions()))*Eigen::Matrix<double,1,Eigen::Dynamic>::Map(&(shapeFunctions[0]),GetNumShapeFunctions());
        for (int count=0; count<GetNumShapeFunctions(); count++)
        {
            for (int count2=0; count2<GetNumShapeFunctions(); count2++)
            {
                rResult.mEigenMatrix(2*count,2*count2) += tmpMatrix(count,count2);
                rResult.mEigenMatrix(2*count+1,2*count2+1) += tmpMatrix(count,count2);
            }
        }
    }
    rSymmetry = true;
    // calculate list of global dofs related to the entries in the element stiffness matrix
    this->CalculateGlobalRowDofs(rGlobalDofsRow);
    rGlobalDofsColumn = rGlobalDofsRow;

    return Error::SUCCESSFUL;
}
*/

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point
//! @param rCoordinates local coordinates (return value)
void NuTo::Plane::GetLocalIntegrationPointCoordinates(int rIpNum, double rCoordinates[2])const
{
    this->mElementData->GetIntegrationType()->GetLocalIntegrationPointCoordinates2D(rIpNum, rCoordinates);
    return;
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point
//! @param rCoordinates local coordinates (return value)
void  NuTo::Plane::GetGlobalIntegrationPointCoordinates(int rIpNum, double rCoordinates[3])const
{
    double naturalCoordinates[2];
    double nodeCoordinates[3];
    this->mElementData->GetIntegrationType()->GetLocalIntegrationPointCoordinates2D(rIpNum, naturalCoordinates);
    std::vector<double> shapeFunctions(GetNumNodes());
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
        if (nodePtr->GetNumCoordinates()==2)
            nodePtr->GetCoordinates2D(nodeCoordinates);
        else
            nodePtr->GetCoordinates3D(nodeCoordinates);
        for (int theCoordinate=0; theCoordinate<nodePtr->GetNumCoordinates(); theCoordinate++)
        {
            rCoordinates[theCoordinate]+=shapeFunctions[theNode]*nodeCoordinates[theCoordinate];
        }
    }
    return;
}

//! @brief calculates the integration point data with the current displacements applied
//! @param rIpDataType data type to be stored for each integration point
//! @param rIpData return value with dimension (dim of data type) x (numIp)
/*NuTo::Error::eError NuTo::Plane::GetIpData(NuTo::IpData::eIpStaticDataType rIpDataType, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rIpData)const
{
	//calculate local coordinates
    std::vector<double> localNodeCoord(GetNumLocalDofs());
    CalculateLocalCoordinates(localNodeCoord);

    //std::cout<< "localNodeCoord " << Eigen::Matrix<double,Eigen::Dynamic,1>::Map(&(localNodeCoord[0]),localNodeCoord.size(),1) << std::endl;

    //calculate local displacements
    std::vector<double> localNodeDisp(GetNumLocalDofs());
    CalculateLocalDisplacements(localNodeDisp);

    //allocate space for ip coordinates in natural coordinate system (-1,1)
    double naturalIPCoord[2];

    //allocate space for derivatives of shape functions in natural coordinate system
    std::vector<double> derivativeShapeFunctionsNatural(2*GetNumShapeFunctions());
    //allocate space for derivatives of shape functions in local coordinate system
    std::vector<double> derivativeShapeFunctionsLocal(2*GetNumShapeFunctions());

    //allocate deformation gradient
    DeformationGradient2D deformationGradient;

    //allocate global engineering stress
    EngineeringStrain3D engineeringStrain;

    //allocate global engineering stress
    EngineeringStress3D engineeringStress;

    //InvJacobian and determinant of Jacobian
    double invJacobian[4], detJac;

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
    case NuTo::IpData::ELASTIC_ENERGY:
    case NuTo::IpData::INTERNAL_ENERGY:
        rIpData.Resize(2,GetNumIntegrationPoints());
    break;
    default:
        throw MechanicsException("[NuTo::Plane::GetIpData] Ip data not implemented.");
    }

    //store the data
    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, naturalIPCoord);

        CalculateDerivativeShapeFunctionsNatural(naturalIPCoord, derivativeShapeFunctionsNatural);
        //std::cout<< "derivativeShapeFunctionsNatural " << Eigen::Matrix<double,Eigen::Dynamic,1>::Map(&(derivativeShapeFunctionsNatural[0]),derivativeShapeFunctionsNatural.size(),1) << std::endl;

        CalculateJacobian(derivativeShapeFunctionsNatural,localNodeCoord, invJacobian, detJac);
        //std::cout<< "invJacobian " << Eigen::Matrix<double,2,2>::Map(&(invJacobian[0]),2,2) << std::endl;

        CalculateDerivativeShapeFunctionsLocal(derivativeShapeFunctionsNatural,invJacobian,
                                                derivativeShapeFunctionsLocal);
        //std::cout<< "derivativeShapeFunctionsLocal " << Eigen::Matrix<double,Eigen::Dynamic,1>::Map(&(derivativeShapeFunctionsLocal[0]),derivativeShapeFunctionsLocal.size(),1) << std::endl;

        // determine deformation gradient from the local Displacements and the derivative of the shape functions
        // this is not included in the AddIpStiffness to avoid reallocation of the deformation gradient for each IP
        CalculateDeformationGradient(derivativeShapeFunctionsLocal, localNodeDisp, deformationGradient);

        //call constitutive law and calculate stress
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
        case NuTo::IpData::ELASTIC_ENERGY:
        {
            error = constitutivePtr->GetElasticEnergy_EngineeringStress_EngineeringStrain(this, theIP, deformationGradient, rIpData.mEigenMatrix(0,theIP));
            assert(mSection->GetThickness()>0);
            double factor(mSection->GetThickness()*detJac*(mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP)));
            rIpData.mEigenMatrix(1,theIP) = factor;
        }
        break;
        case NuTo::IpData::INTERNAL_ENERGY:
        {
            error = constitutivePtr->GetInternalEnergy_EngineeringStress_EngineeringStrain(this, theIP, deformationGradient,rIpData.mEigenMatrix(0,theIP));
            assert(mSection->GetThickness()>0);
            double factor(mSection->GetThickness()*detJac*(mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP)));
            rIpData.mEigenMatrix(1,theIP) = factor;
        }
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

//! @brief Allocates static data for an integration point of an element
//! @param rConstitutiveLaw constitutive law, which is called to allocate the static data object
//! actually, both - the element type and the constitutive law are required to determine the static data object actually required
NuTo::ConstitutiveStaticDataBase* NuTo::Plane::AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw)const
{
    return rConstitutiveLaw->AllocateStaticDataEngineeringStress_EngineeringStrain2D(this);
}

// interpolate geometry
void NuTo::Plane::InterpolateCoordinatesFrom2D(double rNaturalCoordinates[2], double rGlobalCoordinates[3]) const
{
    // calculate shape functions
    std::vector<double> ShapeFunctions(this->GetNumNodesGeometry());
    this->CalculateShapeFunctionsGeometry(rNaturalCoordinates, ShapeFunctions);

    // start interpolation
    rGlobalCoordinates[0] = 0.0;
    rGlobalCoordinates[1] = 0.0;
    rGlobalCoordinates[2] = 0.0;
    for (int NodeCount = 0; NodeCount < this->GetNumNodesGeometry(); NodeCount++)
    {
        // get node coordinate
        double NodeCoordinate[3];
        const NodeBase *nodePtr(GetNodeGeometry(NodeCount));
        if (nodePtr->GetNumCoordinates()==2)
            nodePtr->GetCoordinates2D(NodeCoordinate);
        else
            nodePtr->GetCoordinates3D(NodeCoordinate);

        // add node contribution
        for (int theCoordinate=0; theCoordinate<nodePtr->GetNumCoordinates(); theCoordinate++)
        {
            rGlobalCoordinates[theCoordinate] += ShapeFunctions[NodeCount] *  NodeCoordinate[theCoordinate];
        }
    }
}

// interpolate displacements
void NuTo::Plane::InterpolateDisplacementsFrom2D(int rTimeDerivative, double rNaturalCoordinates[2], double rGlobalDisplacements[3]) const
{
    // calculate shape functions
    std::vector<double> ShapeFunctions(this->GetNumNodesField());
    this->CalculateShapeFunctionsField(rNaturalCoordinates, ShapeFunctions);

    // start interpolation
    rGlobalDisplacements[0] = 0.0;
    rGlobalDisplacements[1] = 0.0;
    rGlobalDisplacements[2] = 0.0;
    for (int NodeCount = 0; NodeCount < this->GetNumNodesField(); NodeCount++)
    {
        // get node displacements
        double NodeDisplacement[3];
        const NodeBase *nodePtr(GetNode(NodeCount));
        if (nodePtr->GetNumDisplacements()==2)
            nodePtr->GetDisplacements2D(rTimeDerivative, NodeDisplacement);
        else
            nodePtr->GetDisplacements3D(rTimeDerivative, NodeDisplacement);

        // add node contribution
        for (int theDisplacement=0; theDisplacement<nodePtr->GetNumDisplacements(); theDisplacement++)
        {
            rGlobalDisplacements[theDisplacement] += ShapeFunctions[NodeCount] *  NodeDisplacement[theDisplacement];
        }
    }
}

// check element definition
void NuTo::Plane::CheckElement()
{
    // check nodes
    for (int nodeCount = 0; nodeCount < this->GetNumNodesGeometry(); nodeCount++)
    {
        int numCoordinates(GetNodeGeometry(nodeCount)->GetNumCoordinates());
        if (numCoordinates<2 || numCoordinates>3)
        {
            throw MechanicsException("[NuTo::Plane::CheckElement] invalid node type (check node definition for coordinates).");
        }
    }

    // check node ordering (element length must be positive) and for changing sign in jacobian determinant
    // calculate coordinates
    std::vector<double> nodeCoord(2*this->GetNumNodesGeometry());
    this->CalculateLocalCoordinates(nodeCoord);
    /*for (int count=0; count<GetNumNodesGeometry(); count++)
    {
        std::cout << "Node " << count+1 << " with coordinates " << nodeCoord[2*count]<<","
                << nodeCoord[2*count+1]<<std::endl;
    }*/

    // check number of integration points
    if (this->GetNumIntegrationPoints() < 1)
    {
        throw MechanicsException("[NuTo::Plane::CheckElement] invalid integration type.");
    }

    // check sign of the jacobian determinant of the first integration point
    double naturalIPCoord[2];
    this->GetLocalIntegrationPointCoordinates(0, naturalIPCoord);

    std::vector<double> derivativeShapeFunctionsNatural(2*GetNumNodesGeometry());
    this->CalculateDerivativeShapeFunctionsGeometryNatural(naturalIPCoord, derivativeShapeFunctionsNatural);

    double invJacobian[4], detJacobian;
    this->CalculateJacobian(derivativeShapeFunctionsNatural,nodeCoord, invJacobian, detJacobian);
    // reorder nodes if determinant is negative
    if (detJacobian < 0.0)
    {
        this->ReorderNodes();

        // recalculate node coordinates after renumbering
        this->CalculateLocalCoordinates(nodeCoord);
    }

    // check jacobian determinant for all integration points for positive sign and calculate element volume
    double volume = 0;
    for (int ipCount = 0; ipCount < this->GetNumIntegrationPoints(); ipCount++)
    {
        // calculate jacobian determinant
        this->GetLocalIntegrationPointCoordinates(ipCount, naturalIPCoord);
        this->CalculateDerivativeShapeFunctionsGeometryNatural(naturalIPCoord, derivativeShapeFunctionsNatural);
        this->CalculateJacobian(derivativeShapeFunctionsNatural,nodeCoord, invJacobian, detJacobian);
        //std::cout << "Jacobian " << detJacobian << std::endl;
        if (detJacobian <= 0)
        {
            std::cout << "jac " << detJacobian << std::endl;
            throw MechanicsException("[NuTo::Plane::CheckElement] element is not properly defined by this nodes (zero or negative jacobian determinant).");
        }
        volume += this->GetIntegrationPointWeight(ipCount) * detJacobian;
    }

    // check element volume
    if (volume < 1e-14)
    {
        throw MechanicsException("[NuTo::Plane::CheckElement] element with zero volume (check nodes).");
    }
}

//! @brief sets the section of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need a section
//! @param rSection pointer to section
//! @return pointer to constitutive law
void NuTo::Plane::SetSection(const SectionBase* rSection)
{
    mSection = rSection;
}

//! @brief returns a pointer to the section of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need a section
//! @return pointer to section
const NuTo::SectionBase* NuTo::Plane::GetSection()const
{
    return mSection;
}

//! @brief calculates the volume of an integration point (weight * detJac)
//! @param rVolume  vector for storage of the ip volumes (area in 2D)
void NuTo::Plane::GetIntegrationPointVolume(std::vector<double>& rVolume)const
{
   //calculate local coordinates
    std::vector<double> localNodeCoord(2.*GetNumNodesGeometry());
    CalculateLocalCoordinates(localNodeCoord);

    //allocate space for ip coordinates in natural coordinate system (-1,1)
    double naturalIPCoord[2];

    //allocate space for derivatives of shape functions in natural coordinate system
    std::vector<double> derivativeShapeFunctionsNatural(2*GetNumNodesGeometry());

    //determinant of Jacobian
    double detJac;

    rVolume.resize(GetNumIntegrationPoints());

    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, naturalIPCoord);

        CalculateDerivativeShapeFunctionsGeometryNatural(naturalIPCoord, derivativeShapeFunctionsNatural);

        CalculateJacobian(derivativeShapeFunctionsNatural,localNodeCoord, 0, detJac);

        //attention in 2D, this is just the area, but that is required for the nonlocal model
        rVolume[theIP] = detJac * mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP);
    }
}

//! @brief cast the base pointer to an ElementPlane, otherwise throws an exception
const NuTo::Plane* NuTo::Plane::AsPlane()const
{
    return this;
}

//! @brief cast the base pointer to an Plane, otherwise throws an exception
NuTo::Plane* NuTo::Plane::AsPlane()
{
    return this;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::Plane::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Plane::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Plane::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Plane::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Plane::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Plane::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Plane::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Plane" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementBase)
       & BOOST_SERIALIZATION_NVP(mSection);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Plane" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Plane)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::Plane)
#endif // ENABLE_SERIALIZATION

