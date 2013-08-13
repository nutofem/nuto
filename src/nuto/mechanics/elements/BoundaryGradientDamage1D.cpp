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
    		bool rEdgeRealBoundaryElement,
    		ElementData::eElementDataType rElementDataType,
    		IntegrationType::eIntegrationType rIntegrationType,
    		IpData::eIpDataType rIpDataType
    		) :
    		NuTo::ElementBase::ElementBase(rStructure, rElementDataType, rIntegrationType, rIpDataType)
{
	mNodes = rNodes;
	mEdgeRealBoundaryElement = rEdgeRealBoundaryElement;
	mRealBoundaryElement = rRealBoundaryElement;
}

//! @brief calculates output data fo the elmement
//! @param eOutput ... coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
//!                    @param updateStaticData (with DummyOutput), IPData, globalrow/column dofs etc.
NuTo::Error::eError NuTo::BoundaryGradientDamage1D::Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput)
{

	if (mStructure->GetHessianConstant(1)==false)
    	throw MechanicsException("[NuTo::BoundaryGradientDamage::Evaluate] only implemented for a constant Hessian for the first derivative (damping).");
    if (mStructure->GetHessianConstant(2)==false)
    	throw MechanicsException("[NuTo::BoundaryGradientDamage::Evaluate] only implemented for a constant Hessian for the second derivative (mass).");

    try
    {
    	//***********************************************************************************************
        //First calculate the relevant informations for the real boundary element on the actual boundary
    	//***********************************************************************************************

    	//this requires the calculation of stresses and the derivatives of stresses with respect to all dofs of the
    	//real boundary element
    	// get section information determining which input on the constitutive level should be used
		const SectionBase* sectionReal(mRealBoundaryElement->GetSection());
		if (sectionReal==0)
			throw MechanicsException("[NuTo::BoundaryGradientDamage::Evaluate] no section allocated for real boundary element.");

		//calculate coordinates
		int numCoordinatesReal(mRealBoundaryElement->GetNumShapeFunctions());
		std::vector<double> localNodeCoordReal(numCoordinatesReal);
		mRealBoundaryElement->CalculateLocalCoordinates(localNodeCoordReal);

		//calculate local displacements, velocities and accelerations
		//the difference between disp and dispdof is a problem where the displacements are fixed, but enter the constitutive equation
		//for example in a two stage problem, first solve mechanics, then thermal and so on
		int numDispReal(mRealBoundaryElement->GetNumShapeFunctions());
		int numDispDofsReal = sectionReal->GetIsDisplacementDof() ? numDispReal : 0;

		int numNonlocalTotalStrainReal(1); //only the value on the boundary is required, in 1D one value, in 2D all nodes on an edge
		int numNonlocalTotalStrainDofsReal(sectionReal->GetIsNonlocalTotalStrainDof() ? numNonlocalTotalStrainReal : 0);

		std::vector<double> localNodeDispReal;
		double nodeNonlocalTotalStrainReal;

		//calculate local displacements, velocities and accelerations
		if (numDispDofsReal>0)
		{
			localNodeDispReal.resize(numDispReal);
			mRealBoundaryElement->CalculateLocalDisplacements(0,localNodeDispReal);
		}
		if (numNonlocalTotalStrainDofsReal>0 || sectionReal->GetInputConstitutiveIsNonlocalTotalStrain())
		{
			if (mEdgeRealBoundaryElement==0)
				nodeNonlocalTotalStrainReal = (mRealBoundaryElement->GetNode(0))->GetNonlocalTotalStrain(0);
			else
				nodeNonlocalTotalStrainReal = (mRealBoundaryElement->GetNode(mRealBoundaryElement->GetNumNodes()-1))->GetNonlocalTotalStrain(0);
		}

		//allocate space for local ip coordinates
		double localIPCoordReal;

		//allocate space for local shape functions
		std::vector<double> derivativeShapeFunctionsNaturalReal(mRealBoundaryElement->GetLocalDimension()*mRealBoundaryElement->GetNumShapeFunctions());  //allocate space for derivatives of shape functions
		std::vector<double> derivativeShapeFunctionsLocalReal(mRealBoundaryElement->GetLocalDimension()*mRealBoundaryElement->GetNumShapeFunctions());    //allocate space for derivatives of shape functions
		std::vector<double> shapeFunctionsReal(mRealBoundaryElement->GetNumShapeFunctions());                                       //allocate space for derivatives of shape functions

		//allocate deformation gradient
		DeformationGradient1D deformationGradientReal;

		//allocate global engineering plastic strain
		EngineeringStrain3D engineeringPlasticStrain3DReal;

		//allocate  damage (output of constitutive relation)
		Damage damageReal;

		//allocate nonlocal eq plastic strain (nodal dof value, input of constitutive relation)
		EngineeringStrain1D nonlocalTotalStrainReal;
		EngineeringStrain1D localTotalStrainReal;
		EngineeringStrain3D localTotalStrain3DReal;

		//allocate global engineering stress
		EngineeringStress1D engineeringStress1DReal;
		EngineeringStress3D engineeringStress3DReal;

		//allocate tangents
		ConstitutiveTangentLocal<1,1> tangentStressRealStrainReal;
		ConstitutiveTangentLocal<1,1> tangentStressRealNonlocalTotalStrainReal;

		//define inputs and outputs
		std::map< NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase* > constitutiveInputListReal;
		std::map< NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase* > constitutiveOutputListReal;

		if (numDispDofsReal>0)
		{
			constitutiveInputListReal[NuTo::Constitutive::Input::DEFORMATION_GRADIENT_1D] = &deformationGradientReal;
		}

		if (numNonlocalTotalStrainDofsReal>0)
		{
			constitutiveInputListReal[NuTo::Constitutive::Input::NONLOCAL_TOTAL_STRAIN_1D] = &nonlocalTotalStrainReal;
		}
		constitutiveOutputListReal[NuTo::Constitutive::Output::ENGINEERING_STRESS_1D] = &engineeringStress1DReal;

    	//*****************************************************************************************************
        //Now calculate the relevant informations for the virtual boundary element within the element boundary (nodal values
    	//*****************************************************************************************************
		//calculate coordinates
		int numCoordinatesVirt(this->GetNumShapeFunctions());
		std::vector<double> localNodeCoordVirt(numCoordinatesVirt);
        this->CalculateLocalCoordinates(localNodeCoordVirt);

		//calculate local displacements, velocities and accelerations
		//the difference between disp and dispdof is a problem where the displacements are fixed, but enter the constitutive equation
		//for example in a two stage problem, first solve mechanics, then thermal and so on
		int numNonlocalTotalStrainVirt(this->GetNumShapeFunctions());
		int numNonlocalTotalStrainDofsVirt(sectionReal->GetIsNonlocalTotalStrainDof() ? numNonlocalTotalStrainVirt : 0);

		std::vector<double> nodeNonlocalTotalStrainVirt;

		if (numNonlocalTotalStrainDofsVirt>0 || sectionReal->GetInputConstitutiveIsNonlocalTotalStrain())
		{
			nodeNonlocalTotalStrainVirt.resize(numNonlocalTotalStrainVirt);
            this->CalculateNodalNonlocalTotalStrain(0,nodeNonlocalTotalStrainVirt);
		}

		//allocate space for local ip coordinates
		double localIPCoordVirt;

		//allocate space for local shape functions
		std::vector<double> derivativeShapeFunctionsNaturalVirt(this->GetLocalDimension()*this->GetNumShapeFunctions());  //allocate space for derivatives of shape functions
		std::vector<double> derivativeShapeFunctionsLocalVirt(this->GetLocalDimension()*this->GetNumShapeFunctions());    //allocate space for derivatives of shape functions
		std::vector<double> shapeFunctionsVirt(this->GetNumShapeFunctions());                                             //allocate space for derivatives of shape functions

		//allocate nonlocal eq plastic strain (nodal dof value, input of constitutive relation)
		EngineeringStrain1D nonlocalTotalStrainVirt;
		EngineeringStrain1D localTotalStrainVirt;

		EngineeringStrain3D localTotalStrain3DVirt;
		EngineeringStrain3D engineeringPlasticStrain3DVirt;

		Damage damageVirt;

		//allocate tangents
		ConstitutiveTangentLocal<1,1> tangentStrainVirtStressReal;
		ConstitutiveTangentLocal<1,1> tangentStrainVirtNonlocalTotalStrainVirt;

		//define inputs and outputs
		std::map< NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase* > constitutiveInputListVirt;
		std::map< NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase* > constitutiveOutputListVirt;

		if (numNonlocalTotalStrainDofsVirt>0)
		{
			constitutiveInputListVirt[NuTo::Constitutive::Input::ENGINEERING_STRESS_1D] = &engineeringStress1DReal;
			constitutiveInputListVirt[NuTo::Constitutive::Input::NONLOCAL_TOTAL_STRAIN_1D] = &nonlocalTotalStrainVirt;
		}

		//define outputs for the element output and the constitutive law
		for (auto it = rElementOutput.begin(); it!=rElementOutput.end(); it++)
		{
			switch(it->first)
			{
			case Element::INTERNAL_GRADIENT:
				it->second->GetFullVectorDouble().Resize(numNonlocalTotalStrainDofsVirt);
				//if the stiffness matrix is constant, the corresponding internal force is calculated via the Kd
				//on the global level
				if (mStructure->GetHessianConstant(0)==false)
				{
					if (numNonlocalTotalStrainDofsVirt>0)
					{
						constitutiveOutputListVirt[NuTo::Constitutive::Output::ENGINEERING_STRAIN_1D] = &localTotalStrainVirt;
					}
				}
			break;
			case Element::HESSIAN_0_TIME_DERIVATIVE:
				{
					it->second->GetFullMatrixDouble().Resize(numNonlocalTotalStrainDofsVirt,numNonlocalTotalStrainDofsVirt + numDispDofsReal);
					it->second->SetSymmetry(false);
					it->second->SetConstant(false);
					if (numDispDofsReal>0)
					{
						constitutiveOutputListReal[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_1D] = &tangentStressRealStrainReal;
						//mixed terms
						if (numNonlocalTotalStrainDofsReal>0)
						{
							constitutiveOutputListReal[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_TOTAL_STRAIN_1D] = &tangentStressRealNonlocalTotalStrainReal;
						}
					}
					if (numNonlocalTotalStrainDofsVirt>0)
					{
						constitutiveOutputListVirt[NuTo::Constitutive::Output::D_ENGINEERING_STRAIN_VIRT_D_STRESS_REAL_1D] = &tangentStrainVirtStressReal;
						constitutiveOutputListVirt[NuTo::Constitutive::Output::D_ENGINEERING_STRAIN_VIRT_D_NONLOCAL_TOTAL_STRAIN_VIRT_1D] = &tangentStrainVirtNonlocalTotalStrainVirt;
					}
				}
			break;
			case Element::HESSIAN_1_TIME_DERIVATIVE:
			{
				it->second->GetFullMatrixDouble().Resize(numNonlocalTotalStrainDofsVirt,numNonlocalTotalStrainDofsVirt + numDispDofsReal+numNonlocalTotalStrainDofsReal);
				it->second->SetSymmetry(true);
				it->second->SetConstant(true);
				// Rayleigh damping should be introduced on the global level
			}
			break;
			case Element::HESSIAN_2_TIME_DERIVATIVE:
			{
				it->second->GetFullMatrixDouble().Resize(numNonlocalTotalStrainDofsVirt,numNonlocalTotalStrainDofsVirt + numDispDofsReal+numNonlocalTotalStrainDofsReal);
				it->second->SetSymmetry(true);
				it->second->SetConstant(true);
				//there is only a constant mass part for the mechanics problem
			}
			break;
			case Element::UPDATE_STATIC_DATA:
				constitutiveOutputListReal[NuTo::Constitutive::Output::UPDATE_STATIC_DATA] = 0;
				constitutiveOutputListVirt[NuTo::Constitutive::Output::UPDATE_STATIC_DATA] = 0;
			break;
			case Element::UPDATE_TMP_STATIC_DATA:
				constitutiveOutputListReal[NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA] = 0;
				constitutiveOutputListVirt[NuTo::Constitutive::Output::UPDATE_STATIC_DATA] = 0;
			break;
			case Element::IP_DATA:
				switch(it->second->GetIpDataType())
				{
				case NuTo::IpData::ENGINEERING_STRAIN:
					it->second->GetFullMatrixDouble().Resize(6,GetNumIntegrationPoints());
					 //define outputs
					constitutiveOutputListReal[NuTo::Constitutive::Output::ENGINEERING_STRAIN_3D] = &localTotalStrain3DReal;
					constitutiveOutputListVirt[NuTo::Constitutive::Output::ENGINEERING_STRAIN_3D] = &localTotalStrain3DVirt;
				break;
				case NuTo::IpData::ENGINEERING_STRESS:
					it->second->GetFullMatrixDouble().Resize(6,GetNumIntegrationPoints());
					 //define outputs
					constitutiveOutputListReal[NuTo::Constitutive::Output::ENGINEERING_STRESS_3D] = &engineeringStress3DReal;
					//constitutiveOutputListVirt[NuTo::Constitutive::Output::ENGINEERING_STRESS_3D] = &engineeringStress3DReal;
				break;
				case NuTo::IpData::ENGINEERING_PLASTIC_STRAIN:
					it->second->GetFullMatrixDouble().Resize(6,GetNumIntegrationPoints());
					 //define outputs
					 constitutiveOutputListReal[NuTo::Constitutive::Output::ENGINEERING_PLASTIC_STRAIN_3D] = &engineeringPlasticStrain3DReal;
					 constitutiveOutputListVirt[NuTo::Constitutive::Output::ENGINEERING_PLASTIC_STRAIN_3D] = &engineeringPlasticStrain3DVirt;
					 break;
				break;
				case NuTo::IpData::DAMAGE:
					it->second->GetFullMatrixDouble().Resize(1,GetNumIntegrationPoints());
					//define outputs
					constitutiveOutputListReal[NuTo::Constitutive::Output::DAMAGE] = &damageReal;
					constitutiveOutputListVirt[NuTo::Constitutive::Output::DAMAGE] = &damageVirt;
				break;
				default:
					throw MechanicsException("[NuTo::BoundaryGradientDamage::Evaluate1D] this ip data type is not implemented.");
				}
			break;
			case Element::GLOBAL_ROW_DOF:
				this->CalculateGlobalRowDofs(it->second->GetVectorInt(),numNonlocalTotalStrainDofsVirt);
			break;
			case Element::GLOBAL_COLUMN_DOF:
				this->CalculateGlobalColumnDofs(it->second->GetVectorInt(),numNonlocalTotalStrainDofsVirt,numDispDofsReal);
			break;
			default:
				throw MechanicsException("[NuTo::BoundaryGradientDamage::Evaluate1D] element output not implemented.");
			}
		}

		// loop over the integration points along the boundary (this is a point in 1D, a line in 2D or a surface in 3D)
		for (int theIP=0; theIP<this->GetNumIntegrationPoints(); theIP++)
		{
			//the integration points are sorted the way that the ones with the zero weight refer to the real boundary element (on their surface)
			//with the stresses used for the next integration points (with nonzero weight) until the next ip with zero weight is obtained
			double weight = mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP);

			if (weight==0.)
			{
				//the ip is located on the boundary of the real element
                this->GetLocalIntegrationPointCoordinatesReal(localIPCoordReal);

				//derivative in natural coordinate system
				mRealBoundaryElement->CalculateDerivativeShapeFunctions(localIPCoordReal, derivativeShapeFunctionsNaturalReal);

				//determinant of the Jacobian
				double detJReal(mRealBoundaryElement->DetJacobian(derivativeShapeFunctionsNaturalReal,localNodeCoordReal));

				//derivative in local coordinate system
				for (unsigned int count=0; count<derivativeShapeFunctionsNaturalReal.size(); count++)
				{
					derivativeShapeFunctionsLocalReal[count] = derivativeShapeFunctionsNaturalReal[count]/detJReal;
				}

				if (numDispDofsReal)
				{
					// determine deformation gradient from the local Displacements and the derivative of the shape functions
					mRealBoundaryElement->CalculateDeformationGradient(derivativeShapeFunctionsLocalReal, localNodeCoordReal, localNodeDispReal, deformationGradientReal);
				}

				if (numNonlocalTotalStrainDofsReal)
				{
					//mRealBoundaryElement->CalculateShapeFunctions(localIPCoordReal, shapeFunctionsReal);
					//mRealBoundaryElement->CalculateNonlocalTotalStrain(shapeFunctionsReal, nodeNonlocalTotalStrainReal, nonlocalTotalStrainReal);
					nonlocalTotalStrainReal(0) = nodeNonlocalTotalStrainReal;
				}

				ConstitutiveBase* constitutivePtr = this->GetConstitutiveLaw(theIP);
				Error::eError error = constitutivePtr->Evaluate1D(this, theIP,
						constitutiveInputListReal, constitutiveOutputListReal);
				if (error!=Error::SUCCESSFUL)
					return error;
				for (auto it = rElementOutput.begin(); it!=rElementOutput.end(); it++)
				{
					switch(it->first)
					{
						case Element::IP_DATA:
							switch (it->second->GetIpDataType())
							{
							case NuTo::IpData::ENGINEERING_STRAIN:
								//error = constitutivePtr->GetEngineeringStrain(this, theIP, deformationGradient, engineeringStrain);
								memcpy(&(it->second->GetFullMatrixDouble().data()[theIP*6]),localTotalStrain3DReal.GetData(),6*sizeof(double));
							break;
							case NuTo::IpData::ENGINEERING_STRESS:
								//error = constitutivePtr->GetEngineeringStressFromEngineeringStrain(this, theIP, deformationGradient, engineeringStress);
								memcpy(&(it->second->GetFullMatrixDouble().data()[theIP*6]),engineeringStress3DReal.GetData(),6*sizeof(double));
							break;
							case NuTo::IpData::ENGINEERING_PLASTIC_STRAIN:
								//error = constitutivePtr->GetEngineeringPlasticStrain(this, theIP, deformationGradient, engineeringStrain);
								memcpy(&(it->second->GetFullMatrixDouble().data()[theIP*6]),engineeringPlasticStrain3DReal.GetData(),6*sizeof(double));
							break;
							case NuTo::IpData::DAMAGE:
								//error = constitutivePtr->GetDamage(this, theIP, deformationGradient, rIpData.mEigenMatrix.data()[theIP]);
								memcpy(&(it->second->GetFullMatrixDouble().data()[theIP]),damageReal.GetData(),sizeof(double));
							break;
							default:
								throw MechanicsException("[NuTo::BoundaryGradientDamage::Evaluate1D] Ip data not implemented.");
							}
						break;
						default:
                        break;
					}
				}
			}
			else
			{
				//the ip is located within the virtual boundary element
				this->GetLocalIntegrationPointCoordinatesVirt(theIP, localIPCoordVirt);

				//derivative in natural coordinate system
				this->CalculateDerivativeShapeFunctions(localIPCoordVirt, derivativeShapeFunctionsNaturalVirt);

				//determinant of the Jacobian
				double detJVirt(this->DetJacobian(derivativeShapeFunctionsNaturalVirt,localNodeCoordVirt));

				//derivative in local coordinate system
				for (unsigned int count=0; count<derivativeShapeFunctionsNaturalVirt.size(); count++)
				{
					derivativeShapeFunctionsLocalVirt[count] = derivativeShapeFunctionsNaturalVirt[count]/detJVirt;
				}

				if (numNonlocalTotalStrainDofsVirt)
				{
					this->CalculateShapeFunctions(localIPCoordVirt, shapeFunctionsVirt);
					this->CalculateNonlocalTotalStrain(shapeFunctionsVirt, nodeNonlocalTotalStrainVirt, nonlocalTotalStrainVirt);
				}

				ConstitutiveBase* constitutivePtr = GetConstitutiveLaw(theIP);
				Error::eError error = constitutivePtr->Evaluate1D(this, theIP,
						constitutiveInputListVirt, constitutiveOutputListVirt);
				if (error!=Error::SUCCESSFUL)
					return error;

				double factorVirt (fabs(detJVirt)*sectionReal->GetArea()*
						       (mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP)));


				FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> KkkVirt;
				if (numNonlocalTotalStrainDofsVirt>0)
				{
					//calculate Kkk detJ*(cBtB+NtN)
					//the nonlocal radius is in a gradient formulation is different from the nonlocal radius in an integral formulation
					CalculateKkk(shapeFunctionsVirt,derivativeShapeFunctionsLocalVirt,constitutivePtr->GetNonlocalRadius(),factorVirt,KkkVirt);
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
							if (numNonlocalTotalStrainDofsVirt>0)
							{
								//add Kkk*nonlocalTotalStrain-detJ*F
								AddDetJRnonlocalTotalStrain(shapeFunctionsVirt,localTotalStrainVirt, KkkVirt, nodeNonlocalTotalStrainVirt, factorVirt, it->second->GetFullVectorDouble());
							}
						}
					}
					break;
					case Element::HESSIAN_0_TIME_DERIVATIVE:
						{
							if (numNonlocalTotalStrainDofsVirt>0)
							{
								//derivative of R(nonlocaltotalStrain) with respect to all unknowns
								it->second->GetFullMatrixDouble().AddBlock(0,0,KkkVirt);
								Sub_K_epsilonNonlocalVirt_NonlocalVirt
								   (factorVirt,shapeFunctionsVirt,tangentStrainVirtNonlocalTotalStrainVirt,
									it->second->GetFullMatrixDouble());

								if (numDispDofsReal>0)
								{
									Sub_K_epsilonNonlocalVirt_DispReal
									   (factorVirt,shapeFunctionsVirt,tangentStrainVirtStressReal,tangentStressRealStrainReal,
										derivativeShapeFunctionsLocalReal,
										numNonlocalTotalStrainDofsVirt,it->second->GetFullMatrixDouble());

								}
								if (numNonlocalTotalStrainDofsReal>0)
								{
									int col(0);

									// subtract detJ Nt dEpsilonVirt_dsigmaReal dSigmaReal_dEpsilonReal Nreal (NReal=1)
									//std::cout << "before \n" << it->second->GetFullMatrixDouble() << std::endl;
									Sub_K_epsilonNonlocalVirt_NonlocalStrainReal
									   (factorVirt,shapeFunctionsVirt,tangentStrainVirtStressReal,
									    tangentStressRealNonlocalTotalStrainReal,
										col,
										it->second->GetFullMatrixDouble());
									//std::cout << "after \n" << it->second->GetFullMatrixDouble() << std::endl;
								}
							}
							//BlowLocalMatrixToGlobal(it->second->GetFullMatrixDouble());
						}
					break;
					case Element::HESSIAN_1_TIME_DERIVATIVE:
					{
/*						if (numDispDofs>0)
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
						//BlowLocalMatrixToGlobal(it->second->GetFullMatrixDouble());
*/
					}
					break;
					case Element::HESSIAN_2_TIME_DERIVATIVE:
					{
/*
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
						*/
					}
					break;
					case Element::UPDATE_STATIC_DATA:
					case Element::UPDATE_TMP_STATIC_DATA:
					{
						//this is just a temporary check remove it afterwards
						localIPCoordVirt = 1;

						//derivative in natural coordinate system
						this->CalculateDerivativeShapeFunctions(localIPCoordVirt, derivativeShapeFunctionsNaturalVirt);

						//determinant of the Jacobian
						double detJVirt(this->DetJacobian(derivativeShapeFunctionsNaturalVirt,localNodeCoordVirt));

						//derivative in local coordinate system
						for (unsigned int count=0; count<derivativeShapeFunctionsNaturalVirt.size(); count++)
						{
							derivativeShapeFunctionsLocalVirt[count] = derivativeShapeFunctionsNaturalVirt[count]/detJVirt;
						}

						if (numNonlocalTotalStrainDofsVirt)
						{
							this->CalculateShapeFunctions(localIPCoordVirt, shapeFunctionsVirt);
							this->CalculateNonlocalTotalStrain(shapeFunctionsVirt, nodeNonlocalTotalStrainVirt, nonlocalTotalStrainVirt);
						}
						//Calculate derivative of nonlocal total strain at the outer boundary
						double derNonlocalTotalStrain(0);
						for (int count=0; count<derivativeShapeFunctionsLocalVirt.size(); count++)
							derNonlocalTotalStrain+= derivativeShapeFunctionsLocalVirt[count] * nodeNonlocalTotalStrainVirt[count];
						std::cout << "derNonlocalTotalStrain = " <<  derNonlocalTotalStrain << std::endl;
					}
					break;
					case Element::IP_DATA:
						switch (it->second->GetIpDataType())
						{
						case NuTo::IpData::ENGINEERING_STRAIN:
							//error = constitutivePtr->GetEngineeringStrain(this, theIP, deformationGradient, engineeringStrain);
							memcpy(&(it->second->GetFullMatrixDouble().data()[theIP*6]),localTotalStrain3DVirt.GetData(),6*sizeof(double));
						break;
						case NuTo::IpData::ENGINEERING_STRESS:
							//error = constitutivePtr->GetEngineeringStressFromEngineeringStrain(this, theIP, deformationGradient, engineeringStress);
							memcpy(&(it->second->GetFullMatrixDouble().data()[theIP*6]),engineeringStress3DReal.GetData(),6*sizeof(double));
						break;
						case NuTo::IpData::ENGINEERING_PLASTIC_STRAIN:
							//error = constitutivePtr->GetEngineeringPlasticStrain(this, theIP, deformationGradient, engineeringStrain);
							memcpy(&(it->second->GetFullMatrixDouble().data()[theIP*6]),engineeringPlasticStrain3DVirt.GetData(),6*sizeof(double));
						break;
						case NuTo::IpData::DAMAGE:
							//error = constitutivePtr->GetDamage(this, theIP, deformationGradient, rIpData.mEigenMatrix.data()[theIP]);
							memcpy(&(it->second->GetFullMatrixDouble().data()[theIP]),damageVirt.GetData(),sizeof(double));
						break;
						default:
							throw MechanicsException("[NuTo::BoundaryGradientDamage::Evaluate1D] Ip data not implemented.");
						}
					break;
					case Element::GLOBAL_ROW_DOF:
					case Element::GLOBAL_COLUMN_DOF:

					break;

					default:
						throw MechanicsException("[NuTo::BoundaryGradientDamage::Evaluate1D] element output not implemented.");
					}
				}//auto it = rElementOutput.begin(); it!=rElementOutput.end(); it++
			}
		}
    }
    catch (NuTo::MechanicsException e)
    {
        std::stringstream ss;
        ss << mStructure->ElementGetId(this);
    	e.AddMessage("[NuTo::BoundaryGradientDamage::Evaluate] Error evaluating element data of element"	+ ss.str() + ".");
        throw e;
    }

    return Error::SUCCESSFUL;
}

// interpolate geometry
void NuTo::BoundaryGradientDamage1D::InterpolateCoordinatesFrom1D(double rLocalCoordinates, double rGlobalCoordinates[3]) const
{
    // calculate shape functions
    std::vector<double> ShapeFunctions(this->GetNumNodes());
    this->CalculateShapeFunctions(rLocalCoordinates, ShapeFunctions);

    // start interpolation
    rGlobalCoordinates[0] = 0.0;
    rGlobalCoordinates[1] = 0.0;
    rGlobalCoordinates[2] = 0.0;
    for (int theNode = 0; theNode < this->GetNumNodes(); theNode++)
    {
        // get node coordinate
        double NodeCoordinate;
        GetNode(theNode)->GetCoordinates1D(&NodeCoordinate);

        // add node contribution
        rGlobalCoordinates[0] += ShapeFunctions[theNode] *  NodeCoordinate;
    }
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point
//! @param rCoordinates local coordinates (return value)
void  NuTo::BoundaryGradientDamage1D::GetGlobalIntegrationPointCoordinates(int rIpNum, double rCoordinates[3])const
{
    double naturalCoordinates;
    double nodeCoordinates[3];
    std::vector<double> shapeFunctions(GetNumNodes());
    GetLocalIntegrationPointCoordinatesVirt(rIpNum, naturalCoordinates);
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

//! @brief calculates the local coordinates of the nodes
//! @param localCoordinates vector with already correct size allocated
//! this can be checked with an assertation
void NuTo::BoundaryGradientDamage1D::CalculateLocalCoordinates(std::vector<double>& rLocalCoordinates)const
{
    assert((int)rLocalCoordinates.size()==GetNumNodes());
    for (int theNode=0; theNode<GetNumNodes(); theNode++)
    {
        this->GetNode(theNode)->GetCoordinates1D(&(rLocalCoordinates[theNode]));
    }
}

//! @brief returns the local coordinates of an integration point on the boundary (in the real boundary element)
//! @param rCoordinates local coordinates (return value)
void NuTo::BoundaryGradientDamage1D::GetLocalIntegrationPointCoordinatesReal(double& rCoordinates)const
{
    if (mEdgeRealBoundaryElement==0)
    	rCoordinates = -1.;
    else
    	rCoordinates = 1.;
    return;
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point
//! @param rCoordinates local coordinates (return value)
void NuTo::BoundaryGradientDamage1D::GetLocalIntegrationPointCoordinatesVirt(int rIpNum, double& rCoordinates)const
{
    this->mElementData->GetIntegrationType()->GetLocalIntegrationPointCoordinates1D(rIpNum, rCoordinates);
    return;
}

//! @brief Allocates static data for an integration point of an element
//! @param rConstitutiveLaw constitutive law, which is called to allocate the static data object
NuTo::ConstitutiveStaticDataBase* NuTo::BoundaryGradientDamage1D::AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw)const
{
	return rConstitutiveLaw->AllocateStaticDataEngineeringStress_EngineeringStrain1D(this);
}


//! @brief stores the nonlocal eq plastic strain of the nodes
//! @param time derivative (0 damage, 1 damage rate, 2 second time derivative of damage)
//! @param nonlocal eq plastic strain vector with already correct size allocated (2*nodes)
//! this can be checked with an assertation
void NuTo::BoundaryGradientDamage1D::CalculateNodalNonlocalEqPlasticStrain(int rTimeDerivative, std::vector<double>& rNodalNonlocalEqPlasticStrain)const
{
    assert(rTimeDerivative>=0);
    assert(rTimeDerivative<3);
	assert((int)rNodalNonlocalEqPlasticStrain.size()==2*GetNumShapeFunctions());
	double nonlocalEqPlasticStrain[2];
    for (int count=0; count<GetNumShapeFunctions(); count++)
    {
        if (GetNode(count)->GetNumNonlocalEqPlasticStrain()!=2)
            throw MechanicsException("[NuTo::BoundaryGradientDamage1D::CalculateNodalDamage] Damage is required as input to the constitutive model, but the node does not have this data.");
        GetNode(count)->GetNonlocalEqPlasticStrain(nonlocalEqPlasticStrain);
        rNodalNonlocalEqPlasticStrain[count] = nonlocalEqPlasticStrain[0];
        rNodalNonlocalEqPlasticStrain[count+GetNumShapeFunctions()] = nonlocalEqPlasticStrain[1];
    }
}

//! @brief stores the nonlocal total strain of the nodes
//! @param time derivative (0 damage, 1 damage rate, 2 second time derivative of damage)
//! @param nonlocal total strain vector with already correct size allocated (1*nodes)
//! this can be checked with an assertation
void NuTo::BoundaryGradientDamage1D::CalculateNodalNonlocalTotalStrain(int rTimeDerivative, std::vector<double>& rNodalNonlocalTotalStrain)const
{
    assert(rTimeDerivative>=0);
    assert(rTimeDerivative<3);
	assert((int)rNodalNonlocalTotalStrain.size()==GetNumShapeFunctions());
	double nonlocalTotalStrain[1];
    for (int count=0; count<GetNumShapeFunctions(); count++)
    {
        if (GetNode(count)->GetNumNonlocalTotalStrain()!=1)
            throw MechanicsException("[NuTo::BoundaryGradientDamage1D::CalculateNodalNonlocalStrain] nonlocal strain is required as input to the constitutive model, but the node does not have this data.");
        GetNode(count)->GetNonlocalTotalStrain1D(nonlocalTotalStrain);
        rNodalNonlocalTotalStrain[count] = nonlocalTotalStrain[0];
    }
}

//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::BoundaryGradientDamage1D::CalculateShapeFunctions(double rLocalCoordinates, std::vector<double>& rShapeFunctions)const
{
	assert(rShapeFunctions.size()==mNodes.size());
	switch (mNodes.size())
	{
	case 3:
		rShapeFunctions[0] = 0.5*(1.-rLocalCoordinates)-0.5*(1.-rLocalCoordinates*rLocalCoordinates);
		rShapeFunctions[1] = 1.-rLocalCoordinates*rLocalCoordinates;
		rShapeFunctions[2] = 0.5*(1.+rLocalCoordinates)-0.5*(1.-rLocalCoordinates*rLocalCoordinates);
		break;
	case 4:
	{
		double x2=rLocalCoordinates * rLocalCoordinates;
		double x3=rLocalCoordinates * x2;
		rShapeFunctions[0] = 1./16.*(-9.*x3+9.*x2+rLocalCoordinates-1.);
		rShapeFunctions[1] = 9./16.*( 3.*x3-x2-3.*rLocalCoordinates+1.);
		rShapeFunctions[2] = 9./16.*(-3.*x3-x2+3.*rLocalCoordinates+1.);
		rShapeFunctions[3] = 1./16.*( 9.*x3+9.*x2-rLocalCoordinates-1.);
	}
	break;
	case 5:
	{
		double x1=rLocalCoordinates;
		double x2=rLocalCoordinates * rLocalCoordinates;
		double x3=rLocalCoordinates * x2;
		double x4=rLocalCoordinates * x2;
		rShapeFunctions[0] =  1./6. * (0. + 1.*x1 -1. *x2 -4.*x3 + 4.*x4);
		rShapeFunctions[1] =  1./6. * (0. - 8.*x1 +16.*x2 +8.*x3 -16.*x4);
		rShapeFunctions[2] =  1./6. * (6. + 0.*x1 -30 *x2 -0.*x3 +24.*x4);
		rShapeFunctions[3] =  1./6. * (0. + 8.*x1 +16.*x2 -8.*x3 -16.*x4);
		rShapeFunctions[4] =  1./6. * (0. - 1.*x1 -1. *x2 +4.*x3 + 4.*x4);
	}
	break;
	default:
		throw MechanicsException("[NuTo::BoundaryGradientDamage1D::CalculateShapeFunctions] only implemented for 3 or 4 nodes.");
	}

}

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
void NuTo::BoundaryGradientDamage1D::CalculateDerivativeShapeFunctions(double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)const
{
	assert(rDerivativeShapeFunctions.size()==mNodes.size());
	switch (mNodes.size())
	{
	case 3:
		rDerivativeShapeFunctions[0] = -0.5 + rLocalCoordinates;
		rDerivativeShapeFunctions[1] = -2.0 * rLocalCoordinates;
		rDerivativeShapeFunctions[2] =  0.5 + rLocalCoordinates;
		break;
	case 4:
	{
		double x2=rLocalCoordinates * rLocalCoordinates;
		rDerivativeShapeFunctions[0] = 1./16.*(-27.*x2+18.*rLocalCoordinates+1.);
		rDerivativeShapeFunctions[1] = 9./16.*(  9.*x2- 2.*rLocalCoordinates-3.);
		rDerivativeShapeFunctions[2] = 9./16.*( -9.*x2- 2.*rLocalCoordinates+3.);
		rDerivativeShapeFunctions[3] = 1./16.*( 27.*x2+18.*rLocalCoordinates-1.);
	}
		break;
	case 5:
	{
		double x1=rLocalCoordinates;
		double x2=rLocalCoordinates * rLocalCoordinates;
		double x3=rLocalCoordinates * x2;
		rDerivativeShapeFunctions[0] = 1./6.*( 1. - 2.*x1 -12.*x2 +16.*x3);
		rDerivativeShapeFunctions[1] = 1./6.*(-8. +32.*x1 +24.*x2 -64.*x3);
		rDerivativeShapeFunctions[2] = 1./6.*( 0. -60.*x1 + 0.*x2 +96.*x3);
		rDerivativeShapeFunctions[3] = 1./6.*( 8. +32.*x1 -24.*x2 -64.*x3);
		rDerivativeShapeFunctions[4] = 1./6.*(-1. - 2.*x1 +12.*x2 +16.*x3);
	}
		break;
	default:
		throw MechanicsException("[NuTo::BoundaryGradientDamage1D::CalculateShapeFunctions] only implemented for 3 or 4 nodes.");
	}

}

//! @brief returns determinant of the Jacobian
//! @param derivativeShapeFunctions derivatives of the shape functions
//! @param localCoord local coordinates
//! @return determinant of the Jacobian
double NuTo::BoundaryGradientDamage1D::DetJacobian(const std::vector<double>& derivativeShapeFunctions,const std::vector<double>& localCoord)const
{
    assert((int)localCoord.size()==GetNumNodes() && (int)derivativeShapeFunctions.size()==GetNumNodes());
    double detJ(0);
    for (int count=0; count<GetNumNodes(); count++)
        detJ+=derivativeShapeFunctions[count]*localCoord[count];
    return detJ;
}

//! @brief returns the nonlocal eq plastic strain interpolated from the nodal values
//! @param shapeFunctionsGlobal shape functions
//! @param rNodeDamage nonlocal eq plastic strain values of the nodes
//! @param rNonlocalEqentPlasticStrain return value
void NuTo::BoundaryGradientDamage1D::CalculateNonlocalEqPlasticStrain(const std::vector<double>& shapeFunctions,
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
void NuTo::BoundaryGradientDamage1D::CalculateNonlocalTotalStrain(const std::vector<double>& shapeFunctions,
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

// build global row dofs
void NuTo::BoundaryGradientDamage1D::CalculateGlobalRowDofs(std::vector<int>& rGlobalRowDofs, int rNumNonlocalTotalStrainDofsVirt) const
{
    rGlobalRowDofs.resize(rNumNonlocalTotalStrainDofsVirt);
    int index(0);
    int numNodes(this->GetNumNodes());
    for (int nodeCount = 0; nodeCount < numNodes; nodeCount++)
    {
        const NodeBase * nodePtr(GetNode(nodeCount));

		for (int countNonlocalTotalStrain=0; countNonlocalTotalStrain<nodePtr->GetNumNonlocalTotalStrain(); countNonlocalTotalStrain++)
		{
			rGlobalRowDofs[index] = nodePtr->GetDofNonlocalTotalStrain(countNonlocalTotalStrain);
			index++;
        }
    }
}

// build global col dofs
void NuTo::BoundaryGradientDamage1D::CalculateGlobalColumnDofs(std::vector<int>& rGlobalColumnDofs, int rNumNonlocalTotalStrainDofsVirt, int rNumDispDofsReal) const
{
	rGlobalColumnDofs.resize(rNumNonlocalTotalStrainDofsVirt+rNumDispDofsReal);
    int indexNonlocalTotalStrainVirt(0);
    int numNodesVirt(this->GetNumNodes());
    for (int nodeCount = 0; nodeCount < numNodesVirt; nodeCount++)
    {
        const NodeBase * nodePtr(GetNode(nodeCount));

        for (int countNonlocalTotalStrain=0; countNonlocalTotalStrain<nodePtr->GetNumNonlocalTotalStrain(); countNonlocalTotalStrain++)
		{
			rGlobalColumnDofs[indexNonlocalTotalStrainVirt] = nodePtr->GetDofNonlocalTotalStrain(countNonlocalTotalStrain);
			indexNonlocalTotalStrainVirt++;
        }
    }

    int indexDispReal(0);
    int numNodesReal(mRealBoundaryElement->GetNumNodes());
    for (int nodeCount = 0; nodeCount < numNodesReal; nodeCount++)
    {
    	const NodeBase* theRealNode(mRealBoundaryElement->GetNode(nodeCount));

    	for (int countDisp=0; countDisp<theRealNode->GetNumDisplacements(); countDisp++)
		{
    		rGlobalColumnDofs[rNumNonlocalTotalStrainDofsVirt + indexDispReal] = theRealNode->GetDofDisplacement(countDisp);
			indexDispReal++;
		}
    }
}



//! @brief add Kkk*omega+detJ*F (detJ is already included in Koo)
//! @param rShapeFunctions of the ip for all shape functions
//! @param rLocalTotalStrain local total strain values
//! @param rKkk stiffness matrix Kkk
//! @param rNodeNonlocalTotalStrain nodal nonlocal total strain values
//! @param rFactor factor including detJ and area
//! @param rResult result
void NuTo::BoundaryGradientDamage1D::AddDetJRnonlocalTotalStrain(const std::vector<double>& rShapeFunctions,const EngineeringStrain1D& rLocalTotalStrain,
		const FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rKkk,
		const std::vector<double>& rNodeNonlocalTotalStrain, double rFactor, FullVector<double,Eigen::Dynamic>& rResult)
{
	assert(rResult.GetNumRows()>=(int)(rShapeFunctions.size()));
	assert(rShapeFunctions.size()==rNodeNonlocalTotalStrain.size());

	//first component
	FullVector<double,Eigen::Dynamic> tmpResult = rKkk*Eigen::Map<Eigen::VectorXd>(const_cast<double*>(&(rNodeNonlocalTotalStrain[0])),rShapeFunctions.size());
	double localTotalStrainMulFactor = rLocalTotalStrain(0)*rFactor;
	for (unsigned int count=0; count<rShapeFunctions.size(); count++)
	{
		tmpResult(count)-=rShapeFunctions[count]*localTotalStrainMulFactor;
	}
	rResult.segment(0,rShapeFunctions.size()) +=tmpResult;

}

void NuTo::BoundaryGradientDamage1D::Sub_K_epsilonNonlocalVirt_NonlocalVirt(
		double rFactorVirt,
		const std::vector<double>& rShapeFunctionsVirt,
		const FullMatrix<double,1,1>& tangentStrainVirtNonlocalTotalStrainVirt,
		FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult)
{
	assert(rResult.GetNumRows()==(int)(rShapeFunctionsVirt.size()));
	assert(rResult.GetNumColumns()>=(int)(rShapeFunctionsVirt.size()));

	double factor = rFactorVirt * tangentStrainVirtNonlocalTotalStrainVirt(0,0);

	for (unsigned int countRow=0; countRow<rShapeFunctionsVirt.size(); countRow++)
	{
		for (unsigned int countCol=0; countCol<rShapeFunctionsVirt.size(); countCol++)
		{
    		rResult(countRow,countCol)-=rShapeFunctionsVirt[countRow]*factor * rShapeFunctionsVirt[countCol];
		}
	}

}

//! @brief add Nt depsilonVirtdSigmaReal dSigmaRealdepsilonReal B
//! @param factorVirt factor
//! @param shapeFunctionsVirt
//! @param tangentStrainVirtStressReal
//! @param tangentStressRealStrainReal
//! @param derivativeShapeFunctionsLocalReal
//! @param rCol column where to add the submatrix
//! @param rResult return value (added)
void NuTo::BoundaryGradientDamage1D::Sub_K_epsilonNonlocalVirt_DispReal
  ( double factorVirt,
    const std::vector<double>& rShapeFunctionsVirt,
	const FullMatrix<double,1,1>& rTangentStrainVirtStressReal,
	const FullMatrix<double,1,1>& rTangentStressRealStrainReal,
	const std::vector<double>& rDerivativeShapeFunctionsLocalReal,
	int rCol,
	FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult)
{
	assert(rResult.GetNumRows()==(int)(rShapeFunctionsVirt.size()));
	assert(rResult.GetNumColumns()>=rCol+(int)(rDerivativeShapeFunctionsLocalReal.size()));

	double factor = factorVirt * rTangentStrainVirtStressReal * rTangentStressRealStrainReal;

	for (unsigned int countRow=0; countRow<rShapeFunctionsVirt.size(); countRow++)
	{
		for (unsigned int countCol=0; countCol<rDerivativeShapeFunctionsLocalReal.size(); countCol++)
		{
    		rResult(countRow,rCol+countCol)-=rShapeFunctionsVirt[countRow]*factor * rDerivativeShapeFunctionsLocalReal[countCol];
		}
	}
}

//! @brief add Nt depsilonVirtdSigmaReal dSigmaRealdepsilonReal B
//! @param factorVirt factor
//! @param shapeFunctionsVirt
//! @param tangentStrainVirtStressReal
//! @param rTangentStressRealNonlocalStrainReal
//! @param rShapeFunctionsReal
//! @param rCol column where to add the submatrix
//! @param rResult return value (added)
void NuTo::BoundaryGradientDamage1D::Sub_K_epsilonNonlocalVirt_NonlocalStrainReal
  ( double factorVirt,
    const std::vector<double>& rShapeFunctionsVirt,
	const FullMatrix<double,1,1>& rTangentStrainVirtStressReal,
	const FullMatrix<double,1,1>& rTangentStressRealNonlocalStrainReal,
	int rCol,
	FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult)
{
	assert(rResult.GetNumRows()==(int)(rShapeFunctionsVirt.size()));
	assert(rResult.GetNumColumns()>rCol);

	double factor = factorVirt * rTangentStrainVirtStressReal * rTangentStressRealNonlocalStrainReal;

	for (unsigned int countRow=0; countRow<rShapeFunctionsVirt.size(); countRow++)
	{
		rResult(countRow,rCol)-=rShapeFunctionsVirt[countRow]*factor;
	}
}

//! @brief calculates the Kkk matrix
//! @param shapeFunctions of the ip for all shape functions
//! @param derivativeShapeFunctions of the ip for all shape functions
//! @param c nonlocal gradient radius
//! @param factor multiplication factor (detJ area..)
//! @param Kkk return matrix with detJ * NtT+cBtB
void NuTo::BoundaryGradientDamage1D::CalculateKkk(const std::vector<double>& shapeFunctions,const std::vector<double>& derivativeShapeFunctions,double nonlocalGradientRadius,double factor,
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
