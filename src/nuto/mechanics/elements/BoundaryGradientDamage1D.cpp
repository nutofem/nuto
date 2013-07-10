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
    		)
{
	mEdgeRealBoundaryElement = rEdgeRealBoundaryElement;
	mRealBoundaryElement = rRealBoundaryElement;
}

//! @brief calculates output data fo the elmement
//! @param eOutput ... coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
//!                    @param updateStaticData (with DummyOutput), IPData, globalrow/column dofs etc.
NuTo::Error::eError NuTo::BoundaryGradientDamage1D::Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput)
{

	if (mStructure->GetHessianConstant(1)==false)
    	throw MechanicsException("[NuTo::Truss::Evaluate] only implemented for a constant Hessian for the first derivative (damping).");
    if (mStructure->GetHessianConstant(2)==false)
    	throw MechanicsException("[NuTo::Truss::Evaluate] only implemented for a constant Hessian for the second derivative (mass).");

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
			throw MechanicsException("[NuTo::Truss::Evaluate] no section allocated for real boundary element.");

		//calculate coordinates
		int numCoordinatesReal(mRealBoundaryElement->GetNumShapeFunctions());
		std::vector<double> localNodeCoordReal(numCoordinatesReal);
		mRealBoundaryElement->CalculateLocalCoordinates(localNodeCoordReal);

		//calculate local displacements, velocities and accelerations
		//the difference between disp and dispdof is a problem where the displacements are fixed, but enter the constitutive equation
		//for example in a two stage problem, first solve mechanics, then thermal and so on
		int numDispReal(mRealBoundaryElement->GetNumShapeFunctions());
		int numDispDofsReal = sectionReal->GetIsDisplacementDof() ? numDispReal : 0;

		int numTempReal(mRealBoundaryElement->GetNumShapeFunctions());
		int numTempDofsReal(sectionReal->GetIsTemperatureDof() ? numTempReal : 0);

		int numNonlocalEqPlasticStrainReal(2*mRealBoundaryElement->GetNumShapeFunctions());
		int numNonlocalEqPlasticStrainDofsReal(sectionReal->GetIsNonlocalEqPlasticStrainDof() ? numNonlocalEqPlasticStrainReal : 0);

		int numNonlocalTotalStrainReal(mRealBoundaryElement->GetNumShapeFunctions());
		int numNonlocalTotalStrainDofsReal(sectionReal->GetIsNonlocalTotalStrainDof() ? numNonlocalTotalStrainReal : 0);

		std::vector<double> localNodeDispReal,nodeTempReal,nodeNonlocalEqPlasticStrainReal,nodeNonlocalTotalStrainReal;

		//calculate local displacements, velocities and accelerations
		if (numDispDofsReal>0)
		{
			localNodeDispReal.resize(numDispReal);
			mRealBoundaryElement->CalculateLocalDisplacements(0,localNodeDispReal);
		}
		if (numTempDofsReal>0)
		{
			nodeTempReal.resize(numTempReal);
			mRealBoundaryElement->CalculateNodalTemperatures(0,nodeTempReal);
		}
		if (numNonlocalEqPlasticStrainDofsReal>0 || sectionReal->GetInputConstitutiveIsNonlocalEqPlasticStrain())
		{
			nodeNonlocalEqPlasticStrainReal.resize(numNonlocalEqPlasticStrainReal);
			mRealBoundaryElement->CalculateNodalNonlocalEqPlasticStrain(0,nodeNonlocalEqPlasticStrainReal);
		}
		if (numNonlocalTotalStrainDofsReal>0 || sectionReal->GetInputConstitutiveIsNonlocalTotalStrain())
		{
			nodeNonlocalTotalStrainReal.resize(numNonlocalTotalStrainReal);
			mRealBoundaryElement->CalculateNodalNonlocalTotalStrain(0,nodeNonlocalTotalStrainReal);
		}

		//allocate space for local ip coordinates
		double localIPCoordReal;

		//allocate space for local shape functions
		std::vector<double> derivativeShapeFunctionsNaturalReal(mRealBoundaryElement->GetLocalDimension()*mRealBoundaryElement->GetNumShapeFunctions());  //allocate space for derivatives of shape functions
		std::vector<double> derivativeShapeFunctionsLocalReal(mRealBoundaryElement->GetLocalDimension()*mRealBoundaryElement->GetNumShapeFunctions());    //allocate space for derivatives of shape functions
		std::vector<double> shapeFunctionsReal(mRealBoundaryElement->GetNumShapeFunctions());                                       //allocate space for derivatives of shape functions

		//allocate deformation gradient
		DeformationGradient1D deformationGradientReal;

		EngineeringStrain3D engineeringStrain3DReal;

		//allocate global engineering plastic strain
		EngineeringStrain3D engineeringPlasticStrain3DReal;

		//allocate  damage (output of constitutive relation)
		Damage damageReal;

		//allocate nonlocal eq plastic strain (nodal dof value, input of constitutive relation)
		NonlocalEqPlasticStrain nonLocalEqPlasticStrainReal;

		EngineeringStrain1D nonlocalTotalStrainReal;

		EngineeringStrain1D localTotalStrainReal;

		//allocate temperature
		Temperature temperatureReal;

		//allocate global engineering stress
		EngineeringStress1D engineeringStress1DReal;
		EngineeringStress3D engineeringStress3DReal;

		//allocate tangents
		ConstitutiveTangentLocal<1,1> tangentStressStrainReal;
		ConstitutiveTangentLocal<1,1> tangentStressTemperatureReal;
		ConstitutiveTangentLocal<1,2> tangentStressNonlocalEqPlasticStrainReal;
		ConstitutiveTangentLocal<1,1> tangentStressNonlocalTotalStrainReal;

		//define inputs and outputs
		std::map< NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase* > constitutiveInputListReal;
		std::map< NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase* > constitutiveOutputListReal;

		if (numDispDofsReal>0)
		{
			constitutiveInputListReal[NuTo::Constitutive::Input::DEFORMATION_GRADIENT_1D] = &deformationGradientReal;
		}

		if (sectionReal->GetInputConstitutiveIsTemperature())
		{
			constitutiveInputListReal[NuTo::Constitutive::Input::TEMPERATURE] = &temperatureReal;
		}

		if (numNonlocalEqPlasticStrainDofsReal>0)
		{
			constitutiveInputListReal[NuTo::Constitutive::Input::NONLOCAL_EQ_PLASTIC_STRAIN] = &nonLocalEqPlasticStrainReal;
		}

		if (numNonlocalTotalStrainDofsReal>0)
		{
			constitutiveInputListReal[NuTo::Constitutive::Input::NONLOCAL_TOTAL_STRAIN_1D] = &nonlocalTotalStrainReal;
		}

    	//*****************************************************************************************************
        //Now calculate the relevant informations for the virtual boundary element within the element boundary (nodal values
    	//*****************************************************************************************************
		//calculate coordinates
		int numCoordinatesVirt(this->GetNumShapeFunctions());
		std::vector<double> localNodeCoordVirt(numCoordinatesVirt);
//this->CalculateLocalCoordinates(localNodeCoordVirt);

		//calculate local displacements, velocities and accelerations
		//the difference between disp and dispdof is a problem where the displacements are fixed, but enter the constitutive equation
		//for example in a two stage problem, first solve mechanics, then thermal and so on
		int numNonlocalEqPlasticStrainVirt(2*this->GetNumShapeFunctions());
		int numNonlocalEqPlasticStrainDofsVirt(sectionReal->GetIsNonlocalEqPlasticStrainDof() ? numNonlocalEqPlasticStrainVirt : 0);

		int numNonlocalTotalStrainVirt(this->GetNumShapeFunctions());
		int numNonlocalTotalStrainDofsVirt(sectionReal->GetIsNonlocalTotalStrainDof() ? numNonlocalTotalStrainVirt : 0);

		std::vector<double> nodeNonlocalEqPlasticStrainVirt,nodeNonlocalTotalStrainVirt;

		if (numNonlocalEqPlasticStrainDofsVirt>0 || sectionReal->GetInputConstitutiveIsNonlocalEqPlasticStrain())
		{
			nodeNonlocalEqPlasticStrainVirt.resize(numNonlocalEqPlasticStrainVirt);
//this->CalculateNodalNonlocalEqPlasticStrain(0,nodeNonlocalEqPlasticStrainVirt);
		}
		if (numNonlocalTotalStrainDofsVirt>0 || sectionReal->GetInputConstitutiveIsNonlocalTotalStrain())
		{
			nodeNonlocalTotalStrainVirt.resize(numNonlocalTotalStrainVirt);
//this->CalculateNodalNonlocalTotalStrain(0,nodeNonlocalTotalStrainVirt);
		}

		//allocate space for local ip coordinates
		double localIPCoordVirt;

		//allocate space for local shape functions
		std::vector<double> derivativeShapeFunctionsNaturalVirt(this->GetLocalDimension()*this->GetNumShapeFunctions());  //allocate space for derivatives of shape functions
		std::vector<double> derivativeShapeFunctionsLocalVirt(this->GetLocalDimension()*this->GetNumShapeFunctions());    //allocate space for derivatives of shape functions
		std::vector<double> shapeFunctionsVirt(this->GetNumShapeFunctions());                                       //allocate space for derivatives of shape functions

		//allocate nonlocal eq plastic strain (nodal dof value, input of constitutive relation)
		NonlocalEqPlasticStrain nonLocalEqPlasticStrainVirt;

		EngineeringStrain1D nonlocalTotalStrainVirt;

		EngineeringStrain1D localTotalStrainVirt;

		//allocate tangents
		ConstitutiveTangentLocal<1,1> tangentBoundaryStrainBoundaryStress;
		ConstitutiveTangentLocal<1,1> tangentBoundaryStrainNonlocalTotalStrainVirt;
		ConstitutiveTangentLocal<1,1> tangentBoundaryStrainNonlocalEqPlasticStrainVirt;

		//define inputs and outputs
		std::map< NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase* > constitutiveInputListVirt;
		std::map< NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase* > constitutiveOutputListVirt;

		if (numNonlocalEqPlasticStrainDofsVirt>0)
		{
			constitutiveInputListVirt[NuTo::Constitutive::Input::ENGINEERING_STRESS_1D] = &engineeringStress1DReal;
			constitutiveInputListVirt[NuTo::Constitutive::Input::NONLOCAL_EQ_PLASTIC_STRAIN] = &nonLocalEqPlasticStrainVirt;
		}

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
				it->second->GetFullVectorDouble().Resize(numNonlocalEqPlasticStrainDofsVirt+numNonlocalTotalStrainDofsVirt);
				//if the stiffness matrix is constant, the corresponding internal force is calculated via the Kd
				//on the global level
				if (mStructure->GetHessianConstant(0)==false)
				{
					if (numDispDofsReal>0)
					{
						constitutiveOutputListReal[NuTo::Constitutive::Output::ENGINEERING_STRESS_1D] = &engineeringStress1DReal;
					}
					if (numNonlocalTotalStrainDofsVirt>0)
					{
						constitutiveOutputListVirt[NuTo::Constitutive::Output::ENGINEERING_STRAIN_FROM_BOUNDARY_1D] = &localTotalStrainVirt;
					}
					if (numNonlocalEqPlasticStrainDofsVirt>0)
					{
						constitutiveOutputListVirt[NuTo::Constitutive::Output::ENGINEERING_STRAIN_FROM_BOUNDARY_1D] = &localTotalStrainVirt;
					}
				}
			break;
			case Element::HESSIAN_0_TIME_DERIVATIVE:
				{
					it->second->GetFullMatrixDouble().Resize(numNonlocalEqPlasticStrainDofsVirt+numNonlocalTotalStrainDofsVirt,numNonlocalEqPlasticStrainDofsVirt+numNonlocalTotalStrainDofsVirt + numDispDofsReal+numTempDofsReal+numNonlocalEqPlasticStrainDofsReal+numNonlocalTotalStrainDofsReal);
					it->second->SetSymmetry(false);
					it->second->SetConstant(false);
					if (numDispDofsReal>0)
					{
						constitutiveOutputListReal[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_1D] = &tangentStressStrainReal;
						//mixed terms
						if (numTempDofsReal>0)
							constitutiveOutputListReal[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_TEMPERATURE_1D] = &tangentStressTemperatureReal;
						if (numNonlocalEqPlasticStrainDofsReal>0)
						{
							constitutiveOutputListReal[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_PLASTIC_STRAIN_1D] = &tangentStressNonlocalEqPlasticStrainReal;
						}
						if (numNonlocalTotalStrainDofsReal>0)
						{
							constitutiveOutputListReal[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_TOTAL_STRAIN_1D] = &tangentStressNonlocalTotalStrainReal;
						}
					}
					if (numNonlocalTotalStrainDofsVirt>0)
					{
						constitutiveOutputListVirt[NuTo::Constitutive::Output::D_ENGINEERING_STRAIN_FROM_BOUNDARY_D_STRESS_1D] = &tangentBoundaryStrainBoundaryStress;
						constitutiveOutputListVirt[NuTo::Constitutive::Output::D_ENGINEERING_STRAIN_FROM_BOUNDARY_D_NONLOCAL_TOTAL_STRAIN_1D] = &tangentBoundaryStrainNonlocalTotalStrainVirt;
					}
					if (numNonlocalTotalStrainDofsVirt>0)
					{
						constitutiveOutputListVirt[NuTo::Constitutive::Output::D_ENGINEERING_STRAIN_FROM_BOUNDARY_D_STRESS_1D] = &tangentBoundaryStrainBoundaryStress;
						constitutiveOutputListVirt[NuTo::Constitutive::Output::D_ENGINEERING_STRAIN_FROM_BOUNDARY_D_NONLOCAL_TOTAL_STRAIN_1D] = &tangentBoundaryStrainNonlocalTotalStrainVirt;
					}
				}
			break;
			case Element::HESSIAN_1_TIME_DERIVATIVE:
			{
				it->second->GetFullMatrixDouble().Resize(numNonlocalEqPlasticStrainDofsVirt+numNonlocalTotalStrainDofsVirt,numNonlocalEqPlasticStrainDofsVirt+numNonlocalTotalStrainDofsVirt + numDispDofsReal+numTempDofsReal+numNonlocalEqPlasticStrainDofsReal+numNonlocalTotalStrainDofsReal);
				it->second->SetSymmetry(true);
				it->second->SetConstant(true);
				// Rayleigh damping should be introduced on the global level
			}
			break;
			case Element::HESSIAN_2_TIME_DERIVATIVE:
			{
				it->second->GetFullMatrixDouble().Resize(numNonlocalEqPlasticStrainDofsVirt+numNonlocalTotalStrainDofsVirt,numNonlocalEqPlasticStrainDofsVirt+numNonlocalTotalStrainDofsVirt + numDispDofsReal+numTempDofsReal+numNonlocalEqPlasticStrainDofsReal+numNonlocalTotalStrainDofsReal);
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
			break;
			case Element::GLOBAL_ROW_DOF:
				//this->CalculateGlobalRowDofs(it->second->GetVectorInt(),numDispDofs,numTempDofs,numNonlocalEqPlasticStrainDofs,numNonlocalTotalStrainDofs);
			break;
			case Element::GLOBAL_COLUMN_DOF:
				//this->CalculateGlobalColumnDofs(it->second->GetVectorInt(),numDispDofs,numTempDofs,numNonlocalEqPlasticStrainDofs,numNonlocalTotalStrainDofs);
			break;
			default:
				throw MechanicsException("[NuTo::Truss::Evaluate] element output not implemented.");
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
/*
				//the ip is located on the boundary of the real element
                this->GetLocalIntegrationPointCoordinatesReal(theIP, localIPCoordReal);

				//derivative in natural coordinate system
				mRealBoundaryElement->CalculateDerivativeShapeFunctions(localIPCoordReal, derivativeShapeFunctionsNaturalReal);

				//determinant of the Jacobian
				double detJReal(mRealBoundaryElement->DetJacobian(derivativeShapeFunctionsNaturalReal,localNodeCoordReal));

				//derivative in local coordinate system
				for (unsigned int count=0; count<derivativeShapeFunctionsNaturalReal.size(); count++)
				{
					derivativeShapeFunctionsLocalReal[count] = derivativeShapeFunctionsNaturalReal[count]/detJ;
				}

				if (numDispDofsReal)
				{
					// determine deformation gradient from the local Displacements and the derivative of the shape functions
					CalculateDeformationGradient(derivativeShapeFunctionsLocalReal, localNodeCoordReal, localNodeDispReal, deformationGradientReal);
				}

				if (section->GetInputConstitutiveIsTemperature())
				{
					// determine temperature
					throw MechanicsException("[Nuto::Truss::Evaluate] temperature not yet implemented.");
				}

				if (numNonlocalEqPlasticStrainDofsReal)
				{
					mRealBoundaryElement->CalculateShapeFunctions(localIPCoordReal, shapeFunctionsReal);
					mRealBoundaryElement->CalculateNonlocalEqPlasticStrain(shapeFunctionsReal, nodeNonlocalEqPlasticStrainReal, nonLocalEqPlasticStrainReal);
				}

				if (numNonlocalTotalStrainDofsReal)
				{
					CalculateShapeFunctions(localIPCoordReal, shapeFunctionsReal);
					CalculateNonlocalTotalStrain(shapeFunctionsReal, nodeNonlocalTotalStrainReal, nonlocalTotalStrainReal);
				}

				ConstitutiveBase* constitutivePtr = this->GetConstitutiveLaw(theIP);
				Error::eError error = constitutivePtr->Evaluate1D(this, theIP,
						constitutiveInputListReal, constitutiveOutputListReal);
				if (error!=Error::SUCCESSFUL)
					return error;
					*/
			}
			else
			{
				//the ip is located within the virtual boundary element
			}


			//double factor (detJ*mSection->GetArea()*
			//		       (mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP)));
		}

/*

     	//*********************************************************************************************************************
		// Perform the integration in the virtual boundary element taking into account the data from the real boundary element
		//*********************************************************************************************************************



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
			constitutiveInputList[NuTo::Constitutive::Input::NONLOCAL_EQ_PLASTIC_STRAIN] = &nonLocalEqPlasticStrain;
		}

		if (numNonlocalTotalStrainDofs>0)
		{
			constitutiveInputList[NuTo::Constitutive::Input::NONLOCAL_TOTAL_STRAIN_1D] = &nonlocalTotalStrain;
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
*/
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
