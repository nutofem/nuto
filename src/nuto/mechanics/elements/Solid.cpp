// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include <boost/assign/list_of.hpp>

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"
#include "nuto/mechanics/constitutive/mechanics/Damage.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/thermal/Temperature.h"
#include "nuto/mechanics/constitutive/thermal/HeatFlux3D.h"
#include "nuto/mechanics/constitutive/thermal/LinearHeatFlux.h"
#include "nuto/mechanics/constitutive/thermal/Temperature.h"
#include "nuto/mechanics/constitutive/thermal/TemperatureGradient3D.h"
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixInt.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixDouble.h"
#include "nuto/mechanics/elements/ElementOutputVectorInt.h"
#include "nuto/mechanics/elements/Solid.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/structures/StructureBase.h"

//! @brief constructor
NuTo::Solid::Solid(const StructureBase* rStructure, ElementData::eElementDataType rElementDataType,
        IntegrationType::eIntegrationType rIntegrationType, IpData::eIpDataType rIpDataType) :
        NuTo::ElementBase::ElementBase(rStructure, rElementDataType, rIntegrationType, rIpDataType)
{
    mSection = 0;
}

//! @brief calculates output data of the element
//! @param eOutput ... coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
//! @param updateStaticData (with DummyOutput), IPData, globalrow/column dofs etc.
NuTo::Error::eError NuTo::Solid::Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput)
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
			throw MechanicsException("[NuTo::Solid::Evaluate] no section allocated for element.");

		//calculate coordinates
		int numCoordinates(3*GetNumNodesGeometry());
		std::vector<double> nodeCoord(numCoordinates);
		CalculateCoordinates(nodeCoord);

		//calculate local displacements
		int numDisp(3*GetNumNodesField());
		std::vector<double> nodeDisp(numDisp);
		if (section->GetInputConstitutiveIsDeformationGradient())
		{
			CalculateDisplacements(nodeDisp);
		}
		int numDispDofs(0);
		if (section->GetIsDisplacementDof())
			numDispDofs = numDisp;

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

		//allocate space for local ip coordinates
		double localIPCoord[3];

		//allocate space for derivatives of shape functions
		std::vector<double> derivativeShapeFunctionsGeometryNatural(3*GetNumNodesGeometry());
		std::vector<double> derivativeShapeFunctionsFieldNatural(3*GetNumNodesField());
		std::vector<double> derivativeShapeFunctionsGeometryGlobal(3*GetNumNodesGeometry());
		std::vector<double> derivativeShapeFunctionsFieldGlobal(3*GetNumNodesField());
		std::vector<double> shapeFunctionsGeometry(GetNumNodesGeometry());    //allocate space for shape functions
		std::vector<double> shapeFunctionsField(GetNumNodesGeometry());    //allocate space for shape functions

		//allocate deformation gradient
		DeformationGradient3D deformationGradient;

		EngineeringStrain3D engineeringStrain;

		//allocate global engineering plastic strain
		EngineeringStrain3D engineeringPlasticStrain;

		//allocate damage
		Damage damage;

		//allocate temperature
         Temperature temperature;

		//allocate temperature gradient
		TemperatureGradient3D temperatureGradient;

		//allocate global engineering stress
		EngineeringStress3D engineeringStress;

		//allocate global heat flux
		HeatFlux3D heatFlux;

		//allocate tangents
		ConstitutiveTangentLocal<6,6> tangentStressStrain;
		ConstitutiveTangentLocal<6,1> tangentStressTemperature;
		ConstitutiveTangentLocal<3,3> tangentHeatFluxTemperatureGradient;

		//InvJacobian and determinant of Jacobian
		double invJacobian[9], detJac;

		//for the lumped mass calculation
		double total_mass(0.);

		//define inputs and outputs
		std::map< NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase* > constitutiveInputList;
		std::map< NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase* > constitutiveOutputList;

		if (numDispDofs>0)
		{
			constitutiveInputList[NuTo::Constitutive::Input::DEFORMATION_GRADIENT_3D] = &deformationGradient;
		}

		if (numTempDofs)
		{
			constitutiveInputList[NuTo::Constitutive::Input::TEMPERATURE_GRADIENT_3D] = &temperatureGradient;
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
				it->second->GetFullVectorDouble().Resize(numDispDofs+numTempDofs);
				//if the stiffness matrix is constant, the corresponding internal force is calculated via the Kd
				//on the global level
				if (mStructure->GetHessianConstant(0)==false)
				{
					if (numDispDofs>0)
					{
						constitutiveOutputList[NuTo::Constitutive::Output::ENGINEERING_STRESS_3D] = &engineeringStress;
					}
					if (numTempDofs>0)
					{
						constitutiveOutputList[NuTo::Constitutive::Output::HEAT_FLUX_3D] = &heatFlux;
					}
				}
			break;
			case Element::HESSIAN_0_TIME_DERIVATIVE:
				{
					it->second->GetFullMatrixDouble().Resize(numDispDofs+numTempDofs,numDispDofs+numTempDofs);
					it->second->SetSymmetry(true);
					it->second->SetConstant(true);
					if (numDispDofs>0)
					{
						constitutiveOutputList[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_3D] = &tangentStressStrain;
						//mixed term
						if (numTempDofs)
							constitutiveOutputList[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_TEMPERATURE_3D] = &tangentStressTemperature;
					}
					if (numTempDofs>0)
					{
						constitutiveOutputList[NuTo::Constitutive::Output::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT_3D] = &tangentHeatFluxTemperatureGradient;
						//mixed term
						//if (numDispDofs)
						//    constitutiveOutputList.insert(std::pair<NuTo::Constitutive::eOutput, ConstitutiveOutputBase*>(NuTo::Constitutive::eOutput::D_HEAT_FLUX_D_ENGINEERING_STRAIN_3D, &tangentHeatFluxEngineeringStrain[timeDerivative]));
					}
				}
			break;
			case Element::HESSIAN_1_TIME_DERIVATIVE:
				break;
			case Element::HESSIAN_2_TIME_DERIVATIVE:
				it->second->GetFullMatrixDouble().Resize(numDispDofs,numDispDofs);
				it->second->SetSymmetry(true);
				it->second->SetConstant(true);
			break;
			case Element::LUMPED_HESSIAN_2_TIME_DERIVATIVE:
				it->second->GetFullVectorDouble().Resize(numDispDofs);
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
				switch(it->second->GetIpDataType())
				{
				case NuTo::IpData::ENGINEERING_STRAIN:
					it->second->GetFullMatrixDouble().Resize(6,GetNumIntegrationPoints());
					 //define outputs
					constitutiveOutputList[NuTo::Constitutive::Output::ENGINEERING_STRAIN_3D] = &engineeringStrain;
				break;
				case NuTo::IpData::ENGINEERING_STRESS:
					it->second->GetFullMatrixDouble().Resize(6,GetNumIntegrationPoints());
					 //define outputs
					 constitutiveOutputList[NuTo::Constitutive::Output::ENGINEERING_STRESS_3D] = &engineeringStress;
				break;
				case NuTo::IpData::ENGINEERING_PLASTIC_STRAIN:
					it->second->GetFullMatrixDouble().Resize(6,GetNumIntegrationPoints());
					 //define outputs
					 constitutiveOutputList[NuTo::Constitutive::Output::ENGINEERING_PLASTIC_STRAIN_3D] = &engineeringPlasticStrain;
					 break;
				break;
				case NuTo::IpData::DAMAGE:
					it->second->GetFullMatrixDouble().Resize(1,GetNumIntegrationPoints());
					//define outputs
					constitutiveOutputList[NuTo::Constitutive::Output::DAMAGE] = &damage;
				break;
				default:
					throw MechanicsException("[NuTo::Solid::Evaluate] this ip data type is not implemented.");
				}
			break;
			case Element::GLOBAL_ROW_DOF:
				this->CalculateGlobalRowDofs(it->second->GetVectorInt(),numDispDofs,numTempDofs);
			break;
			case Element::GLOBAL_COLUMN_DOF:
				this->CalculateGlobalColumnDofs(it->second->GetVectorInt(),numDispDofs,numTempDofs);
			break;
			default:
				throw MechanicsException("[NuTo::Solid::Evaluate] element output not implemented.");
			}
		}

		// loop over the integration points
		for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
		{
			GetLocalIntegrationPointCoordinates(theIP, localIPCoord);

			CalculateDerivativeShapeFunctionsGeometryNatural(localIPCoord, derivativeShapeFunctionsGeometryNatural);

			CalculateJacobian(derivativeShapeFunctionsGeometryNatural,nodeCoord, invJacobian, detJac);

			if (GetNumNodesGeometry()==GetNumNodesField())
			{
				//isoparametric
				CalculateDerivativeShapeFunctionsGlobal(derivativeShapeFunctionsGeometryNatural,invJacobian,
					derivativeShapeFunctionsFieldGlobal);
			}
			else
			{
				//not isoparametric
				CalculateDerivativeShapeFunctionsFieldNatural(localIPCoord, derivativeShapeFunctionsFieldNatural);
				CalculateDerivativeShapeFunctionsGlobal(derivativeShapeFunctionsFieldNatural,invJacobian,
					derivativeShapeFunctionsFieldGlobal);

			}

			if (section->GetInputConstitutiveIsDeformationGradient())
			{
				// determine deformation gradient from the local Displacements and the derivative of the shape functions
				// this is not included in the AddIpStiffness to avoid reallocation of the deformation gradient for each IP
				CalculateDeformationGradient(derivativeShapeFunctionsFieldGlobal, nodeDisp, deformationGradient);
			}
			if (section->GetInputConstitutiveIsTemperatureGradient())
			{
				// determine deformation gradient from the local Displacements and the derivative of the shape functions
				// this is not included in the AddIpStiffness to avoid reallocation of the deformation gradient for each IP
				CalculateTemperatureGradient(derivativeShapeFunctionsFieldGlobal, nodeTemp, temperatureGradient);
			}

	        ConstitutiveBase* constitutivePtr = GetConstitutiveLaw(theIP);
			try
			{
				Error::eError error = constitutivePtr->Evaluate3D(this, theIP,
						constitutiveInputList, constitutiveOutputList);
				if (error!=Error::SUCCESSFUL)
					return error;
			}
			catch (NuTo::MechanicsException &e)
			{
				e.AddMessage("[NuTo::Solid::Evaluate] error evaluating the constitutive model.");
				throw e;
			}

			//calculate output
			for (auto it = rElementOutput.begin(); it!=rElementOutput.end(); it++)
			{
				switch(it->first)
				{
				case Element::INTERNAL_GRADIENT:
				{
					// Jacobian
					double factor(fabs(detJac*(mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP))));
					if (numDispDofs>0)
					{
						AddDetJBtSigma(derivativeShapeFunctionsFieldGlobal,engineeringStress, factor, 0, it->second->GetFullVectorDouble());
					}
					if (numTempDofs>0)
					{
						AddDetJBtHeatFlux(derivativeShapeFunctionsFieldGlobal,heatFlux, factor, numDispDofs, it->second->GetFullVectorDouble());
					}
				}
				break;
				case Element::HESSIAN_0_TIME_DERIVATIVE:
					{
						double factor(fabs(detJac*(mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP))));

						if (numDispDofs>0)
						{
							AddDetJBtCB(derivativeShapeFunctionsFieldGlobal, tangentStressStrain, factor, 0,0, it->second->GetFullMatrixDouble());
							if (tangentStressStrain.GetSymmetry()==false)
								it->second->SetSymmetry(false);
							if (tangentStressStrain.GetConstant()==false)
								it->second->SetConstant(false);
							if (numTempDofs>0)
								throw MechanicsException("[NuTo::Solid::Evaluate] mixed terms not yet implemented.");
						}
						if (numTempDofs>0)
						{
							AddDetJBtCB(derivativeShapeFunctionsFieldGlobal, tangentHeatFluxTemperatureGradient, factor, numDispDofs,numDispDofs, it->second->GetFullMatrixDouble());
							if (tangentHeatFluxTemperatureGradient.GetSymmetry()==false)
								it->second->SetSymmetry(false);
							if (tangentHeatFluxTemperatureGradient.GetConstant()==false)
								it->second->SetConstant(false);
							if (numDispDofs>0)
								throw MechanicsException("[NuTo::Solid::Evaluate] mixed terms not yet implemented.");
						}
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
					 }
                 }
                 break;
				 case Element::HESSIAN_2_TIME_DERIVATIVE:
                    if (numDispDofs>0)
                    {
						this->CalculateShapeFunctionsField(localIPCoord, shapeFunctionsField);
						// calculate local mass matrix (the nonlocal terms are zero)
						// don't forget to include determinant of the Jacobian and area
						// detJ * area * density * HtH, :
				        double factor(detJac*(mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP))*
				        		constitutivePtr->GetDensity());
				        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> tmpMatrix;
				        tmpMatrix = (factor*Eigen::Matrix<double,Eigen::Dynamic,1>::Map(&(shapeFunctionsField[0]),GetNumNodesField()))*
				        		Eigen::Matrix<double,1,Eigen::Dynamic>::Map(&(shapeFunctionsField[0]),GetNumNodesField());
				        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& result(it->second->GetFullMatrixDouble());
				        for (int count=0; count<GetNumNodesField(); count++)
				        {
				            for (int count2=0; count2<GetNumNodesField(); count2++)
				            {
				            	result(3*count,3*count2)     += tmpMatrix(count,count2);
				            	result(3*count+1,3*count2+1) += tmpMatrix(count,count2);
				            	result(3*count+2,3*count2+2) += tmpMatrix(count,count2);
				            }
				        }
                    }
                    if (numTempDofs>0)
                    {
                            //no temperature terms
                    }
                break;
				case Element::LUMPED_HESSIAN_2_TIME_DERIVATIVE:
                    if (numDispDofs>0)
                    {
						this->CalculateShapeFunctionsField(localIPCoord, shapeFunctionsField);
						// calculate local mass matrix (the nonlocal terms are zero)
						// don't forget to include determinant of the Jacobian and area
						// detJ * area * density * HtH, :
				        double factor(detJac*(mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP))*
				        		constitutivePtr->GetDensity());
						FullVector<double,Eigen::Dynamic>& result(it->second->GetFullVectorDouble());
						total_mass+=factor;
            			//calculate for the translational dofs the diagonal entries
						for (int count=0; count<GetNumNodesField(); count++)
						{
							result(3*count)+=shapeFunctionsField[count]*shapeFunctionsField[count]*factor;
						}

						if (theIP+1==GetNumIntegrationPoints())
						{
							//calculate sum of diagonal entries (is identical for all directions, that's why only x direction is calculated
							double sum_diagonal(0);
							for (int count=0; count<GetNumNodesField(); count++)
							{
								sum_diagonal+= result(3*count);
							}

							//scale so that the sum of the diagonals represents the full mass
							double scaleFactor = total_mass/sum_diagonal;
							for (int count=0; count<GetNumNodesField(); count++)
							{
								result(3*count) *= scaleFactor;
								result(3*count+1) = result(3*count);
								result(3*count+2) = result(3*count);
							}
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
						memcpy(&(it->second->GetFullMatrixDouble().data()[theIP*6]),engineeringStrain.GetData(),6*sizeof(double));
					break;
					case NuTo::IpData::ENGINEERING_STRESS:
						//error = constitutivePtr->GetEngineeringStressFromEngineeringStrain(this, theIP, deformationGradient, engineeringStress);
						memcpy(&(it->second->GetFullMatrixDouble().data()[theIP*6]),engineeringStress.GetData(),6*sizeof(double));
					break;
					case NuTo::IpData::ENGINEERING_PLASTIC_STRAIN:
						//error = constitutivePtr->GetEngineeringPlasticStrain(this, theIP, deformationGradient, engineeringStrain);
						memcpy(&(it->second->GetFullMatrixDouble().data()[theIP*6]),engineeringPlasticStrain.GetData(),6*sizeof(double));
					break;
					case NuTo::IpData::DAMAGE:
						//error = constitutivePtr->GetDamage(this, theIP, deformationGradient, rIpData.mEigenMatrix.data()[theIP]);
						memcpy(&(it->second->GetFullMatrixDouble().data()[theIP]),damage.GetData(),sizeof(double));
					break;
					default:
						throw MechanicsException("[NuTo::Solid::Evaluate] Ip data not implemented.");
					}
				break;
				case Element::GLOBAL_ROW_DOF:
				case Element::GLOBAL_COLUMN_DOF:
				break;
				default:
					throw MechanicsException("[NuTo::Solid::Evaluate] element output not implemented.");
				}
			}
		}
    }
    catch (NuTo::MechanicsException e)
    {
        std::stringstream ss;
        ss << mStructure->ElementGetId(this);
    	e.AddMessage("[NuTo::Solid::Evaluate] Error evaluating element data of element"	+ ss.str() + ".");
        throw e;
    }

    return Error::SUCCESSFUL;
}

//! @brief adds to a matrix the product B^tCB, where B contains the derivatives of the shape functions and C is the constitutive tangent
//! eventually include also area/width of an element (that's the mechanics solution)
//! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param ConstitutiveTangentBase constitutive tangent matrix
//! @param rFactor factor including determinant of Jacobian and IP weight
//! @param rRow row, where to start to add the submatrix
//! @param rCoefficientMatrix to be added to
void NuTo::Solid::AddDetJBtCB(const std::vector<double>& rDerivativeShapeFunctionsGlobal,
                              const ConstitutiveTangentLocal<6,6>& rConstitutiveTangent, double rFactor,
                              int rRow, int rCol,
                              FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoefficientMatrix)const
{
    const double *C = rConstitutiveTangent.data();
    double x1,x2,y1,y2,z1,z2,x2x1,y2x1,z2x1,x2y1,y2y1,z2y1,x2z1,y2z1,z2z1;
    for (int theNode1=0; theNode1<GetNumNodes(); theNode1++)
    {
        int node1mul3 = 3*theNode1;
        int node1mul3plus1 = node1mul3+1;
        int node1mul3plus2 = node1mul3plus1+1;

        assert((int)rDerivativeShapeFunctionsGlobal.size()>node1mul3plus2);
        x1 = rFactor * rDerivativeShapeFunctionsGlobal[node1mul3];
        y1 = rFactor * rDerivativeShapeFunctionsGlobal[node1mul3plus1];
        z1 = rFactor * rDerivativeShapeFunctionsGlobal[node1mul3plus2];
        node1mul3 +=rRow;
        node1mul3plus1 +=rRow;
        node1mul3plus2 +=rRow;
        for (int theNode2=0; theNode2<GetNumNodes(); theNode2++)
        {
            int node2mul3 = 3*theNode2;
            int node2mul3plus1 = node2mul3+1;
            int node2mul3plus2 = node2mul3plus1+1;
            node2mul3 +=rCol;
            node2mul3plus1 +=rCol;
            node2mul3plus2 +=rCol;

            assert((int)rDerivativeShapeFunctionsGlobal.size()>node2mul3plus2);
            x2 = rDerivativeShapeFunctionsGlobal[node2mul3];
            y2 = rDerivativeShapeFunctionsGlobal[node2mul3plus1];
            z2 = rDerivativeShapeFunctionsGlobal[node2mul3plus2];

            x2x1 = x2*x1;
            y2x1 = y2*x1;
            z2x1 = z2*x1;
            x2y1 = x2*y1;
            y2y1 = y2*y1;
            z2y1 = z2*y1;
            x2z1 = x2*z1;
            y2z1 = y2*z1;
            z2z1 = z2*z1;
            assert(rCoefficientMatrix.GetNumRows()>node1mul3plus2 && rCoefficientMatrix.GetNumColumns()>node1mul3plus2);
            assert(rCoefficientMatrix.GetNumRows()>node2mul3plus2 && rCoefficientMatrix.GetNumColumns()>node2mul3plus2);

            rCoefficientMatrix(node1mul3,node2mul3)          +=x2x1*C[0] +x2y1*C[3] +x2z1*C[5] +y2x1*C[18]+y2y1*C[21]+y2z1*C[23]+z2x1*C[30]+z2y1*C[33]+z2z1*C[35];
            rCoefficientMatrix(node1mul3,node2mul3plus1)     +=y2x1*C[6] +y2y1*C[9] +y2z1*C[11]+x2x1*C[18]+x2y1*C[21]+x2z1*C[23]+z2x1*C[24]+z2y1*C[27]+z2z1*C[29];
            rCoefficientMatrix(node1mul3,node2mul3plus2)     +=z2x1*C[12]+z2y1*C[15]+z2z1*C[17]+y2x1*C[24]+y2y1*C[27]+y2z1*C[29]+x2x1*C[30]+x2y1*C[33]+x2z1*C[35];
            rCoefficientMatrix(node1mul3plus1,node2mul3)     +=x2y1*C[1] +x2x1*C[3] +x2z1*C[4] +y2y1*C[19]+y2x1*C[21]+y2z1*C[22]+z2y1*C[31]+z2x1*C[33]+z2z1*C[34];
            rCoefficientMatrix(node1mul3plus1,node2mul3plus1)+=y2y1*C[7] +y2x1*C[9] +y2z1*C[10]+x2y1*C[19]+x2x1*C[21]+x2z1*C[22]+z2y1*C[25]+z2x1*C[27]+z2z1*C[28];
            rCoefficientMatrix(node1mul3plus1,node2mul3plus2)+=z2y1*C[13]+z2x1*C[15]+z2z1*C[16]+y2y1*C[25]+y2x1*C[27]+y2z1*C[28]+x2y1*C[31]+x2x1*C[33]+x2z1*C[34];
            rCoefficientMatrix(node1mul3plus2,node2mul3)     +=x2z1*C[2] +x2y1*C[4] +x2x1*C[5] +y2z1*C[20]+y2y1*C[22]+y2x1*C[23]+z2z1*C[32]+z2y1*C[34]+z2x1*C[35];
            rCoefficientMatrix(node1mul3plus2,node2mul3plus1)+=y2z1*C[8] +y2y1*C[10]+y2x1*C[11]+x2z1*C[20]+x2y1*C[22]+x2x1*C[23]+z2z1*C[26]+z2y1*C[28]+z2x1*C[29];
            rCoefficientMatrix(node1mul3plus2,node2mul3plus2)+=z2z1*C[14]+z2y1*C[16]+z2x1*C[17]+y2z1*C[26]+y2y1*C[28]+y2x1*C[29]+x2z1*C[32]+x2y1*C[34]+x2x1*C[35];
        }
    }
}

//! @brief adds to a matrix the product B^tCB, where B contains the derivatives of the shape functions and C is the constitutive tangent
//! eventually include also area/width of an element (that's the thermal solution)
//! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param ConstitutiveTangentBase constitutive tangent matrix
//! @param rFactor factor including determinant of Jacobian and IP weight
//! @param rRow row, where to start to add the submatrix
//! @param rCol col, where to start to add the submatrix
//! @param rCoefficientMatrix to be added to
void NuTo::Solid::AddDetJBtCB(const std::vector<double>& rDerivativeShapeFunctionsGlobal,
                              const ConstitutiveTangentLocal<3,3>& rConstitutiveTangent, double rFactor,
                              int rRow, int rCol,
                              FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCoefficientMatrix)const
{
    const double *C = rConstitutiveTangent.data();
    double x1,x2,y1,y2,z1,z2;
    for (int theNode1=0; theNode1<GetNumNodes(); theNode1++)
    {
        int node1mul3 = 3*theNode1;
        int node1mul3plus1 = node1mul3+1;
        int node1mul3plus2 = node1mul3plus1+1;

        assert((int)rDerivativeShapeFunctionsGlobal.size()>node1mul3plus2);
        x1 = rFactor * rDerivativeShapeFunctionsGlobal[node1mul3];
        y1 = rFactor * rDerivativeShapeFunctionsGlobal[node1mul3plus1];
        z1 = rFactor * rDerivativeShapeFunctionsGlobal[node1mul3plus2];
        node1mul3 +=rRow;
        node1mul3plus1 +=rRow;
        node1mul3plus2 +=rRow;
        for (int theNode2=0; theNode2<GetNumNodes(); theNode2++)
        {
            int node2mul3 = 3*theNode2;
            int node2mul3plus1 = node2mul3+1;
            int node2mul3plus2 = node2mul3plus1+1;
            node2mul3 +=rCol;
            node2mul3plus1 +=rCol;
            node2mul3plus2 +=rCol;

            assert((int)rDerivativeShapeFunctionsGlobal.size()>node2mul3plus2);
            x2 = rDerivativeShapeFunctionsGlobal[node2mul3];
            y2 = rDerivativeShapeFunctionsGlobal[node2mul3plus1];
            z2 = rDerivativeShapeFunctionsGlobal[node2mul3plus2];

            assert(rCoefficientMatrix.GetNumRows()>node1mul3plus2 && rCoefficientMatrix.GetNumColumns()>node1mul3plus2);
            assert(rCoefficientMatrix.GetNumRows()>node2mul3plus2 && rCoefficientMatrix.GetNumColumns()>node2mul3plus2);

            rCoefficientMatrix(node1mul3,node2mul3)          +=x1*C[0]*x2;
            rCoefficientMatrix(node1mul3,node2mul3plus1)     +=x1*C[3]*y2;
            rCoefficientMatrix(node1mul3,node2mul3plus2)     +=x1*C[6]*z2;
            rCoefficientMatrix(node1mul3plus1,node2mul3)     +=y1*C[1]*x2;
            rCoefficientMatrix(node1mul3plus1,node2mul3plus1)+=y1*C[4]*y2;
            rCoefficientMatrix(node1mul3plus1,node2mul3plus2)+=y1*C[7]*z2;
            rCoefficientMatrix(node1mul3plus2,node2mul3)     +=z1*C[2]*x2;
            rCoefficientMatrix(node1mul3plus2,node2mul3plus1)+=z1*C[5]*y2;
            rCoefficientMatrix(node1mul3plus2,node2mul3plus2)+=z1*C[8]*z2;
        }
    }
}


//! @brief adds up the internal force vector
//! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param rEngineeringStress stress
//! @param factor factor including det Jacobian area and integration point weight
//! @param rRow start row (in case of a multifield problem)
//! @param rResult resforce vector
void NuTo::Solid::AddDetJBtSigma(const std::vector<double>& rDerivativeShapeFunctionsGlobal,
                                 const EngineeringStress3D& rEngineeringStress,
                                 double rFactor,
                                 int rRow,
                                 FullVector<double,Eigen::Dynamic>& rResult)const
{
    const double *s = rEngineeringStress.GetData();

    double x1,y1,z1;
    for (int theNode1=0; theNode1<GetNumNodes(); theNode1++)
    {
        int node1mul3 = 3*theNode1;
        int node1mul3plus1 = node1mul3+1;
        int node1mul3plus2 = node1mul3plus1+1;

        assert(rResult.GetNumRows()>node1mul3plus2);
        x1 = rFactor * rDerivativeShapeFunctionsGlobal[node1mul3];
        y1 = rFactor * rDerivativeShapeFunctionsGlobal[node1mul3plus1];
        z1 = rFactor * rDerivativeShapeFunctionsGlobal[node1mul3plus2];

        rResult(rRow + node1mul3)     +=x1*s[0]+y1*s[3]+z1*s[5];
        rResult(rRow + node1mul3plus1)+=y1*s[1]+x1*s[3]+z1*s[4];
        rResult(rRow + node1mul3plus2)+=z1*s[2]+y1*s[4]+x1*s[5];
    }
}

//! @brief adds up the internal force vector
//! @param rDerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param rHeatFlux stress
//! @param factor factor including det Jacobian area and integration point weight
//! @param rRow start row (in case of a multifield problem)
//! @param rResult resforce vector
void NuTo::Solid::AddDetJBtHeatFlux(const std::vector<double>& rDerivativeShapeFunctionsGlobal,
                                 const HeatFlux3D& rHeatFlux,
                                 double rFactor,
                                 int rRow,
                                 FullVector<double,Eigen::Dynamic>& rResult)const
{
    const double *s = rHeatFlux.GetData();
    double x1,y1,z1;
    for (int theNode1=0; theNode1<GetNumNodes(); theNode1++)
    {
        int node1mul3 = 3*theNode1;
        int node1mul3plus1 = node1mul3+1;
        int node1mul3plus2 = node1mul3plus1+1;

        assert(rResult.GetNumRows()>node1mul3plus2);
        x1 = rFactor * rDerivativeShapeFunctionsGlobal[node1mul3];
        y1 = rFactor * rDerivativeShapeFunctionsGlobal[node1mul3plus1];
        z1 = rFactor * rDerivativeShapeFunctionsGlobal[node1mul3plus2];

        rResult(rRow + node1mul3)     +=x1*s[0];
        rResult(rRow + node1mul3plus1)+=y1*s[1];
        rResult(rRow + node1mul3plus2)+=z1*s[2];
    }
}

//! @brief Calculates the the inverse of the Jacobian and its determinant
//! @param rDerivativeShapeFunctions Derivatives of the shape functions (dN1dx, dN1dy, dN1dz, dN2dx, ..
//! @param rNodeCoordinates Node coordinates (X1,Y1,Z1,X2,Y2,Z2,...
//! @param rInvJacobian inverse Jacobian matrix (return value)
//! @param rDetJac determinant of the Jacobian (return value)
void NuTo::Solid::CalculateJacobian(const std::vector<double>& rDerivativeShapeFunctions,
                                    const std::vector<double>& rNodeCoordinates,
                                    double rInvJacobian[9],
                                    double& rDetJac)const
{
    /*       jacobian
           j0, j1, j2
           j3, j4, j5
           j6, j7, j8*/

    assert((int)rDerivativeShapeFunctions.size()==3*GetNumNodes() && (int)rNodeCoordinates.size()==3*GetNumNodes());
    double  j0(0.),j1(0.),j2(0.),j3(0.),j4(0.),j5(0.),j6(0.),j7(0.),j8(0.),
    j48_57,j27_18,j15_24,x,y,z;

    int theDerivative(0);
    for (int count = 0; count < GetNumNodes(); count++)
    {
        x = rNodeCoordinates[theDerivative];
        y = rNodeCoordinates[theDerivative+1];
        z = rNodeCoordinates[theDerivative+2];

        j0 += rDerivativeShapeFunctions[theDerivative] * x;
        j3 += rDerivativeShapeFunctions[theDerivative] * y;
        j6 += rDerivativeShapeFunctions[theDerivative] * z;
        theDerivative++;

        j1 += rDerivativeShapeFunctions[theDerivative] * x;
        j4 += rDerivativeShapeFunctions[theDerivative] * y;
        j7 += rDerivativeShapeFunctions[theDerivative] * z;
        theDerivative++;

        j2 += rDerivativeShapeFunctions[theDerivative] * x;
        j5 += rDerivativeShapeFunctions[theDerivative] * y;
        j8 += rDerivativeShapeFunctions[theDerivative] * z;
        theDerivative++;
    }

    j48_57 = j4 * j8 - j5 * j7;
    j27_18 = j2 * j7 - j1 * j8;
    j15_24 = j1 * j5 - j2 * j4;
    /**********************************/
    rDetJac = j0 * j48_57 +j3 * j27_18 + j6 * j15_24;

    if (rDetJac==0)
        throw MechanicsException("[NuTo::Solid::CalculateJacobian] Determinant of the Jacobian is zero, no inversion possible.");

    if (rInvJacobian!=0)
    {
        double invDeterminant(1./rDetJac);
        rInvJacobian[0]=j48_57*invDeterminant;
        rInvJacobian[1]=j27_18*invDeterminant;
        rInvJacobian[2]=j15_24*invDeterminant;
        rInvJacobian[3]=(j5*j6-j3*j8)*invDeterminant;
        rInvJacobian[4]=(j0*j8-j2*j6)*invDeterminant;
        rInvJacobian[5]=(j2*j3-j0*j5)*invDeterminant;
        rInvJacobian[6]=(j3*j7-j4*j6)*invDeterminant;
        rInvJacobian[7]=(j1*j6-j0*j7)*invDeterminant;
        rInvJacobian[8]=(j0*j4-j1*j3)*invDeterminant;
    }
}

//! @brief calculates the derivative of the shape functions with respect to global coordinates
//! @param std::vector<double>& rDerivativeShapeFunctions derivatives of the shape functions
//! @param rJacInv inverse of the Jacobian
//! @param rDerivativeShapeFunctionsGlobal derivaties of the shape functions with respect to global coordinates
//! size should already be correct, but can be checked with an assert
//! first all the directions for a single node, and then for the next node
void NuTo::Solid::CalculateDerivativeShapeFunctionsGlobal(const std::vector<double>& rDerivativeShapeFunctionsLocal, const double rJacInv[9], std::vector<double>& rDerivativeShapeFunctionsGlobal)const
{
    assert(rDerivativeShapeFunctionsLocal.size()==rDerivativeShapeFunctionsGlobal.size());
    for (int count=0; count<GetNumNodes(); count++)
    {
        int mul3count = 3*count;
        int mul3countplus1 = mul3count+1;
        int mul3countplus2 = mul3countplus1+1;
        rDerivativeShapeFunctionsGlobal[mul3count] =
            rDerivativeShapeFunctionsLocal[mul3count]     *rJacInv[0]+
            rDerivativeShapeFunctionsLocal[mul3countplus1]*rJacInv[3]+
            rDerivativeShapeFunctionsLocal[mul3countplus2]*rJacInv[6];
        rDerivativeShapeFunctionsGlobal[mul3countplus1] =
            rDerivativeShapeFunctionsLocal[mul3count]     *rJacInv[1]+
            rDerivativeShapeFunctionsLocal[mul3countplus1]*rJacInv[4]+
            rDerivativeShapeFunctionsLocal[mul3countplus2]*rJacInv[7];
        rDerivativeShapeFunctionsGlobal[mul3countplus2] =
            rDerivativeShapeFunctionsLocal[mul3count]     *rJacInv[2]+
            rDerivativeShapeFunctionsLocal[mul3countplus1]*rJacInv[5]+
            rDerivativeShapeFunctionsLocal[mul3countplus2]*rJacInv[8];
    }
}

//! @brief sets the section of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need a section
//! @param rSection pointer to section
//! @return pointer to constitutive law
void NuTo::Solid::SetSection(const SectionBase* rSection)
{
    //check the properties of the section
    for (int nodeCount = 0; nodeCount < this->GetNumNodes(); nodeCount++)
    {
        if (rSection->GetIsDisplacementDof() || rSection->GetInputConstitutiveIsDeformationGradient())
        {
            //check if all nodes have displacements
            if (GetNode(nodeCount)->GetNumDisplacements()!=3)
            {
                throw MechanicsException("[NuTo::Solid::CheckElement] displacements/strains are ,\
                        defined as input to the constitutive model (at the sectin level), but the nodes don't have 3 displacement dofs.");
            }
        }
        if (rSection->GetIsRotationDof() )
        {
            //check if all nodes have displacements
            //element has no rotational dofs, nothing will be done anyway
            //this is only important, when shells are coupled to bricks
        }
        if (rSection->GetIsTemperatureDof() || rSection->GetInputConstitutiveIsTemperature() || rSection->GetInputConstitutiveIsTemperatureGradient())
        {
            //check if all nodes have displacements
            if (GetNode(nodeCount)->GetNumTemperatures()!=1)
            {
                throw MechanicsException("[NuTo::Solid::CheckElement] displacements/strains are ,\
                        defined as input to the constitutive model (at the section level), but the nodes don't have a temperature dof.");
            }
        }
    }

    mSection = rSection;
}

//! @brief returns a pointer to the section of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need a section
//! @return pointer to section
const NuTo::SectionBase* NuTo::Solid::GetSection()const
{
    return mSection;
}

//! @brief calculates the deformation gradient in 3D
//! @param rRerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param rLocalDisp local displacements
//! @param rDeformationGradient (return value)
void NuTo::Solid::CalculateDeformationGradient(const std::vector<double>& rDerivativeShapeFunctionsGlobal,
        const std::vector<double>& rLocalDisp,
        DeformationGradient3D& rDeformationGradient)const
{
    assert((int)rLocalDisp.size()==3*GetNumNodes() && (int)rDerivativeShapeFunctionsGlobal.size()==3*GetNumNodes());

    rDeformationGradient.mDeformationGradient[0] = 1.;
    rDeformationGradient.mDeformationGradient[1] = 0.;
    rDeformationGradient.mDeformationGradient[2] = 0.;
    rDeformationGradient.mDeformationGradient[3] = 0.;
    rDeformationGradient.mDeformationGradient[4] = 1.;
    rDeformationGradient.mDeformationGradient[5] = 0.;
    rDeformationGradient.mDeformationGradient[6] = 0.;
    rDeformationGradient.mDeformationGradient[7] = 0.;
    rDeformationGradient.mDeformationGradient[8] = 1.;

    int theDisp(0);
    double dNdX,dNdY,dNdZ;
    for (int count=0; count<GetNumNodes(); count++)
    {
        dNdX = rDerivativeShapeFunctionsGlobal[theDisp];
        dNdY = rDerivativeShapeFunctionsGlobal[theDisp+1];
        dNdZ = rDerivativeShapeFunctionsGlobal[theDisp+2];

        rDeformationGradient.mDeformationGradient[0]+=rLocalDisp[theDisp]* dNdX;
        rDeformationGradient.mDeformationGradient[1]+=rLocalDisp[theDisp]* dNdY;
        rDeformationGradient.mDeformationGradient[2]+=rLocalDisp[theDisp]* dNdZ;
        theDisp++;
        rDeformationGradient.mDeformationGradient[3]+=rLocalDisp[theDisp]* dNdX;
        rDeformationGradient.mDeformationGradient[4]+=rLocalDisp[theDisp]* dNdY;
        rDeformationGradient.mDeformationGradient[5]+=rLocalDisp[theDisp]* dNdZ;
        theDisp++;
        rDeformationGradient.mDeformationGradient[6]+=rLocalDisp[theDisp]* dNdX;
        rDeformationGradient.mDeformationGradient[7]+=rLocalDisp[theDisp]* dNdY;
        rDeformationGradient.mDeformationGradient[8]+=rLocalDisp[theDisp]* dNdZ;
        theDisp++;
    }
}

//! @brief calculates the temperature gradient in 3D
//! @param rRerivativeShapeFunctions derivatives of the shape functions with respect to global coordinates
//! @param rTemp nodal temperatures
//! @param rTemperatureGradient (return value)
void NuTo::Solid::CalculateTemperatureGradient(const std::vector<double>& rDerivativeShapeFunctionsGlobal,
        const std::vector<double>& rTemp,
        TemperatureGradient3D& rTemperatureGradient)const
{
    assert((int)rTemp.size()==GetNumNodes() && (int)rDerivativeShapeFunctionsGlobal.size()==3*GetNumNodes());

    rTemperatureGradient.mTemperatureGradient[0] = 0.;
    rTemperatureGradient.mTemperatureGradient[1] = 0.;
    rTemperatureGradient.mTemperatureGradient[2] = 0.;

    int theTemp(0);
    double dNdX,dNdY,dNdZ;
    for (int theNode=0; theNode<GetNumNodes(); theNode++)
    {
        dNdX = rDerivativeShapeFunctionsGlobal[theTemp];
        theTemp++;
        dNdY = rDerivativeShapeFunctionsGlobal[theTemp];
        theTemp++;
        dNdZ = rDerivativeShapeFunctionsGlobal[theTemp];
        theTemp++;

        rTemperatureGradient.mTemperatureGradient[0]+=rTemp[theNode]* dNdX;
        rTemperatureGradient.mTemperatureGradient[1]+=rTemp[theNode]* dNdY;
        rTemperatureGradient.mTemperatureGradient[2]+=rTemp[theNode]* dNdZ;
    }
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point
//! @param rCoordinates local coordinates (return value)
void NuTo::Solid::GetLocalIntegrationPointCoordinates(int rIpNum, double rCoordinates[3])const
{
    this->mElementData->GetIntegrationType()->GetLocalIntegrationPointCoordinates3D(rIpNum, rCoordinates);
    return;
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point
//! @param rCoordinates local coordinates (return value)
void  NuTo::Solid::GetGlobalIntegrationPointCoordinates(int rIpNum, double rCoordinates[3])const
{
    double localCoordinates[3];
    double nodeCoordinates[3];
    this->mElementData->GetIntegrationType()->GetLocalIntegrationPointCoordinates3D(rIpNum, localCoordinates);
    std::vector<double> shapeFunctions(GetNumNodesGeometry());
    CalculateShapeFunctionsGeometry(localCoordinates, shapeFunctions);
    rCoordinates[0] = 0.;
    rCoordinates[1] = 0.;
    rCoordinates[2] = 0.;

    nodeCoordinates[0] = 0;
    nodeCoordinates[1] = 0;
    nodeCoordinates[2] = 0;
    for (int theNode=0; theNode<GetNumNodesGeometry(); theNode++)
    {
        GetNodeGeometry(theNode)->GetCoordinates3D(nodeCoordinates);
        rCoordinates[0]+=shapeFunctions[theNode]*nodeCoordinates[0];
        rCoordinates[1]+=shapeFunctions[theNode]*nodeCoordinates[1];
        rCoordinates[2]+=shapeFunctions[theNode]*nodeCoordinates[2];
    }
    return;
}

//! @brief Allocates static data for an integration point of an element
//! @param rConstitutiveLaw constitutive law, which is called to allocate the static data object
//! actually, both - the element type and the constitutive law are required to determine the static data object actually required
NuTo::ConstitutiveStaticDataBase* NuTo::Solid::AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw)const
{
    return rConstitutiveLaw->AllocateStaticDataEngineeringStress_EngineeringStrain3D(this);
}

//! @brief stores the coordinates of the nodes
//! @param localCoordinates vector with already correct size allocated
//! this can be checked with an assertation
void NuTo::Solid::CalculateCoordinates(std::vector<double>& rCoordinates)const
{
    assert((int)rCoordinates.size()==3*GetNumNodes());
    for (int count=0; count<GetNumNodesGeometry(); count++)
    {
        GetNodeGeometry(count)->GetCoordinates3D(&(rCoordinates[3*count]));
    }

}

//! @brief stores the displacements of the nodes
//! @param localDisplacements vector with already correct size allocated
//! this can be checked with an assertation
void NuTo::Solid::CalculateDisplacements(std::vector<double>& rDisplacements)const
{
    assert((int)rDisplacements.size()==3*GetNumNodesField());
    for (int count=0; count<GetNumNodesField(); count++)
    {
        if (GetNodeField(count)->GetNumDisplacements()!=3)
            throw MechanicsException("[NuTo::Solid::CalculateDisplacements] Displacement is required as input to the constitutive model, but the node does not have this data.");
        GetNodeField(count)->GetDisplacements3D(&(rDisplacements[3*count]));
    }
}

//! @brief stores the temperatures of the nodes
//! @param temperature vector with already correct size allocated
//! this can be checked with an assertation
void NuTo::Solid::CalculateTemperatures(std::vector<double>& rTemperatures)const
{
    assert((int)rTemperatures.size()==GetNumNodesField());
    for (int count=0; count<GetNumNodesField(); count++)
    {
        if (GetNodeField(count)->GetNumTemperatures()!=1)
            throw MechanicsException("[NuTo::Solid::CalculateTemperatures] Temperature is required as input to the constitutive model, but the node does not have this data.");
        rTemperatures[count] = GetNodeField(count)->GetTemperature();
    }
}


// interpolate geometry
void NuTo::Solid::InterpolateCoordinatesFrom3D(double rLocalCoordinates[3], double rGlobalCoordinates[3]) const
{
    // calculate shape functions
    std::vector<double> ShapeFunctions(this->GetNumNodesGeometry());
    this->CalculateShapeFunctionsGeometry(rLocalCoordinates, ShapeFunctions);

    // start interpolation
    rGlobalCoordinates[0] = 0.0;
    rGlobalCoordinates[1] = 0.0;
    rGlobalCoordinates[2] = 0.0;
    for (int NodeCount = 0; NodeCount < this->GetNumNodesGeometry(); NodeCount++)
    {
        // get node coordinate
        double NodeCoordinate[3];
        GetNodeGeometry(NodeCount)->GetCoordinates3D(NodeCoordinate);

        // add node contribution
        rGlobalCoordinates[0] += ShapeFunctions[NodeCount] *  NodeCoordinate[0];
        rGlobalCoordinates[1] += ShapeFunctions[NodeCount] *  NodeCoordinate[1];
        rGlobalCoordinates[2] += ShapeFunctions[NodeCount] *  NodeCoordinate[2];
    }
}

// interpolate displacements
void NuTo::Solid::InterpolateDisplacementsFrom3D(int rTimeDerivative, double rLocalCoordinates[3], double rGlobalDisplacements[3]) const
{
    // calculate shape functions
    std::vector<double> ShapeFunctions(this->GetNumNodesField());
    this->CalculateShapeFunctionsField(rLocalCoordinates, ShapeFunctions);

    // start interpolation
    rGlobalDisplacements[0] = 0.0;
    rGlobalDisplacements[1] = 0.0;
    rGlobalDisplacements[2] = 0.0;
    for (int NodeCount = 0; NodeCount < this->GetNumNodesField(); NodeCount++)
    {
        // get node displacements
        double NodeDisplacement[3];
        GetNodeField(NodeCount)->GetDisplacements3D(rTimeDerivative, NodeDisplacement);

        // add node contribution
        rGlobalDisplacements[0] += ShapeFunctions[NodeCount] *  NodeDisplacement[0];
        rGlobalDisplacements[1] += ShapeFunctions[NodeCount] *  NodeDisplacement[1];
        rGlobalDisplacements[2] += ShapeFunctions[NodeCount] *  NodeDisplacement[2];
    }
}

// interpolate temperature
void NuTo::Solid::InterpolateTemperatureFrom3D(double rLocalCoordinates[3], double& rTemperature) const
{
    // calculate shape functions
    std::vector<double> ShapeFunctions(this->GetNumNodesField());
    this->CalculateShapeFunctionsField(rLocalCoordinates, ShapeFunctions);

    // start interpolation
    rTemperature = 0.0;
    for (int NodeCount = 0; NodeCount < this->GetNumNodesField(); NodeCount++)
    {
        // get node coordinate
        double nodeTemperature;
        nodeTemperature = GetNodeField(NodeCount)->GetTemperature();

        rTemperature +=ShapeFunctions[NodeCount] * nodeTemperature;
    }
}


//! @brief ... extract global dofs from nodes (mapping of local column ordering of the element matrices to the global dof ordering)
//! @param rGlobalColumnDofs ... vector of global column dofs
//! @param rNumDisp ... number of displacement dofs
//! @param rNumTemp ... number of temperature dofs
void NuTo::Solid::CalculateGlobalRowDofs(std::vector<int>& rGlobalRowDofs, int numDispDofs, int numTempDofs) const
{
    rGlobalRowDofs.resize(numDispDofs+numTempDofs);
    for (int nodeCount = 0; nodeCount < this->GetNumNodes(); nodeCount++)
    {
        const NodeBase * nodePtr(GetNode(nodeCount));
        if (nodePtr->GetNumDisplacements()>0 && numDispDofs>0)
        {
            rGlobalRowDofs[3 * nodeCount    ] = nodePtr->GetDofDisplacement(0);
            rGlobalRowDofs[3 * nodeCount + 1] = nodePtr->GetDofDisplacement(1);
            rGlobalRowDofs[3 * nodeCount + 2] = nodePtr->GetDofDisplacement(2);
        }
        if (nodePtr->GetNumTemperatures()>0 && numTempDofs>0)
        {
            rGlobalRowDofs[numDispDofs + nodeCount ] = nodePtr->GetDofTemperature();
        }
    }
}

//! @brief ... extract global dofs from nodes (mapping of local column ordering of the element matrices to the global dof ordering)
//! @param rGlobalColumnDofs ... vector of global column dofs
//! @param rNumDisp ... number of displacement dofs
//! @param rNumTemp ... number of temperature dofs
void NuTo::Solid::CalculateGlobalColumnDofs(std::vector<int>& rGlobalColumnDofs, int rNumDisp, int rNumTemp) const
{
    this->CalculateGlobalRowDofs(rGlobalColumnDofs,rNumDisp,rNumTemp);
}


// check element definition
void NuTo::Solid::CheckElement()
{
    // check nodes
    for (int nodeCount = 0; nodeCount < this->GetNumNodes(); nodeCount++)
    {
        if (GetNode(nodeCount)->GetNumCoordinates()!=3)
        {
            throw MechanicsException("[NuTo::Solid::CheckElement] invalid node type (check node definition for coordinates).");
        }
    }

    // check node ordering (element length must be positive) and for changing sign in jacobian determinant
    // calculate coordinates
    std::vector<double> nodeCoord(3*GetNumNodesGeometry());
    this->CalculateCoordinates(nodeCoord);

    // check number of integration points
    if (this->GetNumIntegrationPoints() < 1)
    {
        throw MechanicsException("[NuTo::Solid::CheckElement] invalid integration type.");
    }

    // check sign of the jacobian determinant of the first integration point
    double localIPCoord[3];
    this->GetLocalIntegrationPointCoordinates(0, localIPCoord);

    std::vector<double> derivativeShapeFunctionsLocal(3*GetNumNodesGeometry());
    this->CalculateDerivativeShapeFunctionsGeometryNatural(localIPCoord, derivativeShapeFunctionsLocal);

    double invJacobian[9], detJacobian;
    this->CalculateJacobian(derivativeShapeFunctionsLocal,nodeCoord, invJacobian, detJacobian);
    // reorder nodes if determinant is negative
    if (detJacobian < 0.0)
    {
        std::cout << "coordinates \n" << localIPCoord[0] << " " << localIPCoord[1] << " "<< localIPCoord[2]<< std::endl;
        std::cout << "invJacobian \n" << invJacobian[0] << " " << invJacobian[3] << " "<< invJacobian[6]<< std::endl;
        std::cout <<  invJacobian[1] << " " << invJacobian[4] << " "<< invJacobian[7]<< std::endl;
        std::cout <<  invJacobian[2] << " " << invJacobian[5] << " "<< invJacobian[8] << std::endl;
    	this->ReorderNodes();
        // recalculate node coordinates after renumbering
        this->CalculateCoordinates(nodeCoord);
        std::cout << "reorder element due to detJac=" << detJacobian << std::endl;
    }

    // check jacobian determinant for all integration points for positive sign and calculate element volume
    double volume = 0;
    for (int ipCount = 0; ipCount < this->GetNumIntegrationPoints(); ipCount++)
    {
        // calculate jacobian determinant
        this->GetLocalIntegrationPointCoordinates(ipCount, localIPCoord);
        this->CalculateDerivativeShapeFunctionsGeometryNatural(localIPCoord, derivativeShapeFunctionsLocal);
        this->CalculateJacobian(derivativeShapeFunctionsLocal,nodeCoord, invJacobian, detJacobian);
        if (detJacobian <= 0)
        {
            std::cout << "error due to detJac=" << detJacobian << std::endl;
            for (int count=0; count<GetNumNodes(); count++)
            {
            	std::cout << "node " << count+1<< " coordinates : " << nodeCoord[3*count] << " " << nodeCoord[3*count+1] << " " << nodeCoord[3*count+2] << std::endl;
            }
            throw MechanicsException("[NuTo::Solid::CheckElement] element is not properly defined by this nodes (zero or negative jacobian determinant).");
        }
        volume += this->GetIntegrationPointWeight(ipCount) * detJacobian;
    }

    // check element volume
    if (volume < 1e-14)
    {
        throw MechanicsException("[NuTo::Solid::CheckElement] element with zero volume (check nodes).");
    }
}

//! @brief calculates the volume of an integration point (weight * detJac)
//! @param rVolume  vector for storage of the ip volumes (area in 2D)
void NuTo::Solid::GetIntegrationPointVolume(std::vector<double>& rVolume)const
{
    //calculate coordinates
    std::vector<double> nodeCoord(3*GetNumNodesGeometry());
    CalculateCoordinates(nodeCoord);

    //allocate space for local ip coordinates
    double localIPCoord[3];

    //allocate space for derivatives of shape functions
    std::vector<double> derivativeShapeFunctionsLocal(3*GetNumNodesGeometry());

    //determinant of Jacobian
    double detJac;

    rVolume.resize(GetNumIntegrationPoints());

    for (int theIP=0; theIP<GetNumIntegrationPoints(); theIP++)
    {
        GetLocalIntegrationPointCoordinates(theIP, localIPCoord);

        CalculateDerivativeShapeFunctionsGeometryNatural(localIPCoord, derivativeShapeFunctionsLocal);

        CalculateJacobian(derivativeShapeFunctionsLocal,nodeCoord, 0, detJac);

        rVolume[theIP] = detJac * mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP);
    }
}
//! @brief cast the base pointer to a Solid, otherwise throws an exception
const NuTo::Solid* NuTo::Solid::AsSolid()const
{
    return this;
}

//! @brief cast the base pointer to a Solid, otherwise throws an exception
NuTo::Solid* NuTo::Solid::AsSolid()
{
    return this;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::Solid::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Solid::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Solid::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Solid::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Solid::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Solid::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Solid::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Solid" << std::endl;
#endif
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementBase)
           & BOOST_SERIALIZATION_NVP(mSection);
    }
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Solid" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Solid)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::Solid)
#endif // ENABLE_SERIALIZATION
