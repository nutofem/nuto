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
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixInt.h"
#include "nuto/mechanics/elements/ElementOutputFullMatrixDouble.h"
#include "nuto/mechanics/elements/ElementOutputVectorInt.h"
#include "nuto/mechanics/elements/BoundaryGradientDamage1D.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain1D.h"
#include "nuto/mechanics/constitutive/mechanics/Damage.h"

#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/constitutive/mechanics/LocalEqStrain.h"
#include "nuto/mechanics/constitutive/mechanics/NonlocalEqStrain.h"
#include "nuto/mechanics/sections/SectionTruss.h"

#include "nuto/mechanics/structures/StructureBase.h"

//! @brief constructor
NuTo::BoundaryGradientDamage1D::BoundaryGradientDamage1D(const StructureBase* rStructure,
    		Truss* rRealBoundaryElement,
    		int rSurfaceEdge,
            BoundaryCondition::eType rBoundaryConditionType,
    		ElementData::eElementDataType rElementDataType,
    		IntegrationType::eIntegrationType rIntegrationType,
    		IpData::eIpDataType rIpDataType
    		) :
    		NuTo::ElementBase::ElementBase(rStructure, rElementDataType, rIntegrationType, rIpDataType)
{
    mSurfaceEdge = rSurfaceEdge;
	mRealBoundaryElement = rRealBoundaryElement;
	mBoundaryConditionType = rBoundaryConditionType;
}

NuTo::BoundaryGradientDamage1D::BoundaryGradientDamage1D(const StructureBase* rStructure,
            Truss* rRealBoundaryElement,
            NodeBase* rSurfaceNode,
            BoundaryCondition::eType rBoundaryConditionType,
            ElementData::eElementDataType rElementDataType,
            IntegrationType::eIntegrationType rIntegrationType,
            IpData::eIpDataType rIpDataType
            ) :
            NuTo::ElementBase::ElementBase(rStructure, rElementDataType, rIntegrationType, rIpDataType)
{
    mRealBoundaryElement = rRealBoundaryElement;
    mSurfaceEdge = CalculateSurfaceEdge(rSurfaceNode);
    mBoundaryConditionType = rBoundaryConditionType;
}

int NuTo::BoundaryGradientDamage1D::CalculateSurfaceEdge(const NodeBase* rNode) const
{
    // check edge 0:
    std::vector<const NodeBase*> node(1);
    mRealBoundaryElement->GetSurfaceNodes(0, node);
    if (node[0] == rNode)
        return 0;

    // check edge 1:
    mRealBoundaryElement->GetSurfaceNodes(1, node);
    if (node[0] == rNode)
        return 1;

    throw MechanicsException("[NuTo::BoundaryGradientDamage1D::CalculateSurfaceEdge] The element does not contain this boundary node");
}


//! @brief calculates output data fo the elmement
//! @param eOutput ... coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
//!                    @param updateStaticData (with DummyOutput), IPData, globalrow/column dofs etc.
NuTo::Error::eError NuTo::BoundaryGradientDamage1D::Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput)
{


	if (mStructure->GetHessianConstant(1)==false)
    	throw MechanicsException("[NuTo::BoundaryGradientDamage1D::Evaluate] only implemented for a constant Hessian for the first derivative (damping).");
    if (mStructure->GetHessianConstant(2)==false)
    	throw MechanicsException("[NuTo::BoundaryGradientDamage1D::Evaluate] only implemented for a constant Hessian for the second derivative (mass).");

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
			throw MechanicsException("[NuTo::BoundaryGradientDamage1D::Evaluate] no section allocated for real boundary element.");


		std::vector<double> localNodeCoordReal(mRealBoundaryElement->GetNumNodesGeometry());
        std::vector<double> nodalDisp(mRealBoundaryElement->GetNumNodesField());
        std::vector<double> nodalNonlocalEqStrain(mRealBoundaryElement->GetNumShapeFunctionsNonlocalEqStrain());

        mRealBoundaryElement->CalculateLocalCoordinates(localNodeCoordReal);
		mRealBoundaryElement->CalculateLocalDisplacements(0,nodalDisp);
        mRealBoundaryElement->CalculateNodalNonlocalEqStrain(0, nodalNonlocalEqStrain);

		int numDispReal = mRealBoundaryElement->GetNumNodesField();
		int numDispDofsReal = sectionReal->GetIsDisplacementDof() ? numDispReal : 0;

		int numNonlocalEqStrainReal = mRealBoundaryElement->GetNumShapeFunctionsNonlocalEqStrain();
		int numNonlocalEqStrainDofsReal = sectionReal->GetIsNonlocalEqStrainDof() ? numNonlocalEqStrainReal : 0;


        //allocate space for local shape functions
        std::vector<double> shapeFunctionsNonlocalEqStrain(mRealBoundaryElement->GetNumShapeFunctionsNonlocalEqStrain());
        std::vector<double> shapeFunctionsGeometry(mRealBoundaryElement->GetNumNodesGeometry());



        //define inputs ...
        DeformationGradient1D deformationGradient;
        NonlocalEqStrain nonlocalEqStrain;

        std::map< NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase* > constitutiveInputList;
        constitutiveInputList[NuTo::Constitutive::Input::DEFORMATION_GRADIENT_1D] = &deformationGradient;
        constitutiveInputList[NuTo::Constitutive::Input::NONLOCAL_EQ_STRAIN] = &nonlocalEqStrain;

        // ... and outputs
        LocalEqStrain localEqStrain;
        ConstitutiveTangentLocal<1,1> tangentLocalEqStrainStrain;
        ConstitutiveTangentLocal<1,1> tangentStressNonlocalStrain;
        ConstitutiveTangentLocal<1,1> nonlocalParameterXi;
        Damage omega;

        std::map< NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase* > constitutiveOutputList;
        constitutiveOutputList[NuTo::Constitutive::Output::LOCAL_EQ_STRAIN] = &localEqStrain;
        constitutiveOutputList[NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_D_STRAIN_1D] = &tangentLocalEqStrainStrain;
        constitutiveOutputList[NuTo::Constitutive::Output::DAMAGE] = &omega;
        constitutiveOutputList[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_1D] = &tangentStressNonlocalStrain;
        constitutiveOutputList[NuTo::Constitutive::Output::NONLOCAL_PARAMETER_XI] = &nonlocalParameterXi;





		// for the B matrix in the equation
        std::vector<double> derivativeShapeFunctionsGeometryNatural(mRealBoundaryElement->GetLocalDimension()*mRealBoundaryElement->GetNumNodesGeometry());

        std::vector<double> derivativeShapeFunctionsFieldNatural(mRealBoundaryElement->GetLocalDimension()*mRealBoundaryElement->GetNumNodesField());
        std::vector<double> derivativeShapeFunctionsFieldLocal(mRealBoundaryElement->GetLocalDimension()*mRealBoundaryElement->GetNumNodesField());


        // Loop over IP element has only one IP on the boundary

        //the ip is located on the boundary of the real element
        std::vector<double> naturalSurfaceCoordinates(1);
        mRealBoundaryElement->CalculateNaturalSurfaceCoordinates(mSurfaceEdge, 0, naturalSurfaceCoordinates);
        //allocate space for local ip coordinates
        double naturalIPCoordinate = naturalSurfaceCoordinates[0];

        //derivative in natural coordinate system
        mRealBoundaryElement->CalculateDerivativeShapeFunctionsGeometry(naturalIPCoordinate, derivativeShapeFunctionsGeometryNatural);
        mRealBoundaryElement->CalculateDerivativeShapeFunctionsField(naturalIPCoordinate, derivativeShapeFunctionsFieldNatural);

        //determinant of the Jacobian
        double detJReal = mRealBoundaryElement->DetJacobian(derivativeShapeFunctionsGeometryNatural,localNodeCoordReal);

        //derivative in local coordinate system
        for (unsigned int count=0; count<derivativeShapeFunctionsFieldNatural.size(); count++)
            derivativeShapeFunctionsFieldLocal[count] = derivativeShapeFunctionsFieldNatural[count]/detJReal;


        // deformation gradient
        mRealBoundaryElement->CalculateDeformationGradient(derivativeShapeFunctionsFieldLocal, nodalDisp, deformationGradient);
        // nonlocal eq strain
        mRealBoundaryElement->CalculateNonlocalEqStrain(shapeFunctionsNonlocalEqStrain, nodalNonlocalEqStrain, nonlocalEqStrain);

        ConstitutiveBase* constitutivePtr = GetConstitutiveLaw(0);
        Error::eError error = constitutivePtr->Evaluate1D(this, 0,
                constitutiveInputList, constitutiveOutputList);
        if (error!=Error::SUCCESSFUL)
            return error;


        // calculate normal vector

        // calculate local IP coordinates
        double localIPcoordinate = 0;
        for (unsigned int iNode = 0; iNode < localNodeCoordReal.size(); ++iNode)
        {
            localIPcoordinate += localNodeCoordReal[iNode] * shapeFunctionsGeometry[iNode];
        }
        double area = sectionReal->GetArea() * sectionReal->AsSectionTruss()->GetAreaFactor(localIPcoordinate);


        // calculate the KkkMod matrix
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> KkkMod;


        // quick and dirty solution for mcauley combined BCs:
        BoundaryCondition::eType currentBCType = mBoundaryConditionType;

        double alpha = std::sqrt(nonlocalParameterXi.GetValue(0));

        double e0 = constitutivePtr->GetTensileStrength() / constitutivePtr->GetYoungsModulus();
        if (mBoundaryConditionType == BoundaryCondition::MACAULAY)
        {
            // determine state ...
            bool switchToNeumann = (localEqStrain[0] > nonlocalEqStrain[0]) and (localEqStrain[0] > e0);

            // ... and use existing implementations
            if(switchToNeumann)
            {
                currentBCType = BoundaryCondition::NEUMANN_HOMOGENEOUS;
                std::cout << "Macaulay Culcin helps out in element " <<
                        mRealBoundaryElement->GetStructure()->ElementGetId(this) << std::endl;
            }
            else
                currentBCType = BoundaryCondition::ROBIN_INHOMOGENEOUS;
        }



        switch (currentBCType)
        {
            case BoundaryCondition::NOT_SET:
            {
                throw MechanicsException("[NuTo::BoundaryGradientDamage1D::Evaluate] Boundary condition type not set! ");
            }
            break;
            case BoundaryCondition::NEUMANN_HOMOGENEOUS:
            {
                KkkMod.resize(shapeFunctionsNonlocalEqStrain.size(), shapeFunctionsNonlocalEqStrain.size());
                KkkMod.setZero();
            }
            break;
            case BoundaryCondition::ROBIN_INHOMOGENEOUS:
            {
                CalculateKkkMod(
                        shapeFunctionsNonlocalEqStrain,
                        shapeFunctionsNonlocalEqStrain,
                        1, 1./alpha, area, KkkMod);
            }
            break;
            default:
                break;
        }

        int numDofs = numDispDofsReal+numNonlocalEqStrainDofsReal;
        for (auto it = rElementOutput.begin(); it!=rElementOutput.end(); it++)
        {
            switch(it->first)
            {
            case Element::INTERNAL_GRADIENT:
            {
                NuTo::FullVector<double, Eigen::Dynamic>& f = it->second->GetFullVectorDouble();
                f.Resize(numDofs);

                if(currentBCType == BoundaryCondition::ROBIN_INHOMOGENEOUS)
                {
                    // add Nt c eqStrain
//                    std::cout << localEqStrain(0) << std::endl;
                    for (int iRow = 0; iRow < numNonlocalEqStrainDofsReal; ++iRow)
                        f(numDispDofsReal+iRow) +=
                                shapeFunctionsNonlocalEqStrain[iRow] *
                                area/alpha*
                                (nonlocalEqStrain[0] - localEqStrain[0]);

                } else {
                    // add Kkk * nonlocalEqStrain
                    for (int iRow = 0; iRow < numNonlocalEqStrainDofsReal; ++iRow)
                    {
                        for (int iCol = 0; iCol < numNonlocalEqStrainDofsReal; ++iCol)
                        {
                            f(numDispDofsReal+iRow) += KkkMod(iRow,iCol) * nodalNonlocalEqStrain[iCol];
                        }
                    }
                }
            }
            break;
            case Element::HESSIAN_0_TIME_DERIVATIVE:
            {
                it->second->GetFullMatrixDouble().Resize(numDofs, numDofs);
                it->second->SetSymmetry(false);
                it->second->SetConstant(false);
                if (numDispDofsReal>0 and numNonlocalEqStrainDofsReal>0)
                {

                    // add the modified Kkk to the element output
                    it->second->GetFullMatrixDouble().AddBlock(numDispDofsReal, numDispDofsReal, KkkMod);
                    if(currentBCType == BoundaryCondition::ROBIN_INHOMOGENEOUS)
                    {
                        mRealBoundaryElement->AddDetJNtdLocalEqStraindEpsilonB(
                                shapeFunctionsNonlocalEqStrain,
                                tangentLocalEqStrainStrain,
                                derivativeShapeFunctionsFieldLocal,
                                area/alpha, numDispDofsReal,0,
                                it->second->GetFullMatrixDouble());
                    }

                }
            }
            break;
            case Element::GLOBAL_ROW_DOF:
                CalculateGlobalRowDofs(it->second->GetVectorInt(),numDispDofsReal, numNonlocalEqStrainDofsReal);
                break;
            case Element::GLOBAL_COLUMN_DOF:
                CalculateGlobalRowDofs(it->second->GetVectorInt(),numDispDofsReal, numNonlocalEqStrainDofsReal);
                break;
            default:

                break;
          }
        }

    }
    catch (NuTo::MechanicsException& e)
    {
        std::stringstream ss;
        ss << mStructure->ElementGetId(this);
    	e.AddMessage("[NuTo::BoundaryGradientDamage1D::Evaluate] Error evaluating element data of element"	+ ss.str() + ".");
        throw e;
    }

    return Error::SUCCESSFUL;
}

//! @brief Allocates static data for an integration point of an element
//! @param rConstitutiveLaw constitutive law, which is called to allocate the static data object
NuTo::ConstitutiveStaticDataBase* NuTo::BoundaryGradientDamage1D::AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw)const
{
	return rConstitutiveLaw->AllocateStaticDataEngineeringStress_EngineeringStrain1D(this);
}



// build global row dofs
void NuTo::BoundaryGradientDamage1D::CalculateGlobalRowDofs(
        std::vector<int>& rGlobalRowDofs,
        int rNumDispDofs,
        int rNumNonlocalEqStrainDofs) const
{
    rGlobalRowDofs.resize(rNumDispDofs+rNumNonlocalEqStrainDofs);



    // indices for the current dof type
    int iGlobalDisp                    = 0;
    int iGlobalNonlocalEqStrain        = 0;

    int numNodes(mRealBoundaryElement->GetNumNodes());

    for (int iNode = 0; iNode < numNodes; iNode++)
    {
        const NodeBase * node = mRealBoundaryElement->GetNode(iNode);

        if (rNumDispDofs>0)
        {
            for (int iLocalLocal = 0; iLocalLocal < node->GetNumDisplacements(); iLocalLocal++)
            {
                rGlobalRowDofs[iGlobalDisp] = node->GetDofDisplacement(iLocalLocal);
                iGlobalDisp++;
            }
        }

        int iStartNonlocalEqStrain = rNumDispDofs;
        if (node->GetNumNonlocalEqStrain()>0 && rNumNonlocalEqStrainDofs>0)
        {
            rGlobalRowDofs[iStartNonlocalEqStrain + iGlobalNonlocalEqStrain] = node->GetDofNonlocalEqStrain();
            iGlobalNonlocalEqStrain ++;
        }
    }

}

//! @brief calculates the boundary integral of Nt * c * n * B
//! @param shapeFunctions of the ip for all shape functions
//! @param derivativeShapeFunctions of the ip for all shape functions
//! @param c nonlocal gradient radius
//! @param factor multiplication factor (detJ area..)
//! @param Kkkmod return matrix with detJ * Nt * c * n * B
void NuTo::BoundaryGradientDamage1D::CalculateKkkMod(
        const std::vector<double>& rShapeFunctions,
        const std::vector<double>& rDerivativeShapeFunctions,
        double rNonlocalGradientRadius, double rNormalVector, double rfactor,
		FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rKkkMod) const
{
	//resize and set to zero
	rKkkMod.Resize(rShapeFunctions.size(),rShapeFunctions.size());
	rKkkMod.setZero();
    //add NtN
	for (unsigned int iRow=0; iRow<rShapeFunctions.size();iRow++)
	{
		for (unsigned int iCol=0; iCol<rShapeFunctions.size();iCol++)
		{
			rKkkMod(iRow,iCol)=+rfactor*
					      rShapeFunctions[iRow]*
					      rNonlocalGradientRadius*
					      rNormalVector*rDerivativeShapeFunctions[iCol];

		}
	}
}

void NuTo::BoundaryGradientDamage1D::CalculateKedMod(
        const std::vector<double>& rShapeFunctions,
        const std::vector<double>& rDerivativeShapeFunctions,
        double rNonlocalGradientRadius, double rfactor,
        FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rKedMod) const
{
    //resize and set to zero
    rKedMod.Resize(rShapeFunctions.size(),rDerivativeShapeFunctions.size());
    rKedMod.setZero();
    //add NtN
    for (unsigned int iRow=0; iRow<rShapeFunctions.size();iRow++)
    {
        for (unsigned int iCol=0; iCol<rShapeFunctions.size();iCol++)
        {
            rKedMod(iRow,iCol)=-rfactor*
                          rShapeFunctions[iRow]*
                          rNonlocalGradientRadius*
                          rDerivativeShapeFunctions[iCol];

        }
    }
}

//! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
//! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
void NuTo::BoundaryGradientDamage1D::ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
{
    NodeBase* node = GetNode(0);
    if (node == rOldPtr)
    {
        node = rNewPtr;
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
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementBase)
           & BOOST_SERIALIZATION_NVP(mRealBoundaryElement)
           & BOOST_SERIALIZATION_NVP(mSurfaceEdge)
           & BOOST_SERIALIZATION_NVP(mBoundaryConditionType);

    }
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize BoundaryGradientDamage1D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::BoundaryGradientDamage1D)
#endif // ENABLE_SERIALIZATION
