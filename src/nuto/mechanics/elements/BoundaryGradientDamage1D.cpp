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

void NuTo::BoundaryGradientDamage1D::ApplyConstraints(BoundaryCondition::eType rType, StructureBase* rStructure)
{
    mBoundaryConditionType = rType;


    // the boundary conditions are fulfilled in a weak sense, see Evaluate()
    return;




    // get integration point
    std::vector<double> naturalSurfaceCoordinates(1);
    mRealBoundaryElement->CalculateNaturalSurfaceCoordinates(mSurfaceEdge,0,naturalSurfaceCoordinates);
    double ip = naturalSurfaceCoordinates[0];

    // allocate misc stuff
    std::vector<double> localNodeCoord(mRealBoundaryElement->GetNumNodesGeometry());
    std::vector<double> nodeDisp(mRealBoundaryElement->GetNumNodesField());
    std::vector<double> nodeNonlocalEqStrain(mRealBoundaryElement->GetNumShapeFunctionsNonlocalEqStrain());

    std::vector<double> derivativeShapeFunctionsGeometryNatural(mRealBoundaryElement->GetLocalDimension()*mRealBoundaryElement->GetNumNodesGeometry());

    std::vector<double> derivativeShapeFunctionsFieldNatural(mRealBoundaryElement->GetLocalDimension()*mRealBoundaryElement->GetNumNodesField());
    std::vector<double> derivativeShapeFunctionsFieldLocal(mRealBoundaryElement->GetLocalDimension()*mRealBoundaryElement->GetNumNodesField());

    std::vector<double> derivativeShapeFunctionsNonlocalEqStrainNatural(mRealBoundaryElement->GetLocalDimension()*mRealBoundaryElement->GetNumShapeFunctionsNonlocalEqStrain());
    std::vector<double> derivativeShapeFunctionsNonlocalEqStrainLocal(mRealBoundaryElement->GetLocalDimension()*mRealBoundaryElement->GetNumShapeFunctionsNonlocalEqStrain());

    // get nodal values
    mRealBoundaryElement->CalculateLocalCoordinates(localNodeCoord);
    mRealBoundaryElement->CalculateLocalDisplacements(0, nodeDisp);
    mRealBoundaryElement->CalculateNodalNonlocalEqStrain(0, nodeNonlocalEqStrain);


    // get natural derivative shape functions and calculate determinant of the Jacobian
    mRealBoundaryElement->CalculateDerivativeShapeFunctionsGeometry(ip, derivativeShapeFunctionsGeometryNatural);
    mRealBoundaryElement->CalculateDerivativeShapeFunctionsField(ip, derivativeShapeFunctionsFieldNatural);
    mRealBoundaryElement->CalculateDerivativeShapeFunctionsNonlocalEqStrain(ip, derivativeShapeFunctionsNonlocalEqStrainNatural);

    double detJ = mRealBoundaryElement->DetJacobian(derivativeShapeFunctionsGeometryNatural, localNodeCoord);


    // get surface node
    std::vector<const NodeBase*> surfaceNodePtr(1);
    mRealBoundaryElement->GetSurfaceNodes(mSurfaceEdge, surfaceNodePtr);
    int surfaceNode = mStructure->NodeGetId(surfaceNodePtr[0]);

    double sqrtNonlocalParameterC = std::sqrt(GetConstitutiveLaw(0)->GetNonlocalRadius());


    for (unsigned int i = 0; i < derivativeShapeFunctionsFieldNatural.size(); ++i)
        derivativeShapeFunctionsFieldLocal[i] = derivativeShapeFunctionsFieldNatural[i] / detJ;

    for (unsigned int i = 0; i < derivativeShapeFunctionsNonlocalEqStrainNatural.size(); ++i)
        derivativeShapeFunctionsNonlocalEqStrainLocal[i] = derivativeShapeFunctionsNonlocalEqStrainNatural[i] / detJ;


    switch (mBoundaryConditionType)
    {
        case BoundaryCondition::NEUMANN_HOMOGENEOUS:
        {
            // BC: grad nonlocalEqStrain * n = 0
            // B*nonlocalEqStrain = 0

            // create constraint with first term:
            int nodeId = mStructure->NodeGetId(mRealBoundaryElement->GetNodeNonlocalEqStrain(0));
            int constraint = rStructure->ConstraintLinearEquationCreate(nodeId, "nonlocalEqStrain", derivativeShapeFunctionsNonlocalEqStrainLocal[0], 0.);

            // create other terms via loop:
            for (int iEqStrainNode = 1; iEqStrainNode < mRealBoundaryElement->GetNumShapeFunctionsNonlocalEqStrain(); ++iEqStrainNode)
            {
                nodeId = mStructure->NodeGetId(mRealBoundaryElement->GetNodeNonlocalEqStrain(iEqStrainNode));
                rStructure->ConstraintLinearEquationAddTerm(constraint, nodeId, "nonlocalEqStrain", derivativeShapeFunctionsNonlocalEqStrainLocal[iEqStrainNode]);
            }
        }
        break;

        case BoundaryCondition::DIRICHLET_INHOMOGENEOUS:
        {
            // BC: local eq strains - nonlocal eq strains = 0
            // get surface node index


             // add constraints
            int constraint = rStructure->ConstraintLinearEquationCreate(surfaceNode, "nonlocalEqStrain", -1., 0);

            for (int iNodeDisp = 0; iNodeDisp < mRealBoundaryElement->GetNumNodesField(); ++iNodeDisp)
            {
                const NodeBase* nodeDisp = mRealBoundaryElement->GetNodeField(iNodeDisp);
                int nodeIndex = mStructure->NodeGetId(nodeDisp);
                rStructure->ConstraintLinearEquationAddTerm(constraint, nodeIndex, "X_DISPLACEMENT", derivativeShapeFunctionsFieldLocal[iNodeDisp]);
            }
        }
        break;

        case BoundaryCondition::ROBIN_INHOMOGENEOUS:
        {
            // BC: l*grad nonlocal eq strain * n + nonlocal eq strain - eq strain = 0
            // basically a combination of type 0 and 1


            double normalVector = mSurfaceEdge == 0 ? -1 : 1;

            // 0 = nonlocal eq strain ...
            int constraint = rStructure->ConstraintLinearEquationCreate(surfaceNode, "nonlocalEqStrain", 1., 0);

            // + n * B * nonlocal strain
            for (int iEqStrainNode = 0; iEqStrainNode < mRealBoundaryElement->GetNumShapeFunctionsNonlocalEqStrain(); ++iEqStrainNode)
            {
                int nodeId = mStructure->NodeGetId(mRealBoundaryElement->GetNodeNonlocalEqStrain(iEqStrainNode));
                rStructure->ConstraintLinearEquationAddTerm(constraint, nodeId, "nonlocalEqStrain", sqrtNonlocalParameterC*derivativeShapeFunctionsNonlocalEqStrainLocal[iEqStrainNode]*normalVector);
            }

            // - B * displacement
            for (int iNodeDisp = 0; iNodeDisp < mRealBoundaryElement->GetNumNodesField(); ++iNodeDisp)
            {
                const NodeBase* nodeDisp = mRealBoundaryElement->GetNodeField(iNodeDisp);
                int nodeIndex = mStructure->NodeGetId(nodeDisp);
                rStructure->ConstraintLinearEquationAddTerm(constraint, nodeIndex, "X_DISPLACEMENT", -derivativeShapeFunctionsFieldLocal[iNodeDisp]);
            }
        }
        break;
        case BoundaryCondition::ROBIN_HOMOGENEOUS:
        {
            // BC: l*grad nonlocal eq strain * n + nonlocal eq strain  = 0

            double normalVector = mSurfaceEdge == 0 ? -1 : 1;

            // 0 = nonlocal eq strain ...
            int constraint = rStructure->ConstraintLinearEquationCreate(surfaceNode, "nonlocalEqStrain", 1., 0);

            // + n * B * nonlocal strain
            for (int iEqStrainNode = 0; iEqStrainNode < mRealBoundaryElement->GetNumShapeFunctionsNonlocalEqStrain(); ++iEqStrainNode)
            {
                int nodeId = mStructure->NodeGetId(mRealBoundaryElement->GetNodeNonlocalEqStrain(iEqStrainNode));
                rStructure->ConstraintLinearEquationAddTerm(constraint, nodeId, "nonlocalEqStrain", sqrtNonlocalParameterC*derivativeShapeFunctionsNonlocalEqStrainLocal[iEqStrainNode]*normalVector);
            }

        }
        break;
        case BoundaryCondition::MACAULAY:
            // the boundary conditions are fulfilled by the system. TODO understand exactly why...
            break;
        default:
            break;
    }



//    if(rType == 3)
//    {
//        // grad eq strain - grad nonlocal eq strain = 0
//        std::vector<double> d2ShapeFunctions(3); // hardcode!
//        d2ShapeFunctions[0] =  1 / detJ / detJ;
//        d2ShapeFunctions[1] = -2. / detJ / detJ;
//        d2ShapeFunctions[2] =  1 / detJ / detJ;
//
//        // create constraint with first term:
//        int nodeId = mStructure->NodeGetId(mRealBoundaryElement->GetNodeNonlocalEqStrain(0));
//        int constraint = rStructure->ConstraintLinearEquationCreate(nodeId, "nonlocalEqStrain", derivativeShapeFunctionsNonlocalEqStrainLocal[0], 0.);
//
//        // create other terms via loop:
//        for (int iEqStrainNode = 1; iEqStrainNode < mRealBoundaryElement->GetNumShapeFunctionsNonlocalEqStrain(); ++iEqStrainNode)
//        {
//            nodeId = mStructure->NodeGetId(mRealBoundaryElement->GetNodeNonlocalEqStrain(iEqStrainNode));
//            rStructure->ConstraintLinearEquationAddTerm(constraint, nodeId, "nonlocalEqStrain", derivativeShapeFunctionsNonlocalEqStrainLocal[iEqStrainNode]);
//        }
//
//        for (int iNodeDisp = 0; iNodeDisp < mRealBoundaryElement->GetNumNodesField(); ++iNodeDisp)
//        {
//            const NodeBase* nodeDisp = mRealBoundaryElement->GetNodeField(iNodeDisp);
//            int nodeIndex = mStructure->NodeGetId(nodeDisp);
//            rStructure->ConstraintLinearEquationAddTerm(constraint, nodeIndex, "X_DISPLACEMENT", -d2ShapeFunctions[iNodeDisp]);
//        }
//
//    }


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
        ConstitutiveTangentLocal<1,1> variableNonlocalParameter;
        Damage omega;

        std::map< NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase* > constitutiveOutputList;
        constitutiveOutputList[NuTo::Constitutive::Output::LOCAL_EQ_STRAIN] = &localEqStrain;
        constitutiveOutputList[NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_D_STRAIN_1D] = &tangentLocalEqStrainStrain;
        constitutiveOutputList[NuTo::Constitutive::Output::DAMAGE] = &omega;
        constitutiveOutputList[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN_1D] = &tangentStressNonlocalStrain;
        constitutiveOutputList[NuTo::Constitutive::Output::VARIABLE_NONLOCAL_RADIUS] = &variableNonlocalParameter;





		// for the B matrix in the equation
		std::vector<double> derivativeShapeFunctionsNonlocalEqStrainNatural(mRealBoundaryElement->GetLocalDimension()*mRealBoundaryElement->GetNumShapeFunctionsNonlocalEqStrain());
        std::vector<double> derivativeShapeFunctionsNonlocalEqStrainLocal(mRealBoundaryElement->GetLocalDimension()*mRealBoundaryElement->GetNumShapeFunctionsNonlocalEqStrain());

        std::vector<double> derivativeShapeFunctionsGeometryNatural(mRealBoundaryElement->GetLocalDimension()*mRealBoundaryElement->GetNumNodesGeometry());
        std::vector<double> derivativeShapeFunctionsGeometryLocal(mRealBoundaryElement->GetLocalDimension()*mRealBoundaryElement->GetNumNodesGeometry());

        std::vector<double> derivativeShapeFunctionsFieldNatural(mRealBoundaryElement->GetLocalDimension()*mRealBoundaryElement->GetNumNodesField());
        std::vector<double> derivativeShapeFunctionsFieldLocal(mRealBoundaryElement->GetLocalDimension()*mRealBoundaryElement->GetNumNodesField());


        // element has only one IP on the boundary

        //the ip is located on the boundary of the real element
        std::vector<double> naturalSurfaceCoordinates(1);

        mRealBoundaryElement->CalculateNaturalSurfaceCoordinates(mSurfaceEdge, 0, naturalSurfaceCoordinates);
        //allocate space for local ip coordinates
        double naturalIPCoordinate = naturalSurfaceCoordinates[0];

        mRealBoundaryElement->CalculateShapeFunctionsGeometry(naturalIPCoordinate, shapeFunctionsGeometry);
        mRealBoundaryElement->CalculateShapeFunctionsNonlocalEqStrain(naturalIPCoordinate, shapeFunctionsNonlocalEqStrain);
        //derivative in natural coordinate system
        mRealBoundaryElement->CalculateDerivativeShapeFunctionsNonlocalEqStrain(naturalIPCoordinate, derivativeShapeFunctionsNonlocalEqStrainNatural);
        mRealBoundaryElement->CalculateDerivativeShapeFunctionsGeometry(naturalIPCoordinate, derivativeShapeFunctionsGeometryNatural);
        mRealBoundaryElement->CalculateDerivativeShapeFunctionsField(naturalIPCoordinate, derivativeShapeFunctionsFieldNatural);

        //determinant of the Jacobian
        double detJReal = mRealBoundaryElement->DetJacobian(derivativeShapeFunctionsGeometryNatural,localNodeCoordReal);

        //derivative in local coordinate system
        for (unsigned int count=0; count<derivativeShapeFunctionsNonlocalEqStrainNatural.size(); count++)
            derivativeShapeFunctionsNonlocalEqStrainLocal[count] = derivativeShapeFunctionsNonlocalEqStrainNatural[count]/detJReal;

        for (unsigned int count=0; count<derivativeShapeFunctionsGeometryNatural.size(); count++)
            derivativeShapeFunctionsGeometryLocal[count] = derivativeShapeFunctionsGeometryNatural[count]/detJReal;

        for (unsigned int count=0; count<derivativeShapeFunctionsFieldNatural.size(); count++)
            derivativeShapeFunctionsFieldLocal[count] = derivativeShapeFunctionsFieldNatural[count]/detJReal;


        // deformation gradient
        mRealBoundaryElement->CalculateDeformationGradient(derivativeShapeFunctionsGeometryLocal, nodalDisp, deformationGradient);
        // nonlocal eq strain
        mRealBoundaryElement->CalculateNonlocalEqStrain(shapeFunctionsNonlocalEqStrain, nodalNonlocalEqStrain, nonlocalEqStrain);

        ConstitutiveBase* constitutivePtr = GetConstitutiveLaw(0);
        Error::eError error = constitutivePtr->Evaluate1D(this, 0,
                constitutiveInputList, constitutiveOutputList);
        if (error!=Error::SUCCESSFUL)
            return error;


        double E = constitutivePtr->GetYoungsModulus();
        EngineeringStrain1D strain1D;
        deformationGradient.GetEngineeringStrain(strain1D);


        double dOmega = - tangentStressNonlocalStrain[0] / (E * strain1D[0]);

        // calculate normal vector
        double normalVector = mSurfaceEdge == 0 ? -1 : 1;

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

        double alpha = std::sqrt(variableNonlocalParameter(0,0));
//        std::cout << variableNonlocalParameter[0] << std::endl;
//        std::cout << variableNonlocalParameter(0,0) << std::endl;

//        double alpha = 1;//std::sqrt(constitutivePtr->GetNonlocalRadius());


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
//                std::cout << localEqStrain[0] << "\t" << nonlocalEqStrain[0] << std::endl;
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
            case BoundaryCondition::DIRICHLET_INHOMOGENEOUS:
            {
                CalculateKkkMod(
                        shapeFunctionsNonlocalEqStrain,
                        derivativeShapeFunctionsNonlocalEqStrainLocal,
                        1, normalVector, -area, KkkMod);
            }
            break;
            case BoundaryCondition::ROBIN_INHOMOGENEOUS:
            case BoundaryCondition::ROBIN_HOMOGENEOUS:
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
                        tangentLocalEqStrainStrain[0] = 1; // TODO fix that!
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

//                std::cout << "I want to calculate other stuff as well! " << std::endl;
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
