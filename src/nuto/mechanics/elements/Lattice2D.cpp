// $Id:$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/mechanics/DeformationGradient2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/LatticeStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/LatticeStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/LatticeStress2D.h"
#include "nuto/mechanics/constitutive/mechanics/LatticeStress3D.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal2x2.h"
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/elements/Lattice2D.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveEngineeringStressStrain.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveLatticeStressStrain.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/structures/StructureBase.h"

//just for test purpose
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataLatticeConcrete2D.h"

#include "nuto/math/FullMatrix.h"
#include <eigen3/Eigen/Dense>

//! @brief constructor
NuTo::Lattice2D::Lattice2D(const NuTo::StructureBase* rStructure,
		std::vector<NuTo::NodeBase* >& rNodes,
		ElementData::eElementDataType rElementDataType,
		IpData::eIpDataType rIpDataType) :
		NuTo::ElementBase::ElementBase(rStructure, rElementDataType, 3, rIpDataType)
{
    mSection = 0;

    mNodes[0] = rNodes[0];
    mNodes[1] = rNodes[1];
    mNodes[2] = rNodes[2];

    //calculate weights and coordinates of the integration points
    //get coordinates
    boost::array<boost::array<double, 2>,3 > coordinatesNodes;
    for (int count=0; count<3; count++)
    {
    	mNodes[count]->GetCoordinates2D(&(coordinatesNodes[count][0]));
    }

    //calculate global coordinates of the edge points and face points
    boost::array<boost::array<double, 2>,3 > globalCoordinatesEdgePoints;
    boost::array<double, 2> globalCoordinatesFacePoint;
    CalculateGlobalCoordinatesEdgeFacePoints(globalCoordinatesEdgePoints, globalCoordinatesFacePoint);

    //calculate help values for the transformation global local
    Eigen::Matrix<double,2,2> Ainv;
    boost::array<double, 2> coordinatesNode0;
    CalculateHelpTransformationGlobalLocal(Ainv, coordinatesNode0);

    for (int theIp=0; theIp<3; theIp++)
    {
        boost::array<double, 2> coordinatesIpGlobal;
    	coordinatesIpGlobal[0]= 0.5*(globalCoordinatesFacePoint[0]+globalCoordinatesEdgePoints[theIp][0]);
    	coordinatesIpGlobal[1]= 0.5*(globalCoordinatesFacePoint[1]+globalCoordinatesEdgePoints[theIp][1]);

        boost::array<double, 2> coordinatesIpLocal;

        //std::cout << "coordinatesIpGlobal " << coordinatesIpGlobal[0] << " " << coordinatesIpGlobal[1] << "\n";

        CalculateLocalCoordinateFromGlobal(Ainv, coordinatesNode0 ,
        		coordinatesIpGlobal, coordinatesIpLocal);

        //std::cout << "coordinatesIpLocal " << coordinatesIpLocal[0] << " " << coordinatesIpLocal[1] << "\n";

    	mElementData->SetLocalIntegrationPointCoordinates2D(theIp,coordinatesIpLocal);

//boost::array<double, 3> coordinatesIpGlobalTest;
//InterpolateCoordinatesFrom2D(coordinatesIpLocal,coordinatesIpGlobalTest);
//std::cout << "coordinatesIpGlobalTest " << coordinatesIpGlobalTest[0] << " " << coordinatesIpGlobalTest[1] << "\n";
    	//project the facet (connection edge point face point onto the normal to the edge
    	boost::array<double, 2> deltaEdge;
    	deltaEdge[0] = coordinatesNodes[(theIp+1)%3][0]-coordinatesNodes[theIp][0];
    	deltaEdge[1] = coordinatesNodes[(theIp+1)%3][1]-coordinatesNodes[theIp][1];

    	boost::array<double, 2> deltaFacet;
    	deltaFacet[0] = globalCoordinatesFacePoint[0]-globalCoordinatesEdgePoints[theIp][0];
    	deltaFacet[1] = globalCoordinatesFacePoint[1]-globalCoordinatesEdgePoints[theIp][1];

    	double l2 = (deltaEdge[0]*deltaEdge[0]+deltaEdge[1]*deltaEdge[1]);
    	double p = (deltaEdge[0]*deltaFacet[0]+deltaEdge[1]*deltaFacet[1])/l2;
    	deltaFacet[0] -= p * deltaEdge[0];
    	deltaFacet[1] -= p * deltaEdge[1];

    	double w = sqrt(deltaFacet[0]*deltaFacet[0] + deltaFacet[1]*deltaFacet[1])*sqrt(l2);
    	mElementData->SetIntegrationPointWeight(theIp,w);
    }
}

void NuTo::Lattice2D::CalculateHelpTransformationGlobalLocal(
		Eigen::Matrix<double,2,2>& rAinv,
		boost::array<double, 2>& rCoordinatesNode0)const
{
    //get coordinates
    boost::array<boost::array<double, 2>,3 > coordinatesNodes;
    for (int count=0; count<3; count++)
    {
    	mNodes[count]->GetCoordinates2D(&(coordinatesNodes[count][0]));
    }

    //calculate help matrix A inverse
    Eigen::Matrix<double,2,2> A;
	A(0,0) = -coordinatesNodes[0][0] +coordinatesNodes[1][0];
	A(0,1) = -coordinatesNodes[0][0] +coordinatesNodes[2][0];
	A(1,0) = -coordinatesNodes[0][1] +coordinatesNodes[1][1];
	A(1,1) = -coordinatesNodes[0][1] +coordinatesNodes[2][1];
    rAinv = A.inverse();
    rCoordinatesNode0[0] = coordinatesNodes[0][0];
    rCoordinatesNode0[1] = coordinatesNodes[0][1];
}


void NuTo::Lattice2D::CalculateLocalCoordinateFromGlobal(
		const Eigen::Matrix<double,2,2>& Ainv,
		const boost::array<double, 2>& rCoordinatesNode0,
		const boost::array<double, 2>& rGlobalCoordinates,
		boost::array<double, 2>& rLocalCoordinates)const
{
    //calculate local point coordinates
	Eigen::Matrix<double,2,1> b;

   	b(0,0) = rGlobalCoordinates[0] -rCoordinatesNode0[0];
   	b(1,0) = rGlobalCoordinates[1] -rCoordinatesNode0[1];
	Eigen::Matrix<double,2,1> sol(Ainv*b);
	rLocalCoordinates[0] = sol(0,0);
	rLocalCoordinates[1] = sol(1,0);
}

void NuTo::Lattice2D::CalculateGlobalCoordinatesEdgeFacePoints(
		boost::array<boost::array<double, 2>,3 >& rCoordinatesEdgePoints,
		boost::array<double, 2>& rCoordinateFacePoint)const
{
	//calculate weights and coordinates of the integration points
	//get coordinates and radii
	boost::array<boost::array<double, 2>,3 > coordinatesNodes;
	std::array<double, 3 > radius;
	for (int count=0; count<3; count++)
	{
		mNodes[count]->GetCoordinates2D(&(coordinatesNodes[count][0]));
		mNodes[count]->GetRadius(&(radius[count]));
	}

	//calculate edge points Eij
	for (int count=0; count<3; count++)
	{
		boost::array<double, 2> delta;
		delta[0] = coordinatesNodes[(count+1)%3][0]-coordinatesNodes[count][0];
		delta[1] = coordinatesNodes[(count+1)%3][1]-coordinatesNodes[count][1];
		double lEdge = sqrt(delta[0]*delta[0]+delta[1]*delta[1]);
		double w1 = (radius[count]+0.5*(lEdge - radius[count] - radius[(count+1)%3]))/lEdge;
		rCoordinatesEdgePoints[count][0]=coordinatesNodes[count][0] + w1*delta[0];
		rCoordinatesEdgePoints[count][1]=coordinatesNodes[count][1] + w1*delta[1];
	}

	//calculate tmp face points Fij
	boost::array<boost::array<double, 2>,3 > coordinatesTmpFacePoints;
	for (int count=0; count<3; count++)
	{
		boost::array<double, 2> delta;
		delta[0] = coordinatesNodes[(count+2)%3][0]-rCoordinatesEdgePoints[count][0];
		delta[1] = coordinatesNodes[(count+2)%3][1]-rCoordinatesEdgePoints[count][1];
		double lEdge = sqrt(delta[0]*delta[0]+delta[1]*delta[1]);
		double w1 = (0.5*(lEdge - radius[(count+2)%3]))/lEdge;
		coordinatesTmpFacePoints[count][0]=rCoordinatesEdgePoints[count][0] + w1*delta[0];
		coordinatesTmpFacePoints[count][1]=rCoordinatesEdgePoints[count][1] + w1*delta[1];
	}

	//calculate face points Fi
	rCoordinateFacePoint[0]=0.;
	rCoordinateFacePoint[1]=0.;
	for (int count=0; count<3; count++)
	{
		rCoordinateFacePoint[0] += coordinatesTmpFacePoints[count][0];
		rCoordinateFacePoint[1] += coordinatesTmpFacePoints[count][1];
	}
	rCoordinateFacePoint[0] /= 3.;
	rCoordinateFacePoint[1] /= 3.;
}


// old version left for comparison in here (to see if there is a significant speed advantage by not calling a subroutine
// calculating the strain directly without the multiplication with the bmatrix
//! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the stiffness matrix
/*NuTo::Error::eError NuTo::Lattice2D::CalculateCoefficientMatrix_0(NuTo::FullMatrix<double>& rCoefficientMatrix,
        std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn, bool& rSymmetry)const
{
for (int count=0; count<100000; count++)
{
    //allocate deformation gradient
    LatticeStrain2D latticeStrain;
    ConstitutiveTangentLocal2x2 stiffness;

    bool areAllIpsSymmetric(true);
    rCoefficientMatrix.Resize(9,9);
    for (int theIP=0; theIP<3; theIP++)
    {
    	//get the edge
    	NodeBase *node1,*node2;
   		node1 = mNodes[theIP];
   		node2 = mNodes[(theIP+1)%3];

   	    //get coordinates, displacements and radius
   		boost::array<double,2> coordNode1, dispNode1, coordNode2, dispNode2;
   		double rotNode1,rotNode2;
		node1->GetCoordinates2D(&(coordNode1[0]));
		node1->GetDisplacements2D(&(dispNode1[0]));
		node1->GetRotations2D(&rotNode1);
		node2->GetCoordinates2D(&(coordNode2[0]));
		node2->GetDisplacements2D(&(dispNode2[0]));
		node2->GetRotations2D(&rotNode2);

   		//global coordinates of center of facet
   		boost::array<double,3> globalCoordinatesFacet;
   		boost::array<double,2> localCoordinatesFacet;
   		this->mElementData->GetLocalIntegrationPointCoordinates2D(theIP, localCoordinatesFacet);
   		//Get the Weight of an IP from the ip data
   		double weight = this->mElementData->GetIntegrationPointWeight(theIP);

   		InterpolateCoordinatesFrom2D(localCoordinatesFacet, globalCoordinatesFacet);

   		//shape functions for the rotations of the two nodes of the corresponding edge
   		boost::array<double,2> shapeFunctionsRotations1,shapeFunctionsRotations2;
   		shapeFunctionsRotations1 = CalculateShapeFunctionsRotations(globalCoordinatesFacet,coordNode1);
   		shapeFunctionsRotations2 = CalculateShapeFunctionsRotations(globalCoordinatesFacet,coordNode2);

   		//calculate the gap at the facet center
   		boost::array<double,2> deltaCrack;
   		deltaCrack[0] = dispNode2[0] + shapeFunctionsRotations2[0]*rotNode2-
   		                dispNode1[0] - shapeFunctionsRotations1[0]*rotNode1;
   		deltaCrack[1] = dispNode2[1] + shapeFunctionsRotations2[1]*rotNode2-
   		                dispNode1[1] - shapeFunctionsRotations1[1]*rotNode1;

    	//edge vector
    	double edge[2];
    	edge[0] = coordNode2[0]-coordNode1[0];
    	edge[1] = coordNode2[1]-coordNode1[1];

    	//edge normal
    	double edgeNormal[2];
    	double lEdge=sqrt(edge[0]*edge[0]+edge[1]*edge[1]);
    	double lEdgeInv = 1./lEdge;
    	edgeNormal[0] = edge[0]*lEdgeInv;
    	edgeNormal[1] = edge[1]*lEdgeInv;

    	//calculate lattice strain vector (epsilonN, epsilonT : normal and tangent strain
    	latticeStrain.mLatticeStrain[0] = (deltaCrack[0]*edgeNormal[0]+deltaCrack[1]*edgeNormal[1])*lEdgeInv;
    	latticeStrain.mLatticeStrain[1] = (deltaCrack[1]*edgeNormal[0]-deltaCrack[0]*edgeNormal[1])*lEdgeInv;

    	//calculate constitutive law
        //material pointer
        const ConstitutiveLatticeStressStrain *constitutivePtr;
        constitutivePtr = GetConstitutiveLaw(theIP)->AsConstitutiveLatticeStressStrain();
		Error::eError error = constitutivePtr->GetTangent_LatticeStress_LatticeStrain(this, theIP,
				latticeStrain, &stiffness);
		if (error!=Error::SUCCESSFUL)
			return error;

		//calculate bMatrix for first node
		Eigen::Matrix<double,2,6> bMatrix;
		bMatrix(0,0)=-edgeNormal[0];
		bMatrix(0,1)=-edgeNormal[1];
		bMatrix(0,2)=-edgeNormal[0]*shapeFunctionsRotations1[0]-edgeNormal[1]*shapeFunctionsRotations1[1];
		bMatrix(1,0)=+edgeNormal[1];
		bMatrix(1,1)=-edgeNormal[0];
		bMatrix(1,2)=+edgeNormal[1]*shapeFunctionsRotations1[0]-edgeNormal[0]*shapeFunctionsRotations1[1];

		//calculate bMatrix for second node
		bMatrix(0,3)=edgeNormal[0];
		bMatrix(0,4)=edgeNormal[1];
		bMatrix(0,5)=edgeNormal[0]*shapeFunctionsRotations2[0]+edgeNormal[1]*shapeFunctionsRotations2[1];
		bMatrix(1,3)=-edgeNormal[1];
		bMatrix(1,4)=+edgeNormal[0];
		bMatrix(1,5)=-edgeNormal[1]*shapeFunctionsRotations2[0]+edgeNormal[0]*shapeFunctionsRotations2[1];

		bMatrix*=lEdgeInv;

		//check bmatrix
		//Eigen::Matrix<double,6,1> nodalDofs;
		//nodalDofs << dispNode1[0], dispNode1[1], rotNode1, dispNode2[0], dispNode2[1], rotNode2;
		//Eigen::Matrix<double,2,1> deltaCrackLocal = bMatrix*nodalDofs;
		//std::cout << "delta crack using bMatrix         : "       << deltaCrackLocal.transpose() << "\n";
		//std::cout << "delta crack using shape functions : " << latticeStrain.mLatticeStrain[0] << " " << latticeStrain.mLatticeStrain[1] << "\n";

    	//add BtCB with l * projected area/length of facet add thickness
        //factor for the numerical integration
        assert(mSection->GetThickness()>0);
        double factor(mSection->GetThickness()*weight);

		areAllIpsSymmetric &= stiffness.GetSymmetry();

		// calculate element stiffness matrix
		Eigen::Matrix<double,6,6> latticeStiffness;

		latticeStiffness = bMatrix.transpose()*
				           (factor*Eigen::Map<const Eigen::Matrix<double,2,2 > >(stiffness.GetData(),2,2))*
				           bMatrix;

		AddLatticeStiffnessToElementMatrix(latticeStiffness,theIP, (theIP+1)%3, rCoefficientMatrix);
    }
    rSymmetry = areAllIpsSymmetric;

    // calculate list of global dofs related to the entries in the element stiffness matrix
    this->CalculateGlobalRowDofs(rGlobalDofsRow);
    rGlobalDofsColumn = rGlobalDofsRow;

    assert((int)rGlobalDofsRow.size()==rCoefficientMatrix.GetNumRows());
    assert((int)rGlobalDofsColumn.size()==rCoefficientMatrix.GetNumColumns());
}
    return Error::SUCCESSFUL;
}
*/
// old version left for comparison in here (to see if there is a significant speed advantage by not calling a subroutine
// calculating the strain directly without the multiplication with the bmatrix
//! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the stiffness matrix
NuTo::Error::eError NuTo::Lattice2D::CalculateCoefficientMatrix_0(NuTo::FullMatrix<double>& rCoefficientMatrix,
        std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn, bool& rSymmetry)const
{
	//allocate lattice strain and stiffness
    LatticeStrain2D latticeStrain;
    ConstitutiveTangentLocal2x2 stiffness;

    bool areAllIpsSymmetric(true);
    rCoefficientMatrix.Resize(9,9);
    Eigen::Matrix<double,2,6> bMatrix;
    Eigen::Matrix<double,6,1> nodalDofs;
    for (int theIP=0; theIP<3; theIP++)
    {
        //calculate Bmatrix
    	CalculateBMatrixAndLatticeStrain(theIP, bMatrix, latticeStrain);

    	//std::cout << "check in stiffness routine of lattice 2d" << "\n";
//latticeStrain.mLatticeStrain[0] = -7.05825e-05;
//latticeStrain.mLatticeStrain[1] = -0.000371756;


    	//calculate constitutive law
        //material pointer
        const ConstitutiveLatticeStressStrain *constitutivePtr;
        constitutivePtr = GetConstitutiveLaw(theIP)->AsConstitutiveLatticeStressStrain();

//const_cast<ConstitutiveLatticeStressStrain *>(constitutivePtr)->UpdateStaticData_LatticeStress_LatticeStrain(const_cast<Lattice2D*>(this), theIP, latticeStrain);
//latticeStrain.mLatticeStrain[0] *= -0.5;
//latticeStrain.mLatticeStrain[1] *= -1.5;

        Error::eError error = constitutivePtr->GetTangent_LatticeStress_LatticeStrain(this, theIP,
				latticeStrain, &stiffness);
		if (error!=Error::SUCCESSFUL)
			return error;

/*		{
			//central difference
			FullMatrix<double> stiffnessCDF(2,2);
			LatticeStress2D latticeStress2;
			double deltaStrain(1e-14);
			LatticeStress2D latticeStress;
			error = constitutivePtr->GetLatticeStressFromLatticeStrain(this, theIP,
					latticeStrain, latticeStress);
	    	//std::cout << "stress " <<latticeStress.mLatticeStress[0] << " " <<  latticeStress.mLatticeStress[1] << "\n";

			for (int count=0; count<2; count++)
			{
				latticeStrain.mLatticeStrain[count]+=deltaStrain;
				error = constitutivePtr->GetLatticeStressFromLatticeStrain(this, theIP,
						latticeStrain, latticeStress2);
				latticeStrain.mLatticeStrain[count]-=deltaStrain;
				stiffnessCDF(0,count) = (latticeStress2.mLatticeStress[0] - latticeStress.mLatticeStress[0])/deltaStrain;
				stiffnessCDF(1,count) = (latticeStress2.mLatticeStress[1] - latticeStress.mLatticeStress[1])/deltaStrain;
			}
			if (fabs(stiffness.GetData()[0]-stiffnessCDF(0,0))+fabs(stiffness.GetData()[2]-stiffnessCDF(0,1))+
				fabs(stiffness.GetData()[1]-stiffnessCDF(1,0))+fabs(stiffness.GetData()[3]-stiffnessCDF(1,1))>1)
			{
		    	std::cout << "strain " <<latticeStrain.mLatticeStrain[0] << " " <<  latticeStrain.mLatticeStrain[1] << "\n";
		    	std::cout << "stress " <<latticeStress.mLatticeStress[0] << " " <<  latticeStress.mLatticeStress[1] << "\n";
		    	std::cout << "max epsilon " << (this->GetStaticData(theIP))->AsConstitutiveStaticDataLatticeConcrete2D()->mEpsilonMax << "\n";
		    	std::cout << "edge length " << GetIpEdgeLength(theIP) << "\n";

		    	std::cout << "stiffness analytic " <<"\n";
				std::cout << stiffness.GetData()[0] << " " <<  stiffness.GetData()[2] << "\n";
				std::cout << stiffness.GetData()[1] << " " <<  stiffness.GetData()[3] << "\n";

				std::cout << "stiffness CDF " << "\n" << stiffnessCDF << "\n";
				std::cout <<" error in stiffness calculation " << "\n";
				//exit(-1);
			}
		}
*/
    	//add BtCB with l * projected area/length of facet add thickness
        //factor for the numerical integration
        assert(mSection->GetThickness()>0);
        double factor(mSection->GetThickness()*this->mElementData->GetIntegrationPointWeight(theIP));

		areAllIpsSymmetric &= stiffness.GetSymmetry();

		// calculate element stiffness matrix
		Eigen::Matrix<double,6,6> latticeStiffness;

		latticeStiffness = bMatrix.transpose()*
				           (factor*Eigen::Map<const Eigen::Matrix<double,2,2 > >(stiffness.GetData(),2,2))*
				           bMatrix;

		AddLatticeStiffnessToElementMatrix(latticeStiffness,theIP, (theIP+1)%3, rCoefficientMatrix);
    }
    rSymmetry = areAllIpsSymmetric;

    // calculate list of global dofs related to the entries in the element stiffness matrix
    this->CalculateGlobalRowDofs(rGlobalDofsRow);
    rGlobalDofsColumn = rGlobalDofsRow;

    assert((int)rGlobalDofsRow.size()==rCoefficientMatrix.GetNumRows());
    assert((int)rGlobalDofsColumn.size()==rCoefficientMatrix.GetNumColumns());

    return Error::SUCCESSFUL;
}

//! @brief calculates the local B-matrix for an edge connecting two nodes (with one ip)
void NuTo::Lattice2D::CalculateBMatrixAndLatticeStrain(int rIP, Eigen::Matrix<double,2,6>& rBMatrix, LatticeStrain2D& rLatticeStrain)const
{
	//get the edge
	NodeBase *node1,*node2;
	node1 = mNodes[rIP];
	node2 = mNodes[(rIP+1)%3];

	//get coordinates, displacements and radius
	boost::array<double,2> coordNode1, dispNode1, coordNode2, dispNode2;
	double rotNode1,rotNode2;
	node1->GetCoordinates2D(&(coordNode1[0]));
	node1->GetDisplacements2D(&(dispNode1[0]));
	node1->GetRotations2D(&rotNode1);
	node2->GetCoordinates2D(&(coordNode2[0]));
	node2->GetDisplacements2D(&(dispNode2[0]));
	node2->GetRotations2D(&rotNode2);

	//global coordinates of center of facet
	boost::array<double,3> globalCoordinatesFacet;
	boost::array<double,2> localCoordinatesFacet;
	this->mElementData->GetLocalIntegrationPointCoordinates2D(rIP, localCoordinatesFacet);

	InterpolateCoordinatesFrom2D(localCoordinatesFacet, globalCoordinatesFacet);

	//shape functions for the rotations of the two nodes of the corresponding edge
	boost::array<double,2> shapeFunctionsRotations1,shapeFunctionsRotations2;
	shapeFunctionsRotations1 = CalculateShapeFunctionsRotations(globalCoordinatesFacet,coordNode1);
	shapeFunctionsRotations2 = CalculateShapeFunctionsRotations(globalCoordinatesFacet,coordNode2);

	//calculate the gap at the facet center
	boost::array<double,2> deltaCrack;
	deltaCrack[0] = dispNode2[0] + shapeFunctionsRotations2[0]*rotNode2-
					dispNode1[0] - shapeFunctionsRotations1[0]*rotNode1;
	deltaCrack[1] = dispNode2[1] + shapeFunctionsRotations2[1]*rotNode2-
					dispNode1[1] - shapeFunctionsRotations1[1]*rotNode1;

	//edge vector
	double edge[2];
	edge[0] = coordNode2[0]-coordNode1[0];
	edge[1] = coordNode2[1]-coordNode1[1];

	//edge normal
	double edgeNormal[2];
	double lEdge=sqrt(edge[0]*edge[0]+edge[1]*edge[1]);
	double lEdgeInv = 1./lEdge;
	edgeNormal[0] = edge[0]*lEdgeInv;
	edgeNormal[1] = edge[1]*lEdgeInv;

	//calculate bMatrix for first node
	rBMatrix(0,0)=-edgeNormal[0];
	rBMatrix(0,1)=-edgeNormal[1];
	rBMatrix(0,2)=-edgeNormal[0]*shapeFunctionsRotations1[0]-edgeNormal[1]*shapeFunctionsRotations1[1];
	rBMatrix(1,0)=+edgeNormal[1];
	rBMatrix(1,1)=-edgeNormal[0];
	rBMatrix(1,2)=+edgeNormal[1]*shapeFunctionsRotations1[0]-edgeNormal[0]*shapeFunctionsRotations1[1];

	//calculate bMatrix for second node
	rBMatrix(0,3)=edgeNormal[0];
	rBMatrix(0,4)=edgeNormal[1];
	rBMatrix(0,5)=edgeNormal[0]*shapeFunctionsRotations2[0]+edgeNormal[1]*shapeFunctionsRotations2[1];
	rBMatrix(1,3)=-edgeNormal[1];
	rBMatrix(1,4)=+edgeNormal[0];
	rBMatrix(1,5)=-edgeNormal[1]*shapeFunctionsRotations2[0]+edgeNormal[0]*shapeFunctionsRotations2[1];

	rBMatrix*=lEdgeInv;

	//calculate lattice strain vector (epsilonN, epsilonT : normal and tangent strain
	rLatticeStrain.mLatticeStrain[0] = (deltaCrack[0]*edgeNormal[0]+deltaCrack[1]*edgeNormal[1])*lEdgeInv;
	rLatticeStrain.mLatticeStrain[1] = (deltaCrack[1]*edgeNormal[0]-deltaCrack[0]*edgeNormal[1])*lEdgeInv;

	//calculate lattice strain vector (epsilonN, epsilonT : normal and tangent strain
	//Eigen::Matrix<double,6,1> nodalDofs;
	//nodalDofs << dispNode1[0], dispNode1[1], rotNode1, dispNode2[0], dispNode2[1], rotNode2;
	//Eigen::Map<Eigen::Vector2d>(rLatticeStrain.mLatticeStrain) = rBMatrix*nodalDofs;
}


//! @brief calculate the length of an edge (belonging to an integration point
//! @param rIp integration point
//! @return edge length
double NuTo::Lattice2D::GetIpEdgeLength(int rIp)const
{
    //get coordinates
    double coord1[2],coord2[2],delta[2];
   	mNodes[rIp]->GetCoordinates2D(coord1);
   	mNodes[(rIp+1)%3]->GetCoordinates2D(coord2);

   	delta[0] = coord2[0]-coord1[0];
   	delta[1] = coord2[1]-coord1[1];

   	return sqrt(delta[0]*delta[0]+delta[1]*delta[1]);
}

//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::Lattice2D::CalculateShapeFunctionsCoordinates(const boost::array<double,2>& rNaturalCoordinates, boost::array<double,3>& rShapeFunctions)const
{
	rShapeFunctions[0] = 1-rNaturalCoordinates[0]-rNaturalCoordinates[1];
    rShapeFunctions[1] = rNaturalCoordinates[0];
    rShapeFunctions[2] = rNaturalCoordinates[1];
}

//! @brief calculates the shape functions for rotations
//! @param rGlobalCoordinatesFacet coordinates of the integration point
//! @param rCoordNode coordinate of the node
//! @return shape functions for the node
boost::array<double,2> NuTo::Lattice2D::CalculateShapeFunctionsRotations(const boost::array<double,3>& rGlobalCoordinatesFacet,
		const boost::array<double,2>& rCoordNode)const
{
	boost::array<double,2> shapeFunctions;
	shapeFunctions[0] = rCoordNode[1] - rGlobalCoordinatesFacet[1];
	shapeFunctions[1] = rGlobalCoordinatesFacet[0] - rCoordNode[0];
	return shapeFunctions;
}

//! @brief adds the stiffness of one lattice (one edge of the triangle to the global matrix
void NuTo::Lattice2D::AddLatticeStiffnessToElementMatrix(Eigen::Matrix<double,6,6>& latticeStiffness,
		int rLocalNode1, int rLocalNode2,
		NuTo::FullMatrix<double>& rCoefficientMatrix)const
{
	rCoefficientMatrix.mEigenMatrix.block<3,3>(3*rLocalNode1,3*rLocalNode1) += latticeStiffness.block<3,3>(0,0);
	rCoefficientMatrix.mEigenMatrix.block<3,3>(3*rLocalNode1,3*rLocalNode2) += latticeStiffness.block<3,3>(0,3);
	rCoefficientMatrix.mEigenMatrix.block<3,3>(3*rLocalNode2,3*rLocalNode1) += latticeStiffness.block<3,3>(3,0);
	rCoefficientMatrix.mEigenMatrix.block<3,3>(3*rLocalNode2,3*rLocalNode2) += latticeStiffness.block<3,3>(3,3);
}


//! @brief calculates the gradient of the internal potential
//! for a mechanical problem, this corresponds to the internal force vector
NuTo::Error::eError NuTo::Lattice2D::CalculateGradientInternalPotential(NuTo::FullMatrix<double>& rResult,
        std::vector<int>& rGlobalDofs)const
{
    //allocate lattice stress and strain
    LatticeStrain2D latticeStrain;
    LatticeStress2D latticeStress;

    rResult.Resize(9,1);
    Eigen::Matrix<double,2,6> bMatrix;
    for (int theIP=0; theIP<3; theIP++)
    {
        //calculate Bmatrix
    	CalculateBMatrixAndLatticeStrain(theIP, bMatrix, latticeStrain);

    	//calculate constitutive law
        //material pointer
        const ConstitutiveLatticeStressStrain *constitutivePtr;

        constitutivePtr = GetConstitutiveLaw(theIP)->AsConstitutiveLatticeStressStrain();
		Error::eError error = constitutivePtr->GetLatticeStressFromLatticeStrain(this, theIP,
				latticeStrain, latticeStress);
		if (error!=Error::SUCCESSFUL)
			return error;

    	//add BtCB with l * projected area/length of facet add thickness
        //factor for the numerical integration
        assert(mSection->GetThickness()>0);
        double factor(mSection->GetThickness()*this->mElementData->GetIntegrationPointWeight(theIP));

		// calculate element stiffness matrix
		Eigen::Matrix<double,6,1> latticeIntForce;

		latticeIntForce = bMatrix.transpose()*
				           (factor*Eigen::Map<const Eigen::Matrix<double,2,1 > >(latticeStress.GetData(),2,1));

		rResult.mEigenMatrix.block<3,1>(3*theIP,0) += latticeIntForce.block<3,1>(0,0);
		rResult.mEigenMatrix.block<3,1>(3*((theIP+1)%3),0) += latticeIntForce.block<3,1>(3,0);
    }

    // calculate list of global dofs related to the entries in the element stiffness matrix
    this->CalculateGlobalRowDofs(rGlobalDofs);

    assert((int)rGlobalDofs.size()==rResult.GetNumRows());

    return Error::SUCCESSFUL;
}

//! @brief Update the static data of an element
NuTo::Error::eError NuTo::Lattice2D::UpdateStaticData(NuTo::Element::eUpdateType rUpdateType)
{
    //allocate lattice strain
    LatticeStrain2D latticeStrain;

    Eigen::Matrix<double,2,6> bMatrix;
    for (int theIP=0; theIP<3; theIP++)
    {
        //calculate Bmatrix
    	CalculateBMatrixAndLatticeStrain(theIP, bMatrix, latticeStrain);

    	//calculate constitutive law
        //material pointer
        const ConstitutiveLatticeStressStrain *constitutivePtr;
        constitutivePtr = GetConstitutiveLaw(theIP)->AsConstitutiveLatticeStressStrain();
        Error::eError error(Error::SUCCESSFUL);
        switch(rUpdateType)
        {
        case NuTo::Element::STATICDATA:
        	error = constitutivePtr->UpdateStaticData_LatticeStress_LatticeStrain(this, theIP, latticeStrain);
        break;
        case NuTo::Element::TMPSTATICDATA:
        	error = constitutivePtr->UpdateTmpStaticData_LatticeStress_LatticeStrain(this, theIP, latticeStrain);
        break;
        case NuTo::Element::SWITCHMULTISCALE2NONLINEAR:
        	error = constitutivePtr->MultiscaleSwitchToNonlinear(this, theIP, latticeStrain);
        break;
        default:
            throw MechanicsException("[NuTo::Lattice2D::UpdateStaticData] update static data flag not known (neither static not tmpstatic data");
        }
        if (error!=Error::SUCCESSFUL)
        	return error;
    }

    return Error::SUCCESSFUL;
}

//! @brief calculates the coefficient matrix for the 1-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the damping matrix
NuTo::Error::eError NuTo::Lattice2D::CalculateCoefficientMatrix_1(NuTo::FullMatrix<double>& rResult,
        std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn, bool& rSymmetry)const
{
    throw MechanicsException("[NuTo::Lattice2D::CalculateCoefficientMatrix_1] to be implemented.");
}

//! @brief calculates the coefficient matrix for the 2-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the Mass matrix
NuTo::Error::eError NuTo::Lattice2D::CalculateCoefficientMatrix_2(NuTo::FullMatrix<double>& rResult,
        std::vector<int>& rGlobalDofsRow, std::vector<int>& rGlobalDofsColumn, bool& rSymmetry)const
{
	rResult.Resize(9,9);

    //calculate global coordinates of the edge points and face points
    boost::array<boost::array<double, 2>,3 > globalCoordinatesEdgePoints;
    boost::array<double, 2> globalCoordinatesFacePoint;
    CalculateGlobalCoordinatesEdgeFacePoints(globalCoordinatesEdgePoints, globalCoordinatesFacePoint);

    for (int theIP=0; theIP<3; theIP++)
    {
    	//calculate constitutive law
        //material pointer
        const ConstitutiveLatticeStressStrain *constitutivePtr;
        constitutivePtr = GetConstitutiveLaw(theIP)->AsConstitutiveLatticeStressStrain();
    	double density(constitutivePtr->GetDensity());
        for (int countNode=0; countNode<2; countNode++)
        {
			//calculate area/volume of the facet with the node
			boost::array<double,3 > coordinatesNodes;
			int theNode((theIP+countNode)%3);
			mNodes[theNode]->GetCoordinates2D(&(coordinatesNodes[0]));

			//calculate area of the triangle (node, edgePoint, Facepoint
			double mass = density * 0.5*(coordinatesNodes[0]*globalCoordinatesEdgePoints[theIP][1]-coordinatesNodes[1]*globalCoordinatesEdgePoints[theIP][0]+
					           globalCoordinatesEdgePoints[theIP][0]*globalCoordinatesFacePoint[1]-globalCoordinatesEdgePoints[theIP][1]*globalCoordinatesFacePoint[0]+
					           globalCoordinatesFacePoint[0]*coordinatesNodes[1]-globalCoordinatesFacePoint[1]*coordinatesNodes[0]);
			//this is the translational mass
			rResult(3*theNode,3*theNode)+=mass;
			rResult(3*theNode+1,3*theNode+1)+=mass;

			//centroid w.r.t. node
			boost::array<double,2 > centroid;
			boost::array<double,2 > delta2;
			delta2[0] = globalCoordinatesEdgePoints[theIP][0] - coordinatesNodes[0];
			delta2[0] = globalCoordinatesEdgePoints[theIP][1] - coordinatesNodes[1];
			boost::array<double,2 > delta3;
			delta3[0] = globalCoordinatesFacePoint[0] - coordinatesNodes[0];
			delta3[0] = globalCoordinatesFacePoint[1] - coordinatesNodes[1];
			centroid[0] = (delta2[0]+delta3[0])/3.;
			centroid[1] = (delta2[1]+delta3[1])/3.;

			boost::array<double,2 > IiMulRho;
			IiMulRho[0]=mass*centroid[0];
			IiMulRho[1]=mass*centroid[1];

			double IiiMulRho = mass/6.*(delta3[0]*delta3[0]+delta3[1]*delta3[1]+delta2[0]*delta2[0]+delta2[1]*delta2[1]+delta2[0]*delta3[0]+delta2[1]*delta3[1]);

			//this is the rotational mass
			rResult(3*theNode,3*theNode+2)-=IiMulRho[1];
			rResult(3*theNode+2,3*theNode)-=IiMulRho[1];
			rResult(3*theNode+1,3*theNode+2)+=IiMulRho[0];
			rResult(3*theNode+2,3*theNode+1)+=IiMulRho[0];
			rResult(3*theNode+2,3*theNode+2)+=IiiMulRho;
        }
    }
    rSymmetry = true;
    return Error::SUCCESSFUL;
}

// these headers are just for test purpose
//#include <eigen2/Eigen/LU>
//#include <eigen2/Eigen/Array>

//! @brief calculates the integration point data with the current displacements applied
//! @param rIpDataType data type to be stored for each integration point
//! @param rIpData return value with dimension (dim of data type) x (numIp)
NuTo::Error::eError NuTo::Lattice2D::GetIpData(NuTo::IpData::eIpStaticDataType rIpDataType, FullMatrix<double>& rIpData)const
{
    //allocate lattice strain
    LatticeStrain2D latticeStrain;
    LatticeStrain3D plasticLatticeStrain3D;
    LatticeStress3D latticeStress3D;

    Eigen::Matrix<double,2,6> bMatrix;

    //allocate and initialize result matrix
    switch (rIpDataType)
    {
    case NuTo::IpData::LATTICE_STRAIN:
    case NuTo::IpData::LATTICE_STRESS:
    case NuTo::IpData::LATTICE_PLASTIC_STRAIN:
           rIpData.Resize(3,GetNumIntegrationPoints());
    break;
    case NuTo::IpData::DAMAGE:
           rIpData.Resize(1,GetNumIntegrationPoints());
    break;
    case NuTo::IpData::ELASTIC_ENERGY:
    case NuTo::IpData::TOTAL_ENERGY:
        rIpData.Resize(2,GetNumIntegrationPoints());
    break;
    default:
        throw MechanicsException("[NuTo::Lattice2D::GetIpData] Ip data not implemented.");
    }

    for (int theIP=0; theIP<3; theIP++)
    {
        //calculate Bmatrix
    	CalculateBMatrixAndLatticeStrain(theIP, bMatrix, latticeStrain);

    	//calculate constitutive law
        //material pointer
        const ConstitutiveLatticeStressStrain *constitutivePtr;
        constitutivePtr = GetConstitutiveLaw(theIP)->AsConstitutiveLatticeStressStrain();
        Error::eError error(Error::SUCCESSFUL);
        switch (rIpDataType)
        {
        case NuTo::IpData::LATTICE_STRAIN:
        	//error = constitutivePtr->GetLatticeStrain(this, theIP, latticeStrain, latticeStrain3D);
            memcpy(&(rIpData.mEigenMatrix.data()[theIP*3]),latticeStrain.GetData(),2*sizeof(double));
        break;
        case NuTo::IpData::LATTICE_STRESS:
        	error = constitutivePtr->GetLatticeStressFromLatticeStrain(this, theIP, latticeStrain, latticeStress3D);
            memcpy(&(rIpData.mEigenMatrix.data()[theIP*3]),latticeStress3D.GetData(),3*sizeof(double));
        break;
        case NuTo::IpData::LATTICE_PLASTIC_STRAIN:
        	error = constitutivePtr->GetLatticePlasticStrain(this, theIP, latticeStrain, plasticLatticeStrain3D);
            memcpy(&(rIpData.mEigenMatrix.data()[theIP*3]),plasticLatticeStrain3D.GetData(),3*sizeof(double));
        break;
        case NuTo::IpData::DAMAGE:
        	error = constitutivePtr->GetDamage(this, theIP, latticeStrain, rIpData.mEigenMatrix.data()[theIP]);
        break;
        case NuTo::IpData::ELASTIC_ENERGY:
        {
        	error = constitutivePtr->GetElasticEnergy_LatticeStress_LatticeStrain(this, theIP, latticeStrain, rIpData.mEigenMatrix(0,theIP));
            assert(mSection->GetThickness()>0);
            double factor(mSection->GetThickness()*this->mElementData->GetIntegrationPointWeight(theIP));
            rIpData.mEigenMatrix(1,theIP) = factor;
        }
        break;
        case NuTo::IpData::TOTAL_ENERGY:
        {
        	error = constitutivePtr->GetTotalEnergy_LatticeStress_LatticeStrain(this, theIP, latticeStrain,rIpData.mEigenMatrix(0,theIP));
            assert(mSection->GetThickness()>0);
            double factor(mSection->GetThickness()*this->mElementData->GetIntegrationPointWeight(theIP));
            rIpData.mEigenMatrix(1,theIP) = factor;
        }
        break;
        default:
            throw MechanicsException("[NuTo::Lattice2D::GetIpData] Ip data not implemented.");
        }
        if (error!=Error::SUCCESSFUL)
        	return error;
    }

    return Error::SUCCESSFUL;
/*    //calculate local coordinates
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
    case NuTo::IpData::TOTAL_ENERGY:
        rIpData.Resize(2,GetNumIntegrationPoints());
    break;
    default:
        throw MechanicsException("[NuTo::Lattice2D::GetIpData] Ip data not implemented.");
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
        case NuTo::IpData::TOTAL_ENERGY:
        {
        	error = constitutivePtr->GetTotalEnergy_EngineeringStress_EngineeringStrain(this, theIP, deformationGradient,rIpData.mEigenMatrix(0,theIP));
            assert(mSection->GetThickness()>0);
            double factor(mSection->GetThickness()*detJac*(mElementData->GetIntegrationType()->GetIntegrationPointWeight(theIP)));
            rIpData.mEigenMatrix(1,theIP) = factor;
        }
        break;
        default:
            throw MechanicsException("[NuTo::Lattice2D::GetIpData] Ip data not implemented.");
        }
        if (error!=Error::SUCCESSFUL)
        	return error;
    }
*/
    return Error::SUCCESSFUL;
}


//! @brief Allocates static data for an integration point of an element
//! @param rConstitutiveLaw constitutive law, which is called to allocate the static data object
//! actually, both - the element type and the constitutive law are required to determine the static data object actually required
NuTo::ConstitutiveStaticDataBase* NuTo::Lattice2D::AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw)const
{
    const ConstitutiveLatticeStressStrain *constitutivePtr = rConstitutiveLaw->AsConstitutiveLatticeStressStrain();
    if (constitutivePtr==0)
        throw MechanicsException("[NuTo::Lattice2D::AllocateStaticData] Constitutive law can not deal with lattice stresses and strains");
    return constitutivePtr->AllocateStaticDataLatticeStress_LatticeStrain2D(this);

}

// interpolate geometry
void NuTo::Lattice2D::InterpolateCoordinatesFrom2D(const boost::array<double, 2> rNaturalCoordinates, boost::array<double, 3>& rGlobalCoordinates) const
{
	boost::array<double,3> shapeFunctions;
	this->CalculateShapeFunctionsCoordinates(rNaturalCoordinates, shapeFunctions);

    // start interpolation
    rGlobalCoordinates[0] = 0.0;
    rGlobalCoordinates[1] = 0.0;
    rGlobalCoordinates[2] = 0.0;
    for (int NodeCount = 0; NodeCount < this->GetNumNodes(); NodeCount++)
    {
        // get node coordinate
        double NodeCoordinate[3];
        const NodeBase *nodePtr(GetNode(NodeCount));
        if (nodePtr->GetNumCoordinates()==2)
            nodePtr->GetCoordinates2D(NodeCoordinate);
        else
            nodePtr->GetCoordinates3D(NodeCoordinate);

        // add node contribution
        for (int theCoordinate=0; theCoordinate<nodePtr->GetNumCoordinates(); theCoordinate++)
        {
            rGlobalCoordinates[theCoordinate] += shapeFunctions[NodeCount] *  NodeCoordinate[theCoordinate];
        }
    }
}


// interpolate geometry
void NuTo::Lattice2D::InterpolateCoordinatesFrom2D(double rNaturalCoordinates[2], double rGlobalCoordinates[3]) const
{
    // calculate shape functions
	boost::array<double, 3> globalCoordinates;
	boost::array<double, 2> naturalCoordinates;
	// this is silly, but we should work on a conversion from doubl[] to boost arrays, until this is done for InterpolateCoordinatesFrom2D as well, do this additional step
	naturalCoordinates[0] = rNaturalCoordinates[0];
	naturalCoordinates[1] = rNaturalCoordinates[1];
	InterpolateCoordinatesFrom2D(naturalCoordinates, globalCoordinates);
	rGlobalCoordinates[0] = globalCoordinates[0];
	rGlobalCoordinates[1] = globalCoordinates[1];
	rGlobalCoordinates[2] = globalCoordinates[2];
}

// interpolate displacements (only for visualize, since the rotations are not taken into account)
void NuTo::Lattice2D::InterpolateDisplacementsFrom2D(const boost::array<double, 2> rNaturalCoordinates, boost::array<double, 3>& rGlobalDisplacements) const
{
	boost::array<double,3> shapeFunctions;
	this->CalculateShapeFunctionsCoordinates(rNaturalCoordinates, shapeFunctions);

    // start interpolation
	rGlobalDisplacements[0] = 0.0;
	rGlobalDisplacements[1] = 0.0;
	rGlobalDisplacements[2] = 0.0;
    for (int NodeCount = 0; NodeCount < this->GetNumNodes(); NodeCount++)
    {
        // get node coordinate
        double NodeDisplacements[3];
        const NodeBase *nodePtr(GetNode(NodeCount));
        if (nodePtr->GetNumDisplacements()==2)
            nodePtr->GetDisplacements2D(NodeDisplacements);
        else
            nodePtr->GetDisplacements3D(NodeDisplacements);

        // add node contribution
        for (int theDisplacement=0; theDisplacement<nodePtr->GetNumDisplacements(); theDisplacement++)
        {
        	rGlobalDisplacements[theDisplacement] += shapeFunctions[NodeCount] *  NodeDisplacements[theDisplacement];
        }
    }
}

// interpolate geometry
void NuTo::Lattice2D::InterpolateDisplacementsFrom2D(double rNaturalCoordinates[2], double rDisplacements[3]) const
{
    // calculate shape functions
	boost::array<double, 3> displacements;
	boost::array<double, 2> naturalCoordinates;
	// this is silly, but we should work on a conversion from doubl[] to boost arrays, until this is done for InterpolateCoordinatesFrom2D as well, do this additional step
	naturalCoordinates[0] = rNaturalCoordinates[0];
	naturalCoordinates[1] = rNaturalCoordinates[1];
	InterpolateDisplacementsFrom2D(naturalCoordinates, displacements);
	rDisplacements[0] = displacements[0];
	rDisplacements[1] = displacements[1];
	rDisplacements[2] = displacements[2];
}

// check element definition
void NuTo::Lattice2D::CheckElement()
{
/*
    // check nodes
    for (int nodeCount = 0; nodeCount < this->GetNumNodes(); nodeCount++)
    {
        int numCoordinates(GetNode(nodeCount)->GetNumCoordinates());
        if (numCoordinates<2 || numCoordinates>3)
        {
            throw MechanicsException("[NuTo::Lattice2D::CheckElement] invalid node type (check node definition for coordinates).");
        }
    }

    // check node ordering (element length must be positive) and for changing sign in jacobian determinant
    // calculate coordinates
    std::vector<double> nodeCoord(2*this->GetNumNodes());
    this->CalculateLocalCoordinates(nodeCoord);

    // check number of integration points
    if (this->GetNumIntegrationPoints() < 1)
    {
        throw MechanicsException("[NuTo::Lattice2D::CheckElement] invalid integration type.");
    }

    // check sign of the jacobian determinant of the first integration point
    double naturalIPCoord[2];
    this->GetLocalIntegrationPointCoordinates(0, naturalIPCoord);

    std::vector<double> derivativeShapeFunctionsNatural(2*GetNumNodes());
    this->CalculateDerivativeShapeFunctionsNatural(naturalIPCoord, derivativeShapeFunctionsNatural);

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
        this->CalculateDerivativeShapeFunctionsNatural(naturalIPCoord, derivativeShapeFunctionsNatural);
        this->CalculateJacobian(derivativeShapeFunctionsNatural,nodeCoord, invJacobian, detJacobian);
        //std::cout << "Jacobian " << detJacobian << std::endl;
        if (detJacobian <= 0)
        {
            std::cout << "jac " << detJacobian << std::endl;
        	throw MechanicsException("[NuTo::Lattice2D::CheckElement] element is not properly defined by this nodes (zero or negative jacobian determinant).");
        }
        volume += this->GetIntegrationPointWeight(ipCount) * detJacobian;
    }

    // check element volume
    if (volume < 1e-14)
    {
        throw MechanicsException("[NuTo::Lattice2D::CheckElement] element with zero volume (check nodes).");
    }
*/
}

//! @brief sets the section of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need a section
//! @param rSection pointer to section
//! @return pointer to constitutive law
void NuTo::Lattice2D::SetSection(const SectionBase* rSection)
{
    mSection = rSection;
}

//! @brief returns a pointer to the section of an element
//! implemented with an exception for all elements, reimplementation required for those elements
//! which actually need a section
//! @return pointer to section
const NuTo::SectionBase* NuTo::Lattice2D::GetSection()const
{
    return mSection;
}

//! @brief returns the number of nodes in this element
//! @return number of nodes
int NuTo::Lattice2D::GetNumNodes()const
{
    return 3;
}

//! @brief returns a pointer to the i-th node of the element
//! @param local node number
//! @return pointer to the node
NuTo::NodeBase* NuTo::Lattice2D::GetNode(int rLocalNodeNumber)
{
	assert(0<=rLocalNodeNumber && rLocalNodeNumber<3);
    return mNodes[rLocalNodeNumber];
}

//! @brief returns a pointer to the i-th node of the element
//! @param local node number
//! @return pointer to the node
const NuTo::NodeBase* NuTo::Lattice2D::GetNode(int rLocalNodeNumber)const
{
	assert(0<=rLocalNodeNumber && rLocalNodeNumber<3);
    return mNodes[rLocalNodeNumber];
}

//! @brief sets the rLocalNodeNumber-th node of the element
//! @param local node number
//! @param pointer to the node
void NuTo::Lattice2D::SetNode(int rLocalNodeNumber, NuTo::NodeBase* rNode)
{
	assert(0<=rLocalNodeNumber && rLocalNodeNumber<3);
    mNodes[rLocalNodeNumber] = rNode;
}

//! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
//! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
void NuTo::Lattice2D::ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr)
{
    for (int count=0; count<3; count++)
    {
        if (this->mNodes[count]==rOldPtr)
        {
            this->mNodes[count]=rNewPtr;
            break;
        }
    }
}

//! @brief calculates the volume of an integration point (weight * detJac)
//! @param rVolume  vector for storage of the ip volumes (area in 2D)
void NuTo::Lattice2D::GetIntegrationPointVolume(std::vector<double>& rVolume)const
{
	throw MechanicsException("[NuTo::Lattice2D::GetIntegrationPointVolume] not implemented.");
}

//! @brief returns the coordinates of an integration point
//! @param rIpNum integration point
//! @param rCoordinates coordinates to be returned
void  NuTo::Lattice2D::GetGlobalIntegrationPointCoordinates(int rIpNum, double rCoordinates[3])const
{
	boost::array<double,2> localCoordinatesFacet;
	boost::array<double,3> coordinates;
	this->mElementData->GetLocalIntegrationPointCoordinates2D(rIpNum, localCoordinatesFacet);
	InterpolateCoordinatesFrom2D(localCoordinatesFacet, coordinates);
	rCoordinates[0] = coordinates[0];
	rCoordinates[1] = coordinates[1];
	rCoordinates[2] = coordinates[2];
}

// reorder nodes such that the sign of the length/area/volume of the element changes
void NuTo::Lattice2D::ReorderNodes()
{
	//std::cout << "reorder element nodes" << std::endl;
    NodeBase* tmp = this->mNodes[1];
    this->mNodes[1] = this->mNodes[2];
    this->mNodes[2] = tmp;
}

//! @brief ... extract global dofs from nodes (mapping of local row ordering of the element matrices to the global dof ordering)
//! @param rGlobalRowDofs ... vector of global row dofs
void NuTo::Lattice2D::CalculateGlobalRowDofs(std::vector<int>& rGlobalRowDofs) const
{
	rGlobalRowDofs.resize(3 * this->GetNumNodes());
    for (int nodeCount = 0; nodeCount < this->GetNumNodes(); nodeCount++)
    {
        const NodeBase *nodePtr = this->GetNode(nodeCount);
        rGlobalRowDofs[3 * nodeCount    ] = nodePtr->GetDofDisplacement(0);
        rGlobalRowDofs[3 * nodeCount + 1] = nodePtr->GetDofDisplacement(1);
        rGlobalRowDofs[3 * nodeCount + 2] = nodePtr->GetDofRotation(0);
    }

}

//! @brief ... extract global dofs from nodes (mapping of local column ordering of the element matrices to the global dof ordering)
//! @param rGlobalColumnDofs ... vector of global column dofs
void NuTo::Lattice2D::CalculateGlobalColumnDofs(std::vector<int>& rGlobalColumnDofs) const
{
	rGlobalColumnDofs.resize(3 * this->GetNumNodes());
    for (int nodeCount = 0; nodeCount < this->GetNumNodes(); nodeCount++)
    {
        const NodeBase *nodePtr = this->GetNode(nodeCount);
        rGlobalColumnDofs[3 * nodeCount    ] = nodePtr->GetDofDisplacement(0);
        rGlobalColumnDofs[3 * nodeCount + 1] = nodePtr->GetDofDisplacement(1);
        rGlobalColumnDofs[3 * nodeCount + 2] = nodePtr->GetDofRotation(0);
    }
}


//! @brief cast the base pointer to an ElementLattice2D, otherwise throws an exception
const NuTo::Lattice2D* NuTo::Lattice2D::AsLattice2D()const
{
    return this;
}

//! @brief cast the base pointer to an Lattice2D, otherwise throws an exception
NuTo::Lattice2D* NuTo::Lattice2D::AsLattice2D()
{
    return this;
}


#ifdef ENABLE_VISUALIZE
void NuTo::Lattice2D::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::CellBase::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    //calculate global coordinates of the edge points and face points
    boost::array<boost::array<double, 2>,3 > globalCoordinatesEdgePoints;
    boost::array<double, 2> globalCoordinatesFacePoint;
    CalculateGlobalCoordinatesEdgeFacePoints(globalCoordinatesEdgePoints, globalCoordinatesFacePoint);

    //calculate help values for the transformation global local
    Eigen::Matrix<double,2,2> Ainv;
    boost::array<double, 2> coordinatesNode0;
    CalculateHelpTransformationGlobalLocal(Ainv, coordinatesNode0);


    //transform global coordinates to local coordinates
    boost::array<boost::array<double, 2>,3 > localCoordinatesEdgePoints;
    boost::array<double, 2> localCoordinatesFacePoint;

    //transform the global coordinates to local coordinates
    for (int count=0; count<3; count++)
    {
    	CalculateLocalCoordinateFromGlobal(Ainv, coordinatesNode0 ,
    			globalCoordinatesEdgePoints[count],localCoordinatesEdgePoints[count]);
    	//std::cout << "local coordinates edge point " << localCoordinatesEdgePoints[count][0] << " " << localCoordinatesEdgePoints[count][1] << "\n";
    }
    CalculateLocalCoordinateFromGlobal(Ainv, coordinatesNode0 ,
    		globalCoordinatesFacePoint,localCoordinatesFacePoint);

    NumVisualizationPoints = 4;
    VisualizationPointLocalCoordinates.push_back(localCoordinatesEdgePoints[0][0]);
    VisualizationPointLocalCoordinates.push_back(localCoordinatesEdgePoints[0][1]);
    VisualizationPointLocalCoordinates.push_back(localCoordinatesEdgePoints[1][0]);
    VisualizationPointLocalCoordinates.push_back(localCoordinatesEdgePoints[1][1]);
    VisualizationPointLocalCoordinates.push_back(localCoordinatesEdgePoints[2][0]);
    VisualizationPointLocalCoordinates.push_back(localCoordinatesEdgePoints[2][1]);
    VisualizationPointLocalCoordinates.push_back(localCoordinatesFacePoint[0]);
    VisualizationPointLocalCoordinates.push_back(localCoordinatesFacePoint[1]);

    NumVisualizationCells = 3;
    VisualizationCellType.push_back(NuTo::CellBase::LINE);
    VisualizationCellType.push_back(NuTo::CellBase::LINE);
    VisualizationCellType.push_back(NuTo::CellBase::LINE);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIP.push_back(0);
    VisualizationCellsIP.push_back(1);
    VisualizationCellsIP.push_back(2);
}
#endif// ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::Lattice2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Lattice2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Lattice2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Lattice2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Lattice2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Lattice2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Lattice2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Lattice2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementBase)
          & BOOST_SERIALIZATION_NVP(mNodes)
          & BOOST_SERIALIZATION_NVP(mSection);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Lattice2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Lattice2D)
#endif // ENABLE_SERIALIZATION

