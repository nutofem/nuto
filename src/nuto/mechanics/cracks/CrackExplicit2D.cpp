// $Id$

#include <iostream>
#include <cmath>
#include <limits>

#include <boost/foreach.hpp>

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/cracks/CrackExplicit2D.h"
#include "nuto/mechanics/nodes/NodeCoordinates.h"

//! @brief constructor
NuTo::CrackExplicit2D::CrackExplicit2D() : CrackBase(), std::list<NuTo::NodeBase*>()
{}

NuTo::CrackExplicit2D::CrackExplicit2D(const std::string rName): CrackBase(rName), std::list<NuTo::NodeBase*>()
{}

NuTo::CrackExplicit2D::~CrackExplicit2D()
{}


//! @brief info for the crack
//! @param rVerboseLevel verbose Level
void NuTo::CrackExplicit2D::Info(int rVerboseLevel)const
{
	std::cout << "    This crack " << mName << std::endl;
	std::cout << "Number of crack points" << this->size() << std::endl;
	/// output of the cracknodes
    if (rVerboseLevel>3)
    {
		int i=0;
        for (std::list<NuTo::NodeBase*>::const_iterator it = this->begin(); it!= this->end(); it++, i++)
        {
			if (rVerboseLevel>4)
			{
					std::cout << "\tcoordinates of node " << i << " :";
					for(unsigned short iDof=0; iDof < (*it)->GetNumCoordinates(); ++iDof)
							std::cout << "\t" << (*it)->GetCoordinate(iDof);
			}
			std::cout << std::endl;
        }

    }
}

//! @brief initiate the crack into the structure
//! @param rElements ... vector of element pointer
//! @param rCrackedElements ... vector of cracked element pointers
void NuTo::CrackExplicit2D::Initiate(std::vector<NuTo::ElementBase*>& rElements, std::vector<NuTo::ElementBase*>& rCrackedElements)
{
	/// go through all given elements
	typedef std::vector<NuTo::ElementBase*> elemVec_t;
	BOOST_FOREACH(NuTo::ElementBase* thisElement,rElements)
	{
		bool isCracked=false;
		bool lastSgn=0;

		const int numElemNodes=thisElement->GetNumNodes ();
		for(int i=0; i<numElemNodes; i++)
		{

			const NuTo::NodeBase* thisNode=thisElement->GetNode(i);
			const double thisDist=NodeDistance(thisNode);
			if(i==0)
			{
				lastSgn=std::signbit(thisDist);
				continue;
			}

			/// if not all distances are on the same side -> element is cracked by this crack
			/// std::signbit(thisDist) -> 1 if negative 0 otherwise
			if(lastSgn != std::signbit(thisDist))// not all nodes on the same side --> element cracked
			{
				/// check if at least 1 element edge is intersected by the crack
				for(int j=0; j<numElemNodes; j++)
				{
					NuTo::NodeBase* nodeI=new NuTo::NodeCoordinates<2>();
					double relCoor=0.0;
					size_t seg;
					isCracked=Intersect(thisElement->GetNode(j),thisElement->GetNode(((j+1)%numElemNodes)),nodeI, relCoor, seg);
					delete nodeI;

					if(isCracked)
					{
						//! @todo delete thisElement and create a new one with different rElementDataType, rIpDataType
						rCrackedElements.push_back(thisElement);
						break;
					}
				}
				//! if cracked element found: goto next element
				if(isCracked) break;
			}
		}
	}

/*
	//! check if the cracknodes are inside the cracked elements
	for(std::list<NuTo::NodeBase*>::const_iterator nodeIt=this->begin(); nodeIt!=this->end(); ++nodeIt )
	{
		double 	thisCoords[2];
		unsigned int	inside=0;
		(*nodeIt)->GetCoordinates2D(thisCoords);
		DBG_PRINT_VEC(thisCoords)
		BOOST_FOREACH(NuTo::ElementBase* thisElement,rCrackedElements)
		{
			inside+=thisElement->CheckPointInside(thisCoords);
		}

		DBG_PRINT_VAL(inside)
		if(inside==0)
			throw NuTo::MechanicsException("[NuTo::CrackExplicit2D::Intersect] Node outside");

	}
*/

}


//! @brief computes the signed distance of a node to the crack
//! @param rNode ... node to be evaluated
//! @todo  move this function to geometry class
double NuTo::CrackExplicit2D::NodeDistance(const NuTo::NodeBase* rNode) const
{
	double nodeCoor[2];
	double signedDist=HUGE_VAL;//LONG_MAX ;
	rNode->GetCoordinates2D(nodeCoor);

	/// go through all crack segments to check the minimal Distance of the node

	std::list<NuTo::NodeBase*>::const_iterator crackPtIt = this->begin();
	double xA[2],xB[2];

	while (crackPtIt!= this->end())
    {
		/// get coordinates as long as the second segment point is the end of the list
    	(*crackPtIt)->GetCoordinates2D(xA);
		if( ++crackPtIt == this->end()) 	break;

		(*crackPtIt)->GetCoordinates2D(xB);

		/// calculate the normal to the cracksegment
 		/*!
		    In the two-dimensional case, the normal $\boldsymbol{n}$ of the crack is defined as
		    \f[\boldsymbol{n} =
				 \dfrac{ \left( \boldsymbol{x}_B - \boldsymbol{x}_A \right) \times \boldsymbol{e}_z}
				       { \left| \left( \boldsymbol{x}_B - \boldsymbol{x}_A \right) \times \boldsymbol{e}_z \right|}
		    \f]
		*/
		double normalVec[2] = { xB[1] - xA[1] , xA[0] - xB[0]};
		double norm = sqrt(normalVec[0]*normalVec[0]+normalVec[1]*normalVec[1]);
		normalVec[0] /= norm;
		normalVec[1] /= norm;

		/// calculate the distance of the node to the cracksegment
		/*!
		    In the two-dimensional case, with the normalized normal $\boldsymbol{n}$ of the crack the distance $d$ is defined as
		    \f[ d = \left( \boldsymbol{x}_B - \boldsymbol{x}_A \right) \cdot \boldsymbol{n}
		    \f]
		*/
		double dist= (nodeCoor[0] - xA[0])*normalVec[0] +(nodeCoor[1] - xA[1])*normalVec[1];

		if(std::abs(dist) < std::abs(signedDist) )
			signedDist = dist;
    }
	return signedDist;
}

//! @brief computes the intersection point of a crack and a line
//! @return a bool if the line segment AB is cracked or not
//! @param rNodeA (Input) ... const reference to first node of the line
//! @param rNodeB (Input) ... const reference to second node of the line
//! @param rNodeI (Output) ... reference to intersection node of the line with the crack
//! @param rDist  (Output) ... relative coordinate of the intersection
//! @param rSeg   (Output) ... id of the crack segment which intersects the line segment AB (max(size_t) if AB is not cracked) (remark: sero-based)
//! @todo  move this function to geometry class
bool NuTo::CrackExplicit2D::Intersect(const NuTo::NodeBase* rNodeA, const NuTo::NodeBase* rNodeB, NuTo::NodeBase* rNodeI, double &rDist, size_t& rSeg)
{
	//! check if input nodes are 2D-nodes
	if(rNodeA->GetNumCoordinates()!=2 || rNodeB->GetNumCoordinates()!=2 )
		throw NuTo::MechanicsException("[NuTo::CrackExplicit2D::Intersect] don't got an 2D-Node!!!");

	bool isCracked=false;
	double dist2=0.0;
	rSeg=-1;

	std::list<NuTo::NodeBase*>::const_iterator crackPtIt = this->begin();
	size_t iSeg=0;
	while (crackPtIt!= this->end())
	{
		/// get coordinates as long as the second segment point is the end of the list
		const NuTo::NodeBase* nodeC(*crackPtIt);
		if( ++crackPtIt == this->end()) 	break;

		const NuTo::NodeBase* nodeD(*crackPtIt);

		if( IntersectSegmentSegment(rNodeA,rNodeB,nodeC,nodeD,rNodeI,rDist,dist2) ){
			isCracked=true;
			break;
		}
		++iSeg;

	}
	if(isCracked) rSeg=iSeg;
	return isCracked;
}

//! @brief computes the intersection point between the end-ray of this crack and a line segment
//! @return an unsigned short if the line segment AB is intersected or not (0=no intersection, 1=front, 2=end)
//! @param rNodeA (Input) ... const reference to first node of the line
//! @param rNodeB (Input) ... const reference to second node of the line
//! @param rNodeI (Output) ... reference to intersection node of the line with the crack
//! @param rDist  (Output) ... relative coordinate of the segment intersection
unsigned short NuTo::CrackExplicit2D::ExtendEnd(const NuTo::NodeBase* rNodeA, const NuTo::NodeBase* rNodeB, NuTo::NodeBase* rNodeI, double& rDist)
{
	//! check if input nodes are 2D-nodes
	if(rNodeA->GetNumCoordinates()!=2 || rNodeB->GetNumCoordinates()!=2 )
		throw NuTo::MechanicsException("[NuTo::CrackExplicit2D::ExtendEnd] don't got an 2D-Node!!!");

	unsigned short isCracked=0;

	//! Variables for the ray subroutine
	bool rayIntersects=false;
	double distRay=0.0;

	std::list<NuTo::NodeBase*>::const_reverse_iterator crackPtIt1 = this->rbegin();	++crackPtIt1;
	std::list<NuTo::NodeBase*>::const_reverse_iterator crackPtIt2 = this->rbegin();

	//! compute the nomalized direction of the end-ray
	double dir[2];
	double coor1[2], coor2[2];
	(*crackPtIt1)->GetCoordinates2D(coor1);
	(*crackPtIt2)->GetCoordinates2D(coor2);
	dir[0]=coor2[0]-coor1[0];
	dir[1]=coor2[1]-coor1[1];
	const double len=sqrt(dir[0]*dir[0]+dir[1]*dir[1]);
	dir[0]/=len;
	dir[1]/=len;

	//! check if there is an intersection with the end-ray of this crack
	const NuTo::NodeBase* nodeC(*crackPtIt2);
	rayIntersects=IntersectSegmentRay(rNodeA,rNodeB,nodeC,dir,rNodeI,rDist, distRay);
	if(rayIntersects && 0<distRay)
	{
		isCracked=2; //!< projection of the end-tip (end of the vector)
	}else{
		//! now do the same stuff with the other crack tip
		std::list<NuTo::NodeBase*>::const_iterator crackPtIt1 = this->begin();	++crackPtIt1;
		std::list<NuTo::NodeBase*>::const_iterator crackPtIt2 = this->begin();

		//! compute the nomalized direction of the end-ray
		double dir[2];
		double coor1[2], coor2[2];
		(*crackPtIt1)->GetCoordinates2D(coor1);
		(*crackPtIt2)->GetCoordinates2D(coor2);
		dir[0]=coor2[0]-coor1[0];
		dir[1]=coor2[1]-coor1[1];
		const double len=sqrt(dir[0]*dir[0]+dir[1]*dir[1]);
		dir[0]/=len;
		dir[1]/=len;

		const NuTo::NodeBase* nodeC(*crackPtIt2);
		rayIntersects=IntersectSegmentRay(rNodeA,rNodeB,nodeC,dir,rNodeI,rDist, distRay);
		if(rayIntersects && 0<distRay)
		{
			isCracked=1; //!< projection of the front-tip (begin of the vector)
		}
	}

	return isCracked;
}

//! @brief computes the intersection point of two straight segments
//! @return bool intersect?
//! @param rNodeA (Input) ... const reference to first node of the first segment
//! @param rNodeB (Input) ... const reference to second node of the first segment
//! @param rNodeC (Input) ... const reference to first node of the second segment
//! @param rNodeD (Input) ... const reference to second node of the second segment
//! @param rNodeI (Output) ... reference to intersection node of the segments
//! @param rDist1 (Output) ... normalized position of the intersection on the first segment ( xI=xA+r*(xB-xA) )
//! @param rDist2 (Output) ... normalized position of the intersection on the sedond segment ( xI=xC+s*(xD-xC) )
//! @todo  move this function to geometry class
bool NuTo::CrackExplicit2D::IntersectSegmentSegment(	const NuTo::NodeBase* rNodeA, const NuTo::NodeBase* rNodeB,
														const NuTo::NodeBase* rNodeC, const NuTo::NodeBase* rNodeD,
														NuTo::NodeBase* rNodeI, double & rDist1, double & rDist2)
{
	//! check if input nodes are 2D-nodes
	if(rNodeA->GetNumCoordinates()!=2 || rNodeB->GetNumCoordinates()!=2 || rNodeC->GetNumCoordinates()!=2 || rNodeD->GetNumCoordinates()!=2 || rNodeI->GetNumCoordinates()!=2 )
		throw NuTo::MechanicsException("[NuTo::CrackExplicit2D::Intersect] don't got an 2D-Node!!!");

	double xA[2],xB[2], xC[2],xD[2], xI[2], dir1[] = {0.0 , 0.0}, dir2[] = {0.0 , 0.0};
	rNodeA->GetCoordinates2D(xA);
	rNodeB->GetCoordinates2D(xB);
	rNodeC->GetCoordinates2D(xC);
	rNodeD->GetCoordinates2D(xD);

	/*!
	 * Building up Vector 1
	 *
	 * \f{align*}{
	 *  \vec{P}_1 = { A_1 \choose A_2 } + r { dir1_1 \choose dir1_2 }
	 *  \qquad ; \qquad
	 *  \vec{dir1} = \vec{B} - \vec{A}
	 * \f}
	 */
	dir1[0] = xB[0] - xA[0];
	dir1[1] = xB[1] - xA[1];

	/*!
	 * Building up Vector 2
	 *
	 * \f{align*}{
	 *  \vec{P}_2 = { C_1 \choose C_2 } - s { dir2_1 \choose dir2_2 }
	 *  \qquad ; \qquad
	 *  \vec{dir2} = \vec{D} - \vec{C}
	 * \f}
	 */
	dir2[0] = xD[0] - xC[0];
	dir2[1] = xD[1] - xC[1];

	/*!
	 * Solving System of equations
	 * \f{equation*}
	 *   \begin{bmatrix}
	 *     dir1_1 & dir2_1 \\
	 * 	dir1_2 & dir2_2
	 *   \end{bmatrix}
	 *   \begin{Bmatrix}
	 *     r \\
	 * 	s
	 *   \end{Bmatrix}
	 *   =
	 *   \begin{Bmatrix}
	 *     C_1 - A_1 \\
	 * 	   C_2 - A_2
	 *   \end{Bmatrix}
	 * \f}
	 *
	 * using determinant of the coeficiantmatrix
	 * \f[
	 *   det = dir1_1 \cdot dir2_2 - dir2_1 \cdot dir1_2
	 * \f]
	 */
	const double det(dir1[0]*dir2[1] - dir2[0]*dir1[1]);
	//! return false if \f$ det=0 \f$
	if( std::abs(det)< std::numeric_limits<double>::epsilon()*10  )
		return false;

	/*!
	 * Solution with respect to \f$r\f$ :
	 * \f$
	 * r = \frac{1}{det} \left[ dir2_2 \left(C_1-A_1\right) - dir2_1 \left(C_2-A_2\right) \right]
	 * \f$
	 */
	rDist1=(dir2[1]*(xC[0]-xA[0])-dir2[0]*(xC[1]-xA[1]))/det;

	/*!
	 * Solution in Respect to $s$
	 * \f[
	 * s = \frac{1}{det} \left[ dir1_2 \left(C_1-A_1\right) - dir1_1 \left(C_2-A_2\right) \right]
	 * \f]
	 */
	rDist2=(dir1[1]*(xC[0]-xA[0])-dir1[0]*(xC[1]-xA[1]))/det;

	/*!
	 *
	 * Intersection point
	 * \f$
	 *   \vec{x}_I = { A_1 \choose A_2 } + r { dir1_1 \choose dir1_2 }
	 * \f$
	 */
	xI[0]=xA[0]+rDist1*dir1[0];
	xI[1]=xA[1]+rDist1*dir1[1];
	rNodeI->SetCoordinates2D(xI);

	if( rDist1<=1 && 0<=rDist1 && rDist2<=1 && 0<=rDist2 ){
		return true;
	}else{
		return false;
	}


}


//! @brief computes the intersection point of a straight segment and a ray
//! @return bool intersect?
//! @param rNodeA (Input) ... const reference to first node of the first segment
//! @param rNodeB (Input) ... const reference to second node of the first segment
//! @param rNodeC (Input) ... const reference to first node of the second segment
//! @param rDdir2 (Input) ... const reference to the direction of the ray
//! @param rNodeI (Output) ... reference to intersection node of the segments
//! @param rDist1 (Output) ... normalized position of the intersection on the first segment ( xI=xA+r*(xB-xA) )
//! @todo  move this function to geometry class
bool NuTo::CrackExplicit2D::IntersectSegmentRay(	const NuTo::NodeBase* rNodeA, const NuTo::NodeBase* rNodeB,
														const NuTo::NodeBase* rNodeC, const double rDir2[2],
														NuTo::NodeBase* rNodeI, double & rDist1, double & rDist2)
{
	//! check if input nodes are 2D-nodes
	if(rNodeA->GetNumCoordinates()!=2 || rNodeB->GetNumCoordinates()!=2 || rNodeC->GetNumCoordinates()!=2 || rNodeI->GetNumCoordinates()!=2 )
		throw NuTo::MechanicsException("[NuTo::CrackExplicit2D::Intersect] don't got an 2D-Node!!!");

	double xA[2],xB[2], xC[2], xI[2], dir1[] = {0.0 , 0.0};
	rNodeA->GetCoordinates2D(xA);
	rNodeB->GetCoordinates2D(xB);
	rNodeC->GetCoordinates2D(xC);

	/*!
	 * Building up Vector 1
	 *
	 * \f{align*}{
	 *  \vec{P}_1 = { A_1 \choose A_2 } + r { dir1_1 \choose dir1_2 }
	 *  \qquad ; \qquad
	 *  \vec{dir1} = \vec{B} - \vec{A}
	 * \f}
	 */
	dir1[0] = xB[0] - xA[0];
	dir1[1] = xB[1] - xA[1];

	/*!
	 * Solving System of equations
	 * \f{equation*}
	 *   \begin{bmatrix}
	 *     dir1_1 & dir2_1 \\
	 * 	dir1_2 & dir2_2
	 *   \end{bmatrix}
	 *   \begin{Bmatrix}
	 *     r \\
	 * 	s
	 *   \end{Bmatrix}
	 *   =
	 *   \begin{Bmatrix}
	 *     C_1 - A_1 \\
	 * 	   C_2 - A_2
	 *   \end{Bmatrix}
	 * \f}
	 *
	 * using determinant of the coeficiantmatrix
	 * \f[
	 *   det = dir1_1 \cdot dir2_2 - dir2_1 \cdot dir1_2
	 * \f]
	 */
	const double det(dir1[0]*rDir2[1] - rDir2[0]*dir1[1]);
	//! return false if \f$ det=0 \f$
	if( std::abs(det)< std::numeric_limits<double>::epsilon()*10  )
		return false;

	/*!
	 * Solution with respect to \f$r\f$ :
	 * \f$
	 * r = \frac{1}{det} \left[ dir2_2 \left(C_1-A_1\right) - dir2_1 \left(C_2-A_2\right) \right]
	 * \f$
	 */
	rDist1=(rDir2[1]*(xC[0]-xA[0])-rDir2[0]*(xC[1]-xA[1]))/det;

	/*!
	 * Solution in Respect to $s$
	 * \f[
	 * s = \frac{1}{det} \left[ dir1_2 \left(C_1-A_1\right) - dir1_1 \left(C_2-A_2\right) \right]
	 * \f]
	 */
	rDist2=(dir1[1]*(xC[0]-xA[0])-dir1[0]*(xC[1]-xA[1]))/det;

	/*!
	 *
	 * Intersection point
	 * \f$
	 *   \vec{x}_I = { A_1 \choose A_2 } + r { dir1_1 \choose dir1_2 }
	 * \f$
	 */
	xI[0]=xA[0]+rDist1*dir1[0];
	xI[1]=xA[1]+rDist1*dir1[1];
	rNodeI->SetCoordinates2D (xI);

	//! return true if \f$ 0 \leq {r,s} \leq 1 \f$
	//! return false otherwise
	if( rDist1<=1 && 0<=rDist1 ){
		return true;
	}else{
		return false;
	}


}

#ifdef ENABLE_VISUALIZE
void NuTo::CrackExplicit2D::Visualize(VisualizeUnstructuredGrid& rVisualize) const
{
    unsigned int NumVisualizationPoints=size();
    if(NumVisualizationPoints < 2)
    	throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] at least 2 nodes necessary for cracks");

    unsigned int NumVisualizationCells=size()-1;



    // calculate global point coordinates and store point at the visualize object
    std::vector<unsigned int> PointIdVec;
    for (std::list<NuTo::NodeBase*>::const_iterator it = this->begin(); it!= this->end(); it++)
    {
        double GlobalPointCoor[3] = {(*it)->GetCoordinate(0) , (*it)->GetCoordinate(1), 0.0};

        unsigned int PointId = rVisualize.AddPoint(GlobalPointCoor);
        PointIdVec.push_back(PointId);
    }


    // store cells at the visualize object
    std::vector<unsigned int> CellIdVec;
    for (unsigned int CellCount = 0; CellCount < NumVisualizationCells; CellCount++)
    {
			unsigned int Points[2];
            Points[0] = PointIdVec[CellCount    ];
            Points[1] = PointIdVec[CellCount + 1];
            unsigned int CellId = rVisualize.AddLineCell(Points);
            CellIdVec.push_back(CellId);
    }

}
#endif // ENABLE_VISUALIZE
