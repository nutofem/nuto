// $Id$
#pragma once


#include "nuto/mechanics/cracks/CrackBase.h"
#include <list>







namespace NuTo
{
class ElementBase;
class NodeBase;
#ifdef ENABLE_VISUALIZE
class VisualizeUnstructuredGrid;
#endif // ENABLE_VISUALIZE


//! @author Daniel Arnold, ISM
//! @date October 2010
//! @brief ... class for 2D explicit cracks
class CrackExplicit2D : public NuTo::CrackBase, public std::list<NuTo::NodeBase*>
{

public:
    //! @brief constructor
    CrackExplicit2D();
    CrackExplicit2D(const std::string);

	//! @brief destructor
	~CrackExplicit2D();

    //! @brief gives the crack informations
    //! @return group type
    void Info(int rVerboseLevel)const;

	//! @brief initiate the crack into the structure
    //! @param rElements ... vector of element pointer
    //! @param rCrackedElements ... vector of cracked element pointers
	void Initiate(std::vector<ElementBase*>& rElements, std::vector<ElementBase*>& rCrackedElements);

	//! @brief gives the number of crackpoints
    //! @return number of group members
    int GetNumPoints()const
    {
        return (int)this->size();
    }

    //! @brief appends a crack point to the end of the crack
    //! @param rMember new member
    void PushBack(NuTo::NodeBase* rMember)
    {
        this->push_back(rMember);
    }

    //! @brief appends a crack point to the begin of the crack
    //! @param rMember new member
    void PushFront(NuTo::NodeBase* rMember)
    {
        this->push_front(rMember);
    }

    //! @brief removes the last crack point
    void PopBack()
    {
        this->pop_back();
    }

    //! @brief removes the last crack point
    void PopFront()
    {
        this->pop_front();
    }

    // @brief computes the signed distance of a node to the crack
    // @return signed distance of the node
    // @param rNode ... const reference to node to be evaluated
    // @todo  move this function to geometry class
    double NodeDistance(const NuTo::NodeBase* rNode) const;

    // @brief computes the intersection point of a crack and a line
    // @return a bool if the line segment AB is cracked or not
    // @param rNodeA (Input) ... const reference to first node of the line
    // @param rNodeB (Input) ... const reference to second node of the line
    // @param rNodeI (Output) ... reference to intersection node of the line with the crack
    // @param rDist  (Output) ... relative coordinate of the intersection
    // @param rSeg   (Output) ... id of the crack segment which intersects the line segment AB (max(size_t) if AB is not cracked)
    // @todo  move this function to geometry class
    bool Intersect(const NuTo::NodeBase* rNodeA, const NuTo::NodeBase* rNodeB, NuTo::NodeBase* rNodeI, double& rDist, size_t& rSeg);

    // @brief computes the intersection point between the end-ray of this crack and a line segment
    // @return an unsigned short if the line segment AB is intersected or not (0=no intersection, 1=front, 2=end)
    // @param rNodeA (Input) ... const reference to first node of the line
    // @param rNodeB (Input) ... const reference to second node of the line
    // @param rNodeI (Output) ... reference to intersection node of the line with the crack
    // @param rDist  (Output) ... relative coordinate of the segment intersection
    unsigned short ExtendEnd(const NuTo::NodeBase* rNodeA, const NuTo::NodeBase* rNodeB, NuTo::NodeBase* rNodeI, double& rDist);

    // @brief computes the intersection point of two straight segments
	// @return bool intersect?
	// @param rNodeA (Input) ... const reference to first node of the first segment
	// @param rNodeB (Input) ... const reference to second node of the first segment
	// @param rNodeC (Input) ... const reference to first node of the second segment
	// @param rNodeD (Input) ... const reference to second node of the second segment
	// @param rNodeI (Output) ... reference to intersection node of the segments
	// @param rDist1 (Output) ... normalized position of the intersection on the first segment ( xI=xA+r*(xB-xA) )
	// @param rDist2 (Output) ... normalized position of the intersection on the sedond segment ( xI=xC+s*(xD-xC) )
	// @todo  move this function to geometry class
    bool IntersectSegmentSegment(	const NuTo::NodeBase* rNodeA, const NuTo::NodeBase* rNodeB,  const NuTo::NodeBase* rNodeC, const NuTo::NodeBase* rNodeD,
									NuTo::NodeBase* rNodeI, double & rDist1, double & rDist2);

    // @brief computes the intersection point of a straight segment and a ray
    // @return bool intersect?
    // @param rNodeA (Input) ... const reference to first node of the first segment
    // @param rNodeB (Input) ... const reference to second node of the first segment
    // @param rNodeC (Input) ... const reference to first node of the second segment
    // @param rDdir2 (Input) ... const reference to the direction of the ray
    // @param rNodeI (Output) ... reference to intersection node of the segments
    // @param rDist1 (Output) ... normalized position of the intersection on the first segment ( xI=xA+r*(xB-xA) )
    // @todo  move this function to geometry class
    bool IntersectSegmentRay(	const NuTo::NodeBase* rNodeA, const NuTo::NodeBase* rNodeB, const NuTo::NodeBase* rNodeC, const double rDir2[2],
    							NuTo::NodeBase* rNodeI, double & rDist1, double & rDist2);
#ifdef ENABLE_VISUALIZE
    virtual void Visualize(VisualizeUnstructuredGrid& rVisualize) const;
#endif // ENABLE_VISUALIZE

};
}//namespace NuTo

