// $Id$

#include <assert.h>

#include <cmath>
#include <sstream>
#include <iostream>
#include <iterator>

#include <boost/foreach.hpp>

#include "nuto/base/Debug.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NModTriangle.h"
#include "nuto/mechanics/integrationtypes/IntegrationPointBase.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss3Ip.h"
#include "nuto/mechanics/MechanicsException.h"

#ifdef HAVE_CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel CGAL_Kernel;
typedef CGAL_Kernel::Point_2 CGAL_Point;
typedef CGAL::Triangulation_vertex_base_2<CGAL_Kernel> CGAL_Vb;
typedef CGAL::Delaunay_mesh_face_base_2<CGAL_Kernel> CGAL_Fb;
typedef CGAL::Triangulation_data_structure_2<CGAL_Vb, CGAL_Fb> CGAL_Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<CGAL_Kernel, CGAL_Tds> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;

typedef CDT::Vertex_handle CGAL_Vertex_handle;
#endif //HAVE_CGAL


//! @brief constructor
//! @author Daniel Arnold, ISM
//! @date February 2011
NuTo::IntegrationType2D4NModTriangle::IntegrationType2D4NModTriangle(const std::string rName)
{
	mName=rName;
	//~ 
//~ DBG_POSITION_INFO("creating triangulated integration scheme")
//~ DBG_PRINT_VAL(mName)


}

//! @brief creates new integration-cells/order/area
//! @author Daniel Arnold, ISM
//! @date May 2012
//! @param rArea (Input) polygonal surface of integration area
//! @param rOrder (Input) integration order (or number of integration points)
void NuTo::IntegrationType2D4NModTriangle::AddIntegrationPoints(std::vector< std::vector<double> > & rArea, const unsigned short rOrder)
{
#ifndef HAVE_CGAL
	throw MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Cannot triangulate the integration area without CGAL library!");
#else //HAVE_CGAL
	//! triangulation of the domain
	/// preparing the domain
	std::vector<CGAL_Point> integrationArea;
	for(std::vector< std::vector<double> >::iterator itArea=rArea.begin();itArea!=rArea.end();itArea++){
		assert((*itArea).size()==2);
		integrationArea.push_back((CGAL_Point((*itArea)[0] , (*itArea)[1])));
	}
	/// building up the constraints for the triangulations
	CDT cdtReg;
	for(size_t i=0; i<integrationArea.size();i++)
		cdtReg.insert_constraint(integrationArea[i],integrationArea[(i+1)%integrationArea.size()]);
	//! refinement of the mesh
	CGAL::refine_Delaunay_mesh_2(cdtReg, Criteria(0.25, 0.5));
	std::vector<double> shapeFunctions(3);
	for (CDT::Face_iterator fiter=cdtReg.faces_begin(); fiter!=cdtReg.faces_end(); fiter++){
		if(!fiter->is_valid()) continue;
		/// get reference integration scheme
		NuTo::IntegrationType2D3NGauss3Ip intTypeRef;
		intTypeRef=IntegrationType2D3NGauss3Ip();
		const int numIp=intTypeRef.GetNumIntegrationPoints();
#ifdef ENABLE_VISUALIZE
		unsigned int NumVisualizationPoints;
		std::vector< double > VisualizationPointLocalCoordinates;
		unsigned int NumVisualizationCells;
		std::vector< NuTo::CellBase::eCellTypes > VisualizationCellType;
		std::vector< unsigned int > VisualizationCellsIncidence;
		std::vector< unsigned int > VisualizationCellsIP ;
		intTypeRef.GetVisualizationCells 	( NumVisualizationPoints, VisualizationPointLocalCoordinates, NumVisualizationCells, VisualizationCellType, VisualizationCellsIncidence, VisualizationCellsIP);
#endif //ENABLE_VISUALIZE
		/// build up new integration scheme
		for (unsigned short ip=0; ip<numIp;++ip){
			double ipCoord[2]={0,0};
			double locCoord[2]={0,0};
			std::vector<double> ipCell;
			intTypeRef.GetLocalIntegrationPointCoordinates2D(ip,locCoord);
			/// compute coordinate within the complete area
			CalculateShapeFunctions(locCoord, shapeFunctions);
			//! get ip point coordinates
			for(unsigned short i=0; i<3;++i){
				ipCoord[0]+=shapeFunctions[i]*fiter->vertex(i)->point()[0];
				ipCoord[1]+=shapeFunctions[i]*fiter->vertex(i)->point()[1];
			}
			//! weighting factor
			/** area of the Integration cell (in natural coordinates):
				\f[ d = \frac{1}{2}\left( \left(x_1-x_0\right)\left(y_2-x_0\right) - \left(x_2-x_0\right)\left(y_1-x_0\right)\right)
				\f]
			**/
			const double faceArea=0.5*std::abs(	 (fiter->vertex(1)->point()[0]-fiter->vertex(0)->point()[0])*(fiter->vertex(2)->point()[1]-fiter->vertex(0)->point()[1])
											-(fiter->vertex(2)->point()[0]-fiter->vertex(0)->point()[0])*(fiter->vertex(1)->point()[1]-fiter->vertex(0)->point()[1]));
			const double ipWeight=2*intTypeRef.GetIntegrationPointWeight(ip)*faceArea;
//~ DBG_PRINT_VAL(ipWeight)
#ifdef ENABLE_VISUALIZE
			//! get corresponding visualization cell
			std::vector< unsigned int > ::iterator iFound;
			iFound = find ( VisualizationCellsIP.begin (), VisualizationCellsIP.end (), ip );
			size_t idFound=distance(VisualizationCellsIP.begin (),iFound);
			size_t startId=0;
			unsigned short numPtCell=0;
			/// find start position
			for(unsigned short iVis=0; iVis<idFound;++iVis){
				switch (VisualizationCellType[iVis]) {
					case NuTo::CellBase::TRIANGLE:
						startId+=3;
						break;
					case NuTo::CellBase::QUAD:
						startId+=4;
						break;
					default:
						throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Unknown visualization cell!");
				}
			}
			switch (VisualizationCellType[idFound]) {
				case NuTo::CellBase::TRIANGLE:
					numPtCell=3;
					break;
				case NuTo::CellBase::QUAD:
					numPtCell=4;
					break;
				default:
					throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Unknown visualization cell!");
			}
			for(unsigned short iPtCell=0;iPtCell<numPtCell;++iPtCell){
				locCoord[0]=VisualizationPointLocalCoordinates[VisualizationCellsIncidence[startId+iPtCell]*2];
				locCoord[1]=VisualizationPointLocalCoordinates[VisualizationCellsIncidence[startId+iPtCell]*2+1];
				CalculateShapeFunctions(locCoord, shapeFunctions);
				//! get ip point coordinates
				double ptCoord[2]={0,0};
				for(unsigned short i=0; i<3;++i){
					ptCoord[0]+=shapeFunctions[i]*fiter->vertex(i)->point()[0];
					ptCoord[1]+=shapeFunctions[i]*fiter->vertex(i)->point()[1];
				}
				ipCell.push_back(ptCoord[0]);
				ipCell.push_back(ptCoord[1]);
			}
#endif //ENABLE_VISUALIZE
			std::vector<double> ipCoords(2);
			ipCoords.assign(ipCoord,ipCoord+2);
			IntegrationPointBase newIP(ipCoords,ipWeight,ipCell);
			this->AddIntegrationPoint(newIP);
		}
		
	}
#endif //HAVE_CGAL
}

//! @brief calculates the shape functions for the triangular area
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
void NuTo::IntegrationType2D4NModTriangle::CalculateShapeFunctions(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)const
{
	assert(rShapeFunctions.size()==3);
    rShapeFunctions[0] = 1-rNaturalCoordinates[0]-rNaturalCoordinates[1];
    rShapeFunctions[1] = rNaturalCoordinates[0];
    rShapeFunctions[2] = rNaturalCoordinates[1];
}
