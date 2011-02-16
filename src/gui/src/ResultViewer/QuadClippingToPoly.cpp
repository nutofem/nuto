/**\file
 * Subclass of vtkQuad with different clipping output.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"

#include "QuadClippingToPoly.h"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkIncrementalPointLocator.h>
#include <vtkPointData.h>
#include <vtkPoints.h>

namespace nutogui
{
  /* 
    QUAD CLIPPING
    
    Pretty much copied from vtkQuad implementation.
    
    Differences:
    - Quad clipping case table. (Commented lines are original entries.)
    - Check for degenerate polys generalized.
  */
  
  /*=========================================================================

    Program:   Visualization Toolkit
    Module:    $RCSfile: vtkQuad.h,v $

    Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
    All rights reserved.
    See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

      This software is distributed WITHOUT ANY WARRANTY; without even
      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
      PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/
  
  typedef int QUAD_EDGE_LIST;
  typedef struct {
	QUAD_EDGE_LIST edges[9];
  } QUAD_CASES;
  
  static const QUAD_CASES quadCases[] = { 
  {{  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1}}, // 0
  {{   3, 100,   0,   3,  -1,  -1,  -1,  -1, -1}}, // 1
  {{   3, 101,   1,   0,  -1,  -1,  -1,  -1, -1}}, // 2
  {{   4, 100, 101,   1,   3,  -1,  -1,  -1, -1}}, // 3
  {{   3, 102,   2,   1,  -1,  -1,  -1,  -1, -1}}, // 4
  //{{   3, 100,   0,   3,   3, 102,   2,   1,   4,   0,   1,   2,   3,  -1}}, // 5
  {{   6, 100,   0,   1, 102,   2,   3,  -1, -1}}, // 5
  {{   4, 101, 102,   2,   0,  -1,  -1,  -1, -1}}, // 6
  //{{   3, 100, 101,   3,   3, 101,   2,   3,   3, 101, 102,   2,  -1,  -1}}, // 7
  {{   5, 100, 101, 102,   2,   3,  -1,  -1, -1}}, // 7
  {{   3, 103,   3,   2,  -1,  -1,  -1,  -1, -1}}, // 8
  {{   4, 100,   0,   2, 103,  -1,  -1,  -1, -1}}, // 9
  //{{   3, 101,   1,   0,   3, 103,   3,   2,   4,   0,   1,   2,   3,  -1}}, // 10
  {{   6, 101,   1,   2, 103,   3,  -1,  -1, -1}}, // 10
  //{{   3, 100, 101,   1,   3, 100,   1,   2,   3, 100,   2, 103,  -1,  -1}}, // 11
  {{   5, 100, 101,   1,   2, 103,  -1,  -1, -1}}, // 11
  {{   4, 102, 103,   3,   1,  -1,  -1,  -1, -1}}, // 12
  //{{   3, 100,   0, 103,   3,   0,   1, 103,   3,   1, 102, 103,  -1,  -1}}, // 13
  {{   5, 100,   0,   1, 102, 103,  -1,  -1, -1}}, // 13
  //{{   3,   0, 101, 102,   3,   0, 102,   3,   3, 102, 103,   3,  -1,  -1}}, // 14
  {{   5,   0, 101, 102, 103,   3,  -1,  -1, -1}}, // 14
  {{   4, 100, 101, 102, 103,  -1,  -1,  -1, -1}}, // 15
  };

  static const QUAD_CASES quadCasesComplement[] = { 
  {{  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, -1}}, // 0
  {{   3, 100,   0,   3,  -1,  -1,  -1,  -1, -1}}, // 1
  {{   3, 101,   1,   0,  -1,  -1,  -1,  -1, -1}}, // 2
  {{   4, 100, 101,   1,   3,  -1,  -1,  -1, -1}}, // 3
  {{   3, 102,   2,   1,  -1,  -1,  -1,  -1, -1}}, // 4
  {{   3, 100,   0,   3,   3, 102,   2,   1, -1}}, // 5
  {{   4, 101, 102,   2,   0,  -1,  -1,  -1, -1}}, // 6
  //{{   3, 100, 101,   3,   3, 101,   2,   3,   3, 101, 102,   2,  -1,  -1}}, // 7
  {{   5, 100, 101, 102,   2,   3,  -1,  -1, -1}}, // 7
  {{   3, 103,   3,   2,  -1,  -1,  -1,  -1, -1}}, // 8
  {{   4, 100,   0,   2, 103,  -1,  -1,  -1, -1}}, // 9
  {{   3, 101,   1,   0,   3, 103,   3,   2, -1}}, // 10
  //{{   3, 100, 101,   1,   3, 100,   1,   2,   3, 100,   2, 103,  -1,  -1}}, // 11
  {{   5, 100, 101,   1,   2, 103,  -1,  -1, -1}}, // 11
  {{   4, 102, 103,   3,   1,  -1,  -1,  -1, -1}}, // 12
  //{{   3, 100,   0, 103,   3,   0,   1, 103,   3,   1, 102, 103,  -1,  -1}}, // 13
  {{   5, 100,   0,   1, 102, 103,  -1,  -1, -1}}, // 13
  //{{   3,   0, 101, 102,   3,   0, 102,   3,   3, 102, 103,   3,  -1,  -1}}, // 14
  {{   5,   0, 101, 102, 103,   3,  -1,  -1, -1}}, // 14
  {{   4, 100, 101, 102, 103,  -1,  -1,  -1, -1}}, // 15
  };

  static const int edges[4][2] = { {0,1}, {1,2}, {3,2}, {0,3} };
  
  //----------------------------------------------------------------------------
  // Clip this quad using scalar value provided. Like contouring, except
  // that it cuts the quad to produce other quads and/or triangles.
  void QuadClippingToPoly::Clip(double value, vtkDataArray *cellScalars, 
				vtkIncrementalPointLocator *locator, vtkCellArray *polys,
				vtkPointData *inPd, vtkPointData *outPd,
				vtkCellData *inCd, vtkIdType cellId, vtkCellData *outCd,
				int insideOut)
  {
    static const int CASE_MASK[4] = {1,2,4,8};
    const QUAD_CASES *quadCase;
    const QUAD_EDGE_LIST  *edge;
    int i, j, index;
    const int* vert;
    int e1, e2;
    int newCellId;
    vtkIdType pts[4];
    int vertexId;
    double t, x1[3], x2[3], x[3], deltaScalar;
    double scalar0, scalar1, e1Scalar;

    // Build the index into the case table
    if ( insideOut )
      {    
      for ( i=0, index = 0; i < 4; i++)
	{
	if (cellScalars->GetComponent(i,0) <= value)
	  {
	  index |= CASE_MASK[i];
	  }
	}
      // Select case based on the index and get the list of edges for this case
      quadCase = quadCases + index;
      }    
    else
      {
      for ( i=0, index = 0; i < 4; i++)
	{
	if (cellScalars->GetComponent(i,0) > value)
	  {
	  index |= CASE_MASK[i];
	  }
	}
      // Select case based on the index and get the list of edges for this case
      quadCase = quadCasesComplement + index;
      }

    edge = quadCase->edges;

    // generate each quad
    for ( ; edge[0] > -1; edge += edge[0]+1 )
      {
      for (i=0; i < edge[0]; i++) // insert quad or triangle
	{
	// vertex exists, and need not be interpolated
	if (edge[i+1] >= 100)
	  {
	  vertexId = edge[i+1] - 100;
	  this->Points->GetPoint(vertexId, x);
	  if ( locator->InsertUniquePoint(x, pts[i]) )
	    {
	    outPd->CopyData(inPd,this->PointIds->GetId(vertexId),pts[i]);
	    }
	  }

	else //new vertex, interpolate
	  {
	  vert = edges[edge[i+1]];

	  // calculate a preferred interpolation direction
	  scalar0 = cellScalars->GetComponent(vert[0],0);
	  scalar1 = cellScalars->GetComponent(vert[1],0);
	  deltaScalar = scalar1 - scalar0;

	  if (deltaScalar > 0)
	    {
	    e1 = vert[0]; e2 = vert[1];
	    e1Scalar = scalar0;
	    }
	  else
	    {
	    e1 = vert[1]; e2 = vert[0];
	    e1Scalar = scalar1;
	    deltaScalar = -deltaScalar;
	    }

	  // linear interpolation
	  if (deltaScalar == 0.0)
	    {
	    t = 0.0;
	    }
	  else
	    {
	    t = (value - e1Scalar) / deltaScalar;
	    }

	  this->Points->GetPoint(e1, x1);
	  this->Points->GetPoint(e2, x2);

	  for (j=0; j<3; j++)
	    {
	    x[j] = x1[j] + t * (x2[j] - x1[j]);
	    }

	  if ( locator->InsertUniquePoint(x, pts[i]) )
	    {
	    vtkIdType p1 = this->PointIds->GetId(e1);
	    vtkIdType p2 = this->PointIds->GetId(e2);
	    outPd->InterpolateEdge(inPd,pts[i],p1,p2,t);
	    }
	  }
	}
      // check for degenerate output
      if ( edge[0] == 3 ) //i.e., a triangle
	{
	if (pts[0] == pts[1] || pts[0] == pts[2] || pts[1] == pts[2] )
	  {
	  continue;
	  }
	}
      else if ( edge[0] == 4 ) // a quad
	{
	if ((pts[0] == pts[3] && pts[1] == pts[2]) || 
	    (pts[0] == pts[1] && pts[3] == pts[2]) )
	  {
	  continue;
	  }
	}
      else
      {
	int numSame = 0;
	for (int p = 1; p < edge[0]; p++)
	{
	  for (int p2 = 0; p2 < p; p2++)
	  {
	    if (pts[p] == pts[p2]) numSame++;
	  }
	}
	if (numSame >= edge[0]-2) continue;
      }

      newCellId = polys->InsertNextCell(edge[0],pts);
      outCd->CopyData(inCd,cellId,newCellId);
      }
  }

} // namespace nutogui
