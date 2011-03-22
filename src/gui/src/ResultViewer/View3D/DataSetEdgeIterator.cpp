/**\file
 * Class to iterate over all edges of an VTK data set.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"

#include "DataSetEdgeIterator.h"

#include <vtkCell.h>
#include <vtkDataSet.h>

namespace nutogui
{
  DataSetEdgeIterator::DataSetEdgeIterator (vtkDataSet* dataset)
   : dataset (dataset),
     currentCellId (0),
     currentEdgeCellId (0),
     currentEdgeCellPointId (0)
  {
    if (dataset->GetNumberOfCells () > 0)
    {
      currentCell = dataset->GetCell (currentCellId);
      
      if (currentCell->GetNumberOfEdges () > 0)
      {
	currentEdgeCell = currentCell->GetEdge (currentEdgeCellId);
	
	if (currentEdgeCell->GetNumberOfPoints() >= 2)
	{
	  lastPoint = currentEdgeCell->GetPointId (currentEdgeCellPointId);
	  currentEdgeCellPointId++;
	}
	else
	  // Need at least two points to make a line :P
	  currentEdgeCell = nullptr;
      }
    }
  }
  
  bool DataSetEdgeIterator::HasNext() const
  {
    return currentEdgeCell;
  }
  
  DataSetEdgeIterator::Edge DataSetEdgeIterator::Next()
  {
    vtkIdType pt = currentEdgeCell->GetPointId (currentEdgeCellPointId);
    Edge edge (std::make_pair (lastPoint, pt));
    lastPoint = pt;
    Advance ();
    return edge;
  }

  DataSetEdgeIterator::Edge DataSetEdgeIterator::Next (vtkIdType& origCell)
  {
    origCell = currentCellId;
    vtkIdType pt = currentEdgeCell->GetPointId (currentEdgeCellPointId);
    Edge edge (std::make_pair (lastPoint, pt));
    lastPoint = pt;
    Advance ();
    return edge;
  }

  void DataSetEdgeIterator::Advance ()
  {
    // Next point in cell ...
    currentEdgeCellPointId++;
    // ... if past last point, ...
    while (currentEdgeCell && (currentEdgeCellPointId >= currentEdgeCell->GetNumberOfPoints ()))
    {
      // ... go to next edge cell ...
      currentEdgeCellId++;
      // ... if past last edge cell, ...
      while (currentCell && (currentEdgeCellId >= currentCell->GetNumberOfEdges ()))
      {
	// ... go to next cell ...
	currentCellId++;
	// ... if past last cell, ...
	if (currentCellId >= dataset->GetNumberOfCells ())
	{
	  // ... we're done.
	  currentCell = nullptr;
	}
	else
	{
	  // Else, grab new cell.
	  currentCell = dataset->GetCell (currentCellId);
	  currentEdgeCellId = 0;
	}
      }
      if (currentCell)
      {
	// Grab first edge cell.
	currentEdgeCell = currentCell->GetEdge (currentEdgeCellId);
      }
      else
      {
	currentEdgeCell = nullptr;
      }
      if (currentEdgeCell && (currentEdgeCell->GetNumberOfPoints() >= 2))
      {
	// Grab + store first point of cell. Will be used by next call of Next().
	currentEdgeCellPointId = 0;
	lastPoint = currentEdgeCell->GetPointId (currentEdgeCellPointId);
	currentEdgeCellPointId++;
      }
      else
	// Need at least two points to make a line :P
	currentEdgeCell = nullptr;
    }
  }
} // namespace nutogui
