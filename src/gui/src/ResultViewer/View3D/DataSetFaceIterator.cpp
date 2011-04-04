/**\file
 * Class to iterate over all faces of an VTK data set.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"

#include "DataSetFaceIterator.h"

#include <vtkCell.h>
#include <vtkDataSet.h>
#include <vtkGenericCell.h>

namespace nutogui
{
  DataSetFaceIterator::DataSetFaceIterator (vtkDataSet* dataset)
   : dataset (dataset),
     currentCellId (0),
     currentFaceCellId (0),
     nextStorage (0)
  {
    if (dataset->GetNumberOfCells () > 0)
    {
      for (size_t i = 0; i < sizeof(cellStorage)/sizeof(cellStorage[0]); i++)
	cellStorage[i] = vtkSmartPointer<vtkGenericCell>::New();
      
      currentCell = GetCell (currentCellId);
     
      if (!HasCellTypeFaces (currentCell->GetCellType()))
      {
	nextFaceCell = currentCell;
	currentFaceCellId = currentCell->GetNumberOfFaces ()-1;
      }
      else if (currentCell->GetNumberOfFaces () > 0)
      {
	nextFaceCell = currentCell->GetFace (currentFaceCellId);
	currentFaceCellId++;
      }
    }
  }
  
  bool DataSetFaceIterator::HasNext() const
  {
    return nextFaceCell;
  }
  
  vtkCell* DataSetFaceIterator::Next()
  {
    vtkCell* cell = nextFaceCell;
    Advance ();
    return cell;
  }

  vtkCell* DataSetFaceIterator::Next (vtkIdType& origCell)
  {
    origCell = currentCellId;
    vtkCell* cell = nextFaceCell;
    Advance ();
    return cell;
  }

  void DataSetFaceIterator::Advance ()
  {
    // Next face in cell ...
    currentFaceCellId++;
    // ... if past last face, ...
    while (currentCell && (currentFaceCellId >= currentCell->GetNumberOfFaces ()))
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
	currentCell = GetCell (currentCellId);
	currentFaceCellId = 0;
	if (!HasCellTypeFaces (currentCell->GetCellType()))
	{
	  nextFaceCell = currentCell;
	  currentFaceCellId = currentCell->GetNumberOfFaces ()-1;
	}
      }
    }
    if (currentCell)
    {
      // Grab first face cell.
      if (HasCellTypeFaces (currentCell->GetCellType()))
	nextFaceCell = currentCell->GetFace (currentFaceCellId);
    }
    else
    {
      nextFaceCell = nullptr;
    }
  }
  
  bool DataSetFaceIterator::HasCellTypeFaces (int cellType)
  {
    switch (cellType)
    {
    case VTK_HEXAHEDRON:
    case VTK_WEDGE:
    case VTK_PYRAMID:
    case VTK_PENTAGONAL_PRISM:
    case VTK_HEXAGONAL_PRISM:
      return true;
    }
    
    return false;
  }
  
  vtkCell* DataSetFaceIterator::GetCell (vtkIdType index)
  {
    vtkGenericCell* cell = cellStorage[nextStorage];
    dataset->GetCell (index, cell);
    nextStorage = (nextStorage + 1) % (sizeof(cellStorage)/sizeof(cellStorage[0]));
    return cell;
  }
} // namespace nutogui
