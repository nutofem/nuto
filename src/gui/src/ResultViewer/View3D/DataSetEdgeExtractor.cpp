/**\file
 * VTK data set 'filter' to extract all edges.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"

#include "DataSetEdgeExtractor.h"

#include <vtkCellData.h>
#include <vtkDataObject.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkLine.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

#include "DataSetEdgeIterator.h"
#include "DataSetHelpers.h"

namespace nutogui
{
  int DataSetEdgeExtractor::RequestData (vtkInformation* request,
					 vtkInformationVector** inputVector,
					 vtkInformationVector* outputVector)
  {
    // get the info objects
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation* outInfo = outputVector->GetInformationObject(0);

    // get the input and output
    vtkDataSet* input = vtkDataSet::SafeDownCast (
      inInfo->Get (vtkDataObject::DATA_OBJECT()));
    vtkUnstructuredGrid* output = vtkUnstructuredGrid::SafeDownCast (
      outInfo->Get (vtkDataObject::DATA_OBJECT()));
   
    output->Initialize();
    output->Allocate();
    DataSetHelpers::CopyPoints (output, input);
    output->GetPointData()->ShallowCopy (input->GetPointData());
    
    std::vector<vtkIdType> origCells;
    DataSetEdgeIterator edges (input);
    while (edges.HasNext ())
    {
      vtkIdType origCell;
      DataSetEdgeIterator::Edge edge (edges.Next (origCell));
      
      vtkIdType linePts[2] = { edge.first, edge.second };
      output->InsertNextCell (VTK_LINE, 2, linePts);
      origCells.push_back (origCell);
    }
    
    vtkIdType numCells = output->GetNumberOfCells ();
    output->GetCellData()->CopyAllocate (input->GetCellData(), numCells);
    for (vtkIdType cell = 0; cell < numCells; cell++)
      output->GetCellData()->CopyData (input->GetCellData(), origCells[cell], cell);
    
    output->Squeeze();
    
    return 1;
  }
} // namespace nutogui
