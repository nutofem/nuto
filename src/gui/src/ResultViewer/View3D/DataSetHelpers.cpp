/**\file
 * Helpers to deal with data sets.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"

#include "DataSetHelpers.h"

#include <vtkDataSet.h>
#include <vtkPointSet.h>
#include <vtkSmartPointer.h>

namespace nutogui
{
  void DataSetHelpers::CopyPoints (vtkPointSet* dest, vtkDataSet* src)
  {
    vtkPointSet* srcPtSet = vtkPointSet::SafeDownCast (src);
    if (srcPtSet)
    {
      // Fast path
      dest->SetPoints (srcPtSet->GetPoints ());
    }
    else
    {
      vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New ();
      newPoints->SetNumberOfPoints (src->GetNumberOfPoints());
      for (vtkIdType point = 0; point < src->GetNumberOfPoints(); point++)
      {
	double pt[3];
	src->GetPoint (point, pt);
	newPoints->SetPoint (point, pt);
      }
      dest->SetPoints (newPoints);
    }
  }
} // namespace nutogui
