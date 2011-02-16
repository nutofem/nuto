/**\file
 * VTK data set 'filter' to extract all edges.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __NUTOGUI_DATASETEDGEEXTRACTOR_H__
#define __NUTOGUI_DATASETEDGEEXTRACTOR_H__

#include <vtkUnstructuredGridAlgorithm.h>

namespace nutogui
{
  class DataSetEdgeExtractor : public vtkUnstructuredGridAlgorithm
  {
  public:
    static DataSetEdgeExtractor* New ()
    { return new DataSetEdgeExtractor; }
    
    int RequestData (vtkInformation* request,
		     vtkInformationVector** inputVector,
		     vtkInformationVector* outputVector);
  };
} // namespace nutogui

#endif // __NUTOGUI_DATASETEDGEEXTRACTOR_H__
