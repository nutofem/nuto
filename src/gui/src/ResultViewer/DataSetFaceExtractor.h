/**\file
 * VTK data set 'filter' to extract all faces.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __NUTOGUI_DATASETFACEEXTRACTOR_H__
#define __NUTOGUI_DATASETFACEEXTRACTOR_H__

#include <vtkUnstructuredGridAlgorithm.h>

namespace nutogui
{
  class DataSetFaceExtractor : public vtkUnstructuredGridAlgorithm
  {
  public:
    static DataSetFaceExtractor* New ()
    { return new DataSetFaceExtractor; }
    
    int RequestData (vtkInformation* request,
		     vtkInformationVector** inputVector,
		     vtkInformationVector* outputVector);
  };
} // namespace nutogui

#endif // __NUTOGUI_DATASETFACEEXTRACTOR_H__
