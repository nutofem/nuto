/**\file
 * ResultDataSourceVTK implementation that reads the data from a file.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __RESULTDATASOURCEVTKFROMFILE_H__
#define __RESULTDATASOURCEVTKFROMFILE_H__

#include "ResultDataSourceVTK.h"

#include <vtkSmartPointer.h>

class ResultDataSourceVTKFromFile : public nutogui::ResultDataSourceVTK
{
  vtkSmartPointer<vtkDataSet> dataset;
  
  void ReadDataFromFile (const char* filename);
public:
  ResultDataSourceVTKFromFile (const char* filename);
  
  vtkDataSet* QueryDataSet ();
};

#endif // __RESULTDATASOURCEVTKFROMFILE_H__
