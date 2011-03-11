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

#include <wx/string.h>

#include <vtkSmartPointer.h>

#include <vector>

class ResultDataSourceVTKFromFile : public nutogui::ResultDataSourceVTK
{
  typedef std::pair<vtkSmartPointer<vtkDataSet>, wxString> DataSetPair;
  std::vector<DataSetPair> datasets;
  
  vtkSmartPointer<vtkDataSet> ReadDataFromFile (const char* filename);
public:
  ResultDataSourceVTKFromFile ();
  
  void AddDataSet (const char* filename, const wxString& name);
  
  size_t GetNumDataSets () const { return datasets.size(); }
  vtkDataSet* GetDataSet (size_t index);
  const wxString& GetDataSetName (size_t index) const;
};

#endif // __RESULTDATASOURCEVTKFROMFILE_H__
