/**\file
 * ResultDataSourceVTK implementation that reads the data from a file.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"
#include "ResultDataSourceVTKFromFile.h"

#include <vtkDataSet.h>
#include <vtkDataSetReader.h>

ResultDataSourceVTKFromFile::ResultDataSourceVTKFromFile ()
{
}

void ResultDataSourceVTKFromFile::AddDataSet (const char* filename, const wxString& name)
{
  vtkSmartPointer<vtkDataSet> dataset (ReadDataFromFile (filename));
  if (dataset)
  {
    datasets.push_back (std::make_pair (dataset, name));
  }
}

vtkDataSet* ResultDataSourceVTKFromFile::GetDataSet (size_t index)
{
  return datasets[index].first;
}

const wxString& ResultDataSourceVTKFromFile::GetDataSetName (size_t index) const
{
  return datasets[index].second;
}

vtkSmartPointer<vtkDataSet> ResultDataSourceVTKFromFile::ReadDataFromFile (const char* filename)
{
  vtkSmartPointer<vtkDataSetReader> reader = vtkSmartPointer<vtkDataSetReader>::New();
  reader->SetFileName (filename);
  
  // Read scalars
  reader->ReadAllScalarsOn();
  // Read vectors
  reader->ReadAllVectorsOn();
  // Read tensors
  reader->ReadAllTensorsOn();
  // Read normals
  reader->ReadAllNormalsOn();
  
  // Ignoring, for now: tcoords and fielddata - no idea what to do with them...
  
  reader->Update ();
  return reader->GetOutput();
}
