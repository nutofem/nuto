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

ResultDataSourceVTKFromFile::ResultDataSourceVTKFromFile (const char* filename)
{
  ReadDataFromFile (filename);
}

vtkDataSet* ResultDataSourceVTKFromFile::QueryDataSet ()
{
  return dataset;
}

void ResultDataSourceVTKFromFile::ReadDataFromFile (const char* filename)
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
  dataset = reader->GetOutput();
}
