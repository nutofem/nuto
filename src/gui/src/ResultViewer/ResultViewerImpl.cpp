/**\file
 * Result viewer tab.
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"
#include "ResultViewerImpl.h"

#include "Data.h"
#include "SplitManager.h"
#include "View.h"

#include <wx/wx.h>
#include <boost/make_shared.hpp>
#include <boost/weak_ptr.hpp>

#include <vtkDataSetReader.h>

namespace nutogui
{
  ResultViewerImpl::ResultViewerImpl (const wxString& filename,
				      const wxString& caption,
				      bool deleteFileWhenDone)
   : dataFile (filename), deleteFileWhenDone (deleteFileWhenDone), caption (caption),
     splitMgr (nullptr)
  {}
  
  ResultViewerImpl::~ResultViewerImpl()
  {
  }

  wxWindow* ResultViewerImpl::CreateContents (const GuiFrame::TabCallbackWeakPtr& callback_,
					      wxWindow* parentWindow)
  {
    nutogui::GuiFrame::TabCallbackPtr callback (callback_);
    callback->SetCaption (caption);
    callback->SetCloseable (true);
    
    ReadDataFromFile (dataFile.fn_str());
    if (deleteFileWhenDone)
    {
      wxRemoveFile (dataFile);
    }
    
    splitMgr = new SplitManager (parentWindow);
    View* firstView = new View (splitMgr, splitMgr);
    firstView->SetData (data);
    splitMgr->SetInitialView (firstView);
    return splitMgr;
  }

  /* TODO: The result reading stuff should probably be abstracted into some interface
     or so - one day, we probably want to pull it the stuff directly from NuTo
     without going through a VTK data reader. */
  void ResultViewerImpl::ReadDataFromFile (const char* filename)
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
    vtkDataSet* dataset = reader->GetOutput();
    
    data = boost::make_shared<Data> (dataset);
  }

} // namespace nutogui
