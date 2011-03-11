/**\file
 * Overlay classes/methods from 'nuto' module
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include <boost/python.hpp>

#include "NutoModuleOverlay.h"

#include "SwigConnect.h"

#include <wx/filename.h>
#include <wx/log.h>

#include <boost/make_shared.hpp>
#include <boost/python/raw_function.hpp>

#include "nuto/visualize/VisualizeComponentBase.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

using namespace boost::python;

NutoModuleOverlay& NutoModuleOverlay::GetOverlayInstance()
{
  // obtain NutoModuleOverlay instance
  object nuto_module = import ("nuto");
  object inst = nuto_module.attr("__nutoGuiOverlay__");
  
  return extract<NutoModuleOverlay&> (inst);
}

void NutoModuleOverlay::BindToPython ()
{
  class_<NutoModuleOverlay> ("_NutoModuleOverlay",
			     no_init)
    .def("ExportVtkDataFile", &NutoModuleOverlay::Overlay_ExportVtkDataFile)
    .staticmethod("ExportVtkDataFile")
  ;
}

void NutoModuleOverlay::Overlay_ExportVtkDataFile (const boost::python::object& self, const std::string& exportFile)
{
  wxString exportFileStr (exportFile.c_str(), wxConvLibc);
  {
    wxLogDebug (wxT ("intercepted ExportVtkDataFile() to %s"), exportFileStr.c_str());
  }
  // 'self' contains nuto.Structure
  
  NuTo::Structure* nutoStruct;
  if (!SwigConnect::SwigExtract (nutoStruct, self.ptr(), "NuTo::Structure *"))
  {
    // @@@ Proper error feedback
    std::cerr << "couldn't extract NuTo::Structure" << std::endl;
    return;
  }
  
  /* Look for identifier of instance on which ExportVtkDataFile was called.
     Will be used to name result display */
  wxString resultsTitle;
  {
    dict locals (handle<> (borrowed (PyEval_GetLocals())) );
    list keys (locals.keys ());
    for (boost::python::ssize_t i = 0; i < len (keys); i++)
    {
      object key = keys[i];
      if (locals[key] == self)
      {
	std::string keyStr = extract<std::string> (key);
	resultsTitle = wxString (keyStr.c_str(), wxConvLibc);
	break;
      }
    }
  }
  
  
  wxFileName newOutput;
  newOutput.AssignTempFileName (wxT ("nutogui"));
  if (!newOutput.IsOk())
  {
    // @@@ Better error feedback?
    wxLogError (wxT ("Couldn't generate temp file name"));
    return;
  }
  /* FIXME: Better would be:
   *  - Grab data directly from VTK, w/o going through a temp file
   *  - Build VTK DataSet in ResultViewer
   * 
   * It would already be better to pass out nutoStruct instead of a temp file,
   * however, that is not as easy as it may seem, mainly because the data extraction
   * has to be done on the script thread, right away - otherwise the object may be
   * destroyed...
   */
  std::string outputPath (newOutput.GetFullPath().fn_str());
  nutoStruct->ExportVtkDataFile (outputPath);
  
  wxString resultName;
  {
    wxFileName origExportName (exportFileStr);
    if (origExportName.GetExt() == wxT("vtk"))
      resultName = origExportName.GetName();
    else
      resultName = origExportName.GetFullName();
  }
  
  GetOverlayInstance().callback->Result (newOutput.GetFullPath(), resultsTitle, resultName);
}
