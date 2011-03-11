/**\file
 * Overlay classes/methods from 'nuto' module
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __NUTOMODULEOVERLAY_H__
#define __NUTOMODULEOVERLAY_H__

#include "ResultDataSourceVTK.h"

#include <wx/string.h>

namespace boost
{
  namespace python
  {
    class object;
  }
}

class NutoModuleOverlay
{
  static NutoModuleOverlay& GetOverlayInstance();
public:
  struct Callback
  {
    virtual ~Callback() {}
    
    virtual void Result (const nutogui::ResultDataSourceVTKPtr& result,
			 const wxString& title) = 0;
  };
  typedef boost::shared_ptr<Callback> CallbackPtr;
  
  NutoModuleOverlay (const CallbackPtr& callback) : callback (callback) {}
  
  static void BindToPython ();
  
  static void Overlay_ExportVtkDataFile (const boost::python::object& self, const std::string& exportFile);
protected:
  CallbackPtr callback;
};

#endif // __NUTOMODULEOVERLAY_H__
