/**\file
 * Tab to display output from a script
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __SCRIPTOUTPUTTAB_H__
#define __SCRIPTOUTPUTTAB_H__

#include "TabCommonImpl.h"

#include <wx/wx.h>

class ScriptOutputTab : public nutogui::TabCommonImpl
{
  wxTextCtrl* textCtrl;
  
  enum
  {
    lineError = 1 << 0,
    lineNote = 1 << 1,
    
    // Mask of all styles that change the font
    fontStyles = lineNote
  };
  unsigned int lastStyle;
  void AppendOutput (unsigned int style, const wxString& output);
public:
  ScriptOutputTab ();
  
  wxWindow* CreateContents (const nutogui::GuiFrame::TabCallbackWeakPtr& callback,
			    wxWindow* parentWindow);

  void AppendOutput (const wxString& output);
  void AppendErrorOutput (const wxString& output);
  void AppendNote (const wxString& output);
};

#endif // __SCRIPTOUTPUTTAB_H__
