/**\file
 * Result viewer tab.
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __NUTOGUI_RESULTVIEWERIMPL_H__
#define __NUTOGUI_RESULTVIEWERIMPL_H__

#include "GuiFrame.h"
#include "ResultDataSourceVTK.h"
#include "TabCommonImpl.h"

#include <wx/wx.h>

#ifndef RESULTVIEWER_EXPORTED
  #define RESULTVIEWER_EXPORTED	NUTOGUI_IMPORT
#endif

namespace nutogui
{
  class RESULTVIEWER_EXPORTED ResultViewerImpl : public TabCommonImpl
  {
    class View3D;
    
    class Data;
    typedef boost::shared_ptr<const Data> DataConstPtr;
    DataConstPtr data;
    
    const wxString caption;
    
    void ReadDataFromFile (const char* filename);
    
    class SplitManager;
    SplitManager* splitMgr;
  public:
    ResultViewerImpl (const ResultDataSourceVTKPtr& resultData,
		      const wxString& caption);
    ~ResultViewerImpl ();
    
    wxWindow* CreateContents (const GuiFrame::TabCallbackWeakPtr& callback,
			      wxWindow* parentWindow);
  };
} // namespace nutogui

#endif // __NUTOGUI_RESULTVIEWERIMPL_H__
