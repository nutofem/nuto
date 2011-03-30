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
#include "ViewPanel.h"

#include <wx/wx.h>
#include <boost/make_shared.hpp>
#include <boost/weak_ptr.hpp>

namespace nutogui
{
  ResultViewerImpl::ResultViewerImpl (const ResultDataSourceVTKPtr& resultData,
				      const wxString& caption)
   : caption (caption),
     splitMgr (nullptr)
  {
    data = boost::make_shared<Data> (resultData);
  }
  
  ResultViewerImpl::~ResultViewerImpl()
  {
  }

  wxWindow* ResultViewerImpl::CreateContents (const GuiFrame::TabCallbackWeakPtr& callback_,
					      wxWindow* parentWindow)
  {
    nutogui::GuiFrame::TabCallbackPtr callback (callback_);
    callback->SetCaption (caption);
    callback->SetCloseable (true);
    
    splitMgr = new SplitManager (parentWindow);
    ViewPanel* firstView = new ViewPanel (splitMgr, splitMgr, data);
    splitMgr->SetInitialView (firstView);
    return splitMgr;
  }

} // namespace nutogui
