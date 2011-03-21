/**\file
 * Result viewer view panel.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __RESULTVIEWERIMPL_VIEWPANEL_H__
#define __RESULTVIEWERIMPL_VIEWPANEL_H__

#include "ResultViewerImpl.h"

class wxAuiToolBar;
class wxAuiToolBarEvent;

namespace nutogui
{
  class ResultViewerImpl::ViewPanel : public wxPanel
  {
  public:
    class Content;
    
    ViewPanel (wxWindow* parent, SplitManager* splitMgr);
    
    void SetData (const DataConstPtr& data);
    
    /// Post an event to content children of all other view panels
    void PostToOthers (wxEvent& event);
    
    DECLARE_EVENT_TABLE()
  protected:
    ViewPanel (wxWindow* parent, ViewPanel* cloneFrom);
    
    class PostToOthersViewsTraverser;
    
    SplitManager* splitMgr;
    DataConstPtr data;

    struct SharedViewData;
    boost::shared_ptr<SharedViewData> sharedData;
    void SetupSharedData ();
    
    wxSizer* contentsSizer;
    wxSizer* topBarSizer;
    wxAuiToolBar* contentViewBar;
    wxWindow* topBarContentTools;
    wxAuiToolBar* closeMaxButtonsBar;
    enum ContentType
    {
      content3D = 0,
      contentPlot,
      
      numContentTypes
    };
    ContentType contentType;
    Content* childPanel;
    void CreateChildren ();
    
    void OnSplitHorizontally (wxCommandEvent& event);
    void OnSplitVertically (wxCommandEvent& event);
    void OnSplitUpdateUI (wxUpdateUIEvent& event);
    void OnUnsplit (wxCommandEvent& event);
    void OnUnsplitUpdateUI (wxUpdateUIEvent& event);
    void OnToggleMaximization (wxCommandEvent& event);
    void OnToggleMaximizationUpdateUI (wxUpdateUIEvent& event);
    
    void OnContentViewPopup (wxAuiToolBarEvent& event);
    void OnContentViewChange (wxCommandEvent& event);
  };
} // namespace nutogui

#endif // __RESULTVIEWERIMPL_VIEWPANEL_H__
