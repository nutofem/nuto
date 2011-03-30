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

#include <boost/unordered_map.hpp>

class wxAuiToolBar;
class wxAuiToolBarEvent;

namespace nutogui
{
  class ResultViewerImpl::ViewPanel : public wxPanel
  {
  public:
    class Content;
    
    ViewPanel (wxWindow* parent, SplitManager* splitMgr,
	       const DataConstPtr& data);
    
    /// Post an event to content children of all other view panels
    void PostToOthers (wxEvent& event);
    
    /**\name View panel shared data management
     * @{ */
    /// Base class for individual view panel shared data structures
    struct SharedDataBase {};
    
    /// Query shared data.
    template<typename T>
    boost::shared_ptr<T> GetSharedData () const
    {
      const void* tag = GetSharedDataTag<T> ();
      PanelsSharedDataMap::const_iterator panelSharedData = sharedDataBase->panelsSharedData.find (tag);
      if (panelSharedData != sharedDataBase->panelsSharedData.end())
	return boost::shared_static_cast<T> (panelSharedData->second);
      else
	return boost::shared_ptr<T> ();
    }
    
    /// Set shared data.
    template<typename T>
    void SetSharedData (const boost::shared_ptr<T>& panelSharedData)
    {
      const void* tag = GetSharedDataTag<T> ();
      sharedDataBase->panelsSharedData[tag] = panelSharedData;
    }
    
    /** @} */
    
    DECLARE_EVENT_TABLE()
  protected:
    ViewPanel (wxWindow* parent, ViewPanel* cloneFrom);
    
    class PostToOthersViewsTraverser;
    
    SplitManager* splitMgr;

    /**\name Generic shared data
     * @{ */
    /// Get a unique "ID" for a shared data type
    template<typename T>
    static const void* GetSharedDataTag ()
    {
      static char tag; // type doesn't matter
      return &tag;
    }
    
    typedef boost::unordered_map<const void*, boost::shared_ptr<SharedDataBase> > PanelsSharedDataMap;
    struct SharedViewDataBase
    {
      /// Shared data structures
      PanelsSharedDataMap panelsSharedData;
    };
    /** @} */
    
    struct SharedViewData;
    boost::shared_ptr<SharedViewDataBase> sharedDataBase;
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
    Content* CreateChild (ContentType contentType);
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
