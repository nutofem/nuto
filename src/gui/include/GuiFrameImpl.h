#ifndef __NUTOGUI_GUIFRAMEIMPL_H__
#define __NUTOGUI_GUIFRAMEIMPL_H__

#include "GuiFrame.h"

#include <wx/wx.h>
#include <wx/aui/framemanager.h>
#include <boost/unordered_map.hpp>
#include <vector>

class wxAuiNotebook;
class wxAuiNotebookEvent;

#ifndef GUIFRAME_EXPORTED
  #define GUIFRAME_EXPORTED	NUTOGUI_IMPORT
#endif

namespace nutogui
{
  class GUIFRAME_EXPORTED GuiFrameImpl : public GuiFrame,
				         public wxFrame
  {
    wxAuiManager auiManager;
    wxAuiNotebook* notebook;
    wxStatusBar* grippy;
    wxLogWindow* logWindow;

    class TabCallbackImpl : public TabCallback
    {
      GuiFrameImpl* frame;

      wxString caption;
      wxWindow* page;
      bool closeable;

      int GetPageIndex () const;
      void UpdateCaption () const;
    public:
      TabCallbackImpl (GuiFrameImpl* frame)
       : frame (frame), page (nullptr), closeable (false) {}

      void SetCaption (const wxString& caption);
      const wxString& GetCaption () { return caption; }
      
      void SetCloseable (bool flag);
      bool GetCloseable () { return closeable; }
      
      void Activate ();

      void SetPage (wxWindow* page);
    };
    typedef boost::shared_ptr<TabCallbackImpl> TabCallbackImplPtr;
    struct TabInfo
    {
      TabPtr tab;
      TabCallbackImplPtr callback;

      TabInfo (const TabPtr& tab, const TabCallbackImplPtr& callback)
       : tab (tab), callback (callback) {}
    };
    std::vector<TabInfo> tabs;
    TabInfo currentTab;
    
    class PaneCallbackImpl : public PaneCallback
    {
      GuiFrameImpl* frame;

      wxString caption;
      wxWindow* contents;
      typedef std::pair<long, int> IDMapping;
      std::vector<IDMapping> paneAccels;

      void UpdateCaption () const;
    public:
      PaneCallbackImpl (GuiFrameImpl* frame) : frame (frame), contents (nullptr) {}

      void SetCaption (const wxString& caption);
      const wxString& GetCaption () { return caption; }
      
      void AddAccelerator (const wxAcceleratorEntry& accel);

      void SetContents (wxWindow* contents);
    };
    struct PaneInfo
    {
      PanePtr pane;
      boost::shared_ptr<PaneCallbackImpl> callback;

      PaneInfo (const PanePtr& pane, const boost::shared_ptr<PaneCallbackImpl>& callback)
       : pane (pane), callback (callback) {}
    };
    std::vector<PaneInfo> panes;
    
    typedef std::pair<wxWindow*, int> ForwardAccel;
    typedef boost::unordered_map<long, ForwardAccel> ForwardAccelMap;
    ForwardAccelMap forwardAccels;
    
    std::vector<wxAcceleratorEntry> allAccels;
    void UpdateAccelerators ();
    
    void ShowTabCloseButton (bool flag);
    
    void OnShowLog (wxCommandEvent& event);
    void OnClose (wxCommandEvent& event);
    void OnTabChanged (wxAuiNotebookEvent& event);
    void OnTabClose (wxAuiNotebookEvent& event);
    void OnTabClosed (wxAuiNotebookEvent& event);
    void OnSize (wxSizeEvent& event);
    
    void OnAnyCommand (wxCommandEvent& event);
  public:
    GuiFrameImpl (const wxChar* caption);
    ~GuiFrameImpl ();
    
    void SetLogWindow (wxLogWindow* logWindow) { this->logWindow = logWindow; }

    void AddTab (const TabPtr& tab, bool activate = false);
    void ActivateTab (Tab* tab);
    void AddPane (const wxString& name, DefaultLocation defaultLoc, const PanePtr& pane);
    
    DECLARE_EVENT_TABLE()
  };
} // namespace nutogui

#endif // __NUTOGUI_GUIFRAMEIMPL_H__
