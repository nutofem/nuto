#include "common.h"

#include "GuiFrameImpl.h"

#include <wx/aui/auibook.h>
#include <wx/aui/framemanager.h>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <boost/weak_ptr.hpp>

namespace nutogui
{
  enum
  {
    ID_NoteBook = 1,
    ID_ShowLog
  };
  
  BEGIN_EVENT_TABLE(GuiFrameImpl, wxFrame)
    EVT_MENU(ID_ShowLog, GuiFrameImpl::OnShowLog)
    EVT_MENU(wxID_CLOSE, GuiFrameImpl::OnClose)
    
    EVT_AUINOTEBOOK_PAGE_CHANGED(ID_NoteBook, GuiFrameImpl::OnTabChanged)
    EVT_AUINOTEBOOK_PAGE_CLOSE(ID_NoteBook, GuiFrameImpl::OnTabClose)
    EVT_AUINOTEBOOK_PAGE_CLOSED(ID_NoteBook, GuiFrameImpl::OnTabClosed)
    
    EVT_SIZE(GuiFrameImpl::OnSize)
    
    // Catch-all entry, for pane accelerator forwarding
    EVT_MENU(wxID_ANY, GuiFrameImpl::OnAnyCommand)
  END_EVENT_TABLE()

  GuiFrameImpl::GuiFrameImpl (const wxChar* caption)
   : wxFrame (nullptr, wxID_ANY, caption, wxDefaultPosition, wxSize (800, 600)),
     auiManager (this),
     logWindow (nullptr),
     currentTab (TabPtr(), TabCallbackImplPtr ())
  {
    notebook = new wxAuiNotebook (this, ID_NoteBook, wxDefaultPosition, wxDefaultSize,
				  wxAUI_NB_SCROLL_BUTTONS | wxAUI_NB_TOP
				  | wxBORDER_NONE);

    auiManager.AddPane (notebook,
			wxAuiPaneInfo().CenterPane().PaneBorder(false));
    
    /* Vision for menu: a single button, on the left side of the tab bar.
       But for that would need changes to wxAuiNotebook to support offsetting
       the tabs ... */
    wxMenuBar *menuBar = new wxMenuBar;
    {
      wxMenu* menuWindow = new wxMenu;
      menuWindow->Append (ID_ShowLog, wxT("Show &log"));
      menuWindow->Append (wxID_CLOSE, wxT("&Close"));
      menuBar->Append (menuWindow, wxT("&Window"));
    }
    SetMenuBar (menuBar);
    
    /* Create a status bar. It's not really intended to show any status...
       it's just there for the resize grip. */
    grippy = new wxStatusBar (this);
    grippy->SetFieldsCount (1);
    static const int styles[] = { wxSB_FLAT };
    grippy->SetStatusStyles (1, styles);
    wxSize sbSize (grippy->GetSize ());
    sbSize.SetWidth (sbSize.GetHeight());
    grippy->SetSize (sbSize);
    grippy->SetBackgroundColour (auiManager.GetArtProvider()->GetColour (wxAUI_DOCKART_BACKGROUND_COLOUR));
  }

  GuiFrameImpl::~GuiFrameImpl ()
  {
    auiManager.UnInit();
  }

  void GuiFrameImpl::UpdateAccelerators ()
  {
    wxAcceleratorTable accelTable (allAccels.size(), allAccels.data());
    SetAcceleratorTable (accelTable);
  }
  
  void GuiFrameImpl::ShowTabCloseButton (bool flag)
  {
    long notebookStyle = notebook->GetWindowStyle();
    if (flag)
      notebookStyle |= wxAUI_NB_CLOSE_ON_ACTIVE_TAB;
    else
      notebookStyle &= ~wxAUI_NB_CLOSE_ON_ACTIVE_TAB;
    notebook->SetWindowStyle (notebookStyle);
  }

  void GuiFrameImpl::OnShowLog (wxCommandEvent& event)
  {
    if (logWindow)
      logWindow->Show (true);
  }
  
  void GuiFrameImpl::OnClose (wxCommandEvent& event)
  {
    Close (true);
  }

  void GuiFrameImpl::OnTabChanged (wxAuiNotebookEvent& event)
  {
    if (currentTab.tab)
      currentTab.tab->TabDeactivate();
    
    currentTab = tabs[event.GetSelection()]; // @@@ Extract from Tab window or so?
    currentTab.tab->TabActivate();
    ShowTabCloseButton (currentTab.callback->GetCloseable ());
  }

  void GuiFrameImpl::OnTabClose (wxAuiNotebookEvent& event)
  {
    // Delete tab info in 'close' event so the Tab is destroyed before it's contents
    tabs.erase (tabs.begin() + event.GetSelection());
    currentTab = TabInfo (TabPtr(), TabCallbackImplPtr ());
  }
  
  void GuiFrameImpl::OnTabClosed (wxAuiNotebookEvent& event)
  {
    currentTab = tabs[notebook->GetSelection()];
  }

  void GuiFrameImpl::OnSize (wxSizeEvent& event)
  {
    wxPoint grippyPos;
    grippyPos.x = GetClientSize().GetWidth() - grippy->GetSize().GetWidth();
    grippyPos.y = GetClientSize().GetHeight() - grippy->GetSize().GetHeight();
    grippy->SetPosition (grippyPos);
  }
  
  void GuiFrameImpl::OnAnyCommand (wxCommandEvent& event)
  {
    ForwardAccelMap::const_iterator fwdAccel = forwardAccels.find (event.GetId());
    if (fwdAccel == forwardAccels.end())
    {
      event.Skip();
    }
    else
    {
      // Change ID and forward event
      event.SetId (fwdAccel->second.second);
      fwdAccel->second.first->GetEventHandler()->ProcessEvent (event);
    }
  }

  void GuiFrameImpl::AddTab (const TabPtr& tab, bool activate)
  {
    boost::shared_ptr<TabCallbackImpl> tabCallback (boost::make_shared<TabCallbackImpl> (this));
    tabs.push_back (TabInfo (tab, tabCallback));
    
    wxWindow* newPage = tab->CreateContents (tabCallback, this);
    notebook->AddPage (newPage, tabCallback->GetCaption(), activate);
    tabCallback->SetPage (newPage);
  }

  void GuiFrameImpl::ActivateTab (Tab* tab)
  {
    for (size_t t = 0; t < tabs.size(); t++)
    {
      if (tabs[t].tab.get() == tab)
      {
	notebook->SetSelection (t);
	return;
      }
    }
  }

  void GuiFrameImpl::AddPane (const wxString& name,
			      DefaultLocation defaultLoc,
			      const PanePtr& pane)
  {
    boost::shared_ptr<PaneCallbackImpl> paneCallback (boost::make_shared<PaneCallbackImpl> (this));
    panes.push_back (PaneInfo (pane, paneCallback));
    
    wxWindow* newContents = pane->CreateContents (paneCallback, this);
    
    wxAuiPaneInfo paneInfo;
    paneInfo.Name (name);
    paneInfo.Caption (paneCallback->GetCaption());
    if (pane->IsToolbar()) paneInfo.ToolbarPane ();
    switch (defaultLoc)
    {
    case Left:
      paneInfo.Left();
      break;
    case Top:
      paneInfo.Top();
      break;
    case Right:
      paneInfo.Right();
      break;
    case Bottom:
      paneInfo.Bottom();
      break;
    }
    paneInfo.CloseButton (false);
    
    auiManager.AddPane (newContents, paneInfo);
    auiManager.Update ();
    
    paneCallback->SetContents (newContents);
  }

  //-------------------------------------------------------------------------

  int GuiFrameImpl::TabCallbackImpl::GetPageIndex () const
  {
    return (this->page != 0) ? frame->notebook->GetPageIndex (this->page) : wxNOT_FOUND;
  }

  void GuiFrameImpl::TabCallbackImpl::UpdateCaption () const
  {
    int pageIdx = GetPageIndex ();
    if (pageIdx != wxNOT_FOUND)
    {
      frame->notebook->SetPageText (pageIdx, caption);
    }
  }

  void GuiFrameImpl::TabCallbackImpl::SetCaption (const wxString& caption)
  {
    this->caption = caption;
    UpdateCaption ();
  }

  void GuiFrameImpl::TabCallbackImpl::SetCloseable (bool flag)
  {
    closeable = flag;
    if (frame->currentTab.callback.get() == this)
      frame->ShowTabCloseButton (closeable);
  }

  void GuiFrameImpl::TabCallbackImpl::Activate ()
  {
    int pageIdx = GetPageIndex ();
    if (pageIdx != wxNOT_FOUND)
    {
      frame->notebook->SetSelection (pageIdx);
    }
  }

  void GuiFrameImpl::TabCallbackImpl::SetPage (wxWindow* page)
  {
    this->page = page;
    UpdateCaption ();
  }

  //-------------------------------------------------------------------------

  void GuiFrameImpl::PaneCallbackImpl::UpdateCaption () const
  {
    if (contents)
    {
      wxAuiPaneInfo& paneInfo = frame->auiManager.GetPane (contents);
      paneInfo.Caption (caption);
    }
  }

  void GuiFrameImpl::PaneCallbackImpl::SetCaption (const wxString& caption)
  {
    this->caption = caption;
    UpdateCaption ();
  }
  
  void GuiFrameImpl::PaneCallbackImpl::AddAccelerator (const wxAcceleratorEntry& accel)
  {
    long newID = wxNewId();
    wxAcceleratorEntry newAccel (accel.GetFlags(), accel.GetKeyCode(), newID);
    paneAccels.push_back (std::make_pair (newID, accel.GetCommand()));
    
    if (contents)
      frame->forwardAccels[newID] = std::make_pair (contents, accel.GetCommand());
    
    frame->allAccels.push_back (newAccel);
    frame->UpdateAccelerators();
  }

  void GuiFrameImpl::PaneCallbackImpl::SetContents (wxWindow* contents)
  {
    if (this->contents)
    {
      // Clear out old accelerators
      for (ForwardAccelMap::iterator fwdAccel = frame->forwardAccels.begin();
	   fwdAccel != frame->forwardAccels.end();
	  )
      {
	if (fwdAccel->second.first == this->contents)
	  frame->forwardAccels.erase (fwdAccel);
	else
	  ++fwdAccel;
      }
    }
    
    this->contents = contents;
    
    if (this->contents && (paneAccels.size() > 0))
    {
      // Add accelerators
      BOOST_FOREACH(const IDMapping& map, paneAccels)
      {
	frame->forwardAccels[map.first] = std::make_pair (contents, map.second);
      }
      frame->UpdateAccelerators();
    }
    
    UpdateCaption ();
  }
} // namespace nutogui
