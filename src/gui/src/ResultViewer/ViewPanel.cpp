/**\file
 * Result viewer view panel.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"

#include "ViewPanel.h"

#include "SplitManager.h"
#include "View3D.h"

#include <wx/artprov.h>
#include <wx/aui/auibar.h>
#include <boost/make_shared.hpp>

namespace nutogui
{
  struct ResultViewerImpl::ViewPanel::SharedViewData
  {
    wxSize smallButtonSize;
    wxBitmap imgSplitH;
    wxBitmap imgSplitV;
    wxBitmap imgUnsplit;
    wxBitmap imgMaximize;
    wxBitmap imgUnmaximize;
  };
  
  enum
  {
    ID_SplitHorz = 12345678, // Arbitrary, large ID to avoid conflicts with child IDs
    ID_SplitVert,
    ID_Unsplit,
    ID_ToggleMaximization
  };
  
  BEGIN_EVENT_TABLE(ResultViewerImpl::ViewPanel, wxPanel)
    EVT_MENU(ID_SplitHorz, ResultViewerImpl::ViewPanel::OnSplitHorizontally)
    EVT_UPDATE_UI(ID_SplitHorz, ResultViewerImpl::ViewPanel::OnSplitUpdateUI)
    EVT_MENU(ID_SplitVert, ResultViewerImpl::ViewPanel::OnSplitVertically)
    EVT_UPDATE_UI(ID_SplitVert, ResultViewerImpl::ViewPanel::OnSplitUpdateUI)
    EVT_MENU(ID_Unsplit, ResultViewerImpl::ViewPanel::OnUnsplit)
    EVT_UPDATE_UI(ID_Unsplit, ResultViewerImpl::ViewPanel::OnUnsplitUpdateUI)
    EVT_MENU(ID_ToggleMaximization, ResultViewerImpl::ViewPanel::OnToggleMaximization)
    EVT_UPDATE_UI(ID_ToggleMaximization, ResultViewerImpl::ViewPanel::OnToggleMaximizationUpdateUI)
  END_EVENT_TABLE()
  
  ResultViewerImpl::ViewPanel::ViewPanel (wxWindow* parent, SplitManager* splitMgr)
   : wxPanel (parent), splitMgr (splitMgr)
  {
    SetupSharedData ();
    childPanel = new View3D (this);
    CreateChildren ();
  }

  ResultViewerImpl::ViewPanel::ViewPanel (wxWindow* parent, ViewPanel* cloneFrom)
   : wxPanel (parent), splitMgr (cloneFrom->splitMgr), sharedData (cloneFrom->sharedData)
  {
    childPanel = cloneFrom->childPanel->Clone (this);
    CreateChildren ();
  }

  void ResultViewerImpl::ViewPanel::SetData (const DataConstPtr& data)
  {
    childPanel->SetData (data);
  }

  class ResultViewerImpl::ViewPanel::PostToOthersViewsTraverser : public SplitManager::ViewsTraverser
  {
    ViewPanel* skipView;
    wxEvent& event;
  public:
    PostToOthersViewsTraverser (ViewPanel* skipView, wxEvent& event)
      : skipView (skipView), event (event) {}
    
    bool operator() (wxWindow* view)
    {
      if (view != skipView)
      {
	Content* contentChild = static_cast<ViewPanel*> (view)->childPanel;
	contentChild->GetEventHandler()->ProcessEvent (event);
      }
      return true;
    }
  };

  void ResultViewerImpl::ViewPanel::PostToOthers (wxEvent& event)
  {
    PostToOthersViewsTraverser ptovt (this, event);
    splitMgr->Traverse (ptovt);
  }

  void ResultViewerImpl::ViewPanel::SetupSharedData ()
  {
    sharedData = boost::make_shared<SharedViewData> ();
    
    sharedData->smallButtonSize = wxSize (11, 11);
    sharedData->imgSplitH = wxArtProvider::GetBitmap (wxART_MAKE_ART_ID(split-horz), wxART_TOOLBAR,
						      sharedData->smallButtonSize);
    sharedData->imgSplitV = wxArtProvider::GetBitmap (wxART_MAKE_ART_ID(split-vert), wxART_TOOLBAR,
						      sharedData->smallButtonSize);
    sharedData->imgUnsplit = wxArtProvider::GetBitmap (wxART_CLOSE, wxART_TOOLBAR,
						       sharedData->smallButtonSize);
    sharedData->imgMaximize = wxArtProvider::GetBitmap (wxART_MAKE_ART_ID(maximize), wxART_TOOLBAR,
							sharedData->smallButtonSize);
    sharedData->imgUnmaximize = wxArtProvider::GetBitmap (wxART_MAKE_ART_ID(unmaximize), wxART_TOOLBAR,
							  sharedData->smallButtonSize);
  }

  void ResultViewerImpl::ViewPanel::CreateChildren ()
  {
    wxSizer* sizer = new wxBoxSizer (wxVERTICAL);
    
    topBarSizer = new wxBoxSizer (wxHORIZONTAL);
    
    wxAuiToolBar* splitButtonsBar = new wxAuiToolBar (this, wxID_ANY,
						      wxDefaultPosition, wxDefaultSize,
						      wxAUI_TB_VERTICAL | wxAUI_TB_NO_AUTORESIZE);
    splitButtonsBar->SetToolBitmapSize (sharedData->smallButtonSize);
    splitButtonsBar->SetToolBorderPadding (2);
    splitButtonsBar->SetMargins (1, -1);
    splitButtonsBar->AddTool (ID_SplitHorz, wxEmptyString,
			      sharedData->imgSplitH,
			      wxT ("Split left/right"));
    splitButtonsBar->AddTool (ID_SplitVert, wxEmptyString,
			      sharedData->imgSplitV,
			      wxT ("Split top/bottom"));
    splitButtonsBar->Realize ();
    topBarSizer->Add (splitButtonsBar, 0, wxEXPAND);
    
    wxWindow* topBarContentTools = childPanel->CreateTopTools (this);
    if (topBarContentTools)
      topBarSizer->Add (topBarContentTools, wxSizerFlags (1).Expand());
    else
      topBarSizer->AddStretchSpacer ();
    
    closeMaxButtonsBar = new wxAuiToolBar (this, wxID_ANY,
					   wxDefaultPosition, wxDefaultSize,
					   wxAUI_TB_VERTICAL | wxAUI_TB_NO_AUTORESIZE);
    closeMaxButtonsBar->SetToolBitmapSize (sharedData->smallButtonSize);
    closeMaxButtonsBar->SetToolBorderPadding (2);
    closeMaxButtonsBar->SetMargins (1, -1);
    closeMaxButtonsBar->AddTool (ID_Unsplit, wxEmptyString,
				 sharedData->imgUnsplit,
				 wxT ("Close this view"));
    closeMaxButtonsBar->AddTool (ID_ToggleMaximization, wxEmptyString,
				 sharedData->imgMaximize,
				 wxT ("Maximize this view"));
    closeMaxButtonsBar->Realize ();
    topBarSizer->Add (closeMaxButtonsBar, 0, wxEXPAND);
    
    sizer->Add (topBarSizer, wxSizerFlags (0).Expand());
    
    sizer->Add (childPanel, wxSizerFlags (1).Expand());
    
    SetSizer (sizer);
  }

  void ResultViewerImpl::ViewPanel::OnSplitHorizontally (wxCommandEvent& event)
  {
    ViewPanel* newView = new ViewPanel (splitMgr, this);
    splitMgr->SplitHorizontally (this, newView);
  }
  
  void ResultViewerImpl::ViewPanel::OnSplitVertically (wxCommandEvent& event)
  {
    ViewPanel* newView = new ViewPanel (splitMgr, this);
    splitMgr->SplitVertically (this, newView);
  }

  void ResultViewerImpl::ViewPanel::OnSplitUpdateUI (wxUpdateUIEvent& event)
  {
    event.Enable (splitMgr->CanSplit (this));
  }
  
  void ResultViewerImpl::ViewPanel::OnUnsplit (wxCommandEvent& event)
  {
    splitMgr->Unsplit (this);
  }
  
  void ResultViewerImpl::ViewPanel::OnUnsplitUpdateUI (wxUpdateUIEvent& event)
  {
    event.Enable (splitMgr->CanUnsplit (this));
  }

  void ResultViewerImpl::ViewPanel::OnToggleMaximization (wxCommandEvent& event)
  {
    bool isMaximized = splitMgr->ToggleMaximization (this);
    if (isMaximized)
    {
      closeMaxButtonsBar->SetToolBitmap (ID_ToggleMaximization,
					 sharedData->imgUnmaximize);
      closeMaxButtonsBar->SetToolShortHelp (ID_ToggleMaximization,
					    wxT ("Reduce this view to normal size"));
    }
    else
    {
      closeMaxButtonsBar->SetToolBitmap (ID_ToggleMaximization,
					 sharedData->imgMaximize);
      closeMaxButtonsBar->SetToolShortHelp (ID_ToggleMaximization,
					    wxT ("Maximize this view"));
    }
  }

  void ResultViewerImpl::ViewPanel::OnToggleMaximizationUpdateUI (wxUpdateUIEvent& event)
  {
    event.Enable (splitMgr->CanToggleMaximization (this));
  }

} // namespace nutogui
