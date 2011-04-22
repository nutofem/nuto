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

#include "AllViewSharedData.h"
#include "SplitManager.h"
#include "View3D.h"
#include "ViewPlot.h"

#include <wx/artprov.h>
#include <wx/aui/auibar.h>
#include <boost/make_shared.hpp>

namespace nutogui
{
  struct ResultViewerImpl::ViewPanel::SharedViewData : public SharedViewDataBase
  {
    wxSize smallButtonSize;
    wxBitmap imgSplitH;
    wxBitmap imgSplitV;
    wxBitmap imgUnsplit;
    wxBitmap imgMaximize;
    wxBitmap imgUnmaximize;
    
    wxBitmap imgContentViewButtonImages[numContentTypes];
    wxMenu contentViewMenu;
  };
  
  enum
  {
    ID_SplitHorz = 12345, // Arbitrary, large ID to avoid conflicts with child IDs
    ID_SplitVert,
    ID_Unsplit,
    ID_ToggleMaximization,
    
    ID_ContentViewPopup,
    ID_ContentViewFirst
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
    
    EVT_AUITOOLBAR_TOOL_DROPDOWN(ID_ContentViewPopup, ResultViewerImpl::ViewPanel::OnContentViewPopup)
    EVT_MENU_RANGE(ID_ContentViewFirst,
		   ID_ContentViewFirst + ResultViewerImpl::ViewPanel::numContentTypes - 1,
		   ResultViewerImpl::ViewPanel::OnContentViewChange)
  END_EVENT_TABLE()
  
  ResultViewerImpl::ViewPanel::ViewPanel (wxWindow* parent, SplitManager* splitMgr,
					  const DataConstPtr& data)
   : wxPanel (parent), splitMgr (splitMgr)
  {
    SetupSharedData ();
    
    AllViewSharedDataPtr allSharedData = AllViewSharedData::Setup (this);
    allSharedData->data = data;
    
    contentType = content3D;
    childPanel = CreateChild (contentType);
    CreateChildren ();
  }

  ResultViewerImpl::ViewPanel::ViewPanel (wxWindow* parent, ViewPanel* cloneFrom)
   : wxPanel (parent), splitMgr (cloneFrom->splitMgr), 
     sharedDataBase (cloneFrom->sharedDataBase)
  {
    contentType = cloneFrom->contentType;
    childPanel = CreateChild (contentType);
    CreateChildren ();
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

// Hack for wx 2.8
#ifndef wxART_CLOSE
#define wxART_CLOSE	wxART_MAKE_ART_ID(wxART_CLOSE)
#endif

  void ResultViewerImpl::ViewPanel::SetupSharedData ()
  {
    boost::shared_ptr<SharedViewData> sharedData (boost::make_shared<SharedViewData> ());
    
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

    static const wxChar* const contentTypeNames[numContentTypes] =
    {
      wxT ("&3D"),
      wxT ("&Plot"),
    };
    static const wxChar* const contentTypeArtNames[numContentTypes] =
    {
      wxART_MAKE_ART_ID(view-3d),
      wxART_MAKE_ART_ID(view-plot)
    };
    for (size_t i = 0; i < numContentTypes; i++)
    {
      wxMenuItem* newItem = new wxMenuItem (&sharedData->contentViewMenu,
					    ID_ContentViewFirst + i,
					    contentTypeNames[i],
					    wxEmptyString,
					    wxITEM_NORMAL);
      newItem->SetBitmap (wxArtProvider::GetBitmap (contentTypeArtNames[i], wxART_MENU));
      sharedData->contentViewMenu.Append (newItem);
      
      sharedData->imgContentViewButtonImages[i] = wxArtProvider::GetBitmap (contentTypeArtNames[i],
									    wxART_TOOLBAR);
    }
    
    this->sharedDataBase = sharedData;
  }

  ResultViewerImpl::ViewPanel::Content* ResultViewerImpl::ViewPanel::CreateChild (ContentType contentType)
  {
    Content* newChild = nullptr;
    switch (contentType)
    {
    case content3D:
      newChild = new View3D (this);
      break;
    case contentPlot:
      newChild = new ViewPlot (this);
      break;
    case numContentTypes:
      break;
    }
    assert (newChild);
    return newChild;
  }

  void ResultViewerImpl::ViewPanel::CreateChildren ()
  {
    boost::shared_ptr<SharedViewData> sharedData (
      boost::shared_static_cast<SharedViewData> (this->sharedDataBase));
    
    contentsSizer = new wxBoxSizer (wxVERTICAL);
    
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
    
    contentViewBar = new wxAuiToolBar (this, wxID_ANY,
				       wxDefaultPosition, wxDefaultSize,
				       wxAUI_TB_HORZ_LAYOUT | wxAUI_TB_NO_AUTORESIZE);
    contentViewBar->AddTool (ID_ContentViewPopup, wxEmptyString,
			     sharedData->imgContentViewButtonImages[contentType],
			     wxT ("Select how to display the data"));
    contentViewBar->SetToolDropDown (ID_ContentViewPopup, true);
    contentViewBar->Realize ();
    topBarSizer->Add (contentViewBar, wxSizerFlags (0).Expand());
    
    topBarContentTools = childPanel->CreateTopTools (this);
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
    
    contentsSizer->Add (topBarSizer, wxSizerFlags (0).Expand());
    
    contentsSizer->Add (childPanel, wxSizerFlags (1).Expand());
    
    SetSizer (contentsSizer);
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
    boost::shared_ptr<SharedViewData> sharedData (
      boost::shared_static_cast<SharedViewData> (this->sharedDataBase));
    
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

  void ResultViewerImpl::ViewPanel::OnContentViewPopup (wxAuiToolBarEvent& event)
  {
    boost::shared_ptr<SharedViewData> sharedData (
      boost::shared_static_cast<SharedViewData> (this->sharedDataBase));
    
    wxPoint popupPos (contentViewBar->GetToolRect (ID_ContentViewPopup).GetBottomLeft());
    popupPos = contentViewBar->ClientToScreen (popupPos);
    popupPos = ScreenToClient (popupPos);
    
    PopupMenu (&sharedData->contentViewMenu, popupPos);
  }

  void ResultViewerImpl::ViewPanel::OnContentViewChange (wxCommandEvent& event)
  {
    ContentType newContentType = ContentType (event.GetId() - ID_ContentViewFirst);
    if (contentType == newContentType) return;
    
    Content* newChild = CreateChild (newContentType);
    if (!newChild) return;
    contentType = newContentType;
    
    contentsSizer->Replace (childPanel, newChild);
    wxWindow* newTopBarContentTools = newChild->CreateTopTools (this);
    
    wxSizerItem* topBarContentSizerItem;
    if (newTopBarContentTools)
      topBarContentSizerItem = new wxSizerItem (newTopBarContentTools, wxSizerFlags (1).Expand());
    else
      // Stretch spacer
      topBarContentSizerItem = new wxSizerItem (0, 0, 1, 0, 0, nullptr);
    // @@@ Hardcoded sizer item index
    topBarSizer->Replace (2, topBarContentSizerItem);
    
    childPanel->DestroyTopTools (topBarContentTools);
    delete topBarContentTools;
    topBarContentTools = newTopBarContentTools;
    
    delete childPanel;
    childPanel = newChild;
    
    boost::shared_ptr<SharedViewData> sharedData (
      boost::shared_static_cast<SharedViewData> (this->sharedDataBase));
    contentViewBar->SetToolBitmap (ID_ContentViewPopup,
				   sharedData->imgContentViewButtonImages[newContentType]);
    
    Layout();
  }
} // namespace nutogui
