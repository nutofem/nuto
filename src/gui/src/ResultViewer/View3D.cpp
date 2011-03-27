/**\file
 * Result viewer 3D view panel.
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"
#include "View3D.h"

#include "Data.h"
#include "View3D/DataSetEdgeExtractor.h"
#include "View3D/DataSetEdgeIterator.h"
#include "View3D/DataSetFaceExtractor.h"
#include "View3D/DataSetHelpers.h"
#include "View3D/DisplacementDirectionSizePanel.h"
#include "View3D/Gradients.h"
#include "View3D/UnstructuredGridWithQuadsClippingToPoly.h"

#include "uicommon/Quote.h"

#include "vtkwx.h"
#include <wx/artprov.h>
#include <wx/aui/auibar.h>
#include <boost/make_shared.hpp>
#include <boost/weak_ptr.hpp>

#include <vtkArrowSource.h>
#include <vtkAxesActor.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellPicker.h>
#include <vtkClipDataSet.h>
#include <vtkCommand.h>
#include <vtkDataSetMapper.h>
#include <vtkFloatArray.h>
#include <vtkGlyph3D.h>
#include <vtkInteractorStyleSwitch.h>
#include <vtkLookupTable.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkPlane.h>
#include <vtkPlaneWidget.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkScalarBarActor.h>
#include <vtkScalarBarWidget.h>

namespace nutogui
{
  
  class ResultViewerImpl::View3D::CameraModifiedCallback : public vtkCommand
  {
  public:
    View3D* view;
    
    static CameraModifiedCallback* New ()
    {
      return new CameraModifiedCallback ();
    }
    
    void Execute (vtkObject* caller, unsigned long eventId, void* callData)
    {
      vtkCamera* cam = static_cast<vtkCamera*> (caller);
      view->CameraChanged (cam);
    }
  };
  
  //-------------------------------------------------------------------
  
  class ResultViewerImpl::View3D::ClipPlaneChangedCallback : public vtkCommand
  {
  public:
    View3D* view;
    
    static ClipPlaneChangedCallback* New ()
    {
      return new ClipPlaneChangedCallback ();
    }
    
    void Execute (vtkObject* caller, unsigned long eventId, void* callData)
    {
      view->ClipPlaneChanged ();
    }
  };
  
  //-------------------------------------------------------------------
  
  class ResultViewerImpl::View3D::RenderViewMouseCallback : public vtkCommand
  {
  public:
    View3D* view;
    
    static RenderViewMouseCallback* New ()
    {
      return new RenderViewMouseCallback ();
    }
    
    void Execute (vtkObject* caller, unsigned long eventId, void* callData)
    {
      vtkRenderWindowInteractor* iren = reinterpret_cast<vtkRenderWindowInteractor*> (caller);
      int x, y;
      iren->GetEventPosition (x, y);
      
      switch (eventId)
      {
      case vtkCommand::MouseMoveEvent:
	view->HandleMouseMove (x, y);
	break;
      }
    }
  };
  
  //-------------------------------------------------------------------
  
  BEGIN_DECLARE_EVENT_TYPES()
    DECLARE_EVENT_TYPE(EVENT_UPDATE_CAMERA, 0)
  END_DECLARE_EVENT_TYPES()
  DEFINE_EVENT_TYPE(EVENT_UPDATE_CAMERA)
  
  class ResultViewerImpl::View3D::UpdateCameraEvent : public wxEvent
  {
    vtkCamera* cam;
  public:
    UpdateCameraEvent (vtkCamera* cam) : wxEvent (0, EVENT_UPDATE_CAMERA),
      cam (cam) {}
      
    wxEvent* Clone() const { return new UpdateCameraEvent (cam); }
    
    vtkCamera* GetCamera() const { return cam; }
  };
  
  #define EVT_UPDATE_CAMERA(fn) \
      DECLARE_EVENT_TABLE_ENTRY( EVENT_UPDATE_CAMERA, wxID_ANY, -1, \
      (wxObjectEventFunction) (wxEventFunction) \
      wxStaticCastEvent( UpdateCameraEventFunction, & fn ), (wxObject *) NULL ),

  
  // ------------------------------------------------------------------------

  BEGIN_DECLARE_EVENT_TYPES()
    DECLARE_EVENT_TYPE(EVENT_DATASET_FRAME_CHANGED, 0)
  END_DECLARE_EVENT_TYPES()
  DEFINE_EVENT_TYPE(EVENT_DATASET_FRAME_CHANGED)
  
  struct ResultViewerImpl::View3D::SharedViewData : public ViewPanel::SharedDataBase
  {
    wxMenu renderModeMenu;
    wxBitmap renderModeButtonImages[numRenderModes];
    wxBitmap imgShowLegend;
    wxBitmap imgClipPlane;
    wxBitmap imgLinkViews;
    
    wxBitmap imgDisplacementOffset;
    wxMenu displacementDirMenu;
    wxBitmap displacementDirButtonImages[numDisplacementDirModes];
    
    vtkSmartPointer<vtkArrowSource> displaceDirectionGlyphSource;
    
    Gradients gradients;
    std::vector<wxBitmap> gradientToolImages;
    wxMenu gradientMenu;
    
    // Linked views
    unsigned int linkedViewCount;
    vtkSmartPointer<vtkCamera> camTemplate;
    // Data set frame linking
    size_t linkedFrame;
    
    SharedViewData() : linkedViewCount (0), linkedFrame (0) {}
  };
  
  enum
  {
    ID_Toolbar = 1,
    ID_DisplayData,
    ID_VisOption,
    ID_ActorRenderMode,
    ID_ShowLegend,
    ID_LegendOptions,
    ID_ClipPlane,
    ID_LinkViews,
    ID_DisplacementOffset,
    ID_DisplacementDirMode,
    
    ID_DataSetSlider,
    
    ID_RenderModeFirst = 100,
    ID_DisplacementDirModeFirst = 110
  };
  
  BEGIN_EVENT_TABLE(ResultViewerImpl::View3D, ViewPanelContentVTK)
    EVT_WINDOW_CREATE(ResultViewerImpl::View3D::OnWindowCreate)
    
    EVT_UPDATE_CAMERA(ResultViewerImpl::View3D::OnUpdateCamera)
    
    EVT_DIRECTION_SCALE_CHANGED(ResultViewerImpl::View3D::OnDisplacementScaleChange)
    
    EVT_COMMAND_SCROLL(ID_DataSetSlider, ResultViewerImpl::View3D::OnDataSetSelectionChanged)
    EVT_COMMAND(wxID_ANY, EVENT_DATASET_FRAME_CHANGED, ResultViewerImpl::View3D::OnLinkedDataSetChanged)
  END_EVENT_TABLE()
  
  ResultViewerImpl::View3D::View3D (ViewPanel* parent, const View3D* cloneFrom)
   : ViewPanelContentVTK (parent),
     renderMode (0),
     displacementDirection (ddNone),
     oldDisplacementDirection (ddColored),
     topBar (nullptr),
     useDisplaceData (false),
     displacementData ((size_t)~0),
     displaceDirectionsDataDS (0),
     displacementDirScale (1),
     useClipper (false),
     updatingCam (false),
     useLinkView (false),
     gradientMenuHandler (this)
  {
    sharedData = parent->GetSharedData<SharedViewData> ();
    if (!sharedData)
    {
      SetupSharedData ();
      parent->SetSharedData (sharedData);
    }
    
    wxSizer* sizer = new wxBoxSizer (wxVERTICAL);
    
    vtkwx::RenderWidget* renderWidget = CreateRenderWidget (this);
    renderWidget->EnableKeyboardHandling (false);
    sizer->Add (renderWidget, 1, wxEXPAND);
    
    dataSetSelectionBar = new wxPanel (this);
    wxSizer* dataSetSelectionBarSizer = new wxBoxSizer (wxHORIZONTAL);
    dataSetSelectionSlider = new uicommon::TextCtrlBuddySlider (dataSetSelectionBar, ID_DataSetSlider, 0, 0, 1);
    dataSetSelectionSlider->SetToolTip (wxString::Format (wxT ("Select which result data %s to display"),
							  uicommon::Quote::Double (wxT ("frame")).c_str()));
    dataSetSelectionBarSizer->Add (dataSetSelectionSlider, wxSizerFlags().Proportion(1).Expand());
    dataSetSelectionBar->SetSizer (dataSetSelectionBarSizer);
    sizer->Add (dataSetSelectionBar, wxSizerFlags().Proportion(0).Expand());
    
    SetSizer (sizer);
    
    displacementSizePanel = new DisplacementDirectionSizePanel (this, wxDefaultPosition);
    // Is positioned later, once layouting has happened.
    displacementSizePanel->Hide();
    
    if (cloneFrom)
    {
      this->data = cloneFrom->data;
    }
    
    cellPicker = vtkSmartPointer<vtkCellPicker>::New();
  }
  
  ResultViewerImpl::View3D::~View3D ()
  {
    // To ensure linkedViewCount is correctly decremented
    SetLinkView (false);
  }

  void ResultViewerImpl::View3D::SetupSharedData ()
  {
    sharedData = boost::make_shared<SharedViewData> ();
    
    sharedData->imgShowLegend = wxArtProvider::GetBitmap (wxART_MAKE_ART_ID(show-legend), wxART_TOOLBAR);
    sharedData->imgClipPlane = wxArtProvider::GetBitmap (wxART_MAKE_ART_ID(clip), wxART_TOOLBAR);
    sharedData->imgLinkViews = wxArtProvider::GetBitmap (wxART_MISSING_IMAGE,
							 wxART_TOOLBAR);
    
    static const wxChar* const renderModeNames[numRenderModes] =
    {
      wxT ("&Shaded"),
      wxT ("S&haded + Wireframe"),
      wxT ("&Flat"),
      wxT ("F&lat + Wireframe"),
      wxT ("&Wireframe")
    };
    static const wxChar* const renderModeArtNames[numRenderModes] =
    {
      wxART_MAKE_ART_ID(render-shaded),
      wxART_MAKE_ART_ID(render-shaded+wf),
      wxART_MAKE_ART_ID(render-flat),
      wxART_MAKE_ART_ID(render-flat+wf),
      wxART_MAKE_ART_ID(render-wf)
    };
    for (int m = 0; m < numRenderModes; m++)
    {
      wxMenuItem* newItem = new wxMenuItem (&sharedData->renderModeMenu,
					    ID_RenderModeFirst + m,
					    renderModeNames[m],
					    wxEmptyString,
					    /*wxITEM_RADIO*/wxITEM_NORMAL);
      newItem->SetBitmap (wxArtProvider::GetBitmap (renderModeArtNames[m], wxART_MENU));
      sharedData->renderModeMenu.Append (newItem);
      
      sharedData->renderModeButtonImages[m] = wxArtProvider::GetBitmap (renderModeArtNames[m],
									wxART_TOOLBAR);
    }

    sharedData->imgDisplacementOffset = wxArtProvider::GetBitmap (wxART_MAKE_ART_ID(displace-offset),
								  wxART_TOOLBAR);
    static const wxChar* const displacementDirNames[numDisplacementDirModes] =
    {
      wxT ("Do&n't show displacement directions"),
      wxT ("Show displacement directions &scaled by magnitude"),
      wxT ("Show displacement directions &colored by magnitude")
    };
    static const wxChar* const displacementDirArtNames[numRenderModes] =
    {
      wxART_MAKE_ART_ID(displace-dir-none),
      wxART_MAKE_ART_ID(displace-dir-scaled),
      wxART_MAKE_ART_ID(displace-dir-colored)
    };
    for (int m = 0; m < numDisplacementDirModes; m++)
    {
      wxMenuItem* newItem = new wxMenuItem (&sharedData->displacementDirMenu,
					    ID_DisplacementDirModeFirst + m,
					    displacementDirNames[m],
					    wxEmptyString,
					    /*wxITEM_RADIO*/wxITEM_NORMAL);
      newItem->SetBitmap (wxArtProvider::GetBitmap (displacementDirArtNames[m],
						    wxART_MENU));
      sharedData->displacementDirMenu.Append (newItem);
      
      sharedData->displacementDirButtonImages[m] =
	wxArtProvider::GetBitmap (displacementDirArtNames[m], wxART_TOOLBAR);
    }
    
    wxSize toolImageSize (wxArtProvider::GetSizeHint (wxART_TOOLBAR));
    // Note: for horizontal gradients a width of toolImageSize.x would be better looking
    wxSize menuImageSize (wxArtProvider::GetSizeHint (wxART_MENU));
    for (size_t i = 0; i < sharedData->gradients.GetGradientsNum(); i++)
    {
      sharedData->gradientToolImages.push_back (sharedData->gradients.RenderGradient (i, toolImageSize));
      
      wxMenuItem* newItem = new wxMenuItem (&sharedData->gradientMenu,
					    i,
					    sharedData->gradients.GetGradientName (i),
					    wxEmptyString,
					    wxITEM_NORMAL);
      newItem->SetBitmap (sharedData->gradients.RenderGradient (i, menuImageSize));
      sharedData->gradientMenu.Append (newItem);
    }
  }
  
  void ResultViewerImpl::View3D::SetData (const DataConstPtr& data)
  {
    this->data = data;
    displacementData = (size_t)~0;
    
    dataSetMapper.InterpolateScalarsBeforeMappingOn ();
    
    origDataSetFaces = vtkSmartPointer<DataSetFaceExtractor>::New();
    vtkSmartPointer<UnstructuredGridWithQuadsClippingToPoly> faceExtractOut =
      vtkSmartPointer<UnstructuredGridWithQuadsClippingToPoly>::New ();
    origDataSetFaces->SetOutput (faceExtractOut);
    origDataSetEdges = vtkSmartPointer<DataSetEdgeExtractor>::New ();
    origDataSetEdges->SetInput (origDataSetFaces->GetOutput());
    origDataSetMapper = vtkSmartPointer<vtkDataSetMapper>::New ();
    origDataSetMapper->SetInput (origDataSetEdges->GetOutput());
    
    /* Set a scale for displacement direction arrows; use the min distance
       along all edges of the cells of the model for that. */
    {
      float minEdge = HUGE_VAL;
      vtkDataSet* dataset0 = data->GetDataSet (0);
      DataSetEdgeIterator edges (dataset0);
      
      while (edges.HasNext())
      {
	DataSetEdgeIterator::Edge edge (edges.Next());
	double pt1[3];
	double pt2[3];
	dataset0->GetPoint (edge.first, pt1);
	dataset0->GetPoint (edge.second, pt2);
	double dist = sqrt (((pt1[0] - pt2[0]) * (pt1[0] - pt2[0]))
			  + ((pt1[1] - pt2[1]) * (pt1[1] - pt2[1]))
			  + ((pt1[2] - pt2[2]) * (pt1[2] - pt2[2])));
	minEdge = std::min (minEdge, float (dist));
      }
      displacementDirScale = minEdge;
    }
    
    displayDataChoice->Clear();
    displayDataChoice->Append (wxT ("Solid color"));
    displayDataChoice->SetSelection (0);
    
    for (size_t i = 0; i < data->GetNumDataArrays(); i++)
    {
      displayDataChoice->Append (data->GetDataArrayName (i));
    }
    
    lastVisOpt.clear ();
    lastVisOpt.resize (data->GetNumDataArrays());
    
    for (size_t i = 0; i < data->GetNumDataArrays(); i++)
    {
      if (data->IsDataArrayDisplacement (i))
      {
	displacementData = i;
	
	// Get maximum displacement value
	double displaceMagRange[2];
	data->GetDataArrayMagnitudeRange (i, displaceMagRange);
	/* Compute a default displacement scale which scales the maximum magnitude
	 * to show up with a length of displacementDirScale */
	float defaultDisplacementScale = displacementDirScale/displaceMagRange[1];
	// But never downscale by default.
	defaultDisplacementScale = std::max (defaultDisplacementScale, 1.0f);
	
	displacementSizePanel->SetDisplacementScale (defaultDisplacementScale);
	
	break;
      }
    }
    actorOptionsTB->EnableTool (ID_DisplacementOffset, displacementData != (size_t)~0);
    actorOptionsTB->EnableTool (ID_DisplacementDirMode, displacementData != (size_t)~0);
    if (renderer && origDataSetActor) renderer->RemoveActor (origDataSetActor);
    
    // Choice box strings have changed, so update size
    UpdateToolbarControlMinSize (displayDataChoice, toolbar);
    
    if (data->GetDataSetNum() > 1)
    {
      dataSetSelectionSlider->SetRange (0, data->GetDataSetNum()-1);
      dataSetSelectionBar->Show();
    }
    else
      dataSetSelectionBar->Hide();
    
    SetDisplayedDataSet (0);
  }

  wxWindow* ResultViewerImpl::View3D::CreateTopTools (wxWindow* parentWindow)
  {
    topBar = new wxPanel (parentWindow);
    
    topBarSizer = new wxBoxSizer (wxHORIZONTAL);
    
    toolbar = new wxAuiToolBar (topBar, ID_Toolbar, wxDefaultPosition, wxDefaultSize,
				wxAUI_TB_HORZ_LAYOUT | wxAUI_TB_NO_AUTORESIZE);
    toolbar->SetToolBitmapSize (wxArtProvider::GetSizeHint (wxART_TOOLBAR));
    displayDataChoice = new wxChoice (toolbar, ID_DisplayData);
    displayDataChoice->SetToolTip (wxT ("Data"));
    toolbar->AddControl (displayDataChoice);
    toolbar->Realize();
    topBarSizer->Add (toolbar, wxSizerFlags(0).Expand());
    
    {
      visOptionEmpty = new wxAuiToolBar (topBar, wxID_ANY, wxDefaultPosition, wxDefaultSize,
					 wxAUI_TB_HORZ_LAYOUT | wxAUI_TB_NO_AUTORESIZE);
      visOptionEmpty->Realize();
    }
    topBarSizer->Add (visOptionEmpty, wxSizerFlags (0).Expand());
    
    {
      visOptionChoice = new wxAuiToolBar (topBar, wxID_ANY, wxDefaultPosition, wxDefaultSize,
					  wxAUI_TB_HORZ_LAYOUT | wxAUI_TB_NO_AUTORESIZE);
      visOptionChoice->Hide ();
      visChoiceCtrl = new wxChoice (visOptionChoice, ID_VisOption);
      visChoiceCtrl->SetToolTip (wxT ("Component"));
      visOptionChoice->AddControl (visChoiceCtrl);
      visOptionChoice->Realize ();
    }
    
    actorOptionsTB = new wxAuiToolBar (topBar, wxID_ANY, wxDefaultPosition, wxDefaultSize,
				       wxAUI_TB_HORZ_LAYOUT | wxAUI_TB_NO_AUTORESIZE);
    actorOptionsTB->AddTool (ID_ActorRenderMode, wxEmptyString,
			     sharedData->renderModeButtonImages[0],
			     wxT ("Rendering mode"));
    actorOptionsTB->SetToolDropDown (ID_ActorRenderMode, true);
    actorOptionsTB->AddTool (ID_ShowLegend, wxEmptyString,
			     sharedData->imgShowLegend,
			     wxT ("Show legend"),
			     wxITEM_CHECK);
    actorOptionsTB->AddTool (ID_LegendOptions, wxEmptyString,
			     sharedData->gradientToolImages[0],
			     wxT ("Select gradient"),
			     wxITEM_NORMAL);
    actorOptionsTB->SetToolDropDown (ID_LegendOptions, true);
    actorOptionsTB->AddSeparator ();
    actorOptionsTB->AddTool (ID_ClipPlane, wxEmptyString,
			     sharedData->imgClipPlane,
			     wxT ("Clipping plane"),
			     wxITEM_CHECK);
    actorOptionsTB->AddTool (ID_LinkViews, wxEmptyString,
			     sharedData->imgLinkViews,
			     wxT ("Link display to other views"),
			     wxITEM_CHECK);
    actorOptionsTB->ToggleTool (ID_LinkViews, true);
    actorOptionsTB->AddSeparator ();
    actorOptionsTB->AddTool (ID_DisplacementOffset, wxEmptyString,
			     sharedData->imgDisplacementOffset,
			     wxT ("Apply displacement to model"),
			     wxITEM_CHECK);
    actorOptionsTB->EnableTool (ID_DisplacementOffset, false);
    actorOptionsTB->AddTool (ID_DisplacementDirMode, wxEmptyString,
			     sharedData->displacementDirButtonImages[0],
			     wxT ("Show displacement directions"));
    actorOptionsTB->SetToolDropDown (ID_DisplacementDirMode, true);
    actorOptionsTB->EnableTool (ID_DisplacementDirMode, false);
    actorOptionsTB->Realize ();
    topBarSizer->Add (actorOptionsTB, wxSizerFlags (1).Expand ());
    
    topBar->SetSizer (topBarSizer);
    
    // Connect events
    topBar->Connect (ID_DisplayData,
		     wxEVT_COMMAND_CHOICE_SELECTED,
		     wxCommandEventHandler (ResultViewerImpl::View3D::OnDisplayDataChanged),
		     nullptr, this);
    topBar->Connect (ID_VisOption,
		     wxEVT_COMMAND_CHOICE_SELECTED,
		     wxCommandEventHandler (ResultViewerImpl::View3D::OnVisOptionChanged),
		     nullptr, this);
    
    topBar->Connect (ID_ActorRenderMode,
		     wxEVT_COMMAND_AUITOOLBAR_TOOL_DROPDOWN,
		     wxAuiToolBarEventHandler (ResultViewerImpl::View3D::OnRenderModeDropDown),
		     nullptr, this);
    topBar->Connect (ID_RenderModeFirst,
		     ID_RenderModeFirst+ResultViewerImpl::View3D::numRenderModes-1,
		     wxEVT_COMMAND_MENU_SELECTED,
		     wxCommandEventHandler (ResultViewerImpl::View3D::OnRenderModeCommand),
		     nullptr, this);
    topBar->Connect (ID_ShowLegend,
		     wxEVT_COMMAND_MENU_SELECTED,
		     wxCommandEventHandler (ResultViewerImpl::View3D::OnShowLegend),
		     nullptr, this);
    topBar->Connect (ID_ShowLegend,
		     wxEVT_UPDATE_UI,
		     wxUpdateUIEventHandler (ResultViewerImpl::View3D::OnShowLegendUpdateUI),
		     nullptr, this);
    topBar->Connect (ID_LegendOptions,
		     wxEVT_COMMAND_AUITOOLBAR_TOOL_DROPDOWN,
		     wxAuiToolBarEventHandler (ResultViewerImpl::View3D::OnLegendOptions),
		     nullptr, this);
    topBar->Connect (ID_LegendOptions,
		     wxEVT_UPDATE_UI,
		     wxUpdateUIEventHandler (ResultViewerImpl::View3D::OnLegendOptionsUpdateUI),
		     nullptr, this);
    topBar->Connect (ID_ClipPlane,
		     wxEVT_COMMAND_MENU_SELECTED,
		     wxCommandEventHandler (ResultViewerImpl::View3D::OnClipPlane),
		     nullptr, this);
    topBar->Connect (ID_LinkViews,
		     wxEVT_COMMAND_MENU_SELECTED,
		     wxCommandEventHandler (ResultViewerImpl::View3D::OnLinkViews),
		     nullptr, this);
    topBar->Connect (ID_LinkViews,
		     wxEVT_UPDATE_UI,
		     wxUpdateUIEventHandler (ResultViewerImpl::View3D::OnLinkViewsUpdateUI),
		     nullptr, this);
    topBar->Connect (ID_DisplacementOffset,
		     wxEVT_COMMAND_MENU_SELECTED,
		     wxCommandEventHandler (ResultViewerImpl::View3D::OnDisplacementOffset),
		     nullptr, this);
    topBar->Connect (ID_DisplacementDirMode,
		     wxEVT_COMMAND_AUITOOLBAR_TOOL_DROPDOWN,
		     wxAuiToolBarEventHandler (ResultViewerImpl::View3D::OnDisplacementDirDropDown),
		     nullptr, this);
    topBar->Connect (ID_DisplacementDirModeFirst,
		     ID_DisplacementDirModeFirst+ResultViewerImpl::View3D::numDisplacementDirModes-1,
		     wxEVT_COMMAND_MENU_SELECTED,
		     wxCommandEventHandler (ResultViewerImpl::View3D::OnDisplacementDirCommand),
		     nullptr, this);
    
    if (data) SetData (data);

    return topBar;
  }

  void ResultViewerImpl::View3D::OnWindowCreate (wxWindowCreateEvent& event)
  {
    if (event.GetWindow() == this)
    {
      // Layouting was done, can place panel.
      wxPoint displacementPanelPos (GetRenderWidget()->GetPosition());
      displacementPanelPos.x += 8;
      displacementPanelPos.y += 8;
      displacementSizePanel->SetPosition (displacementPanelPos);
    }
  }
  
  void ResultViewerImpl::View3D::SetupVTKRenderer ()
  {
    if (!renderer)
    {
      SetupRenderer();
      
      // Default render mode
      ApplyRenderMode ();
      SetUIRenderMode ();
      
      // Default displacement display to 'on'
      if (displacementData != (size_t)~0)
      {
	useDisplaceData = true;
	ComputeDisplacementOffset ();
	ShowDisplacementOffset ();

	CheckDisplacementSizePanelVisibility ();
	
	actorOptionsTB->ToggleTool (ID_DisplacementOffset, useDisplaceData);
      }
    }
  }

  void ResultViewerImpl::View3D::SetupRenderer ()
  {
    renderer = vtkSmartPointer<vtkRenderer>::New ();
    GetRenderWidget()->GetRenderWindow()->AddRenderer (renderer);
    
    vtkSmartPointer<RenderViewMouseCallback> mouseCallback = vtkSmartPointer<RenderViewMouseCallback>::New ();
    mouseCallback->view = this;
    GetRenderWidget()->GetRenderWindow()->GetInteractor()->AddObserver (vtkCommand::MouseMoveEvent, mouseCallback);
    
    vtkSmartPointer<vtkInteractorStyleSwitch> istyle = vtkSmartPointer<vtkInteractorStyleSwitch>::New ();
    istyle->SetCurrentStyleToTrackballCamera ();
    GetRenderWidget()->GetRenderWindow()->GetInteractor()->SetInteractorStyle (istyle);
    
    vtkSmartPointer<vtkActor> mapperActor = dataSetMapper.CreateMapperActor();
    renderer->AddActor (mapperActor);
    
    vtkSmartPointer<vtkActor> edgesActor = dataSetMapper.CreateEdgesActor();
    renderer->AddActor (edgesActor);
    
    origDataSetActor = vtkSmartPointer<vtkActor>::New ();
    origDataSetActor->SetMapper (origDataSetMapper);
    origDataSetActor->GetProperty()->SetOpacity (0.3);
    origDataSetActor->GetProperty()->SetColor (0, 0, 0);
    origDataSetActor->GetProperty()->SetLighting (false);
    origDataSetActor->GetProperty()->SetEdgeVisibility (false);
    
    renderer->SetBackground (0.6, 0.6, 0.6);
    
    scalarBar = vtkSmartPointer<vtkScalarBarWidget>::New ();
    scalarBar->SetInteractor (GetRenderWidget()->GetRenderWindow()->GetInteractor());
    scalarBar->KeyPressActivationOff();
    
    orientationMarker = vtkSmartPointer<vtkOrientationMarkerWidget>::New ();
    orientationMarker->SetInteractor (GetRenderWidget()->GetRenderWindow()->GetInteractor());
    vtkSmartPointer<vtkAxesActor> axes (vtkSmartPointer<vtkAxesActor>::New ());
    orientationMarker->SetOrientationMarker (axes);
    orientationMarker->EnabledOn();
    orientationMarker->KeyPressActivationOff();
    orientationMarker->InteractiveOff();
    
    clipPlaneWidget = vtkSmartPointer<vtkPlaneWidget>::New ();
    clipPlaneWidget->SetInteractor (GetRenderWidget()->GetRenderWindow()->GetInteractor());
    clipPlaneWidget->KeyPressActivationOff();
    clipPlaneWidget->SetInput (data->GetDataSet (currentDataSet));
    clipPlaneWidget->PlaceWidget();
    clipPlane = vtkSmartPointer<vtkPlane>::New ();
    clipPlaneWidget->GetPlane (clipPlane);
    
    // Callback for clip plane widget changed
    vtkSmartPointer<ClipPlaneChangedCallback> clipPlaneCallback (vtkSmartPointer<ClipPlaneChangedCallback>::New ());
    clipPlaneCallback->view = this;
    clipPlaneWidget->AddObserver (vtkCommand::InteractionEvent, clipPlaneCallback);
    
    /* Manually center camera on mapper bounds.
       (Automatic centering would also consider the orientation marker.) */
    double mapperBounds[6];
    dataSetMapper.GetBounds (mapperBounds);
    renderer->ResetCamera (mapperBounds);

    // Cam starts out as linked. Copies camera settings if this isn't the first view.
    SetLinkView (true);
    
    // Install callback to listen for camera modifications
    vtkSmartPointer<CameraModifiedCallback> camCallback (vtkSmartPointer<CameraModifiedCallback>::New ());
    camCallback->view = this;
    renderer->GetActiveCamera ()->AddObserver (vtkCommand::ModifiedEvent,
					       camCallback);

    /* Create a LUT for displacements display as we need the displacement
	data's range on it */
    displacementColors = vtkSmartPointer<vtkLookupTable>::New();
    /* Use a black-to-magenta ramp for displacement direction coloring.
	Slightly ugly, but doesn't clash with default colors ... */
    displacementColors->SetHueRange (5.0/6.0, 5.0/6.0);
    displacementColors->SetValueRange (0, 1);
    displacementColors->SetVectorModeToMagnitude();
    if (displacementData != (size_t)~0)
    {
      double range[2];
      data->GetDataArrayMagnitudeRange (displacementData, range);
      displacementColors->SetRange (range);
    }
  }

  bool ResultViewerImpl::View3D::UpdateCameraPositions (vtkCamera* cam)
  {
    vtkCamera* myCam = renderer->GetActiveCamera ();
    if (cam == myCam) return false;
    /* Change focal point + position but keep distance.
       However, changing the distance in VTK changes the focal point
       position. Changing the eye position works better, so do that
       manually. */
    double oldDistance = myCam->GetDistance ();
    double newProjDir[3];
    cam->GetDirectionOfProjection (newProjDir);
    double newPos[3];
    for (int c = 0; c < 3; c++)
    {
      newPos[c] = cam->GetFocalPoint ()[c] - oldDistance*newProjDir[c];
    }
    
    myCam->SetFocalPoint (cam->GetFocalPoint ());
    myCam->SetPosition (newPos);
    myCam->SetViewUp (cam->GetViewUp ());
    
    renderer->ResetCameraClippingRange ();
    
    return true;
  }

  vtkCamera* ResultViewerImpl::View3D::GetCamera () const
  {
    return renderer->GetActiveCamera ();
  }
  
  void ResultViewerImpl::View3D::CameraChanged (vtkCamera* cam)
  {
    if (updatingCam) return;
    // Unlinked, don't propagate
    if (!useLinkView) return;
    
    UpdateCameraEvent event (cam);
    PostToOthers (event);
  }

  void ResultViewerImpl::View3D::ClipPlaneChanged ()
  {
    renderer->ResetCameraClippingRange ();
    clipPlaneWidget->GetPlane (clipPlane);
  }

  bool ResultViewerImpl::View3D::SetLinkView (bool flag)
  {
    if (useLinkView == flag) return false;
    
    bool ret;
    if (flag)
    {
      sharedData->linkedViewCount++;
      if (sharedData->linkedViewCount == 1)
      {
	/* This is the only view marked as "linked".
	   Set the current camera as the template the next views will
	   copy their settings from upon linking. */
	sharedData->camTemplate = renderer->GetActiveCamera ();
	// Also set current data set frame as linked
	sharedData->linkedFrame = currentDataSet;
	ret = false;
      }
      else
      {
	// Other views are linked, use existing template
	UpdateCameraPositions (sharedData->camTemplate);
	// Also copy existing data set frame
	dataSetSelectionSlider->SetValue (sharedData->linkedFrame);
	SetDisplayedDataSet (sharedData->linkedFrame);
	ret = true;
      }
    }
    else
    {
      sharedData->linkedViewCount--;
      ret = false;
    }
    useLinkView = flag;
    return ret;
  }

  void ResultViewerImpl::View3D::OnDisplayDataChanged (wxCommandEvent& event)
  {
    if (event.GetSelection() < 0) return;
    if (!dataSetMapper) return;
    
    wxWindow* visOptions = nullptr;
    bool useScalarBar = false;
    
    int sel = event.GetSelection();
    if (sel <= 0)
    {
      dataSetMapper.SetColorModeToDefault();
      dataSetMapper.ScalarVisibilityOff();
      visOptions = visOptionEmpty;
    }
    else
    {    
      visOptions = visOptionChoice;
      
      switch (data->GetDataArrayAssociation (sel-1))
      {
      case Data::perCell:
	dataSetMapper.SetScalarModeToUseCellFieldData();
	break;
      case Data::perPoint:
	dataSetMapper.SetScalarModeToUsePointFieldData();
	break;
      }
      
      dataSetMapper.SetColorModeToMapScalars();
      dataSetMapper.ScalarVisibilityOn();
      dataSetMapper.SelectColorArray (data->GetDataArrayIndex (sel-1));
      
      UpdateVisOptionChoice (sel-1, lastVisOpt[sel-1].comp);
      SetVisComponent (sel-1, lastVisOpt[sel-1].comp);
      useScalarBar = lastVisOpt[sel-1].legend;
      actorOptionsTB->ToggleTool (ID_ShowLegend, useScalarBar);
      size_t gradient = lastVisOpt[sel-1].gradient;
      UpdateGradientUI (gradient);
      SetGradient (gradient);
      actorOptionsTB->Refresh ();
    }
    
    wxWindow* currentVisOpt = topBarSizer->GetItem (1)->GetWindow();
    if (currentVisOpt != visOptions)
    {
      visOptions->Show ();
      topBarSizer->Replace (currentVisOpt, visOptions);
      topBarSizer->Layout ();
      currentVisOpt->Hide ();
    }
    
    scalarBar->SetEnabled (useScalarBar);
    dataSetMapper.Update();
    Render();
  }
  
  void ResultViewerImpl::View3D::OnVisOptionChanged (wxCommandEvent& event)
  {
    int displaySel = displayDataChoice->GetSelection();
    if (displaySel <= 0) return;
    
    int visComp = event.GetSelection()-1;
    SetVisComponent (displaySel-1, visComp);
    lastVisOpt[displaySel-1].comp = visComp;

    dataSetMapper.Update();
    Render();
  }

  void ResultViewerImpl::View3D::UpdateVisOptionChoice (size_t arrayIndex, int initialSel)
  {
    visChoiceCtrl->Clear ();
    visChoiceCtrl->Append (wxT ("Magnitude"));
    for (int i = 0; i < data->GetDataArrayComponents (arrayIndex); i++)
    {
      visChoiceCtrl->Append (wxString::Format (wxT ("%d"), i + 1));
    }
    visChoiceCtrl->SetSelection (initialSel+1);
    
    /* WX 2.8 computes a slightly different height for visChoiceCtrl,
       so manually force the desired height */
    UpdateToolbarControlMinSize (visChoiceCtrl, visOptionChoice,
				 displayDataChoice->GetBestSize().GetHeight());
  }

  void ResultViewerImpl::View3D::SetVisComponent (size_t arrayIndex, int visComp)
  {
    vtkSmartPointer<vtkScalarsToColors> lut (dataSetMapper.GetLookupTable ());
    double range[2];
    if (visComp < 0)
    {
      lut->SetVectorModeToMagnitude();
      data->GetDataArrayMagnitudeRange (arrayIndex, range);
    }
    else
    {
      lut->SetVectorComponent (visComp);
      lut->SetVectorModeToComponent();
      data->GetDataArrayCompValueRange (arrayIndex, range);
    }
    dataSetMapper.SetScalarRange (range);
  }

  void ResultViewerImpl::View3D::UpdateGradientUI (size_t gradient)
  {
    actorOptionsTB->SetToolBitmap (ID_LegendOptions,
				   sharedData->gradientToolImages[gradient]);
    actorOptionsTB->Refresh ();
				   
    int displaySel = displayDataChoice->GetSelection();
    if (displaySel <= 0) return;
    
    lastVisOpt[displaySel-1].gradient = gradient;
  }
  
  void ResultViewerImpl::View3D::SetGradient (size_t gradient)
  {
    if (gradient == (size_t)~0) return;

    vtkSmartPointer<vtkScalarsToColors> colors (dataSetMapper.GetLookupTable ());
    vtkSmartPointer<vtkLookupTable> newColors (vtkSmartPointer<vtkLookupTable>::New ());
    newColors->DeepCopy (sharedData->gradients.GetGradientColors (gradient)); // TODO: Look into sharing the table data
    newColors->SetVectorMode (colors->GetVectorMode());
    dataSetMapper.SetLookupTable (newColors);
    scalarBar->GetScalarBarActor()->SetLookupTable (newColors);
  }

  void ResultViewerImpl::View3D::OnRenderModeDropDown (wxAuiToolBarEvent& event)
  {
    if (event.IsDropDownClicked())
    {
      wxPoint popupPos (actorOptionsTB->GetToolRect (ID_ActorRenderMode).GetBottomLeft());
      popupPos = actorOptionsTB->ClientToScreen (popupPos);
      popupPos = ScreenToClient (popupPos);
      
      PopupMenu (&sharedData->renderModeMenu, popupPos);
    }
    else
    {
      renderMode = (renderMode + 1) % numRenderModes;
      ApplyRenderMode ();
      SetUIRenderMode ();

      Render();
    }
  }
  
  void ResultViewerImpl::View3D::OnRenderModeCommand (wxCommandEvent& event)
  {
    renderMode = event.GetId() - ID_RenderModeFirst;
    ApplyRenderMode ();
    SetUIRenderMode ();

    Render();
  }

  namespace
  {
    struct RenderModeSettings
    {
      bool lighting;
      int repr;
      bool edgeVis;
    };
  }
  
  void ResultViewerImpl::View3D::ApplyRenderMode ()
  {
    static const RenderModeSettings rmSettings[numRenderModes] =
    {
      // Shaded
      { true,	VTK_SURFACE,	false },
      // Shaded + Wireframe
      { true,	VTK_SURFACE,	true },
      // Flat
      { false,	VTK_SURFACE,	false },
      // Flat + Wireframe
      { false,	VTK_SURFACE,	true },
      // Wireframe
      { false,	VTK_WIREFRAME,	false }
    };
    
    dataSetMapper.SetLighting (rmSettings[renderMode].lighting);
    dataSetMapper.SetRepresentation (rmSettings[renderMode].repr);
    dataSetMapper.SetEdgeVisibility (rmSettings[renderMode].edgeVis);
  }

  void ResultViewerImpl::View3D::SetUIRenderMode ()
  {
    //renderModeMenu.Check (ID_RenderModeFirst + mode, true);
    actorOptionsTB->SetToolBitmap (ID_ActorRenderMode,
				   sharedData->renderModeButtonImages[renderMode]);
    actorOptionsTB->Refresh();
  }

  void ResultViewerImpl::View3D::OnShowLegend (wxCommandEvent& event)
  {
    int displaySel = displayDataChoice->GetSelection();
    if (displaySel <= 0) return;
    
    scalarBar->SetEnabled (event.IsChecked());
    lastVisOpt[displaySel-1].legend = event.IsChecked();
    Render();
  }

  void ResultViewerImpl::View3D::OnShowLegendUpdateUI (wxUpdateUIEvent& event)
  {
    event.Enable (displayDataChoice->GetSelection() > 0);
  }

  void ResultViewerImpl::View3D::OnLegendOptions (wxAuiToolBarEvent& event)
  {
    int displaySel = displayDataChoice->GetSelection();
    if (displaySel <= 0) return;
    
    sharedData->gradientMenu.SetEventHandler (&gradientMenuHandler);
    DoToolDropDown (actorOptionsTB, ID_LegendOptions, &sharedData->gradientMenu);
  }
  
  void ResultViewerImpl::View3D::OnLegendOptionsUpdateUI (wxUpdateUIEvent& event)
  {
    event.Enable (displayDataChoice->GetSelection() > 0);
  }

  void ResultViewerImpl::View3D::OnGradientSelect (wxCommandEvent& event)
  {
    size_t gradientID = event.GetId();
    SetGradient (gradientID);
    UpdateGradientUI (gradientID);
    Render();
  }
  
  inline vtkSmartPointer<vtkDataSet> ClipWrapDataSet (vtkDataSet* dataset,
						      vtkImplicitFunction* clipFunc)
  {
    vtkSmartPointer<vtkClipDataSet> clipData = vtkSmartPointer<vtkClipDataSet>::New ();
    clipData->SetInput (dataset);
    clipData->SetClipFunction (clipFunc);
    return clipData->GetOutput();
  }

  void ResultViewerImpl::View3D::OnClipPlane (wxCommandEvent& event)
  {
    clipPlaneWidget->SetEnabled (event.IsChecked());
    
    if (event.IsChecked())
    {
      // clip original data, change edge extractor
      vtkSmartPointer<vtkDataSet> clippedOrigData =
	ClipWrapDataSet (origDataSetFaces->GetOutput(), clipPlane);
      origDataSetEdges->SetInput (clippedOrigData);
      // Handy method sets clipping automatically on last input
      dataSetMapper.SetClipFunction (clipPlane);
      
      if (displaceDirectionsData && displaceDirectionsGlyphs)
	displaceDirectionsGlyphs->SetInput (ClipWrapDataSet (displaceDirectionsData, clipPlane));
    
      useClipper = true;
    }
    else
    {
      // change edge extractor to unclipped original data
      origDataSetEdges->SetInput (origDataSetFaces->GetOutput());
      // disable clipping in dataSetMapper
      dataSetMapper.SetClipFunction (nullptr);
      
      if (displaceDirectionsData && displaceDirectionsGlyphs)
	displaceDirectionsGlyphs->SetInput (displaceDirectionsData);
    
      useClipper = false;
    }
    
    Render();
  }

  void ResultViewerImpl::View3D::OnLinkViews (wxCommandEvent& event)
  {
    if (SetLinkView (event.IsChecked ()))
      Render();
  }

  void ResultViewerImpl::View3D::OnLinkViewsUpdateUI (wxUpdateUIEvent& event)
  {
    /* Link views makes only sense with > 1 views */
    event.Enable (!sharedData.unique());
  }

  void ResultViewerImpl::View3D::OnDisplacementOffset (wxCommandEvent& event)
  {
    if (event.IsChecked())
    {
      useDisplaceData = true;
      ComputeDisplacementOffset ();
      ShowDisplacementOffset ();
    }
    else
    {
      useDisplaceData = false;
      HideDisplacementOffset ();
    }
    
    CheckDisplacementSizePanelVisibility ();

    Render();
  }

  void ResultViewerImpl::View3D::OnDisplacementDirDropDown (wxAuiToolBarEvent& event)
  {
    if (event.IsDropDownClicked())
    {
      wxPoint popupPos (actorOptionsTB->GetToolRect (ID_DisplacementDirMode).GetBottomLeft());
      popupPos = actorOptionsTB->ClientToScreen (popupPos);
      popupPos = ScreenToClient (popupPos);
      
      PopupMenu (&sharedData->displacementDirMenu, popupPos);
    }
    else
    {
      /* If the button is clicked, toggle between last displacement direction mode and
	 displacement directions off */
      if (displacementDirection == ddNone)
      {
	assert (oldDisplacementDirection != ddNone);
	displacementDirection = oldDisplacementDirection;
	ComputeDisplacementDirections ();
	ShowDisplacementDirections();
      }
      else
      {
	oldDisplacementDirection = displacementDirection;
	displacementDirection = ddNone;
	HideDisplacementDirections();
      }

      actorOptionsTB->SetToolBitmap (ID_DisplacementDirMode,
				    sharedData->displacementDirButtonImages[displacementDirection]);
      actorOptionsTB->Refresh();
      
      CheckDisplacementSizePanelVisibility ();
    
      Render();
    }
  }

  void ResultViewerImpl::View3D::OnDisplacementDirCommand (wxCommandEvent& event)
  {
    DisplacementDir oldDisplacementDirection = displacementDirection;
    displacementDirection = (DisplacementDir)(event.GetId() - ID_DisplacementDirModeFirst);
    
    switch (displacementDirection)
    {
    case ddNone:
      // If displacement directions were disabled, save last mode
      if (oldDisplacementDirection != ddNone)
	this->oldDisplacementDirection = oldDisplacementDirection;
      HideDisplacementDirections();
      break;
    case ddScaled:
    case ddColored:
      ComputeDisplacementDirections ();
      ShowDisplacementDirections();
      break;
    case numDisplacementDirModes:
      assert (false);
      break;
    }
    
    actorOptionsTB->SetToolBitmap (ID_DisplacementDirMode,
				   sharedData->displacementDirButtonImages[displacementDirection]);
    actorOptionsTB->Refresh();
    
    CheckDisplacementSizePanelVisibility ();

    Render();
  }

  void ResultViewerImpl::View3D::ShowDisplacementOffset ()
  {
    dataSetMapper.SetInput (displacedData);
    renderer->AddActor (origDataSetActor);
  }
  
  void ResultViewerImpl::View3D::HideDisplacementOffset ()
  {
    dataSetMapper.SetInput (data->GetDataSet (currentDataSet));
    renderer->RemoveActor (origDataSetActor);
  }

  void ResultViewerImpl::View3D::ComputeDisplacementOffset ()
  {
    if (displacedData) return;
      
    double scale = displacementSizePanel->GetDisplacementScale();
    vtkDataArray* displaceData = data->GetDataArrayRawData (currentDataSet, displacementData);
    
    vtkDataSet* dataset = data->GetDataSet (currentDataSet);
    vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New ();
    newPoints->SetNumberOfPoints (dataset->GetNumberOfPoints());
    for (vtkIdType point = 0; point < dataset->GetNumberOfPoints(); point++)
    {
      double pt[3];
      dataset->GetPoint (point, pt);
      double* ptDisplace = displaceData->GetTuple3 (point);
      for (int c = 0; c < 3; c++)
      {
	pt[c] += ptDisplace[c] * scale;
      }
      newPoints->SetPoint (point, pt);
    }
    vtkSmartPointer<vtkUnstructuredGrid> newDataSet = vtkSmartPointer<vtkUnstructuredGrid>::New ();
    newDataSet->SetPoints (newPoints);
    newDataSet->Allocate (dataset->GetNumberOfCells());
    for (vtkIdType cell = 0; cell < dataset->GetNumberOfCells(); cell++)
    {
      vtkCell* cellPtr = dataset->GetCell (cell);
      newDataSet->InsertNextCell (cellPtr->GetCellType(),
				  cellPtr->GetPointIds());
    }
    newDataSet->GetPointData()->ShallowCopy (dataset->GetPointData());
    newDataSet->GetCellData()->ShallowCopy (dataset->GetCellData());
    
    displacedData = newDataSet;
  }

  void ResultViewerImpl::View3D::ShowDisplacementDirections ()
  {
    if (!displaceDirectionsActor) return;

    if (useClipper)
      displaceDirectionsGlyphs->SetInput (ClipWrapDataSet (displaceDirectionsData, clipPlane));
    else
      displaceDirectionsGlyphs->SetInput (displaceDirectionsData);

    switch (displacementDirection)
    {
    case ddScaled:
      displaceDirectionsGlyphs->SetScaleFactor (displacementSizePanel->GetDisplacementScale());
      displaceDirectionsGlyphs->SetColorModeToColorByScale();
      displaceDirectionsGlyphs->SetScaleModeToScaleByVector();
      break;
    case ddColored:
      displaceDirectionsGlyphs->SetScaleFactor (displacementDirScale);
      displaceDirectionsGlyphs->SetColorModeToColorByVector();
      displaceDirectionsGlyphs->SetScaleModeToDataScalingOff();
      break;
    case ddNone:
    case numDisplacementDirModes:
      assert (false);
      break;
    }
    
    renderer->AddActor (displaceDirectionsActor);
  }
  
  void ResultViewerImpl::View3D::HideDisplacementDirections ()
  {
    if (!displaceDirectionsActor) return;
    renderer->RemoveActor (displaceDirectionsActor);
  }
  
  void ResultViewerImpl::View3D::ComputeDisplacementDirections ()
  {
    if (!sharedData->displaceDirectionGlyphSource)
    {
      sharedData->displaceDirectionGlyphSource = vtkSmartPointer<vtkArrowSource>::New ();
    }
    if (!displaceDirectionsData || (displaceDirectionsDataDS != currentDataSet))
    {
      vtkDataArray* displaceData = data->GetDataArrayRawData (currentDataSet, displacementData);
    
      displaceDirectionsData = vtkSmartPointer<vtkPolyData>::New ();
      
      vtkDataSet* dataset = data->GetDataSet (currentDataSet);
      // Points: simply copy the points from the original data set
      DataSetHelpers::CopyPoints (displaceDirectionsData, dataset);
      vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New ();
      verts->Allocate (dataset->GetNumberOfPoints());
      for (vtkIdType point = 0; point < dataset->GetNumberOfPoints(); point++)
      {
	verts->InsertNextCell (1, &point);
      }
      displaceDirectionsData->SetVerts (verts);
      
      // Normals: displacement vectors
      vtkSmartPointer<vtkFloatArray> normals;
      normals.TakeReference (vtkFloatArray::New());

      normals->SetNumberOfComponents (3);
      normals->SetNumberOfTuples (dataset->GetNumberOfPoints());
      for (vtkIdType norm = 0; norm < dataset->GetNumberOfPoints(); norm++)
      {
	double* ptDisplace = displaceData->GetTuple3 (norm);
	normals->SetTuple (norm, ptDisplace);
      }
      displaceDirectionsData->GetPointData()->SetNormals (normals);
      
      displaceDirectionsDataDS = currentDataSet;
    }
    if (!displaceDirectionsGlyphs)
    {
      displaceDirectionsGlyphs = vtkSmartPointer<vtkGlyph3D>::New ();
      displaceDirectionsGlyphs->SetSourceConnection (sharedData->displaceDirectionGlyphSource->GetOutputPort());
      displaceDirectionsGlyphs->SetVectorModeToUseNormal ();
    }
    
    if (!displaceDirectionsMapper)
    {
      displaceDirectionsMapper = vtkSmartPointer<vtkPolyDataMapper>::New ();
      displaceDirectionsMapper->SetInput (displaceDirectionsGlyphs->GetOutput ());
      displaceDirectionsMapper->UseLookupTableScalarRangeOn();
      displaceDirectionsMapper->SetLookupTable (displacementColors);
    }
    if (!displaceDirectionsActor)
    {
      displaceDirectionsActor = vtkSmartPointer<vtkActor>::New ();
      
      vtkSmartPointer<vtkProperty> prop = vtkSmartPointer<vtkProperty>::New ();
      prop->SetColor (1, 1, 0);
      displaceDirectionsActor->SetProperty (prop);
    }
    displaceDirectionsActor->SetMapper (displaceDirectionsMapper);
  }
  
  void ResultViewerImpl::View3D::DoToolDropDown (wxAuiToolBar* toolbar, int toolID, wxMenu* menu)
  {
    wxPoint popupPos (toolbar->GetToolRect (toolID).GetBottomLeft());
    popupPos = toolbar->ClientToScreen (popupPos);
    popupPos = ScreenToClient (popupPos);
    
    PopupMenu (menu, popupPos);
  }

  void ResultViewerImpl::View3D::UpdateToolbarControlMinSize (wxWindow* control,
							    wxAuiToolBar* toolbar,
							    int forceHeight)
  {
    wxAuiToolBarItem* ctrl_item = nullptr;
    for (size_t i = 0; i < toolbar->GetToolCount(); i++)
    {
      wxAuiToolBarItem* item = toolbar->FindToolByIndex (i);
      if (item->GetWindow() == control)
      {
	ctrl_item = item;
	break;
      }
    }
    if (ctrl_item)
    {
      wxSize newMinSize (control->GetBestSize());
      if (forceHeight > 0) newMinSize.SetHeight (forceHeight);
      ctrl_item->SetMinSize (newMinSize);
      toolbar->Realize();
    }
  }

  void ResultViewerImpl::View3D::OnUpdateCamera (UpdateCameraEvent& event)
  {
    // Not linked, so ignore event
    if (!useLinkView) return;
    
    // Going to change camera which is going to trigger the callback ... ignore that
    updatingCam = true;
    
    if (UpdateCameraPositions (event.GetCamera ()))
    {
      Render();
    }
    
    updatingCam = false;
  }

  void ResultViewerImpl::View3D::CheckDisplacementSizePanelVisibility ()
  {
    displacementSizePanel->Show (useDisplaceData || (displacementDirection == ddScaled));
  }

  void ResultViewerImpl::View3D::OnDisplacementScaleChange (wxEvent& event)
  {
    switch (displacementDirection)
    {
    case ddScaled:
      displaceDirectionsGlyphs->SetScaleFactor (displacementSizePanel->GetDisplacementScale());
      break;
    case ddColored:
      displaceDirectionsData = nullptr;
      ComputeDisplacementDirections();
      ShowDisplacementDirections();
      break;
    default:
      break;
    }
    
    if (useDisplaceData)
    {
      displacedData = nullptr;
      ComputeDisplacementOffset ();
      ShowDisplacementOffset ();
    }
    
    Render();
  }

  void ResultViewerImpl::View3D::OnDataSetSelectionChanged (wxScrollEvent& event)
  {
    SetDisplayedDataSet (event.GetPosition());
    
    Render();
    
    if (useLinkView)
    {
      sharedData->linkedFrame = event.GetPosition();
      wxCommandEvent changeEvent (EVENT_DATASET_FRAME_CHANGED);
      changeEvent.SetInt (event.GetPosition());
      PostToOthers (changeEvent);
    }
  }

  void ResultViewerImpl::View3D::OnLinkedDataSetChanged (wxCommandEvent& event)
  {
    if (!useLinkView) return;
    
    dataSetSelectionSlider->SetValue (event.GetInt());
    SetDisplayedDataSet (event.GetInt());
    
    Render();
  }

  void ResultViewerImpl::View3D::SetDisplayedDataSet (size_t index)
  {
    currentDataSet = index;
    vtkDataSet* dataset = data->GetDataSet (currentDataSet);
    dataSetMapper.SetInput (dataset);
    dataSetMapper.Update ();
    origDataSetFaces->SetInput (dataset);
    origDataSetMapper->Update ();
    
    if (displacementDirection != ddNone)
    {
      ComputeDisplacementDirections ();
      ShowDisplacementDirections();
    }
    
    // Always invalidate (so it gets updated for new dataset even if _not_ visible)
    displacedData = nullptr;
    if (useDisplaceData)
    {
      ComputeDisplacementOffset ();
      ShowDisplacementOffset ();
    }
  }
  
  void ResultViewerImpl::View3D::HandleMouseMove (int x, int y)
  {
    cellPicker->Pick (x, y, 0, renderer);
    
    std::cerr << "cellId: " << cellPicker->GetCellId() << std::endl;
    std::cerr << "subId: " << cellPicker->GetSubId() << std::endl;
  }

  //-------------------------------------------------------------------------
  
  ResultViewerImpl::View3D::DataSetMapperWithEdges::DataSetMapperWithEdges ()
   : currentRepr (VTK_SURFACE), edgeVisibility (0)
  {
    mapper = vtkSmartPointer<vtkDataSetMapper>::New ();
    
    // @@@ Not really the ideal place
    vtkMapper::SetResolveCoincidentTopologyToPolygonOffset ();
    vtkMapper::SetResolveCoincidentTopologyPolygonOffsetParameters (-0.7, -1);
    vtkMapper::SetResolveCoincidentTopologyPolygonOffsetFaces (0);
    
    faceExtract = vtkSmartPointer<DataSetFaceExtractor>::New ();
    /* We want polygons after clipping, so provide the QuadClippingToPoly-producing 
       output data set up front */
    vtkSmartPointer<UnstructuredGridWithQuadsClippingToPoly> faceExtractOut =
      vtkSmartPointer<UnstructuredGridWithQuadsClippingToPoly>::New ();
    faceExtract->SetOutput (faceExtractOut);
    
    edgesMapper = vtkSmartPointer<vtkDataSetMapper>::New ();
    edgesMapper->SetInput (faceExtract->GetOutput());
  }
  
  vtkActor* ResultViewerImpl::View3D::DataSetMapperWithEdges::CreateMapperActor ()
  {
    mapperActor = vtkSmartPointer<vtkActor>::New ();
    mapperActor->SetMapper (mapper);
    return mapperActor;
  }
  
  vtkActor* ResultViewerImpl::View3D::DataSetMapperWithEdges::CreateEdgesActor ()
  {
    edgesActor = vtkSmartPointer<vtkActor>::New ();
    edgesActor->SetMapper (edgesMapper);
    edgesActor->GetProperty()->SetColor (0, 0, 0);
    edgesActor->GetProperty()->LightingOff();
    edgesActor->GetProperty()->SetRepresentation (VTK_WIREFRAME);
    edgesActor->VisibilityOff();
    return edgesActor;
  }
      
  void ResultViewerImpl::View3D::DataSetMapperWithEdges::SetInput (vtkDataSet* dataset)
  {
    originalInput = dataset;
    faceExtract->SetInput (dataset);
    
    SetClipFunction (lastClipper);
  }
  
  void ResultViewerImpl::View3D::DataSetMapperWithEdges::Update ()
  {
    mapper->Update();
    edgesMapper->Update();
  }

  void ResultViewerImpl::View3D::DataSetMapperWithEdges::GetBounds (double bounds[6])
  {
    mapper->GetBounds (bounds);
  }
  
  void ResultViewerImpl::View3D::DataSetMapperWithEdges::SetClipFunction (vtkImplicitFunction* clipFunc)
  {
    if (clipFunc)
    {
      mapper->SetInput (ClipWrapDataSet (originalInput, clipFunc));
      edgesMapper->SetInput (ClipWrapDataSet (faceExtract->GetOutput (), clipFunc));
    }
    else
    {
      mapper->SetInput (originalInput);
      edgesMapper->SetInput (faceExtract->GetOutput ());
    }
    lastClipper = clipFunc;
  }
      
  void ResultViewerImpl::View3D::DataSetMapperWithEdges::SetColorModeToDefault()
  {
    if (currentRepr == VTK_SURFACE)
      mapper->SetColorModeToDefault();
    else
      edgesMapper->SetColorModeToDefault();
  }
  
  void ResultViewerImpl::View3D::DataSetMapperWithEdges::SetColorModeToMapScalars()
  {
    if (currentRepr == VTK_SURFACE)
      mapper->SetColorModeToMapScalars();
    else
      edgesMapper->SetColorModeToMapScalars();
  }
  
  void ResultViewerImpl::View3D::DataSetMapperWithEdges::SetScalarModeToUseCellFieldData()
  {
    if (currentRepr == VTK_SURFACE)
      mapper->SetScalarModeToUseCellFieldData();
    else
      edgesMapper->SetScalarModeToUseCellFieldData();
  }
  
  void ResultViewerImpl::View3D::DataSetMapperWithEdges::SetScalarModeToUsePointFieldData()
  {
    if (currentRepr == VTK_SURFACE)
      mapper->SetScalarModeToUsePointFieldData();
    else
      edgesMapper->SetScalarModeToUsePointFieldData();
  }
    
  void ResultViewerImpl::View3D::DataSetMapperWithEdges::ScalarVisibilityOff()
  {
    if (currentRepr == VTK_SURFACE)
      mapper->ScalarVisibilityOff();
    else
      edgesMapper->ScalarVisibilityOff();
  }
  
  void ResultViewerImpl::View3D::DataSetMapperWithEdges::ScalarVisibilityOn()
  {
    if (currentRepr == VTK_SURFACE)
      mapper->ScalarVisibilityOn();
    else
      edgesMapper->ScalarVisibilityOn();
  }
    
  void ResultViewerImpl::View3D::DataSetMapperWithEdges::SelectColorArray (int array)
  {
    mapper->SelectColorArray (array);
    edgesMapper->SelectColorArray (array);
  }
  
  void ResultViewerImpl::View3D::DataSetMapperWithEdges::SetLookupTable (vtkScalarsToColors* lut)
  {
    mapper->SetLookupTable (lut);
    edgesMapper->SetLookupTable (lut);
  }
  
  vtkScalarsToColors* ResultViewerImpl::View3D::DataSetMapperWithEdges::GetLookupTable ()
  {
    return mapper->GetLookupTable();
  }
  
  void ResultViewerImpl::View3D::DataSetMapperWithEdges::SetScalarRange (double range[2])
  {
    mapper->SetScalarRange (range);
    edgesMapper->SetScalarRange (range);
  }
  
  void ResultViewerImpl::View3D::DataSetMapperWithEdges::InterpolateScalarsBeforeMappingOn ()
  {
    mapper->InterpolateScalarsBeforeMappingOn();
    edgesMapper->InterpolateScalarsBeforeMappingOn();
  }
  
  void ResultViewerImpl::View3D::DataSetMapperWithEdges::SetRepresentation (int repr)
  {
    if (repr == currentRepr) return;
    
    currentRepr = repr;
    if (repr == VTK_WIREFRAME)
    {
      edgesMapper->SetColorMode (mapper->GetColorMode());
      edgesMapper->SetScalarMode (mapper->GetScalarMode());
      edgesActor->GetProperty()->SetLighting (mapperActor->GetProperty()->GetLighting());
      edgesActor->GetProperty()->SetColor (1, 1, 1);
      mapperActor->VisibilityOff();
      edgesActor->VisibilityOn();
    }
    else
    {
      mapper->SetColorMode (edgesMapper->GetColorMode());
      mapper->SetScalarMode (edgesMapper->GetScalarMode());
      edgesMapper->SetColorModeToDefault();
      mapperActor->GetProperty()->SetLighting (edgesActor->GetProperty()->GetLighting());
      edgesActor->GetProperty()->LightingOff();
      edgesActor->GetProperty()->SetColor (0, 0, 0);
      mapperActor->VisibilityOn();
      edgesActor->SetVisibility (edgeVisibility);
    }
  }
  
  void ResultViewerImpl::View3D::DataSetMapperWithEdges::SetLighting (int lighting)
  {
    if (currentRepr == VTK_SURFACE)
      mapperActor->GetProperty()->SetLighting (lighting);
    else
      edgesActor->GetProperty()->SetLighting (lighting);
  }
  
  void ResultViewerImpl::View3D::DataSetMapperWithEdges::SetEdgeVisibility (int edgeVis)
  {
    if (currentRepr != VTK_SURFACE) return;
    edgesActor->SetVisibility (edgeVis);
    this->edgeVisibility = edgeVis;
  }
  
  ResultViewerImpl::View3D::DataSetMapperWithEdges::operator bool() const
  {
    return mapper;
  }
  
  //-------------------------------------------------------------------------
  
  BEGIN_EVENT_TABLE(ResultViewerImpl::View3D::GradientMenuEventHandler, wxEvtHandler)
    EVT_MENU(wxID_ANY, ResultViewerImpl::View3D::GradientMenuEventHandler::OnGradientSelect)
  END_EVENT_TABLE()
  
  void ResultViewerImpl::View3D::GradientMenuEventHandler::OnGradientSelect (wxCommandEvent& event)
  {
    parent->OnGradientSelect (event);
  }

} // namespace nutogui
