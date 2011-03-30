/**\file
 * Result viewer panel for plots.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"

#include "ViewPlot.h"

#include "AllViewSharedData.h"
#include "Data.h"
#include "SelectedCellsChangedEvent.h"

#include "vtkwx.h"
#include <wx/aui/auibar.h>
#include <boost/foreach.hpp>

#include <vtkAxis.h>
#include <vtkChartXY.h>
#include <vtkColorSeries.h>
#include <vtkContextActor.h>
#include <vtkContextScene.h>
#include <vtkFloatArray.h>
#include <vtkPlot.h>
#include <vtkProperty2D.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkTable.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkUnsignedIntArray.h>

namespace nutogui
{
  enum
  {
    ID_Toolbar = 1,
    ID_DisplayData,
    ID_VisOption
  };
  
  BEGIN_EVENT_TABLE(ResultViewerImpl::ViewPlot, ViewPanelContentVTK)
    EVT_SELECTEDCELLS_CHANGED(ResultViewerImpl::ViewPlot::OnSelectedCellsChanged)
  END_EVENT_TABLE()

  ResultViewerImpl::ViewPlot::ViewPlot (ViewPanel* parent)
    : ViewPanelContentVTK (parent),
      currentDataArray ((size_t)~0),
      currentVisComp (-1)
  {
    wxSizer* sizer = new wxBoxSizer (wxVERTICAL);
    
    vtkwx::RenderWidget* renderWidget = CreateRenderWidget (this);
    renderWidget->EnableKeyboardHandling (false);
    sizer->Add (renderWidget, 1, wxEXPAND);
    
    SetSizer (sizer);
    
    sharedAllData = AllViewSharedData::Setup (parent);
    
    plotColors = vtkSmartPointer<vtkColorSeries>::New ();
    
    chart = vtkSmartPointer<vtkChartXY>::New();
    chart->GetAxis (0)->SetTitle ("");
    chart->GetAxis (1)->SetTitle ("Frame");
    chart->SetShowLegend (true);
    
    contextActor = vtkSmartPointer<vtkContextActor>::New();
    // Set up a 2D scene, add an XY chart to it
    vtkSmartPointer<vtkContextScene> scene = contextActor->GetScene();
    scene->AddItem (chart.GetPointer());
    
    noDataMessage = vtkSmartPointer<vtkTextActor>::New();
    noDataMessage->SetInput ("No cells selected");
    noDataMessage->GetTextProperty()->SetJustificationToCentered();
    noDataMessage->GetTextProperty()->SetVerticalJustificationToCentered();
    noDataMessage->GetPositionCoordinate()->SetCoordinateSystemToNormalizedDisplay();
    noDataMessage->GetPositionCoordinate()->SetValue (0.5, 0.5, 0);
    vtkSmartPointer<vtkProperty2D> textProp2d;
    textProp2d.TakeReference (vtkProperty2D::New ());
    textProp2d->SetColor (0, 0, 0);
    noDataMessage->SetProperty (textProp2d);
  }

  wxWindow* ResultViewerImpl::ViewPlot::CreateTopTools (wxWindow* parentWindow)
  {
    wxPanel* topBar = new wxPanel (parentWindow);
    
    wxSizer* topBarSizer = new wxBoxSizer (wxHORIZONTAL);
    
    toolbar = new wxAuiToolBar (topBar, ID_Toolbar, wxDefaultPosition, wxDefaultSize,
				wxAUI_TB_HORZ_LAYOUT | wxAUI_TB_NO_AUTORESIZE);
    wxChoice* displayDataChoiceCtrl = new wxChoice (toolbar, ID_DisplayData);
    displayDataChoiceCtrl->SetToolTip (wxT ("Data"));
    toolbar->AddControl (displayDataChoiceCtrl);
    displayDataChoice.SetControl (displayDataChoiceCtrl);
    
    visChoiceCtrl = new wxChoice (toolbar, ID_VisOption);
    visChoiceCtrl->SetToolTip (wxT ("Component"));
    toolbar->AddControl (visChoiceCtrl);
    
    topBarSizer->Add (toolbar, wxSizerFlags(1).Expand());
    
    topBar->SetSizer (topBarSizer);
    
    // Connect events
    topBar->Connect (ID_DisplayData,
		     wxEVT_COMMAND_CHOICE_SELECTED,
		     wxCommandEventHandler (ResultViewerImpl::ViewPlot::OnDisplayDataChanged),
		     nullptr, this);
    topBar->Connect (ID_VisOption,
		     wxEVT_COMMAND_CHOICE_SELECTED,
		     wxCommandEventHandler (ResultViewerImpl::ViewPlot::OnVisOptionChanged),
		     nullptr, this);
    
    SetData (sharedAllData->data);

    return topBar;
  }

  void ResultViewerImpl::ViewPlot::SetupVTKRenderer ()
  {
    if (!renderer)
    {
      SetupRenderer();
    }
  }

  void ResultViewerImpl::ViewPlot::SetupRenderer ()
  {
    renderer = vtkSmartPointer<vtkRenderer>::New ();
    GetRenderWidget()->GetRenderWindow()->AddRenderer (renderer);
    renderer->SetBackground (1, 1, 1);
    
    renderer->AddActor (contextActor);
    renderer->AddActor (noDataMessage);
  }

  void ResultViewerImpl::ViewPlot::SetData (const DataConstPtr& data)
  {
    displayDataChoice.Clear();
    
    for (size_t i = 0; i < data->GetNumDataArrays(); i++)
    {
      if (data->GetDataArrayAssociation (i) != Data::perCell) continue;
      displayDataChoice.Append (data->GetDataArrayName (i), i);
    }
    
    lastVisOpt.clear ();
    lastVisOpt.resize (displayDataChoice.GetCount());
    
    if (displayDataChoice.GetCount() == 0)
    {
      displayDataChoice.Append (wxT ("No plottable data"), 0);
      displayDataChoice.GetControl()->Disable ();
      visChoiceCtrl->Hide();
      
      // Display message instead of chart
      contextActor->VisibilityOff();
      noDataMessage->SetInput ("Sorry, no plottable data is available");
      noDataMessage->VisibilityOn();
    }
    else
    {
      // Fill component choice
      currentDataArray = displayDataChoice.GetClientData (0);
      currentVisComp = lastVisOpt[0].comp;
      UpdateVisOptionChoice (currentDataArray, currentVisComp);
    }
    
    displayDataChoice.SetSelection (0);

    // Choice box strings have changed, so update size
    UpdateToolbarControlMinSize (displayDataChoice.GetControl(), toolbar);
    
    SetupChart (currentDataArray, currentVisComp);
  }
  
  static double GetArrayTupleMagnitude (vtkDataArray* array, size_t element)
  {
    double mag = 0;
    double* comps = array->GetTuple (element);
    int numComp = array->GetNumberOfComponents();
    for (int c = 0; c < numComp; c++)
    {
      mag += comps[c]*comps[c];
    }
    return sqrt (mag);
  }

  bool ResultViewerImpl::ViewPlot::SetupChart (size_t dataArray, int viscomp)
  {
    if (dataArray == (size_t)~0) return false;
    
    if (sharedAllData->selectedCellIDs.empty())
    {
      contextActor->VisibilityOff();
      noDataMessage->VisibilityOn();
      return true;
    }
    contextActor->VisibilityOn();
    noDataMessage->VisibilityOff();
    
    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New ();
    
    vtkSmartPointer<vtkUnsignedIntArray> frameNumber = vtkSmartPointer<vtkUnsignedIntArray>::New();
    frameNumber->SetName ("Dataset frame");
    table->AddColumn (frameNumber);
    
    DataConstPtr data (sharedAllData->data);
    
    std::vector<vtkIdType> plotData;
    plotData.insert (plotData.begin(), sharedAllData->selectedCellIDs.begin(), sharedAllData->selectedCellIDs.end ());
    std::sort (plotData.begin(), plotData.end());
    
    BOOST_FOREACH(vtkIdType cellID, plotData)
    {
      vtkSmartPointer<vtkFloatArray> cellData = vtkSmartPointer<vtkFloatArray>::New ();
      wxString dataDescr (wxString::Format (wxT("Cell %lu"), (unsigned long)cellID));
      cellData->SetName (dataDescr.mb_str ());
      table->AddColumn (cellData);
    }
    
    table->SetNumberOfRows (data->GetDataSetNum());
    for (size_t frame = 0; frame < data->GetDataSetNum(); frame++)
    {
      table->SetValue (frame, 0, frame);
      
      vtkDataArray* frameData = data->GetDataArrayRawData (frame, dataArray);
      for (size_t plot = 0; plot < plotData.size(); plot++)
      {
	vtkIdType cellID = plotData[plot];
	float val;
	if (viscomp >= 0)
	  val = frameData->GetComponent (cellID, viscomp);
	else
	  val = GetArrayTupleMagnitude (frameData, cellID);
	table->SetValue (frame, 1+plot, val);
      }
    }

    /* Q: Why not use ClearPlots?
     * A: It's broken (doesn't clear 'PlotCorners' member) */
    while (chart->GetNumberOfPlots() > 0)
    {
      chart->RemovePlot (chart->GetNumberOfPlots()-1);
    }
    for (size_t plot = 0; plot < plotData.size(); plot++)
    {
      vtkPlot* line = chart->AddPlot (vtkChart::LINE);
      line->SetInput (table.GetPointer(), 0, 1+plot);
      line->SetWidth (2.0);
      vtkColor3ub color = plotColors->GetColorRepeating (plotData[plot]);
      line->SetColor (color.Red(), color.Green(), color.Blue(), 255);
    }
    
    return true;
  }
  
  void ResultViewerImpl::ViewPlot::OnSelectedCellsChanged (wxCommandEvent& event)
  {
    if (SetupChart (currentDataArray, currentVisComp))
      Render();
  }
  
  void ResultViewerImpl::ViewPlot::OnDisplayDataChanged (wxCommandEvent& event)
  {
    if (event.GetSelection() < 0) return;
    
    int sel = event.GetSelection();
    currentDataArray = displayDataChoice.GetClientData (sel);
    
    currentVisComp = lastVisOpt[sel].comp;
    UpdateVisOptionChoice (currentDataArray, currentVisComp);
    
    if (SetupChart (currentDataArray, currentVisComp))
      Render();
  }
  
  void ResultViewerImpl::ViewPlot::OnVisOptionChanged (wxCommandEvent& event)
  {
    int displaySel = displayDataChoice.GetSelection();
    
    currentVisComp = event.GetSelection()-1;
    lastVisOpt[displaySel].comp = currentVisComp;
    
    if (SetupChart (currentDataArray, currentVisComp))
      Render();
  }

  void ResultViewerImpl::ViewPlot::UpdateVisOptionChoice (size_t arrayIndex, int initialSel)
  {
    DataConstPtr data (sharedAllData->data);
    visChoiceCtrl->Clear ();
    visChoiceCtrl->Append (data->GetDataArrayComponentDisplayName (arrayIndex, Data::compMagnitude));
    for (int i = 0; i < data->GetDataArrayComponents (arrayIndex); i++)
    {
      visChoiceCtrl->Append (data->GetDataArrayComponentDisplayName (arrayIndex, i));
    }
    visChoiceCtrl->SetSelection (initialSel+1);
    
    /* WX 2.8 computes a slightly different height for visChoiceCtrl,
       so manually force the desired height */
    UpdateToolbarControlMinSize (visChoiceCtrl, toolbar,
				 displayDataChoice.GetControl()->GetBestSize().GetHeight());
  }

} // namespace nutogui
