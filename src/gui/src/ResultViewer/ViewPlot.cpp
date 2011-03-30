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
  BEGIN_EVENT_TABLE(ResultViewerImpl::ViewPlot, ViewPanelContentVTK)
    EVT_SELECTEDCELLS_CHANGED(ResultViewerImpl::ViewPlot::OnSelectedCellsChanged)
  END_EVENT_TABLE()

  ResultViewerImpl::ViewPlot::ViewPlot (ViewPanel* parent)
    : ViewPanelContentVTK (parent)
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

  void ResultViewerImpl::ViewPlot::SetData (const DataConstPtr& data)
  {
    this->data = data;
    SetupChart ();
  }

  wxWindow* ResultViewerImpl::ViewPlot::CreateTopTools (wxWindow* parentWindow)
  {
    return nullptr;
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

  void ResultViewerImpl::ViewPlot::SetupChart ()
  {
    if (sharedAllData->selectedCellIDs.empty())
    {
      contextActor->VisibilityOff();
      noDataMessage->VisibilityOn();
      return;
    }
    contextActor->VisibilityOn();
    noDataMessage->VisibilityOff();
    
    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New ();
    
    vtkSmartPointer<vtkUnsignedIntArray> frameNumber = vtkSmartPointer<vtkUnsignedIntArray>::New();
    frameNumber->SetName ("Dataset frame");
    table->AddColumn (frameNumber);
    
    size_t dataArray = (size_t)~0;
    size_t dataComponent = 0;
    for (size_t i = 0; i < data->GetNumDataArrays(); i++)
    {
      if (data->GetDataArrayAssociation (i) == Data::perCell)
      {
	dataArray = i;
	break;
      }
    }
    
    std::vector<vtkIdType> plotData;
    if (dataArray != (size_t)~0)
    {
      plotData.insert (plotData.begin(), sharedAllData->selectedCellIDs.begin(), sharedAllData->selectedCellIDs.end ());
      std::sort (plotData.begin(), plotData.end());
      
      BOOST_FOREACH(vtkIdType cellID, plotData)
      {
	vtkSmartPointer<vtkFloatArray> cellData = vtkSmartPointer<vtkFloatArray>::New ();
	wxString dataDescr (wxString::Format (wxT("Cell %lu"), (unsigned long)cellID));
	cellData->SetName (dataDescr.mb_str ());
	table->AddColumn (cellData);
      }
    }
    
    table->SetNumberOfRows (data->GetDataSetNum());
    for (size_t frame = 0; frame < data->GetDataSetNum(); frame++)
    {
      table->SetValue (frame, 0, frame);
      
      if (dataArray != (size_t)~0)
      {
	vtkDataArray* frameData = data->GetDataArrayRawData (frame, dataArray);
	for (size_t plot = 0; plot < plotData.size(); plot++)
	{
	  vtkIdType cellID = plotData[plot];
	  float val = frameData->GetComponent (cellID, dataComponent);
	  table->SetValue (frame, 1+plot, val);
	}
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
  }
  
  void ResultViewerImpl::ViewPlot::OnSelectedCellsChanged (wxCommandEvent& event)
  {
    SetupChart();
    Render();
  }
} // namespace nutogui
