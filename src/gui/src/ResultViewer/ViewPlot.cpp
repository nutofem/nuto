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

#include "vtkwx.h"

#include <vtkChartXY.h>
#include <vtkContextActor.h>
#include <vtkContextScene.h>
#include <vtkFloatArray.h>
#include <vtkPlot.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkTable.h>

namespace nutogui
{
  ResultViewerImpl::ViewPlot::ViewPlot (ViewPanel* parent)
    : ViewPanelContentVTK (parent)
  {
    wxSizer* sizer = new wxBoxSizer (wxVERTICAL);
    
    vtkwx::RenderWidget* renderWidget = CreateRenderWidget (this);
    renderWidget->EnableKeyboardHandling (false);
    sizer->Add (renderWidget, 1, wxEXPAND);
    
    SetSizer (sizer);
  }

  void ResultViewerImpl::ViewPlot::SetData (const DataConstPtr& data)
  {
  }

  ResultViewerImpl::ViewPanel::Content* ResultViewerImpl::ViewPlot::Clone (ViewPanel* parent)
  {
    return new ViewPlot (parent);
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
    
    vtkSmartPointer<vtkContextActor> actor = vtkSmartPointer<vtkContextActor>::New();
    renderer->AddActor(actor);
    
    // Code below is from VTK TestLinePlot example
    
    // Set up a 2D scene, add an XY chart to it
    vtkSmartPointer<vtkContextScene> scene = actor->GetScene(); // We keep a pointer to this for convenience
    
    vtkSmartPointer<vtkChartXY> chart = vtkSmartPointer<vtkChartXY>::New();
    scene->AddItem(chart.GetPointer());

    // Create a table with some points in it...
    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New ();
    vtkSmartPointer<vtkFloatArray> arrX = vtkSmartPointer<vtkFloatArray>::New();
    arrX->SetName("X Axis");
    table->AddColumn(arrX.GetPointer());
    vtkSmartPointer<vtkFloatArray> arrC = vtkSmartPointer<vtkFloatArray>::New();
    arrC->SetName("Cosine");
    table->AddColumn(arrC.GetPointer());
    vtkSmartPointer<vtkFloatArray> arrS = vtkSmartPointer<vtkFloatArray>::New();
    arrS->SetName("Sine");
    table->AddColumn(arrS.GetPointer());
    vtkSmartPointer<vtkFloatArray> arrS2 = vtkSmartPointer<vtkFloatArray>::New();
    arrS2->SetName("Sine2");
    table->AddColumn(arrS2.GetPointer());
    // Test charting with a few more points...
    int numPoints = 69;
    float inc = 7.5 / (numPoints-1);
    table->SetNumberOfRows(numPoints);
    for (int i = 0; i < numPoints; ++i)
      {
      table->SetValue(i, 0, i * inc);
      table->SetValue(i, 1, cos(i * inc) + 0.0);
      table->SetValue(i, 2, sin(i * inc) + 0.0);
      table->SetValue(i, 3, sin(i * inc) + 0.5);
      }

    // Add multiple line plots, setting the colors etc
    vtkPlot *line = chart->AddPlot(vtkChart::LINE);
    line->SetInput(table.GetPointer(), 0, 1);
    line->SetColor(0, 255, 0, 255);
    line->SetWidth(1.0);
    line = chart->AddPlot(vtkChart::LINE);
    line->SetInput(table.GetPointer(), 0, 2);
    line->SetColor(255, 0, 0, 255);
    line->SetWidth(5.0);
    line = chart->AddPlot(vtkChart::LINE);
    line->SetInput(table.GetPointer(), 0, 3);
    line->SetColor(0, 0, 255, 255);
    line->SetWidth(4.0);
  }
} // namespace nutogui
