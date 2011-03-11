/**\file
 * Result viewer view panel.
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __NUTOGUI_RESULTVIEWER_VIEW_H__
#define __NUTOGUI_RESULTVIEWER_VIEW_H__

#include "ResultViewerImpl.h"

#include "uicommon/ControlWithItemsClientDataWrapper.h"

#include <vector>

#include <vtkSmartPointer.h>

class vtkAbstractArray;
class vtkActor;
class vtkCamera;
class vtkDataArray;
class vtkDataSet;
class vtkDataSetMapper;
class vtkGlyph3D;
class vtkImplicitFunction;
class vtkLookupTable;
class vtkOrientationMarkerWidget;
class vtkPlane;
class vtkPlaneWidget;
class vtkPolyData;
class vtkPolyDataMapper;
class vtkProperty;
class vtkRenderer;
class vtkScalarBarWidget;
class vtkScalarsToColors;

class wxAuiToolBar;
class wxAuiToolBarEvent;

namespace nutogui
{
  class DataSetEdgeExtractor;
  class DataSetFaceExtractor;
  class DisplacementDirectionSizePanel;
  
  class ResultViewerImpl::View : public wxPanel
  {
    class CameraModifiedCallback;
    
    SplitManager* splitMgr;
    
    wxSizer* topBarSizer;
    wxAuiToolBar* toolbar;
    wxAuiToolBar* actorOptionsTB;
    wxAuiToolBar* closeMaxButtonsBar;
    enum { numRenderModes = 5 };
    struct SharedViewData;
    boost::shared_ptr<SharedViewData> sharedData;
    void SetupSharedData ();
    void SetupGradients ();
    
    wxAuiToolBar* visOptionEmpty;
    wxAuiToolBar* visOptionChoice;
    wxChoice* visChoiceCtrl;
    int renderMode;
    // Keep in sync with menu!
    enum DisplacementDir
    {
      ddNone,
      ddScaled,
      ddColored,
      numDisplacementDirModes
    };
    DisplacementDir displacementDirection;
    DisplacementDir oldDisplacementDirection;
    
    DisplacementDirectionSizePanel* displacementSizePanel;
    
    class RenderWidget;
    RenderWidget* renderWidget;
    wxPanel* dataSetSelectionBar;
    wxSlider* dataSetSelectionSlider;
    vtkSmartPointer<vtkRenderer> renderer;
    struct DataSetMapperWithEdges
    {
      DataSetMapperWithEdges ();
      
      vtkActor* CreateMapperActor ();
      vtkActor* CreateEdgesActor ();
      
      void SetInput (vtkDataSet* dataset);
      void Update ();
      void GetBounds (double bounds[6]);
      
      void SetClipFunction (vtkImplicitFunction* clipFunc);
      
      void SetColorModeToDefault();
      void SetColorModeToMapScalars();
      
      void SetScalarModeToUseCellFieldData();
      void SetScalarModeToUsePointFieldData();
      
      void ScalarVisibilityOff();
      void ScalarVisibilityOn();
      
      void SelectColorArray (int array);
      
      void SetLookupTable (vtkScalarsToColors* lut);
      vtkScalarsToColors* GetLookupTable ();
      void SetScalarRange (double range[2]);
      
      void InterpolateScalarsBeforeMappingOn ();
      
      void SetRepresentation (int repr);
      void SetLighting (int lighting);
      
      void SetEdgeVisibility (int edgeVis);
      
      operator bool() const;
    private:
      vtkSmartPointer<vtkDataSetMapper> mapper;
      vtkSmartPointer<vtkDataSetMapper> edgesMapper;
      vtkSmartPointer<DataSetFaceExtractor> faceExtract;
      
      vtkSmartPointer<vtkActor> mapperActor;
      vtkSmartPointer<vtkActor> edgesActor;
      
      vtkSmartPointer<vtkDataSet> originalInput;
      vtkSmartPointer<vtkImplicitFunction> lastClipper;
      
      int currentRepr;
      int edgeVisibility;
    };
    DataSetMapperWithEdges dataSetMapper;
    vtkSmartPointer<DataSetFaceExtractor> origDataSetFaces;
    vtkSmartPointer<DataSetEdgeExtractor> origDataSetEdges;
    vtkSmartPointer<vtkDataSetMapper> origDataSetMapper;
    vtkSmartPointer<vtkActor> origDataSetActor;
    vtkSmartPointer<vtkLookupTable> displacementColors;
    vtkSmartPointer<vtkDataSet> displacedData;
    bool useDisplaceData;
    size_t displacementData;
    vtkSmartPointer<vtkPolyData> displaceDirectionsData;
    size_t displaceDirectionsDataDS;
    vtkSmartPointer<vtkGlyph3D> displaceDirectionsGlyphs;
    vtkSmartPointer<vtkPolyDataMapper> displaceDirectionsMapper;
    vtkSmartPointer<vtkActor> displaceDirectionsActor;
    float displacementDirScale;
    
    vtkSmartPointer<vtkScalarBarWidget> scalarBar;
    vtkSmartPointer<vtkOrientationMarkerWidget> orientationMarker;
    vtkSmartPointer<vtkPlaneWidget> clipPlaneWidget;
    vtkSmartPointer<vtkPlane> clipPlane;
    bool useClipper;
    
    vtkSmartPointer<vtkCamera> camTemplate;
    
    void OnRenderWidgetRealized (wxCommandEvent& event);
    void OnWindowCreate (wxWindowCreateEvent& event);
    void SetupRenderer ();
    bool UpdateCameraPositions (vtkCamera* cam);
    vtkCamera* GetCamera () const;
    
    bool updatingCam;
    void CameraChanged (vtkCamera* cam);

    class ClipPlaneChangedCallback;
    void ClipPlaneChanged ();
    
    DataConstPtr data;
    size_t currentDataSet;
    wxChoice* displayDataChoice;
    struct VisOpt
    {
      int comp;
      bool legend;
      size_t gradient;
      
      VisOpt() : comp (-1), legend (true), gradient (0) {}
    };
    std::vector<VisOpt> lastVisOpt;

    void OnDisplayDataChanged (wxCommandEvent& event);

    void OnVisOptionChanged (wxCommandEvent& event);
    void UpdateVisOptionChoice (size_t dataIndex, int initialSel);
    void SetVisComponent (size_t dataIndex, int visComp);
    void UpdateGradientUI (size_t gradient);
    void SetGradient (size_t gradient);

    void OnRenderModeDropDown (wxAuiToolBarEvent& event);
    void OnRenderModeCommand (wxCommandEvent& event);
    void ApplyRenderMode ();
    void SetUIRenderMode ();
    void OnShowLegend (wxCommandEvent& event);
    void OnShowLegendUpdateUI (wxUpdateUIEvent& event);
    void OnLegendOptions (wxAuiToolBarEvent& event);
    void OnLegendOptionsUpdateUI (wxUpdateUIEvent& event);
    void OnGradientSelect (wxCommandEvent& event);
    void OnClipPlane (wxCommandEvent& event);
    void OnDisplacementOffset (wxCommandEvent& event);
    void OnDisplacementDirDropDown (wxAuiToolBarEvent& event);
    void OnDisplacementDirCommand (wxCommandEvent& event);
    void ShowDisplacementOffset ();
    void HideDisplacementOffset ();
    void ComputeDisplacementOffset ();
    void ShowDisplacementDirections ();
    void HideDisplacementDirections ();
    void ComputeDisplacementDirections ();
    
    void DoToolDropDown (wxAuiToolBar* toolbar, int toolID, wxMenu* menu);
    void UpdateToolbarControlMinSize (wxWindow* control,
				      wxAuiToolBar* toolbar,
				      int forceHeight = 0);
				      
    void OnSplitHorizontally (wxCommandEvent& event);
    void OnSplitVertically (wxCommandEvent& event);
    void OnSplitUpdateUI (wxUpdateUIEvent& event);
    void OnUnsplit (wxCommandEvent& event);
    void OnUnsplitUpdateUI (wxUpdateUIEvent& event);
    void OnToggleMaximization (wxCommandEvent& event);
    void OnToggleMaximizationUpdateUI (wxUpdateUIEvent& event);
    
    class UpdateCameraEvent;
    typedef void (ResultViewerImpl::View::*UpdateCameraEventFunction)(UpdateCameraEvent&);
    void OnUpdateCamera (UpdateCameraEvent& event);
    
    void CheckDisplacementSizePanelVisibility ();
    void OnDisplacementScaleChange (wxEvent& event);
    
    void OnDataSetSelectionChanged (wxScrollEvent& event);
    void SetDisplayedDataSet (size_t index);
    
    /**
     * Class to handle commands from gradient menu.
     * (Don't use event handling in View since gradient menu entries have
     * "special" IDs, clashing with regular IDs.)
     */
    class GradientMenuEventHandler : public wxEvtHandler
    {
      View* parent;
    public:
      GradientMenuEventHandler (View* parent) : parent (parent) {}
      
      void OnGradientSelect (wxCommandEvent& event);
			      
      DECLARE_EVENT_TABLE()
    };
    GradientMenuEventHandler gradientMenuHandler;
  public:
    View (wxWindow* parent, SplitManager* splitMgr,
	  const View* cloneFrom = nullptr);
    
    void SetData (const DataConstPtr& data);
			      
    DECLARE_EVENT_TABLE()
  };
} // namespace nutogui

#endif // __NUTOGUI_RESULTVIEWER_VIEW_H__
