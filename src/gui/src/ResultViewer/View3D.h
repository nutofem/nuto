/**\file
 * Result viewer 3D view panel.
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __NUTOGUI_RESULTVIEWER_VIEW3D_H__
#define __NUTOGUI_RESULTVIEWER_VIEW3D_H__

#include "ResultViewerImpl.h"
#include "ViewPanelContentVTK.h"

#include "uicommon/ControlWithItemsClientDataWrapper.h"
#include "uicommon/TextCtrlBuddySlider.h"

#include <vector>

#include <vtkSmartPointer.h>

class vtkAbstractArray;
class vtkActor;
class vtkCamera;
class vtkCellPicker;
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
  
  class ResultViewerImpl::View3D : public ViewPanelContentVTK
  {
    class CameraModifiedCallback;
    
    class Gradients;
    
    wxSizer* topBarSizer;
    wxAuiToolBar* toolbar;
    wxAuiToolBar* actorOptionsTB;
    enum { numRenderModes = 5 };
    struct SharedViewData;
    boost::shared_ptr<SharedViewData> sharedData;
    void SetupSharedData ();
    
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
    
    wxPanel* topBar;
    wxPanel* dataSetSelectionBar;
    uicommon::TextCtrlBuddySlider* dataSetSelectionSlider;
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
    
    vtkSmartPointer<vtkCellPicker> cellPicker;
    
    void OnWindowCreate (wxWindowCreateEvent& event);
    void SetupVTKRenderer ();
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
    
    bool useLinkView;
    /**
     * Change whether the current view should be linked to all others.
     * \return Whether the current camera view was changed (usually after linking was enabled).
     */
    bool SetLinkView (bool flag);

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
    void OnLinkViews (wxCommandEvent& event);
    void OnLinkViewsUpdateUI (wxUpdateUIEvent& event);
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
    
    class UpdateCameraEvent;
    typedef void (ResultViewerImpl::View3D::*UpdateCameraEventFunction)(UpdateCameraEvent&);
    void OnUpdateCamera (UpdateCameraEvent& event);
    
    void CheckDisplacementSizePanelVisibility ();
    void OnDisplacementScaleChange (wxEvent& event);
    
    void OnDataSetSelectionChanged (wxScrollEvent& event);
    void OnLinkedDataSetChanged (wxCommandEvent& event);
    void SetDisplayedDataSet (size_t index);
    
    class RenderViewMouseCallback;
    void HandleMouseMove (int x, int y);
    
    /**
     * Class to handle commands from gradient menu.
     * (Don't use event handling in View since gradient menu entries have
     * "special" IDs, clashing with regular IDs.)
     */
    class GradientMenuEventHandler : public wxEvtHandler
    {
      View3D* parent;
    public:
      GradientMenuEventHandler (View3D* parent) : parent (parent) {}
      
      void OnGradientSelect (wxCommandEvent& event);
			      
      DECLARE_EVENT_TABLE()
    };
    GradientMenuEventHandler gradientMenuHandler;
  public:
    View3D (ViewPanel* parent, const View3D* cloneFrom = nullptr);
    ~View3D ();
    
    void SetData (const DataConstPtr& data);
    wxWindow* CreateTopTools (wxWindow* parentWindow);
			      
    DECLARE_EVENT_TABLE()
  };
} // namespace nutogui

#endif // __NUTOGUI_RESULTVIEWER_VIEW3D_H__
