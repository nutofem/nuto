/**\file
 * Result viewer view split management.
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"
#include "SplitManager.h"

#include <wx/splitter.h>

#include <boost/make_shared.hpp>

namespace nutogui
{
  ResultViewerImpl::SplitManager::SplitManager (wxWindow* parent)
   : wxPanel (parent), maximizedView (nullptr)
  {
  }

  void ResultViewerImpl::SplitManager::SetInitialView (wxWindow* firstView)
  {
    assert(!splitTree);
    
    wxSizer* sizer = new wxBoxSizer (wxVERTICAL);
    sizer->Add (firstView, wxSizerFlags (1).Expand ());
    SetSizer (sizer);

    splitTree = boost::make_shared<SplitNode> (SplitNodePtr (), firstView);
  }

  void ResultViewerImpl::SplitManager::SplitCommon (wxWindow* view, SplitDirection direction,
						    wxWindow* newView)
  {
    SplitNodePtr node = FindNode (splitTree, view);
    assert (!node->children[0]);
    assert (!node->children[1]);
    
    SplitNodePtr nodeParent (node->parent.lock());
    wxWindow* splitterParent = nodeParent ? nodeParent->window : this;
    wxSplitterWindow* newSplitter = new wxSplitterWindow (splitterParent);
    newSplitter->SetSashGravity (0.5);
    view->Reparent (newSplitter);
    if (nodeParent)
    {
      static_cast<wxSplitterWindow*> (nodeParent->window)->ReplaceWindow (
	view, newSplitter);
    }
    newView->Reparent (newSplitter);
    if (direction == splitV)
      newSplitter->SplitHorizontally (view, newView);
    else
      newSplitter->SplitVertically (view, newView);
    
    if (splitterParent == this)
    {
      GetSizer()->Replace (view, newSplitter);
      GetSizer()->Layout ();
    }
    
    node->window = newSplitter;
    node->children[0] = boost::make_shared<SplitNode> (node, view);
    node->children[1] = boost::make_shared<SplitNode> (node, newView);
    
    RedistributeSize (node);
  }
  
  void ResultViewerImpl::SplitManager::RedistributeSize (const SplitNodePtr& node)
  {
    assert (node->children[0] && node->children[1]);
    
    wxSplitterWindow* splitter = static_cast<wxSplitterWindow*> (node->window);
    int splitDir = splitter->GetSplitMode ();
    SplitNodePtr topNode (node);
    {
      SplitNodePtr parentNode (node->parent.lock());
      while (parentNode)
      {
	wxSplitterWindow* parentSplitter = static_cast<wxSplitterWindow*> (parentNode->window);
	if (parentSplitter->GetSplitMode() != splitDir)
	  break;
	topNode = parentNode;
	parentNode = topNode->parent.lock();
      }
    }
    std::pair<size_t, int> countAndSize (CountAndSizeInDirection (topNode, splitDir));
    if ((countAndSize.first == 0) || (countAndSize.second == 0)) return;
    int newViewSize = countAndSize.second / countAndSize.first;
    SetSizeInDirection (topNode, splitDir, newViewSize);
  }
  
  std::pair<size_t, int>
  ResultViewerImpl::SplitManager::CountAndSizeInDirection (const SplitNodePtr& node,
							   int splitDir)
  {
    if (node->IsLeaf ())
    {
      if (splitDir == wxSPLIT_VERTICAL)
	return std::make_pair (1, node->window->GetSize().GetWidth());
      else
	return std::make_pair (1, node->window->GetSize().GetHeight());
    }
    
    wxSplitterWindow* splitter = static_cast<wxSplitterWindow*> (node->window);
    if (splitter->GetSplitMode() == splitDir)
    {
      std::pair<size_t, int> countAndSize1 (CountAndSizeInDirection (node->children[0], splitDir));
      std::pair<size_t, int> countAndSize2 (CountAndSizeInDirection (node->children[1], splitDir));
      return std::make_pair (countAndSize1.first + countAndSize2.first,
			     countAndSize1.second + countAndSize2.second);
    }
    else
      return std::make_pair (0, 0);
  }
  
  int ResultViewerImpl::SplitManager::ComputeSizeInDirection (const SplitNodePtr& node,
							      int splitDir,
							      int viewSize)
  {
    if (node->IsLeaf ())
      return viewSize;
    
    wxSplitterWindow* splitter = static_cast<wxSplitterWindow*> (node->window);
    if (splitter->GetSplitMode() != splitDir)
    {
      if (splitDir == wxSPLIT_VERTICAL)
	return splitter->GetSize().GetWidth();
      else
	return splitter->GetSize().GetHeight();
    }
    
    return ComputeSizeInDirection (node->children[0], splitDir, viewSize)
	 + ComputeSizeInDirection (node->children[1], splitDir, viewSize)
	 + splitter->GetSashSize ();
  }

  void ResultViewerImpl::SplitManager::SetSizeInDirection (const SplitNodePtr& node,
							   int splitDir,
							   int viewSize)
  {
    if (node->IsLeaf ())
      return;
    
    wxSplitterWindow* splitter = static_cast<wxSplitterWindow*> (node->window);
    if (splitter->GetSplitMode() != splitDir) return;
    
    splitter->SetSashPosition (ComputeSizeInDirection (node->children[0], splitDir, viewSize));
    
    SetSizeInDirection (node->children[0], splitDir, viewSize);
    SetSizeInDirection (node->children[1], splitDir, viewSize);
  }

  bool ResultViewerImpl::SplitManager::CanSplit (wxWindow* view)
  {
    return !maximizedView;
  }

  void ResultViewerImpl::SplitManager::SplitVertically (wxWindow* view, wxWindow* newView)
  {
    SplitCommon (view, splitV, newView);
  }
  
  void ResultViewerImpl::SplitManager::SplitHorizontally (wxWindow* view, wxWindow* newView)
  {
    SplitCommon (view, splitH, newView);
  }

  bool ResultViewerImpl::SplitManager::CanUnsplit (wxWindow* view)
  {
    SplitNodePtr node = FindNode (splitTree, view);
    SplitNodePtr nodeParent (node->parent.lock());
    return nodeParent;
  }
  
  void ResultViewerImpl::SplitManager::Unsplit (wxWindow* view)
  {
    if (maximizedView == view) ToggleMaximization (view);
    
    SplitNodePtr node = FindNode (splitTree, view);
    assert (node);
    SplitNodePtr nodeParent (node->parent.lock());
    if (!nodeParent) return;
    
    // Determine parent splitter window of view
    wxSplitterWindow* splitter = static_cast<wxSplitterWindow*> (nodeParent->window);
    
    // Determine sibling of view to close
    wxWindow* other;
    if (nodeParent->children[0]->window == view)
      other = nodeParent->children[1]->window;
    else
      other = nodeParent->children[0]->window;
    assert (other);
    splitter->Unsplit (other);
    // Unsplit() hides other
    other->Show ();
    
    SplitNodePtr nodeGrandParent (nodeParent->parent.lock());
    if (nodeGrandParent)
    {
      // Splitter is child of some other splitter
      wxSplitterWindow* splitterParent =
	static_cast<wxSplitterWindow*> (nodeGrandParent->window);
      other->Reparent (splitterParent);
      splitterParent->ReplaceWindow (splitter, other);
    }
    else
    {
      // Splitter is child of panel
      other->Reparent (this);
      GetSizer()->Replace (splitter, other);
      GetSizer()->Layout ();
    }
    
    // Update split node
    if (nodeParent->children[0] == node)
    {
      nodeParent->window = nodeParent->children[1]->window;
      nodeParent->children[0] = nodeParent->children[1]->children[0];
      nodeParent->children[1] = nodeParent->children[1]->children[1];
    }
    else
    {
      nodeParent->window = nodeParent->children[0]->window;
      nodeParent->children[1] = nodeParent->children[0]->children[1];
      nodeParent->children[0] = nodeParent->children[0]->children[0];
    }
    if (nodeParent->children[0]) nodeParent->children[0]->parent = nodeParent;
    if (nodeParent->children[1]) nodeParent->children[1]->parent = nodeParent;
    
    splitter->Destroy (); // will also destroy 'view' since it's a child
  }

  bool ResultViewerImpl::SplitManager::CanToggleMaximization (wxWindow* view)
  {
    return !splitTree->IsLeaf();
  }
  
  bool ResultViewerImpl::SplitManager::ToggleMaximization (wxWindow* view)
  {
    assert ((view == maximizedView) || !maximizedView);
    
    SplitNodePtr node = FindNode (splitTree, view);
    assert (node);
    SplitNodePtr nodeParent (node->parent.lock());
    assert (nodeParent);
    
    wxSplitterWindow* splitter = static_cast<wxSplitterWindow*> (nodeParent->window);
    
    if (maximizedView == view)
    {
      view->Reparent (splitter);
      
      wxWindow* splitChild1;
      wxWindow* splitChild2;
      if (view == nodeParent->children[0]->window)
      {
	splitChild1 = view;
	splitChild2 = splitter->GetWindow1 ();
      }
      else
      {
	splitChild1 = splitter->GetWindow1 ();
	splitChild2 = view;
      }
      if (splitter->GetSplitMode() == wxSPLIT_VERTICAL)
	splitter->SplitVertically (splitChild1, splitChild2,
				   oldSplitSashFrac * splitter->GetSize().GetWidth ());
      else
	splitter->SplitHorizontally (splitChild1, splitChild2,
				   oldSplitSashFrac * splitter->GetSize().GetHeight ());
      // Make split tree root the 'main' child again
      GetSizer()->Replace (view, splitTree->window);
      GetSizer()->Layout ();
      splitTree->window->Show ();
      
      maximizedView = nullptr;
      return false;
    }
    else
    {
      oldSplitSashFrac = splitter->GetSashPosition();
      if (splitter->GetSplitMode() == wxSPLIT_VERTICAL)
	oldSplitSashFrac *= 1.0 / splitter->GetSize().GetWidth();
      else
	oldSplitSashFrac *= 1.0 / splitter->GetSize().GetHeight();
      // Remove view from parent splitter, make top level child
      splitter->Unsplit (view);
      view->Reparent (this);
      
      // Replace split tree root with view as the 'main' child
      GetSizer()->Replace (splitTree->window, view);
      splitTree->window->Hide ();
      view->Show ();
      GetSizer()->Layout ();
      
      maximizedView = view;
      return true;
    }
  }
  
  ResultViewerImpl::SplitManager::SplitNodePtr
  ResultViewerImpl::SplitManager::FindNode (const SplitNodePtr& tree, wxWindow* window)
  {
    if (tree->window == window)
    {
      return tree;
    }
   
    if (tree->children[0])
    {
      assert (tree->children[1]);
      
      SplitNodePtr child1 = FindNode (tree->children[0], window);
      if (child1)
      {
	return child1;
      }
      SplitNodePtr child2 = FindNode (tree->children[1], window);
      if (child2)
      {
	return child2;
      }
    }
    return SplitNodePtr ();
  }
  
  void ResultViewerImpl::SplitManager::PostToOthers (wxEvent& event, wxWindow* source)
  {
    SplitNodePtr node = FindNode (splitTree, source);
    assert (node);
    PostToOthersUp (event, node);
  }
  
  void ResultViewerImpl::SplitManager::PostToOthersUp (wxEvent& event, const SplitNodePtr& node)
  {
    SplitNodePtr nodeParent (node->parent.lock());
    if (!nodeParent) return;

    if (nodeParent->children[0] == node)
    {
      PostToOthersDown (event, nodeParent->children[1]);
    }
    else
    {
      PostToOthersDown (event, nodeParent->children[0]);
    }
    
    PostToOthersUp (event, nodeParent);
  }
  
  void ResultViewerImpl::SplitManager::PostToOthersDown (wxEvent& event, const SplitNodePtr& node)
  {
    if (node->IsLeaf ())
    {
      node->window->GetEventHandler()->ProcessEvent (event);
    }
    else
    {
      PostToOthersDown (event, node->children[0]);
      PostToOthersDown (event, node->children[1]);
    }
  }
  
} // namespace nutogui
