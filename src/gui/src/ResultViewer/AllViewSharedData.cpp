/**\file
 * Data shared by pretty much all views.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"

#include "AllViewSharedData.h"

#include "ViewPanel.h"

#include <boost/make_shared.hpp>

namespace nutogui
{
  ResultViewerImpl::AllViewSharedDataPtr ResultViewerImpl::AllViewSharedData::Setup (ViewPanel* panel)
  {
    AllViewSharedDataPtr data (panel->GetSharedData<AllViewSharedData> ());
    if (!data)
    {
      data = boost::make_shared<AllViewSharedData> ();
      panel->SetSharedData (data);
    }
    return data;
  }
} // namespace nutogui
