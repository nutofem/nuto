/**\file
 * Implementation of a Python script runner, running scripts in a separate thread.
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"
#include <boost/python.hpp>

#include "NutoModuleOverlay.h"
#include "ScriptRunnerThreaded.h"
#include "Timing.h"

#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <wx/log.h>

using namespace boost::python;

ScriptRunnerThreaded::ScriptRunnerThreaded ()
 : startupComplete (false), outputNotifyCallbacks (false)
{
  // wxThreadHelper creates a joinable thread
  wxThreadHelper::Create ();
  GetThread()->Run ();
}

ScriptRunnerThreaded::~ScriptRunnerThreaded ()
{
  ThreadStop();
  GetThread()->Wait();
}

void ScriptRunnerThreaded::Result (const wxString& resultTempFile,
				   const wxString& title,
				   const wxString& resultName)
{
  wxCriticalSectionLocker _lock (feedbackProtect);
  
  BOOST_FOREACH(const FeedbackCallback& callback, feedbackCallbacks)
  {
    callback.wrappedCallback->Result (resultTempFile, title, resultName);
  }
}

void ScriptRunnerThreaded::ScriptRunStart ()
{
  wxCriticalSectionLocker _lock (feedbackProtect);

  BOOST_FOREACH(const FeedbackCallback& callback, feedbackCallbacks)
  {
    callback.wrappedCallback->ScriptRunStart ();
  }
}

void ScriptRunnerThreaded::ScriptRunEnd (bool success,
					 const nutogui::ScriptRunner::TracebackPtr& traceback)
{
  wxCriticalSectionLocker _lock (feedbackProtect);

  BOOST_FOREACH(const FeedbackCallback& callback, feedbackCallbacks)
  {
    callback.wrappedCallback->ScriptRunEnd (success, traceback);
  }
}

void ScriptRunnerThreaded::OutputPending (PythonOutputGrabber* source)
{
  if (!outputNotifyCallbacks) return;
  
  wxString outputStr;
  source->GetPendingOutput (outputStr);
  if (outputStr.IsEmpty()) return;
  
  /* Note: this method will be called from main thread, hence the
     originalCallback calls */
  
  wxCriticalSectionLocker _lock (feedbackProtect);
  
  nutogui::ScriptRunner::FeedbackCallback::OutputTarget target =
    (source == grabberStdout) ? nutogui::ScriptRunner::FeedbackCallback::StdOut 
			      : nutogui::ScriptRunner::FeedbackCallback::StdErr;
  BOOST_FOREACH(const FeedbackCallback& callback, feedbackCallbacks)
  {
    callback.originalCallback->ScriptOutput (target, outputStr);
  }
}

static bool HandlePyErr ()
{
  if (!PyErr_ExceptionMatches(PyExc_SystemExit))
  {
    PyErr_Print(); // That shoots the host app down when the exception is PyExc_SystemExit
    return true;
  }
  else
  {
    // Script requested to exit
    PyErr_Clear();
    return false;
  }
}

static void ExtractTraceback (object py_traceback, nutogui::ScriptRunner::Traceback& traceback)
{
  try
  {
    while (!py_traceback.is_none())
    {
      int tb_lineno = extract<int> (py_traceback.attr ("tb_lineno"));
      object tb_frame = py_traceback.attr ("tb_frame");
      object f_code = tb_frame.attr ("f_code");
      std::string co_name = extract<std::string> (f_code.attr ("co_name"));
      std::string co_filename = extract<std::string> (f_code.attr ("co_filename"));
      
      nutogui::ScriptRunner::TracebackEntry tbEntry;
      tbEntry.module = wxString (co_name.c_str(), wxConvLibc);
      tbEntry.filename = wxString (co_filename.c_str(), wxConvLibc);
      tbEntry.line = tb_lineno;
      traceback.push_back (tbEntry);
      
      object tb_next = py_traceback.attr ("tb_next");
      py_traceback = tb_next;
    }
  }
  catch(const error_already_set& e)
  {
    // Bummer.
    HandlePyErr ();
  }
}

void ScriptRunnerThreaded::ActualStartScript (const wxString& source)
{
  bool success = false;
  boost::shared_ptr<nutogui::ScriptRunner::Traceback> traceback;
  ScriptRunStart();
  try
  {
    object main_module = import("__main__");
    object main_namespace = main_module.attr("__dict__");
    
    /* Use the same namespace for both globals and locals
	(so global functions work as expected), but make a
	copy so to not pollute the global NS for other scripts */
    dict script_ns = dict (main_namespace).copy();

    std::string sourceStr (source.mb_str());
    exec (sourceStr.c_str(), script_ns, script_ns);

    success = true;
  }
  catch(const error_already_set& e)
  {
    // HandlePyErr() will clear exception status
    bool isSyntaxError = PyErr_ExceptionMatches (PyExc_SyntaxError);
    if (HandlePyErr ())
    {
      traceback = boost::make_shared<nutogui::ScriptRunner::Traceback> ();
      if (isSyntaxError)
      {
	/* With syntax errors there is usually not a 'regular' trace back.
	   But to allow navigation to the error, a traceback is faked up from the
	   file/line information contained in the exception object. */
	handle<PyObject> last_value_ptr (borrowed (allow_null (PySys_GetObject ("last_value"))));
	if (last_value_ptr)
	{
	  object last_value (last_value_ptr);
	  nutogui::ScriptRunner::TracebackEntry tbEntry;
	  // module: ? - not contained in SyntaxError exception
	  std::string filenameStr (extract<std::string> (last_value.attr ("filename")));
	  tbEntry.filename = wxString (filenameStr.c_str(), wxConvLibc);
	  tbEntry.line = extract<int> (last_value.attr ("lineno"));
	  traceback->push_back (tbEntry);
	}
      }
      else
      {
	handle<PyObject> last_traceback_ptr (borrowed (allow_null (PySys_GetObject ("last_traceback"))));
	if (last_traceback_ptr)
	{
	  object last_traceback (last_traceback_ptr);
	  ExtractTraceback (last_traceback, *traceback);
	}
      }
    }
    else
      success = true; // Script requested to exit
  }
  ScriptRunEnd (success, traceback);
}

void ScriptRunnerThreaded::ActualThreadStop()
{
  runThread = false;
}

void ScriptRunnerThreaded::AddFeedbackCallback (const nutogui::ScriptRunner::FeedbackCallbackPtr& callback)
{
  {
    wxCriticalSectionLocker _lock (feedbackProtect);
    
    nutogui::ScriptRunner::FeedbackCallbackPtr wrappedCallback (
      boost::make_shared<ScriptRunnerFeedbackCallbackWrapper> (callback));
    feedbackCallbacks.push_back (FeedbackCallback (callback, wrappedCallback));
    
    if (!startupComplete)
      return;
  }
  // This is the 'right' thread, call directly
  callback->StartupComplete (startupState, startupMessage);
}

wxThread::ExitCode ScriptRunnerThreaded::Entry()
{
  ThreadStartup();
  
  runThread = true;
  while (runThread)
  {
    TRPCDispatch();
  }
  
  ThreadCleanup();
  return 0;
}

void ScriptRunnerThreaded::ThreadStartup()
{
  nutogui::Timing _timing ("Python thread startup");
  bool startupSuccess = false;
  
  Py_Initialize ();
  /* NOTE: Py_Initialize() may already produce output. To intercept that, we
     would have to switch out the stdout, stderr as used by the Python library
     with out own files/pipes. However, reliably doing that, on all platforms,
     is _very_ hard. */
  
  // Try to import 'nuto' namespace
  try
  {
    // Set up Python output redirection, so we can display script output in a nice UI
    PythonOutputGrabber::BindToPython();
    
    PythonOutputGrabber::OutputPendingNotifierPtr notifyWrap =
      boost::make_shared<OutputPendingNotifier> (this);
    PythonOutputGrabber::OutputPendingNotifierPtr outputNotify =
      boost::make_shared<PythonOutputGrabberOutputPendingNotifierWrapper> (notifyWrap);
    
    grabberStdout = new PythonOutputGrabber (outputNotify);
    grabberStderr = new PythonOutputGrabber (outputNotify);
    {
      object wrap_grabStdout (ptr (grabberStdout));
      object wrap_grabStderr (ptr (grabberStderr));
      PySys_SetObject ("stdout", wrap_grabStdout.ptr());
      PySys_SetObject ("stderr", wrap_grabStderr.ptr());
    }
    
    object main_module = import("__main__");
    object main_namespace = main_module.attr("__dict__");

    static const char import_nuto_code[] =
      // This imports NUTO
      "import nuto\n";
    exec (import_nuto_code, main_namespace, main_namespace);
    
    NutoModuleOverlay::BindToPython();
    {
      NutoModuleOverlay::CallbackPtr callback (
	boost::make_shared<OverlayCallback> (this));
      object overlay = object (NutoModuleOverlay (callback));
      object overlay_ExportVtkDataFile = overlay.attr("ExportVtkDataFile");
      
      object nuto_module = import ("nuto");
      nuto_module.attr("__nutoGuiOverlay__") = overlay;
      object nuto_StructureBase = nuto_module.attr("StructureBase");
      nuto_StructureBase.attr("ExportVtkDataFile") = overlay_ExportVtkDataFile;
    }
    
    startupSuccess = true;
  }
  catch(error_already_set& e)
  {
    HandlePyErr ();
    
    grabberStderr->GetPendingOutput (startupMessage);
    startupMessage.Trim();
    wxLogError (wxT ("Python initialization error: %s"), startupMessage.c_str());
  }
  outputNotifyCallbacks = true;
  
  {
    wxCriticalSectionLocker _lock (feedbackProtect);
    startupComplete = true;
    startupState = startupSuccess;
    
    BOOST_FOREACH(const FeedbackCallback& callback, feedbackCallbacks)
    {
      callback.wrappedCallback->StartupComplete (startupState, startupMessage);
    }
  }
}

void ScriptRunnerThreaded::ThreadCleanup()
{
}
