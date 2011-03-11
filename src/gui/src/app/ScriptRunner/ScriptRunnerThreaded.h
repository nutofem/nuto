/**\file
 * Implementation of a Python script runner, running scripts in a separate thread.
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __SCRIPTRUNNERTHREADED_H__
#define __SCRIPTRUNNERTHREADED_H__

#include "ScriptRunner.h"

#include "ResultDataSourceVTK.h"
#include "platform/ThreadRPC.h"

#include <wx/string.h>
#include <boost/enable_shared_from_this.hpp>
#include <boost/shared_ptr.hpp>

#include <queue>

#include "NutoModuleOverlay.h"
#include "PythonOutputGrabber.h"

#define SRT_RPC_METHODS					\
  /* Top level: methods */				\
  (TRPC_METHOD (public, StartScript,			\
    /* Arguments (Boost PP array) */			\
    TRPC_METHOD_ARG_LIST(1, (const wxString&))		\
  ))							\
  (TRPC_METHOD (protected, ThreadStop,			\
    TRPC_METHOD_ARG_LIST_EMPTY				\
  ))

class ScriptRunnerThreaded : public wxThreadHelper
{
  // Thread main function
  bool runThread;
  wxThread::ExitCode Entry();
  void ThreadStartup();
  void ThreadCleanup();
  
  TRPC_GENERATE_IMPLEMENTATION_DECLS (Actual, SRT_RPC_METHODS)
  TRPC_GENERATE_DISPATCHER (TRPCDispatch, Actual, SRT_RPC_METHODS)
  
  /**\name NutoModuleOverlay::Callback implementation
   * @{ */
  class OverlayCallback : public NutoModuleOverlay::Callback
  {
    ScriptRunnerThreaded* runner;
  public:
    OverlayCallback (ScriptRunnerThreaded* runner) : runner (runner) {}
  
    void Result (const wxString& resultTempFile,
		 const wxString& title,
		 const wxString& resultName)
    { runner->Result (resultTempFile, title, resultName); }
  };
  void Result (const wxString& resultTempFile,
	       const wxString& title,
	       const wxString& resultName);
  /** @} */
  
  wxCriticalSection feedbackProtect;
  bool startupComplete;
  bool startupState;
  wxString startupMessage;
  struct FeedbackCallback
  {
    nutogui::ScriptRunner::FeedbackCallbackPtr originalCallback;
    nutogui::ScriptRunner::FeedbackCallbackPtr wrappedCallback;
    
    FeedbackCallback (const nutogui::ScriptRunner::FeedbackCallbackPtr& originalCallback,
		      const nutogui::ScriptRunner::FeedbackCallbackPtr& wrappedCallback)
     : originalCallback (originalCallback), wrappedCallback (wrappedCallback) {}
  };
  std::vector<FeedbackCallback> feedbackCallbacks;
  void ScriptRunStart ();
  void ScriptRunEnd (bool success, const nutogui::ScriptRunner::TracebackPtr& traceback);
  
  PythonOutputGrabber* grabberStdout;
  PythonOutputGrabber* grabberStderr;
  
  bool outputNotifyCallbacks;
  
  /**\name PythonOutputGrabber::OutputPendingNotifier implementation
   * @{ */
  class OutputPendingNotifier : public PythonOutputGrabber::OutputPendingNotifier
  {
    ScriptRunnerThreaded* runner;
  public:
    OutputPendingNotifier (ScriptRunnerThreaded* runner) : runner (runner) {}
  
    void OutputPending (PythonOutputGrabber* source)
    { runner->OutputPending (source); }
  };
  void OutputPending (PythonOutputGrabber* source);
  /** @} */
public:
  ScriptRunnerThreaded ();
  ~ScriptRunnerThreaded ();
  
  void AddFeedbackCallback (const nutogui::ScriptRunner::FeedbackCallbackPtr& callback);
  TRPC_GENERATE_STUB_DEFNS (BOOST_PP_EMPTY(), SRT_RPC_METHODS)
};

#undef SRT_RPC_METHODS

#define SRFC_RPC_METHODS					\
  (TRPC_METHOD (public, StartupComplete,			\
    TRPC_METHOD_ARG_LIST(2, (bool, const wxString&))		\
  ))								\
  (TRPC_METHOD (public, Result,					\
    TRPC_METHOD_ARG_LIST(3,					\
      (const wxString&,	const wxString&, const wxString&))	\
  ))								\
  (TRPC_METHOD (public, ScriptOutput,				\
    TRPC_METHOD_ARG_LIST(2, 					\
      (nutogui::ScriptRunner::FeedbackCallback::OutputTarget, 	\
      const wxString&))						\
  ))								\
  (TRPC_METHOD (public, ScriptRunStart,				\
    TRPC_METHOD_ARG_LIST_EMPTY					\
  ))								\
  (TRPC_METHOD (public, ScriptRunEnd,				\
    TRPC_METHOD_ARG_LIST(2,					\
      (bool,							\
      const nutogui::ScriptRunner::TracebackPtr&))		\
  ))

class ScriptRunnerFeedbackCallbackWrapper :
  public nutogui::ScriptRunner::FeedbackCallback,
  public platform::ThreadRPC::MainThreadRunner<ScriptRunnerFeedbackCallbackWrapper>
{
  nutogui::ScriptRunner::FeedbackCallbackPtr callback;
  
  friend class platform::ThreadRPC::MainThreadRunner<ScriptRunnerFeedbackCallbackWrapper>;
  TRPC_GENERATE_DISPATCHER_MAIN_FWD (callback->, BOOST_PP_EMPTY(), SRFC_RPC_METHODS)
public:
  ScriptRunnerFeedbackCallbackWrapper (const nutogui::ScriptRunner::FeedbackCallbackPtr& callback)
   : callback (callback) {}
   
  TRPC_GENERATE_STUB_DEFNS (BOOST_PP_EMPTY(), SRFC_RPC_METHODS)
};

#undef SRFC_RPC_METHODS

#define POGOPN_RPC_METHODS					\
  (TRPC_METHOD (public, OutputPending,				\
    TRPC_METHOD_ARG_LIST(1, (PythonOutputGrabber*))		\
  ))

class PythonOutputGrabberOutputPendingNotifierWrapper :
  public PythonOutputGrabber::OutputPendingNotifier,
  public platform::ThreadRPC::MainThreadRunner<PythonOutputGrabberOutputPendingNotifierWrapper>
{
  PythonOutputGrabber::OutputPendingNotifierPtr callback;
  
  friend class platform::ThreadRPC::MainThreadRunner<PythonOutputGrabberOutputPendingNotifierWrapper>;
  TRPC_GENERATE_DISPATCHER_MAIN_FWD (callback->, BOOST_PP_EMPTY(), POGOPN_RPC_METHODS)
public:
  PythonOutputGrabberOutputPendingNotifierWrapper (const PythonOutputGrabber::OutputPendingNotifierPtr& callback)
   : callback (callback) {}
   
  TRPC_GENERATE_STUB_DEFNS (BOOST_PP_EMPTY(), POGOPN_RPC_METHODS)
};

#undef POGOPN_RPC_METHODS

#endif // __SCRIPTRUNNERTHREADED_H__
