/**\file
 * Collection of macros to automatically generate code for cross-thread calls.
 * Methods are not really executed in the calling thread but asynchronously in
 * some other thread. The macros in this file automatically generate the code
 * necessary to marshal method arguments, to maintain a FIFO of methods to call
 * on the 'remote' thread, and to pop off items of this FIFO and execute the
 * attached actual method implementation.
 *
 * Visible to, and callable from, the outside are the so-called stubs.
 * These methods appear the same as the actual implementation - they have the
 * same signature and take their names. The stubs can be called from any thread;
 * they merely pack ("marshal") the method arguments and signal the worker
 * thread that it needs to run the actual implementation.
 *
 * The actual implementations of methods are usually not visible to the
 * outside. (They also need to have a different name to prevent conflicts
 * with the stubs.) They are always executed on the worker thread.
 *
 * On the worker thread, the "dispatcher" manages to collect the requests
 * for asynchronous execution and call the actual implementations.
 *
 * The stub definitionss, actual implementation declarations and dispatcher code
 * are automatically generated through macros.
 * (The macros are based upon Boost.Preprocessor, please have a look at it's
 * documentation when there is talk about "Boost PP something".)
 *
 * The actual method implementations can be written like regular implementations
 * (the only limitation is that it can't be inline, since a forward declaration
 * is made in the class body).
 * At runtime, the dispatch function needs to be called regularly so the methods
 * are actually called. Typically, in a worker thread, this is a simple loop.
 * (Remember to provide a way to break execution of that loop...)
 *
 * The methods that should be callable are defined as a Boost PP sequence.
 * The sequence must consists of TRPC_METHOD entries, which specify the method's
 * visibility (public, protected...), it's name and it's arguments.
 * The arguments are specified as a Boost PP array of argument types
 * (the TRPC_METHOD_ARG_LIST and TRPC_METHOD_ARG_LIST_EMPTY are helpful here).
 * Note there is no return type of the methods. This would make only limited sense
 * since the methods are to be executed asynchronously, so a return type would
 * not be available anyway. To return data from a function, use a different
 * method (such as a callback).
 *
 * In the class declaration that should contain the methods, you need to invoke these macros,
 * in order:
 * 	TRPC_GENERATE_IMPLEMENTATION_DECLS ([prefix], [methods])
 *		[prefix] is a prefix for the names of the actual implementations.
 *		This is usually not empty, to prevent conflicts with the names
 *		of the stubs.
 *		[methods] is methods sequence, described above.
 *		This is usually invoked in private visibility.
 *	TRPC_GENERATE_DISPATCHER ([dispatch-func], [impl-prefix], [methods])
 *		[dispatch-func] is the name of the generated dispatcher function.
 *		[impl-prefix] is the same prefix as used for the implementation
 *		This is usually invoked in private visibility.
 *		declarations.
 *	TRPC_GENERATE_STUB_DEFNS ([prefix], [methods])
 *		[prefix] is a prefix for the names of the stub methods.
 *		This is usually BOOST_PP_EMPTY().
 *		This is usually invoked in public visibility.
 *
 * Should the method execution happen on the main thread, two things must be done:
 * 1. The class implementing the methods must be derived from platform::ThreadRPC::MainThreadRunner<>.
 *    platform::ThreadRPC::MainThreadRunner<> must also be a friend.
 * 2. The TRPC_GENERATE_DISPATCHER_MAIN[_FWD] variants of the dispatcher function generation
 *    must be used.
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __NUTOGUI_THREADRPC_H__
#define __NUTOGUI_THREADRPC_H__

#include <boost/make_shared.hpp>
#include <boost/preprocessor/array.hpp>
#include <boost/preprocessor/list.hpp>
#include <boost/preprocessor/seq.hpp>
#include <boost/type_traits.hpp>

// RPC methods sequence: method entry
#define TRPC_METHOD(Visibility, Name, Args)	(Name, Args, Visibility)
// RPC methods sequence: argument list
#define TRPC_METHOD_ARG_LIST(N, TypesTuple)	(N, TypesTuple)
// RPC methods sequence: empty list
#define TRPC_METHOD_ARG_LIST_EMPTY		(0, BOOST_PP_EMPTY())

// Extract method name
#define _TRPC_METHOD_NAME(Method)					\
  BOOST_PP_TUPLE_ELEM (3, 0, Method)
// Extract method args
#define _TRPC_METHOD_ARGS(Method)					\
  BOOST_PP_TUPLE_ELEM (3, 1, Method)
// Extract method visibility
#define _TRPC_METHOD_VIS(Method)					\
  BOOST_PP_TUPLE_ELEM (3, 2, Method)
// Generate argument list for method
#define _TRPC_GENERATE_METHOD_ARG(z, n, Args)				\
  BOOST_PP_COMMA_IF(n) BOOST_PP_ARRAY_ELEM(n, Args) BOOST_PP_CAT(A, n)
#define _TRPC_METHOD_GENERATE_ARGS(Method)				\
  BOOST_PP_REPEAT (							\
    BOOST_PP_ARRAY_SIZE (_TRPC_METHOD_ARGS(Method)),			\
    _TRPC_GENERATE_METHOD_ARG, _TRPC_METHOD_ARGS(Method))
// Type name for marshalled method args
#define _TRPC_MARSHAL_TYPE_NAME(Method)					\
  BOOST_PP_CAT (_TRPCMArgs, _TRPC_METHOD_NAME (Method))

/*
 * Method stubs (callable from any thread)
 */
// Generate list of method arguments
#define _TRPC_GENERATE_STUB_MARSHAL_ARGS(z, n, unused)			\
  BOOST_PP_COMMA_IF(n) BOOST_PP_CAT(A, n)
// Generate RPC method ID 
#define _TRPC_GENERATE_METHOD_ID(r, Unused, i, Method)				\
  BOOST_PP_COMMA_IF(i) BOOST_PP_CAT (_trpcCall, _TRPC_METHOD_NAME (Method))
// Generate one stub
#define _TRPC_GENERATE_METHOD_STUB(r, Prefix, Method)					\
  _TRPC_METHOD_VIS(Method):								\
  void BOOST_PP_CAT (Prefix, _TRPC_METHOD_NAME (Method))				\
    (_TRPC_METHOD_GENERATE_ARGS (Method) )						\
    {											\
      _TRPCMBasePtr marshalledArgs (							\
	boost::make_shared<_TRPC_MARSHAL_TYPE_NAME(Method)> (				\
	  BOOST_PP_REPEAT (								\
	    BOOST_PP_ARRAY_SIZE (_TRPC_METHOD_ARGS(Method)),				\
	    _TRPC_GENERATE_STUB_MARSHAL_ARGS, ~) ) );					\
      _TRPCPostCommand (_TRPC_GENERATE_METHOD_ID(~, ~, 0, Method), marshalledArgs);	\
    }
// Generate stubs for all methods
#define TRPC_GENERATE_STUB_DEFNS(Prefix, Methods)				\
  BOOST_PP_SEQ_FOR_EACH (_TRPC_GENERATE_METHOD_STUB, Prefix, Methods)

/*
 * Dispatcher (runs in worker thread, extract pending calls from queue, execute)
 */
// Marshalling: pack method args
#define _TRPC_MARSHAL_BASE_TYPE						\
  struct _TRPCMBase							\
  {									\
    virtual ~_TRPCMBase() {}						\
  };									\
  typedef boost::shared_ptr<_TRPCMBase> _TRPCMBasePtr;
// Generate a member of a marshal struct
#define _TRPC_GENERATE_MARSHAL_MEMBER(z, n, Args)		\
  boost::remove_reference<BOOST_PP_ARRAY_ELEM(n, Args) >::type BOOST_PP_CAT(v, n);
// Generate all members of a marshal struct
#define _TRPC_MARSHAL_MEMBERS(Method)					\
  BOOST_PP_REPEAT (							\
    BOOST_PP_ARRAY_SIZE (_TRPC_METHOD_ARGS(Method)),			\
    _TRPC_GENERATE_MARSHAL_MEMBER, _TRPC_METHOD_ARGS(Method))
// Generate constructor initializers for a marshal struct
#define _TRPC_GENERATE_MARSHAL_TYPE_CTOR_INIT(z, n, unused)	\
  BOOST_PP_COMMA_IF(n) BOOST_PP_CAT(v, n) ( BOOST_PP_CAT(A, n) )
// Generate a marshal struct
#define _TRPC_MARSHAL_CTOR_INIT(Method)					\
  BOOST_PP_REPEAT (							\
    BOOST_PP_ARRAY_SIZE (_TRPC_METHOD_ARGS(Method)),			\
    _TRPC_GENERATE_MARSHAL_TYPE_CTOR_INIT, ~)
#define _TRPC_NEED_INITIALIZERS(Method)					\
  BOOST_PP_ARRAY_SIZE (_TRPC_METHOD_ARGS(Method))
#define _TRPC_GENERATE_MARSHAL_TYPE(r, Unused, Method)				\
  struct _TRPC_MARSHAL_TYPE_NAME (Method) : public _TRPCMBase			\
  {										\
    _TRPC_MARSHAL_MEMBERS (Method)						\
    _TRPC_MARSHAL_TYPE_NAME (Method) (_TRPC_METHOD_GENERATE_ARGS (Method)) 	\
      BOOST_PP_IF(								\
	_TRPC_NEED_INITIALIZERS(Method),					\
	BOOST_PP_IDENTITY(:),			\
	BOOST_PP_EMPTY)() _TRPC_MARSHAL_CTOR_INIT(Method) {}							\
  };
// Generate enum with method IDs
#define _TRPC_GENERATE_METHOD_IDS(Methods)				\
  enum _TRPCCommandID							\
  {									\
    BOOST_PP_SEQ_FOR_EACH_I (_TRPC_GENERATE_METHOD_ID, ~, Methods)	\
  };
// Generate method call queue member
#define _TRPC_COMMAND_QUEUE_MEMBER(X)	_trpcQueue.X
#define _TRPC_DECLARE_COMMAND_QUEUE(CondType)	\
  struct _TRPCQueue				\
  {						\
    wxMutex Lock;				\
    std::queue<_TRPCCommand> Queue;		\
    CondType QueueChangedCond;			\
    _TRPCQueue() : QueueChangedCond (Lock) {}	\
  } _trpcQueue;
// Generate list of all members od a marshal struct
#define _TRPC_GENERATE_MARSHAL_MEMBER_LIST_ENTRY(z, i, Prefix)	\
  BOOST_PP_COMMA_IF(i) Prefix BOOST_PP_CAT(v, i)
#define _TRPC_MARSHAL_MEMBERS_LIST(Prefix, Method)			\
  BOOST_PP_REPEAT (							\
    BOOST_PP_ARRAY_SIZE (_TRPC_METHOD_ARGS(Method)),			\
    _TRPC_GENERATE_MARSHAL_MEMBER_LIST_ENTRY, Prefix)
// Generate a case in the dispatcher switch
#define _TRPC_DISPATCHER_CASE(r, ObjAndImplPrefix, Method)			\
  case _TRPC_GENERATE_METHOD_ID(~, ~, 0, Method):				\
    {										\
      BOOST_PP_IF(								\
	BOOST_PP_ARRAY_SIZE (_TRPC_METHOD_ARGS(Method)),			\
	BOOST_PP_IDENTITY(							\
	  _TRPC_MARSHAL_TYPE_NAME(Method)* marshalInfo = 			\
	    static_cast<_TRPC_MARSHAL_TYPE_NAME(Method)*> (cmd.second.get());),	\
	BOOST_PP_EMPTY)()							\
      BOOST_PP_TUPLE_ELEM (2, 0, ObjAndImplPrefix)				\
      BOOST_PP_CAT(BOOST_PP_TUPLE_ELEM (2, 1, ObjAndImplPrefix),		\
	_TRPC_METHOD_NAME (Method)) (						\
	_TRPC_MARSHAL_MEMBERS_LIST (marshalInfo->, Method));			\
    }										\
    break;
 // Generate main dispatcher method
#define _TRPC_DISPATCHER_METHOD(Name, Obj, ImplPrefix, Methods)			\
  bool Name ()									\
  {										\
    _TRPCCommand cmd;								\
    {										\
      wxMutexLocker _lock (_TRPC_COMMAND_QUEUE_MEMBER(Lock));			\
      while (true)								\
      {										\
	if (_TRPC_COMMAND_QUEUE_MEMBER(Queue).size() > 0)			\
	{									\
	  cmd = _TRPC_COMMAND_QUEUE_MEMBER(Queue).front();			\
	  _TRPC_COMMAND_QUEUE_MEMBER(Queue).pop();				\
	  break;								\
	}									\
	else									\
	{									\
	  if (_TRPC_COMMAND_QUEUE_MEMBER(QueueChangedCond).Wait() 		\
	      != wxCOND_NO_ERROR)						\
	    return false;							\
	}									\
      }										\
    }										\
    switch (cmd.first)								\
    {										\
    BOOST_PP_SEQ_FOR_EACH (_TRPC_DISPATCHER_CASE, (Obj, ImplPrefix), Methods)	\
    }										\
    return true;								\
  }
// Generate all members, declarations etc. for the dispatching system
#define TRPC_GENERATE_DISPATCHER_COND_FWD(Dispatcher, Cond, Obj, ImplPrefix, Methods)	\
  _TRPC_MARSHAL_BASE_TYPE								\
  BOOST_PP_SEQ_FOR_EACH (_TRPC_GENERATE_MARSHAL_TYPE, ~, Methods)			\
  _TRPC_GENERATE_METHOD_IDS (Methods)							\
  typedef std::pair<_TRPCCommandID, _TRPCMBasePtr> _TRPCCommand;			\
  _TRPC_DECLARE_COMMAND_QUEUE(Cond)							\
  void _TRPCPostCommand (_TRPCCommandID cmd, const _TRPCMBasePtr& args)			\
  {											\
    wxMutexLocker _lock (_TRPC_COMMAND_QUEUE_MEMBER(Lock));				\
    _TRPC_COMMAND_QUEUE_MEMBER(Queue).push (std::make_pair (cmd, args));		\
    /*  The signal is lost if 1. the thread is currently processing a command or	\
	2. the thread is waiting for the queue lock.					\
	However, we don't care as 1. it'll check the queue again after processing ended	\
	2. it'll check the queue next thing after getting the lock anyway.		\
      */										\
    _TRPC_COMMAND_QUEUE_MEMBER(QueueChangedCond).Signal();				\
  }											\
  _TRPC_DISPATCHER_METHOD (Dispatcher, Obj, ImplPrefix, Methods)
#define TRPC_GENERATE_DISPATCHER_FWD(Dispatcher, Obj, ImplPrefix, Methods)	\
  TRPC_GENERATE_DISPATCHER_COND_FWD (Dispatcher, wxCondition, Obj, ImplPrefix, Methods)
#define TRPC_GENERATE_DISPATCHER(Dispatcher, ImplPrefix, Methods)		\
  TRPC_GENERATE_DISPATCHER_FWD (Dispatcher, BOOST_PP_EMPTY(), ImplPrefix, Methods)
  
#define TRPC_GENERATE_DISPATCHER_MAIN_FWD(Obj, ImplPrefix, Methods)	\
  TRPC_GENERATE_DISPATCHER_COND_FWD (Dispatch, platform::ThreadRPC::ConditionSignalMain, Obj, ImplPrefix, Methods)
#define TRPC_GENERATE_DISPATCHER_MAIN(ImplPrefix, Methods)		\
  TRPC_GENERATE_DISPATCHER_MAIN_FWD (DBOOST_PP_EMPTY(), ImplPrefix, Methods)

/*
 * Method implementations (called in worker thread)
 */
#define _TRPC_GENERATE_IMPL_DECL(r, Prefix, Method)			\
  void BOOST_PP_CAT (Prefix, _TRPC_METHOD_NAME (Method))		\
    (_TRPC_METHOD_GENERATE_ARGS (Method) );
#define TRPC_GENERATE_IMPLEMENTATION_DECLS(Prefix, Methods)		\
  BOOST_PP_SEQ_FOR_EACH (_TRPC_GENERATE_IMPL_DECL, Prefix, Methods)

#include "export.h"

#include <wx/event.h>
#include <wx/thread.h>

#include <boost/enable_shared_from_this.hpp>

namespace platform
{

DECLARE_EXPORTED_EVENT_TYPE(PLATFORM_EXPORTED, trpcEVT_METHOD_PENDING, -1)

struct ThreadRPC
{
  class MethodPendingEvent;
  
  struct MainThreadRunnerBase;
  typedef boost::shared_ptr<MainThreadRunnerBase> MainThreadRunnerBasePtr;
  struct MainThreadRunnerBase :
    public boost::enable_shared_from_this<MainThreadRunnerBase>
  {
    virtual bool Dispatch() = 0;
    
    MainThreadRunnerBasePtr GetShared() { return shared_from_this(); }
  };
  
  static PLATFORM_EXPORTED void InstallMainThreadEventHandler ();
  
  template<typename Class>
  class MainThreadRunner : public MainThreadRunnerBase
  {
  public:
    MainThreadRunner()
    {
      Class* real_me = static_cast<Class*> (this);
      /* *QUITE* hacky, as this runs before the QueueChangedCond ctor */
      real_me->_trpcQueue.QueueChangedCond.mtr = this;
      
      InstallMainThreadEventHandler ();
    }
  };
  
  /**
   * Condition look alike in case the 'worker' thread is the main thread.
   * If this is to be used, the stubbed class must derive from MainThreadRunner!
   */
  class PLATFORM_EXPORTED ConditionSignalMain
  {
  public:
    // Pointer to running class
    MainThreadRunnerBase* mtr;
    
    ConditionSignalMain (wxMutex&) {}
    
    inline wxCondError Wait()
    {
      // Don't wait, instead, pass control back to main event loop
      return wxCOND_MISC_ERROR;
    }
    void Signal(); // will post event to event loop to trigger queue execution
  };
};

} // namespace platform

#endif // __NUTOGUI_THREADRPC_H__
