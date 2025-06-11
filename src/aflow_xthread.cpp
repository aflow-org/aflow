// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Marco Esters
//
// Contains thread manager class xThread, which executes functions in parallel
// and handles all index/iterator management and progress bar updates.
//
// It supports two modes: on-the-fly task distribution (default) or with
// pre-distributed task lists. The latter is useful for many small tasks or
// if the tasks involve expensive pre-computation steps.
//
// This class uses variadic templates. For more information, see here:
// https://kubasejdak.com/variadic-templates
//
// ----------
//
// Usage notes:
//
// Requirements for the functions that can be run:
// For on-the-fly redistribution, the function must have either an index
// (integer type) or an iterable to parallelize over as its first parameter.
//
// For pre-distributed execution, the function must have two integer types
// (start and end points) as first parameters. The pre-distributed scheme
// cannot be used with iterators and currently does not support progress bars.
//
// Additionally, function inputs cannot be prvalues (see here for a distinction
// between value types: https://accu.org/journals/overload/27/150/knatten_2641/).
//
// ----------
//
// Running functions with xThread:
//
// A function can be run in parallel with the on-the-fly scheme using the
// run() function.
//
// For functions with an index: run(max_index, function, args...)
// For functions over an iterable: run(iterable, function, args...)
//
// For the pre-distributed scheme: runPredistributed(max_index, function, args...)
//
// Every function handled by this class needs to be instantiated here to avoid
// linker errors. There is a section on how to do this at the end of this file.
//
// Static functions are called differently than non-static member functions.
// The following examples show the difference.
//
// Example functions:
// void f1(uint i, const vector<int>& vint, vector<double>& vdbl)
//   Parallelize over vdbl
// void f2(vector<int>::iterator& i, const xmatrix<double>& mdbl)
//   Parallelize over vector<int> v that i will iterate over
// void f3(int i)
//
// If functions are static functions:
// Static member functions can directly be plugged into run()
//
// f1:
//   uint ntasks = vdbl.size();
//   xThread xt;
//   xt.run(ntasks, f1, vint, vdbl);
//
// f2:
//   xThread xt;
//   xt.run(v, mdbl);
//
// f3:
//   xThread xt;
//   xt.run(ntasks, f3);
//
//   With ntasks = number of tasks.
//
// Non-static member functions with xThread:
//
// Member functions of a class cannot be directly plugged into run because they
// have to be bound to an instance of the class (this is also the reason why
// "this" had to be added when calling a member function in std::thread). This
// can be done using std::bind. For more information and an example, see:
// https://en.cppreference.com/w/cpp/utility/functional/placeholders
//
// The following examples show how std::bind works with the functions f1 - f3.
// Let all functions belong to class C instantiated as cls. Let f1
// and f2 be called inside another class function of C and let f3 be called
// outside.
//
// f1:
//   uint ntasks = vdbl.size();
//   xThread xt;
//   std::function<void(uint, const vector<int>&, const vector<double>&)> fn1 =
//     std::bind(&C::f1, this,
//       std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
//   xt.run(ntasks, fn1, vint, vdbl);
//
// f2:
//   xThread xt;
//   std::function<void(vector<int>::iterator&, const xmatrix<double>&)> fn2 =
//     std::bind(&C::f2, this, std::placeholders::_1, std::placeholders::_2);
//   xt.run(v, fn2, mdbl);
//
// f3:
//   xThread xt;
//   std::function<void(int)> fn3 = std::bind(&C::f3, &cls, std::placeholders::_1);
//   xt.run(ntasks, fn3);
//
// The template parameter in std::function takes the function return type and
// the type of the function inputs exactly as written in the function definition.
// std::bind takes an address to the function and a pointer to the class
// instance (always "this" when called inside the class). If the function has
// arguments, one std::placeholders::_N for each argument of the function needs
// to be added. This tells std::bind how many arguments to expect, which is why
// f1 needs three and f3 only one placeholder. Unfortunately, the number of
// arguments cannot simply be indicated using an integer because std::function
// needs to expand all the templates.
// Using namespace std::placeholders can make std::bind more legible because it
// allows writing _1 instead of std::placeholders::_1.
//
//
// Lambda functions with xThread:
//
// Lambda functions should be assigned to a std::function. However, they do not
// need to be bound like non-static member functions.
//
// ----------
//
// Setting the number of CPUs:
// The default constructor uses as many threads as possible, but the number of
// CPUs can be passed into the constructor or changed via setCPUs(). A minimum
// number of threads can be set as well. In that case, xThread will wait until
// that minimum number of threads is available.
//
// WARNING: If the minimum is set to a value that is unattainable, the function
// will never run. It is best to leave the minimum number of CPUs to 1, which
// will always start as it doesn't spawn a new thread.
//
// ----------
//
// Progress bars:
// To use a progress bar, pass an ostream into the setProgressBar() function
// before calling run(). unsetProgressBar() removes the progress bar.
//
// Do not set the progress bar if the xThread is called in another threaded
// run to avoid garbled output.
//
// ----------
//
// Thread safety:
// xThread guarantees only that no two instances of the called function write
// to the same index or iterator using mutexes. It cannot guarantee that the
// passed function itself is thread-safe. Functions that perform actions that
// are not thread-safe should pass their own mutex as an input parameter.
//
// ----------

#ifdef AFLOW_MULTITHREADS_ENABLE

#include "aflow.h"
#include "aflow_pocc.h"

// Global mutex that prevents two xThread instances from checking the number
// of available CPUs at the same time.
static std::mutex xthread_cpu_check;

namespace xthread {

  /// @brief Constructor for xThread
  ///
  /// @param nmax Maximum number of CPUs used by xThread (default: 0 for all available CPUs)
  /// @param nmin Mininum number of CPUs required to spawn thread workers default: 1)
  /// @param oss  stream for the progress bar output
  xThread::xThread(int nmax, int nmin) {
    free();
    setCPUs(nmax, nmin);
  }

  xThread::xThread(ostream& oss, int nmax, int nmin) {
    free();
    setCPUs(nmax, nmin);
    setProgressBar(oss);
  }

  xThread::xThread(const xThread& that) {
    copy(that);
  }

  const xThread& xThread::operator=(const xThread& that) {
    copy(that);
    return *this;
  }

  void xThread::copy(const xThread& that) {
    if (this == &that) return;
    // std::mutex should not be copied because
    // it needs to stay immutable
    ncpus_max = that.ncpus_max;
    ncpus_max = that.ncpus_min;
    progress_bar = that.progress_bar;
    progress_bar_counter = that.progress_bar_counter;
    progress_bar_set = that.progress_bar_set;
  }

  xThread::~xThread() {
    free();
  }

  void xThread::free() {
    ncpus_max = 0;
    ncpus_min = 0;
    progress_bar = nullptr;
    progress_bar_set = false;
    progress_bar_counter = 0;
  }

  // Progress bar -------------------------------------------------------------

  /// @brief Tells xThread that it has to update a progress bar when running
  ///
  /// @param oss The ostream for the progress bar
  void xThread::setProgressBar(ostream& oss) {
    progress_bar = &oss;
    progress_bar_set = true;
    progress_bar_counter = 0;
  }

  void xThread::unsetProgressBar() {
    progress_bar = nullptr;
    progress_bar_set = false;
    progress_bar_counter = 0;
  }

  void xThread::initializeProgressBar(unsigned long long int ntasks) {
    progress_bar_counter = 0;
    pflow::updateProgressBar(0, ntasks, *progress_bar);
  }

  // CPU functions ------------------------------------------------------------

  /// @brief Sets the minimum and maximum number of CPUs used for threading
  ///
  /// @param nmax Maximum number of CPUs used by xThread (default: 0 for all available CPUs)
  /// @param nmin Mininum number of CPUs required to spawn thread workers (default: 1);
  //              set to 0 for nmin = nmax
  void xThread::setCPUs(int nmax, int nmin) {
    if (nmax < nmin) std::swap(nmax, nmin);
    ncpus_max = (nmax > 0)?nmax:(KBIN::get_NCPUS());
    ncpus_min = (nmin > 0)?nmin:nmax;
  }

  /// @brief Checks if enough threads are available and "reserves" them
  /// using XHOST.CPU_active
  ///
  /// @param ntasks Number of tasks
  ///
  /// @return Number of threads to use for the function
  ///
  /// This prevents threaded functions that spawn other multi-threaded
  /// processes from allocating more threads than the machine can afford.
  /// This only works within a single AFLOW run
  ///
  /// There are two ways to allocate CPUs for competing processes:
  ///   1) always maximize CPU usage: as soon as enough threads are available
  ///      for a process, allocate and run. This will reduce idle time, but
  ///      will prioritize smaller jobs. This is the default.
  ///   2) first come, first serve: allocation is halted until the first
  ///      process to call this function gets the threads it needs. This will
  ///      make sure that large jobs are not stuck until all small tasks are
  ///      finished, but can lead more idle time. This process can be set with
  ///      the compiler flag -DXTHREAD_USE_QUEUE.
  int xThread::reserveThreads(unsigned long long int ntasks) {
    uint sleep_second = 10;
    int ncpus_max_available = KBIN::get_NCPUS();
    int ncpus = 0, ncpus_available = 0, nmax = 0, nmin = 0;
    // Adjust max. and min. number of CPUs to the number of tasks.
    // Decrement numbers by 1 because the main thread (which is idle
    // during execution) should not count towards the total number
    // of active CPUs
    if (ntasks <= (unsigned long long int) ncpus_max_available) {
      // Downcasting is safe because ntasks has a value in int range
      nmin = std::min(ncpus_min, (int) ntasks) - 1;
      nmax = std::min(ncpus_max, (int) ntasks) - 1;
    } else {
      nmin = ncpus_min - 1;
      nmax = ncpus_max - 1;
    }
#ifdef XTHREAD_USE_QUEUE
    std::lock_guard<std::mutex> lk(xthread_cpu_check);
#endif
    do {
#ifndef XTHREAD_USE_QUEUE
      xthread_cpu_check.lock();
#endif
      ncpus_available = ncpus_max_available - XHOST.CPU_active;
      if (ncpus_available >= nmin) {
        ncpus = std::min(ncpus_available, nmax);
        // "reserve" threads globally
        XHOST.CPU_active += ncpus;
        // Increase by 1 - the main thread does not count towards the
        // total, but should still be used for a worker
        ncpus++;
#ifndef XTHREAD_USE_QUEUE
        xthread_cpu_check.unlock();
#endif
      } else {
#ifndef XTHREAD_USE_QUEUE
        xthread_cpu_check.unlock();
#endif
        aurostd::Sleep(sleep_second);
      }
    } while (ncpus_available < nmin);
    return ncpus;
  }

  /// @brief "Frees" threads in XHOST.CPU_active
  ///
  /// @param ncpus Number of threads to free
  void xThread::freeThreads(int ncpus) {
    std::lock_guard<std::mutex> lk(xthread_cpu_check);
    // Decrease by 1 to account for the main thread not counting
    // towards the number of active CPUs
    XHOST.CPU_active -= (ncpus - 1);
  }

  // On-the-fly scheme --------------------------------------------------------

  /// @brief Overload for running threads using indices to keep track of tasks
  ///
  /// @param ntasks The number of tasks over which the function is parallelized
  /// @param func The function to be called by the worker
  /// @param args The arguments passed into func, if any
  template <typename F, typename... A>
  void xThread::run(int ntasks, F& func, A&... args) {
    if (ntasks <= 0) return;
    int task_index = 0;
    run<int, F, A...>(task_index, ntasks, (unsigned long long int) ntasks, func, args...);
  }

  template <typename F, typename... A>
  void xThread::run(uint ntasks, F& func, A&... args) {
    if (ntasks == 0) return;
    uint task_index = 0;
    run<uint, F, A...>(task_index, ntasks, (unsigned long long int) ntasks, func, args...);
  }

  template <typename F, typename... A>
  void xThread::run(unsigned long long int ntasks, F& func, A&... args) {
    if (ntasks == 0) return;
    unsigned long long int task_index = 0;
    run<unsigned long long int, F, A...>(task_index, ntasks, ntasks, func, args...);
  }

  /// @brief Overload for running threads using iterables
  ///
  /// @param it the iterable over which to parallelize
  /// @param func The function to be called by the worker
  /// @param args The arguments passed into func, if any
  ///
  /// ntasks will be converted to unsigned long long int because
  /// the progress bar counter in pflow is of that type.
  template <typename IT, typename F, typename... A>
  void xThread::run(const IT& it, F& func, A&... args) {
    int dist = (int) std::distance(it.begin(), it.end());
    if (dist <= 0) return; // Cannot iterate backwards (yet)
    unsigned long long int ntasks = (unsigned long long int) dist;
    typename IT::const_iterator start = it.begin();
    typename IT::const_iterator end = it.end();
    run(start, end, ntasks, func, args...);
  }

  template <typename IT, typename F, typename... A>
  void xThread::run(IT& it, F& func, A&... args) {
    int dist = (int) std::distance(it.begin(), it.end());
    if (dist <= 0) return; // Cannot iterate backwards (yet)
    unsigned long long int ntasks = (unsigned long long int) dist;
    typename IT::iterator start = it.begin();
    typename IT::iterator end = it.end();
    run(start, end, ntasks, func, args...);
  }

  /// @brief Executes a function over at least ncpus_min and at most ncus_max threads
  ///
  /// @param it The iterator or task index
  /// @param end The end() iterator or the last task index
  /// @param ntasks Number of tasks
  /// @param func The function to be called by the worker
  /// @param args The arguments passed into func, if any
  template <typename I, typename F, typename... A>
  void xThread::run(I& it, I& end, unsigned long long int ntasks, F& func, A&... args) {
    if (ntasks == 0) return;
    int ncpus = reserveThreads(ntasks);
    if (progress_bar_set) initializeProgressBar(ntasks);

    if (ncpus > 1) {
      vector<std::thread*> threads;
      for (int i = 0; i < ncpus; i++) {
        threads.push_back(new std::thread(&xThread::spawnWorker<I, F, A...>, this,
                                          i, std::ref(it), std::ref(end), ntasks,
                                          std::ref(func), std::ref(args)...)
        );
      }
      for (uint t = 0; t < threads.size(); t++) {
        threads[t]->join();
        delete threads[t];
      }
    } else {
      // No need for thread overhead when ncpus == 1
      for (; it != end; ++it) {
        func(it, args...);
        if (progress_bar_set) pflow::updateProgressBar(++progress_bar_counter, ntasks, *progress_bar);
      }
    }

    freeThreads(ncpus);
  }

  /// @brief Worker called by the threads
  ///
  /// @param it The iterator or task index
  /// @param end The end() iterator or the last task index
  /// @param ntasks Number of tasks (for the progress bar)
  /// @param func The function to be called by the worker
  /// @param args The arguments passed into func, if any
  template <typename I, typename F, typename... A>
  void xThread::spawnWorker(int ithread, I& it, I& end,
                            unsigned long long int ntasks,
                            F& func, A&... args) {
    I icurr = advance(it, end, ntasks);

    while (icurr != end) {
      try {
        func(icurr, args...); // Call function
      } catch (aurostd::xerror& e) {
        string message = "Error in thread " + aurostd::utype2string<int>(ithread) + ": " + e.what();
        throw aurostd::xerror(e.whereFileName(), e.whereFunction(), message, e.whatCode());
      }
      icurr = advance(it, end, ntasks, progress_bar_set);
    }
  }

  /// @brief Advances to the next task
  ///
  /// @param it The iterator or task index
  /// @param end The end() iterator or the last task index
  /// @param ntasks Number of tasks (for the progress bar)
  /// @param update_progress_bar Update progress bar if true (default: false)
  ///
  /// @return The next task index or iterator position
  template <typename I>
  I xThread::advance(I& it, I& end,
                     unsigned long long int ntasks,
                     bool update_progress_bar) {
    std::lock_guard<std::mutex> lk(mtx);
    if (update_progress_bar && (progress_bar_counter <= ntasks)) {
      pflow::updateProgressBar(++progress_bar_counter, ntasks, *progress_bar);
    }
    return (it == end)?end:(it++);
  }

  // Pre-distributed scheme ---------------------------------------------------

  /// @brief Overload for running threads using indices to keep track of tasks
  ///
  /// @param ntasks The number of tasks over which the function is parallelized
  /// @param func The function to be called by the worker
  /// @param args The arguments passed into func, if any
  ///
  /// The overloads are not necessary, but they allow for the same instantiation
  /// scheme as the on-the-fly methods.
  template <typename F, typename... A>
  void xThread::runPredistributed(int ntasks, F& func, A&... args) {
    if (ntasks <= 0) return;
    runPredistributed<int, F, A...>(ntasks, func, args...);
  }

  template <typename F, typename... A>
  void xThread::runPredistributed(uint ntasks, F& func, A&... args) {
    if (ntasks == 0) return;
    runPredistributed<uint, F, A...>(ntasks, func, args...);
  }

  template <typename F, typename... A>
  void xThread::runPredistributed(unsigned long long int ntasks, F& func, A&... args) {
    if (ntasks == 0) return;
    runPredistributed<unsigned long long int, F, A...>(ntasks, func, args...);
  }

  template <typename I, typename F, typename... A>
  void xThread::runPredistributed(I ntasks, F& func, A&... args) {
    if (ntasks <= 0) return;
    int ncpus = reserveThreads(ntasks);

    if (ncpus > 1) {
      vector<std::thread*> threads;
      I n = (I) ncpus; // convert ncpus to same type
      I tasks_per_thread = ntasks/n;
      I remainder = ntasks % n;
      I startIndex = 0, endIndex = 0;
      for (I t = 0; t < n; t++) {
        if (t < remainder) {
          startIndex = (tasks_per_thread + 1) * t;
          endIndex = startIndex + tasks_per_thread + 1;
        } else {
          startIndex = tasks_per_thread * t + remainder;
          endIndex = startIndex + tasks_per_thread;
        }
        threads.push_back(new std::thread(&xThread::spawnWorkerPredistributed<I, F, A...>, this,
                                          (int) t, startIndex, endIndex,
                                          std::ref(func), std::ref(args)...));
      }
      for (uint t = 0; t < threads.size(); t++) {
        threads[t]->join();
        delete threads[t];
      }
    } else {
      // No need for thread overhead when ncpus == 1
      func(0, ntasks, args...);
    }

    freeThreads(ncpus);
  }

  template <typename I, typename F, typename... A>
  void xThread::spawnWorkerPredistributed(int ithread, I startIndex, I endIndex, F& func, A&... args) {
    try {
      func(startIndex, endIndex, args...);
    } catch (aurostd::xerror& e) {
      string message = "Error in thread " + aurostd::utype2string<int>(ithread) + ": " + e.what();
      throw aurostd::xerror(e.whereFileName(), e.whereFunction(), message, e.whatCode());
    }
  }

}

// Template Instantiation
//
// Functions run with xThread need to be instantiated to avoid linker errors.
// The way the instantiation needs to be done depends on whether an index or
// an iterator is used, on whether the function is static or a std::function,
// and on the origin of the arguments passed into run().
//
// On-the-fly and pre-distributed schemes are instantiated the same way though.
//
// ----------
//
// Instantiating functions using an index:
//
// It is best to learn by example. Take the following function f1:
//
// void f1(int, const vector<int>&, const vector<double>&, vector<int>&, vector<double>&);
//
// Let it be called inside function f2:
//
// bool f2(const vector<int>& vint1, vector<double>& vdbl2) {
//   vector<int> vint2;
//   vector<double> vdbl1;
//   uint ntasks = vdbl2.size();
//   xthread::xThread xt;
//   xt.run(ntasks, f1, vint1, vdbl1, vint2, vdbl2);
// }
//
// The instantiation is:
//
// template void xThread::run<
//   void(int, const vector<int>&, const vector<double>&, vector<int>&, vector<double>&),
//   const vector<int>,
//   vector<double>,
//   vector<int>,
//   vector<double>
// >(uint, void(&) (int, const vector<int>&, const vector<double>&, vector<int>&, vector<double>&),
//   const vector<int>&,
//   vector<double>&,
//   vector<int>&,
//   vector<double>&
// );
//
// The first parameter inside <> is function type with the argument types in
// parentheses. The remaining parameters are the argument types but without &.
// The usage of const is different (see explanation later).
//
// The first parameter inside () is the index type, then the function type with
// the argument types. This time, the function type has an additional (&). The
// parentheses are mandatory. The remaining parameters are the argument types,
// but with & this time.
//
// The remaining parameters do not follow the function definition since the
// first vector<double> does not get const, even though one is const & in the
// function definition. On the other hand, the first vector<int> is const &.
// The important part is how it passed into xt.run(), not how it is passed into
// f1. vdbl1 is created inside f2 as vector<double> and is thus passed into
// run() as vector<double>&. vint1 on the other hand is passed into f2 as
// const vector<int>& and will thus be passed into run() as const vector<int>&.
//
// So, if an argument is created in the same function that calls run(), do not
// use const inside (). If it is passed down by another function as const &,
// use const &.
//
//
// If f2 was converted to a std::function, e.g. via std::bind for non-static
// member function, the instantiation is:
//
// template void xThread::run<
//   uint, std::function<void(int, const vector<int>&, const vector<double>&, vector<int>&, vector<double>&)>,
//   const vector<int>,
//   vector<double>,
//   vector<int>,
//   vector<double>
// >(uint, std::function<(int, const vector<int>&, const vector<double>&, vector<int>&, vector<double>&)>&,
//   const vector<int>&,
//   vector<double>&,
//   vector<int>&,
//   vector<double>&
// );
//
// Note: "template void" will always be "template void" since it is the type of
// xThread::run(). It does not depend on the types of f1 and f2.
//
// ----------
//
// Instantiating functions using an iterator:
//
// Take the function f1:
//
// void f1(vector<int>::iterator& it, const vector<double>&);
//
// Let it be called as:
//
//   vector<int> v1;
//   vector<double> v2;
//   xthread::xThread xt;
//   xt.run(v1, f1, v2);
//
// The instantiation is then:
//
//   template void xThread::run<
//     vector<int>,
//     void(vector<int>::iterator&, const vector<double>&),
//     vector<double>
//   >(
//     vector<int>&,
//     void (&) (vector<int>::iterator&, const vector<double>&),
//     vector<double>&
//   );
//
// Using const_iterators work the same way. The first template parameter of run
// cannot be const though.
//
// Note that the iterable has to appear as the first template parameter.
// The rules and caveats for member functions and const references are the same
// as for functions using an index.
//
// ----------
//
// Common mistakes
//
// When encountering linker errors, check if:
//
// - the index variable is exactly the same type (no implicit conversions)
// - all arguments inside () have an & (except for the index variable)
// - no parameter inside <> has an &
// - this includes the function inside () - double-check
// - the parameters inside the function are identical to the function definition
// - only use const when passed in by parent function or when defined as const
// - const inside () and <> are used the same way
// - the iterable is added before functions inside <>
// - on the other hand, the index variable type is not to be added inside <>
//
// ----------
//
// All instantiations should be added below with a comment on which function is
// instantiated. One instantiation can cover multiple functions.
//

#include "APL/aflow_apl.h"

namespace xthread {

  // run ----------------------------------------------------------------------

  //apl::AtomicDisplacements::calculateEigenvectorsInThread
  //apl::DOSCalculator::calculateInOneThread
  //apl::PhononDispersionCalculator::calculateInOneThread
  //apl::TCONDCalculator::calculateTransitionProbabilitiesIsotope
  template void xThread::run<
    std::function<void(int)>
  >(int, std::function<void(int)>&
  );

  //qca::QuasiChemApproxCalculator::calculateProbabilityCluster1D
  //qca::QuasiChemApproxCalculator::calculateProbabilityClusterND
  template void xThread::run<
    std::function<void(int, const int&)>,
    int
  >(int, std::function<void(int, const int&)>&,
    int&
  );

  //lambda function inside aurostd::multithread_execute
  template void xThread::run<
    deque<string>,
    std::function<void(const deque<string>::const_iterator&)>
  >(const deque<string>&, std::function<void(const deque<string>::const_iterator&)>&);

  //apl::PhononCalculator::calculateGroupVelocitiesThread
  template void xThread::run<
    std::function<void(int, vector<vector<double> >&, vector<xmatrix<xcomplex<double> > >&, vector<vector<xvector<double> > >&)>,
    vector<vector<double> >,
    vector<xmatrix<xcomplex<double> > >,
    vector<vector<xvector<double> > >
  >(int, std::function<void(int, vector<vector<double> >&, vector<xmatrix<xcomplex<double> > >&, vector<vector<xvector<double> > >&)>&,
    vector<vector<double> >&,
    vector<xmatrix<xcomplex<double> > >&,
    vector<vector<xvector<double> > >&
  );

  //apl::TCONDCalculator::calculateTransitionProbabilitiesPhonon
  template void xThread::run<
    std::function<void(int, vector<vector<vector<vector<double> > > >&, const vector<vector<vector<xcomplex<double> > > >&)>,
    vector<vector<vector<vector<double> > > >,
    vector<vector<vector<xcomplex<double> > > >
  >(int, std::function<void(int, vector<vector<vector<vector<double> > > >&, const vector<vector<vector<xcomplex<double> > > >&)>&,
    vector<vector<vector<vector<double> > > >&,
    vector<vector<vector<xcomplex<double> > > >&
  );

  //POccCalculator::calculatePhononDOSThread
  template void xThread::run<
    std::function<void(uint, const vector<uint>&, const aurostd::xoption&, vector<apl::DOSCalculator>&, vector<xDOSCAR>&, std::mutex&)>,
    vector<uint>,
    aurostd::xoption,
    vector<apl::DOSCalculator>,
    vector<xDOSCAR>,
    std::mutex
  >(uint, std::function<void(uint, const vector<uint>&, const aurostd::xoption&, vector<apl::DOSCalculator>&, vector<xDOSCAR>&, std::mutex&)>&,
    vector<uint>&,
    aurostd::xoption&,
    vector<apl::DOSCalculator>&,
    vector<xDOSCAR>&,
    std::mutex&
  );

  //POccCalculator::calculatePOccSuperCellUFF
  template void xThread::run<
    std::function<void(int, vector<pocc::POccSuperCell>&, const vector<pocc::POccUFFEnergyAnalyzer>&, const vector<vector<vector<int>>>&, size_t&, std::mutex&, std::mutex&)>,
    vector<pocc::POccSuperCell>,
    vector<pocc::POccUFFEnergyAnalyzer>,
    vector<vector<vector<int>>>,
    size_t,
    std::mutex,
    std::mutex
  >(int, std::function<void(int, vector<pocc::POccSuperCell>&, const vector<pocc::POccUFFEnergyAnalyzer>&, const vector<vector<vector<int>>>&, size_t&, std::mutex&, std::mutex&)>&,
    vector<pocc::POccSuperCell>&,
    vector<pocc::POccUFFEnergyAnalyzer>&,
    vector<vector<vector<int>>>&,
    size_t&,
    std::mutex&,
    std::mutex&
  );

  //POccCalculator::countUniquePOccSuperCellUFF
  template void xThread::run<
    std::function<void(int, std::map<unsigned long long int, std::unordered_map<unsigned long int, unsigned long int>>&, const vector<pocc::POccSuperCell>&, const vector<pocc::POccUFFEnergyAnalyzer>&, const vector<vector<vector<int>>>&, size_t&, std::mutex&, std::mutex&)>,
    std::map<unsigned long long int, std::unordered_map<unsigned long int, unsigned long int>>,
    vector<pocc::POccSuperCell>,
    vector<pocc::POccUFFEnergyAnalyzer>,
    vector<vector<vector<int>>>,
    size_t,
    std::mutex,
    std::mutex
  >(int, std::function<void(int, std::map<unsigned long long int, std::unordered_map<unsigned long int, unsigned long int>>&, const vector<pocc::POccSuperCell>&, const vector<pocc::POccUFFEnergyAnalyzer>&, const vector<vector<vector<int>>>&, size_t&, std::mutex&, std::mutex&)>&,
    std::map<unsigned long long int, std::unordered_map<unsigned long int, unsigned long int>>&,
    vector<pocc::POccSuperCell>&,
    vector<pocc::POccUFFEnergyAnalyzer>&,
    vector<vector<vector<int>>>&,
    size_t&,
    std::mutex&,
    std::mutex&
  );


  //apl::TCONDCalculator::calculateAnharmonicRates
  template void xThread::run<std::function<void(int, const vector<vector<double> >&, vector<vector<double> >&)>,
    const vector<vector<double> >,
    vector<vector<double> >
  >(int, std::function<void(int, const vector<vector<double> >&, vector<vector<double> >&)>&,
    const vector<vector<double> >&,
    vector<vector<double> >&
  );

  //apl::TCONDCalculator::calculateDelta
  template void xThread::run<
    std::function<void(int, const vector<vector<double> >&, const vector<vector<xvector<double> > >&, vector<vector<xvector<double> > >&)>,
    const vector<vector<double> >,
    vector<vector<xvector<double> > >,
    vector<vector<xvector<double> > >
  >(int, std::function<void(int, const vector<vector<double> >&, const vector<vector<xvector<double> > >&, vector<vector<xvector<double> > >&)>&,
    const vector<vector<double> >&,
    vector<vector<xvector<double> > >&,
    vector<vector<xvector<double> > >&
  );

  //XtalFinderCalculator::performStructureConversions
  template void xThread::run<
    std::function<void(uint, const vector<bool>&, const vector<bool>&, const vector<bool>&)>,
    vector<bool>,
    vector<bool>,
    vector<bool>
  >(uint, std::function<void(uint, const vector<bool>&, const vector<bool>&, const vector<bool>&)>&,
    vector<bool>&,
    vector<bool>&,
    vector<bool>&
  );

  //XtalFinderCalculator::runComparisonThreads
  template void xThread::run<
    std::function<void(uint,
      vector<StructurePrototype>&,
      const vector<std::pair<uint, uint> >&,
      const vector<std::pair<uint, uint> >&,
      bool, bool, bool)>,
    vector<StructurePrototype>,
    vector<std::pair<uint, uint> >,
    vector<std::pair<uint, uint> >,
    bool, bool, bool
  >(
    uint,
    std::function<void(uint,
      vector<StructurePrototype>&,
      const vector<std::pair<uint, uint> >&,
      const vector<std::pair<uint, uint> >&,
      bool, bool, bool)>&,
    vector<StructurePrototype>&,
    vector<std::pair<uint, uint> >&,
    vector<std::pair<uint, uint> >&,
    bool&, bool&, bool&
  );

  //XtalFinderCalculator::getPrototypeDesignations
  template void xThread::run<
    vector<StructurePrototype>,
    std::function<void(vector<StructurePrototype>::iterator&)>
  >(
    vector<StructurePrototype>&,
    std::function<void(vector<StructurePrototype>::iterator&)>&
  );

  //XtalFinderCalculator::getMatchingAFLOWPrototypes
  template void xThread::run<
    std::function<void(uint, vector<StructurePrototype>&, const aurostd::xoption&)>,
    vector<StructurePrototype>,
    aurostd::xoption
  >(uint, std::function<void(uint, vector<StructurePrototype>&, const aurostd::xoption&)>&,
    vector<StructurePrototype>&,
    aurostd::xoption&
  );

  //aflowlib::AflowDB::createTable
  template void xThread::run<
    std::function<void(int, const vector<string>&, const vector<string>&)>,
    vector<string>,
    vector<string>
  >(int, std::function<void(int, const vector<string>&, const vector<string>&)>&,
    vector<string>&,
    vector<string>&
  );

  //aflowlib::XPLUG_Directories_ok
  template void xThread::run<
    void(uint, uint, const deque<string>&, deque<string>&, deque<bool>&, std::mutex&),
    uint,
    deque<string>,
    deque<string>,
    deque<bool>,
    std::mutex
  >(uint, void (&) (uint, uint, const deque<string>&, deque<string>&, deque<bool>&, std::mutex&),
    uint&,
    deque<string>&,
    deque<string>&,
    deque<bool>&,
    std::mutex&
  );

  //UnitTest::runUnitTest
  template void xThread::run<
    vector<string>,
    std::function<void(vector<string>::iterator&, const vector<string>&)>,
    vector<string>
  >(
    vector<string>&,
    std::function<void(vector<string>::iterator&, const vector<string>&)>&,
    vector<string>&
  );

  //_testPrototype
  template void xThread::run<
      void(uint, const vector<string>&, vector<uint>&, vector<string>&, std::mutex&),
      vector<string>,
      vector<uint>,
      vector<string>,
      std::mutex
    >(uint,
      void(&) (uint, const vector<string>&, vector<uint>&, vector<string>&, std::mutex&),
      vector<string>&,
      vector<uint>&,
      vector<string>&,
      std::mutex&
  );

  //convexHull  //ME+CO20220630
  template void xthread::xThread::run<
    std::function<void (
        uint,
        const vector<string>&,
        const aurostd::xoption&,
        const _aflags&,
        vector<uint>&,
        std::ostream&
        )>,
    vector<string>,
    aurostd::xoption,
    _aflags,
    vector<uint>,
    std::ostream
      >(uint,
          std::function<void (
            uint,
            const vector<string>&,
            const aurostd::xoption&,
            const _aflags&,
            vector<uint>&,
            std::ostream&
            )>&,
          vector<string>&,
          aurostd::xoption&, _aflags&, vector<uint>&, std::ostream&);

  // runPredistributed --------------------------------------------------------

  //XTalFinderCalculator::calculateSpaceGroups
  template void xThread::runPredistributed<
    std::function<void(uint, uint, uint)>,
    uint
  >(
    uint, std::function<void(uint, uint, uint)>&, uint&
  );

  //XtalFinderCalculator::computeLFAEnvironments
  //XtalFinderCalculator::calculateNearestNeighbors
  template void xThread::runPredistributed<
    std::function<void(uint, uint)>
  >(
    uint, std::function<void(uint, uint)>&
  );

  //XtalFinderCalculator::searchAtomMappings
  template void xThread::runPredistributed<
    std::function<bool(
      uint, uint, const xstructure&,
      const vector<double>&, const xstructure&, const string&,
      vector<xmatrix<double> >&, vector<structure_mapping_info>&, bool, bool
    )>,
    xstructure,
    vector<double>,
    xstructure,
    string,
    vector<xmatrix<double> >,
    vector<structure_mapping_info>,
    bool, bool
  >(
    uint, std::function<bool(uint, uint, const xstructure&,
      const vector<double>&, const xstructure&, const string&,
      vector<xmatrix<double> >&, vector<structure_mapping_info>&, bool, bool)>&,
    xstructure&,
    vector<double>&,
    xstructure&,
    string&,
    vector<xmatrix<double> >&,
    vector<structure_mapping_info>&,
    bool&, bool&
  );

  //aflowlib::AflowDB::getColStats
  template void xThread::runPredistributed<
    std::function<void(
      uint, int, const vector<string>&, vector<aflowlib::DBStats>&
    )>,
    const vector<string>,
    vector<aflowlib::DBStats>
  >(uint,
    std::function<void(
      uint, int, const vector<string>&, vector<aflowlib::DBStats>&
    )>&,
    const vector<string>&,
    vector<aflowlib::DBStats>&
  );

}

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
