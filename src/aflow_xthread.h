
#ifndef AFLOW_XTHREAD_H
#define AFLOW_XTHREAD_H
#ifdef AFLOW_MULTITHREADS_ENABLE

#include <iostream>
#include <mutex>

#include <sys/types.h>

// ME20220130
namespace xthread {
  class xThread {
  public:
    xThread(int nmax = 0, int nmin = 1);
    xThread(std::ostream& oss, int nmax = 0, int nmin = 1);
    xThread(const xThread& xt);
    const xThread& operator=(const xThread& xt);
    ~xThread();

    void clear();

    void setCPUs(int nmax, int nmin = 1);
    void setProgressBar(std::ostream& oss);
    void unsetProgressBar();

    template <typename F, typename... A> void run(uint ntasks, F& func, A&... args);
    template <typename F, typename... A> void run(int ntasks, F& func, A&... args);
    template <typename F, typename... A> void run(unsigned long long int ntasks, F& func, A&... args);
    template <typename IT, typename F, typename... A> void run(const IT& it, F& func, A&... args);
    template <typename IT, typename F, typename... A> void run(IT& it, F& func, A&... args);

    template <typename F, typename... A> void runPredistributed(int ntasks, F& func, A&... args);
    template <typename F, typename... A> void runPredistributed(uint ntasks, F& func, A&... args);
    template <typename F, typename... A> void runPredistributed(unsigned long long int ntasks, F& func, A&... args);

  private:
    void free();
    void copy(const xThread&);

    int ncpus_max;
    int ncpus_min;
#ifdef AFLOW_MULTITHREADS_ENABLE
    std::mutex mtx;
#endif
    std::ostream* progress_bar;
    unsigned long long int progress_bar_counter;
    bool progress_bar_set;

    void initializeProgressBar(unsigned long long int ntasks);

    [[nodiscard]] int reserveThreads(unsigned long long int ntasks) const;
    void freeThreads(int ncpus);

    template <typename I, typename F, typename... A> void run(I& it, I& end, unsigned long long int ntasks, F& func, A&... args);
    template <typename I, typename F, typename... A> void spawnWorker(int ithread, I& it, I& end, unsigned long long int ntasks, F& func, A&... args);
    template <typename I> I advance(I& it, I& end, unsigned long long int ntasks, bool update_progress_bar = false);

    template <typename I, typename F, typename... A> void runPredistributed(I ntasks, F& func, A&... args);
    template <typename I, typename F, typename... A> void spawnWorkerPredistributed(int ithread, I startIndex, I endIndex, F& func, A&... args);
  };
} // namespace xthread
#endif

#endif // AFLOW_XTHREAD_H
