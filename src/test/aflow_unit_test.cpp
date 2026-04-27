// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Written by Marco Esters, 2022
// Based on prior work by Hagen Eckert
//
// Unit test class to conduct parallelized unit tests with unified output.
//
// Tests are grouped into categories so that entire suites can be tested.
// For example, the test "aurostd" calls xscalar, xmatrix, etc. The individual
// tests can be called as well.

#include "aflow_unit_test.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_hash.h"
#include "AUROSTD/aurostd_time.h"
#include "AUROSTD/aurostd_xcomplex.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xfile.h"
#include "AUROSTD/aurostd_xfit.h"
#include "AUROSTD/aurostd_xhttp.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xparser.h"
#include "AUROSTD/aurostd_xparser_json.h"
#include "AUROSTD/aurostd_xrandom.h"
#include "AUROSTD/aurostd_xscalar.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow.h"
#include "aflow_aflowrc.h"
#include "aflow_defs.h"
#include "aflow_mock_json.h"
#include "aflow_xhost.h"
#include "aflow_xthread.h"
#include "aflowlib/aflowlib_entry_loader.h"
#include "aflowlib/aflowlib_libraries.h"
#include "flow/aflow_pflow.h"
#include "flow/aflow_xclasses.h"
#include "modules/COMPARE/aflow_compare_structure.h"
#include "modules/CUMULANTS/aflow_cumulants.h"
#include "modules/HULL/aflow_nhull.h"
#include "modules/HULL/aflow_nhull_entry.h"
#include "modules/HULL/aflow_nhull_facet.h"
#include "modules/HULL/aflow_nhull_util.h"
#include "modules/PROTOTYPES/aflow_anrl.h"  //DX20201104
#include "structure/aflow_lattice.h"
#include "structure/aflow_surface.h"
#include "structure/aflow_xatom.h"
#include "structure/aflow_xstructure.h"

using namespace std::placeholders;

using std::deque;
using std::function;
using std::ifstream;
using std::iostream;
using std::istringstream;
using std::map;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::pair;
using std::string;
using std::stringstream;
using std::tuple;
using std::unordered_map;
using std::vector;

using aurostd::xcomplex;
using aurostd::xmatrix;
using aurostd::xvector;

namespace unittest {

  UnitTest::UnitTest(ostream& oss) : xStream(oss) {
    initialize();
  }

  UnitTest::UnitTest(ofstream& mf, ostream& oss) : xStream(mf, oss) {
    initialize();
  }

  UnitTest::UnitTest(const UnitTest& ut) : xStream(*ut.getOFStream(), *ut.getOSS()) {
    copy(ut);
  }

  const UnitTest& UnitTest::operator=(const UnitTest& ut) {
    copy(ut);
    return *this;
  }

  UnitTest::~UnitTest() {
    free();
  }

  void UnitTest::clear() {
    free();
  }

  void UnitTest::free() {
    aflags.clear();
    test_functions.clear();
    test_groups.clear();
  }

  void UnitTest::copy(const UnitTest& ut) {
    if (this == &ut) {
      return;
    }
    aflags = ut.aflags;
    test_functions = ut.test_functions;
    test_groups = ut.test_groups;
  }

  void UnitTest::initialize() {
    free();
    utd = aurostd::EmbData::get_unit_test();
    const string dir = XHOST.vflag_control.getattachedscheme("DIRECTORY");
    aflags.Directory = (!dir.empty() ? dir : aurostd::getPWD());

    initializeTestFunctions();
    initializeTestGroups();
  }

  /// @brief Initialize unit tests and add them to map of test functions.
  void UnitTest::initializeTestFunctions() {
    // aurostd
    register_test(&UnitTest::xscalarTest, "xscalar", "xscalarTest():", "xscalar functions", "aurostd");
    register_test(&UnitTest::xvectorTest, "xvector", "xvectorTest():", "xvector functions", "aurostd");
    register_test(&UnitTest::xmatrixTest, "xmatrix", "xmatrixTest():", "xmatrix functions", "aurostd");
    register_test(&UnitTest::xfitTest, "xfit", "xfitTest():", "xfit functions", "aurostd");
    register_test(&UnitTest::aurostdMainTest, "aurostd_main", "aurostdMainTest():", "aurostd_main functions", "aurostd");
    register_test(&UnitTest::xfileTest, "xfile", "xfileTest():", "xfile functions", "aurostd");
    register_test(&UnitTest::xparserTest, "xparser", "xparserTest():", "xparser functions", "aurostd");

    // json serialization
    register_test(&UnitTest::jsonSerialization, "json_serialization", "jsonSerialization():", "json serialization", "json");

    // database
    register_test(&UnitTest::schemaTest, "schema", "schemaTest():", "AFLOW schema", "database");

    // xstructure
    register_test(&UnitTest::atomicEnvironmentTest, "atomic_environment", "atomicEnvironmentTest():", "Creating atomic environments", "xstructure");
    register_test(&UnitTest::xstructureParserTest, "xstructure_parser", "xstructureParserTest():", "xstructure parsers", "xstructure");
    register_test(&UnitTest::xstructureTest, "xstructure", "xstructureTest():", "xstructure functions", "xstructure");
    register_test(&UnitTest::kpathTest, "xstructure_kpath", "kpathTest():", "check lattice", "xstructure");

    // structure generation
    register_test(&UnitTest::ceramgenTest, "ceramgen", "ceramgenTest():", "pflow::GENERATE_CERAMICS()", "structure_gen");
    register_test(&UnitTest::prototypeGeneratorTest, "proto", "prototypeGeneratorTest():", "Generate all prototypes and test symmetry", "structure_gen");

    // EntryLoader
    register_test(&UnitTest::entryLoaderTest, "entry_loader", "entryLoaderTest():", "entryLoaderTest():", "entry_loader");

    // Qhull
    register_test(&UnitTest::qhullTest, "new_qhull", "qhullTest():", "qhull functions", "hull");

    // ovasp
    register_test(&UnitTest::ovaspTest, "ovasp", "ovaspTest():", "xOVASP constructions", "ovasp");

    //dielectric
    register_test(&UnitTest::dielectricTest, "dielectric", "dielectricTest():", "dielectric functions", "dielectric");

    // cumulants
    register_test(&UnitTest::cumulantsTest, "cumulants", "cumulantsTest():", "CUMULANTS module", "cumulants");

    // lib2raw local
    register_test(&UnitTest::lib2rawTest, "lib2raw", "lib2rawTest():", "local lib2raw functionality", "lib2raw");

    // surface tests
    register_test(&UnitTest::surface_TriangleArea, "surface_TriangleArea", "surface_TriangleArea():", "surface_TriangleArea", "surface");
    register_test(&UnitTest::surface_PlaneGetABCD, "surface_PlaneGetABCD", "surface_PlaneGetABCD():", "surface_PlaneGetABCD", "surface");
    register_test(&UnitTest::surface_PlaneDistance, "surface_PlaneDistance", "surface_PlaneDistance():", "surface_PlaneDistance", "surface");
    register_test(&UnitTest::surface_PlaneGetProjection, "surface_PlaneGetProjection", "surface_PlaneGetProjection():", "surface_PlaneGetProjection", "surface");
    register_test(&UnitTest::surface_PlaneGetHKL, "surface_PlaneGetHKL", "surface_PlaneGetHKL():", "surface_PlaneGetHKL", "surface");
    register_test(&UnitTest::surface_PlaneGetVVV, "surface_PlaneGetVVV", "surface_PlaneGetVVV():", "surface_PlaneGetVVV", "surface");
    register_test(&UnitTest::surface_PlaneGetVVV_V2, "surface_PlaneGetVVV_V2", "surface_PlaneGetVVV_V2():", "surface_PlaneGetVVV_V2", "surface");
    register_test(&UnitTest::surface_GetPlaneDensityAtoms, "surface_GetPlaneDensityAtoms", "surface_GetPlaneDensityAtoms():", "surface_GetPlaneDensityAtoms", "surface");

    // data
    register_test(&UnitTest::xpseudopotentialDataTest, "pspot_data", "get_pseudopotential_data():", "get_pseudopotential_data()", "data");
  }

  /// @brief Initialize xcheck struct to an empty object.
  xcheck UnitTest::initializeXCheck() {
    xcheck xt;
    resetUnitTest(xt);
    xt.func = nullptr;
    xt.function_name = "";
    xt.task_description = "";
    xt.test_group = "";
    return xt;
  }

  /// @brief Reset a unit test to its pre-run state based on test name.
  void UnitTest::resetUnitTest(const string& test_name) {
    if (test_functions.count(test_name)) {
      resetUnitTest(test_functions[test_name]);
    }
  }

  /// @brief Reset an xcheck object to its pre-run state.
  void UnitTest::resetUnitTest(xcheck& test) {
    test.errors.clear();
    test.finished = false;
    test.passed_checks = 0;
    test.results.clear();
  }

  /// @brief Create unit test groups.
  void UnitTest::initializeTestGroups() {
    test_groups.clear();

    test_groups["aurostd"] = {"xscalar", "xvector", "xmatrix", "xfit", "aurostd_main", "xfile", "xparser"};
    test_groups["database"] = {"schema", "entry_loader"};
    test_groups["structure"] = {"atomic_environment", "xstructure", "xstructure_parser", "xstructure_kpath"}; // AZ20231030 added xstructure_kpath
    test_groups["structure_gen"] = {"ceramgen", "proto"};
    // test_groups["ovasp"] = {"outcar"};

    for (std::map<string, vector<string>>::iterator it = test_groups.begin(); it != test_groups.end(); ++it) {
      const vector<string>& members = (*it).second;
      for (size_t i = 0; i < members.size(); i++) {
        test_functions[members[i]].test_group = (*it).first;
      }
    }
  }

  /// @brief Registers a test to the unit test system
  /// @param test_function The test function to register
  /// @param key the key to register the test under
  /// @param name the name for the test
  /// @param desc the description for the test
  /// @param group the testing group to add to
  /// @authors
  /// @mod{ST,20241208,created function}
  template <typename UTFunc> void UnitTest::register_test(const UTFunc&& test_function, const string& key, const string& name, const string& desc, const string& group) {
    xcheck xchk = initializeXCheck();
    xchk.func = std::bind(test_function, this, _1, _2, _3);
    xchk.function_name = name;
    xchk.task_description = desc;
    xchk.test_group = group;
    test_functions[key] = xchk;
  }

} // namespace unittest

// redirection imlpementation
namespace unittest {

  /// @class RedirectStream
  /// @brief A class to facilitate redirecting and capturing contents of an output stream.
  ///
  /// Upon construction of this class with a buffer and output stream, the buffer of the output stream
  /// is replaced with the supplied buffer. Any subsequent usage of the output stream will use the supplied buffer.
  /// This is used to effectively redirect and capture the output of a stream such as @c std::cout,
  /// the default stream for this class. The stream can be restored to its original buffer by calling @c release().
  /// The stream is also restored upon calling the destructor so the redirection will release upon leaving scope.
  ///
  /// The simplest usage of this class is to call @c run_with_capture_cout() with a function or lambda which takes
  /// no arguments. The captured output will be returned as a @c stringstream.
  ///
  /// Some example usage:
  ///
  /// @code
  /// std::stringstream stream;
  /// RedirectStream redirection(stream.rdbuf(), &std::cout);
  /// do_stuff()
  /// redirect.release();
  /// do_more_stuff(stream.str());
  /// // OR
  /// std::stringstream stream = RedirectStream::run_with_capture_cout(do_stuff());
  /// @endcode
  ///
  /// @authors
  /// @mod{ST,20241120,created class}

  /// @brief Starts capturing upon construction.
  /// @param buffer The new buffer to use to capture from the stream.
  /// @param captured_stream The stream that should be captured.
  RedirectStream::RedirectStream(std::streambuf* buffer, std::ostream* captured_stream) :
      original_stream(captured_stream), original_buffer(captured_stream->rdbuf(buffer)), capture_buffer(buffer), capturing(true) {}

  RedirectStream::~RedirectStream() {
    this->release();
  }

  /// @brief Start capturing. Does nothing if already capturing.
  /// @param buffer The new buffer to use to capture from the stream.
  /// @param captured_stream The stream that should be captured.
  void RedirectStream::capture(std::streambuf* buffer, std::ostream* captured_stream) {
    if (!capturing) {
      original_stream = captured_stream;
      capture_buffer = buffer;
      original_buffer = original_stream->rdbuf(capture_buffer);
      capturing = true;
    }
  }

  /// @brief Release the capturing and restore the stream. Does nothing if not capturing.
  void RedirectStream::release() {
    if (capturing) {
      original_stream->rdbuf(original_buffer);
      original_stream = nullptr;
      original_buffer = nullptr;
      capture_buffer = nullptr;
      capturing = false;
    }
  }

  /// @brief Convenience method to capture std::cout while func is running.
  /// @param func The function to be called while capturing. Must not take any parameters.
  /// @return The captured output.
  template <typename Function> stringstream RedirectStream::run_with_capture_cout(Function&& func) {
    stringstream stream;
    RedirectStream redirect(stream.rdbuf(), &std::cout);
    func();
    redirect.release();
    return stream;
  }
} // namespace unittest

// Run functions
namespace unittest {

  /// @brief Run single unit test.
  ///
  /// @param unit_test Unit test name
  ///
  /// @return Whether all tests were successful.
  bool UnitTest::runTestSuites(const string& unit_test) {
    const vector<string> vunit_test(1, unit_test);
    return runTestSuites(vunit_test);
  }

  /// @brief Run a set of unit tests.
  ///
  /// @param unit_tests_in Set of unit test names.
  ///
  /// @return Whether all tests were successful.
  ///
  /// The function consists of three steps:
  ///   1) Collect the set of tasks and expand test groups into individual tests.
  ///   2) Run all requested test functions, outputting results as test groups finish.
  ///   3) Print final summary.
  bool UnitTest::runTestSuites(const vector<string>& unit_tests_in) {
    long double start = aurostd::get_seconds();
    stringstream message;
    // Create task lists (groups or individual tests).
    // unit_test is the individual small tests over which to parallelize
    vector<string> unit_tests;
    vector<string> tasks;
    for (const auto& test : unit_tests_in) {
      const bool isgroup = (test_groups.find(test) != test_groups.end());
      if (test == "all") {
        tasks.clear();
        for (auto& test_group : test_groups) {
          tasks.push_back(test_group.first);
          for (const auto& m : test_group.second) {
            unit_tests.push_back(m);
          }
        }
        break;
      } else if (!isgroup && (test_functions.find(test) == test_functions.end())) {
        message << "Skipping unrecognized test name " << test << ".";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
      } else if (isgroup && !aurostd::WithinList(tasks, test)) {
        tasks.push_back(test);
        const vector<string>& members = test_groups[test];
        for (const auto& member : members) {
          unit_tests.push_back(member);
        }
      } else if (!aurostd::WithinList(tasks, test_functions[test].test_group)) {
        unit_tests.push_back(test);
        tasks.push_back(test);
      }
    }
    const uint ntasks = tasks.size();
    if (ntasks == 0) {
      message << "No unit tests to run.";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, *p_FileMESSAGE, *p_oss, _LOGGER_NOTICE_);
      return true;
    }

    // Run tests
    // Many AFLOW functions produce output to screen without the opportunity
    // to silence it, which makes unit test output harder to read and can
    // lead to garbled output when run in parallel. To get clean output,
    // silence output globally except for the functions that produce desired
    // unit test output, unless --quiet is requested or --debug is run.
    const bool quiet_copy = XHOST.QUIET;
    vector<string> whitelist;
    if (!XHOST.QUIET && !XHOST.DEBUG) {
      whitelist.emplace_back("unittest::UnitTest::runUnitTest()");
      // Add function names to whitelist for displayResults
      for (const auto& unit_test : unit_tests) {
        whitelist.push_back(test_functions[unit_test].function_name);
      }
      XHOST.QUIET = true;
      for (const auto& i : whitelist) {
        XHOST.LOGGER_WHITELIST.push_back(i);
      }
    }
#ifdef AFLOW_MULTITHREADS_ENABLE
    std::mutex const mtx;
    xthread::xThread xt(KBIN::get_NCPUS());
    std::function<void(vector<string>::iterator&, const vector<string>&)> fn = std::bind(&UnitTest::runUnitTest, this, _1, _2);
    xt.run(unit_tests, fn, tasks);
#else
    for (vector<string>::iterator it = unit_tests.begin(); it != unit_tests.end(); ++it) runUnitTest(it, tasks);
#endif
    XHOST.QUIET = quiet_copy;
    if (!XHOST.QUIET && !XHOST.DEBUG) {
      XHOST.QUIET = quiet_copy;
      for (size_t i = 0; i < whitelist.size(); i++) {
        XHOST.LOGGER_WHITELIST.pop_back();
      }
    }

    // Print final summary
    uint nsuccess = 0;
    vector<vector<string>> summary(tasks.size(), vector<string>(2));
    for (size_t t = 0; t < tasks.size(); t++) {
      const bool success = taskSuccessful(tasks[t]);
      if (success) {
        nsuccess++;
      }
      summary[t][0] = tasks[t];
      summary[t][1] = (success ? "pass" : "fail");
    }
    const long double duration = aurostd::get_delta_seconds(start);
    if (nsuccess == ntasks) {
      message << "Unit tests passed successfully (passing " << ntasks << " task" + string((ntasks == 1) ? "" : "s") + " in " << std::setprecision(3) << duration << " s).";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
    } else {
      message << "Some unit tests failed (" << (ntasks - nsuccess) << " of " << ntasks << " failed).";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
    }
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, formatResultsTable(summary), aflags, *p_FileMESSAGE, *p_oss, _LOGGER_RAW_);
    return (nsuccess == ntasks);
  }

  /// @brief Run unit test function inside a thread.
  ///
  /// @param it    Iterator pointing to the unit test name.
  /// @param tasks Set of task names.
  void UnitTest::runUnitTest(vector<string>::iterator& it, const vector<string>& tasks) {
    long double start = aurostd::get_seconds();
    const string& test_name = (*it);
    xcheck& test = test_functions[test_name];
    resetUnitTest(test);
    test.func(test.passed_checks, test.results, test.errors);
    // Output results
#ifdef AFLOW_MULTITHREADS_ENABLE
    std::lock_guard<std::mutex> const lk(mtx);
#endif
    test.finished = true;
    if (aurostd::WithinList(tasks, test_name)) {
      // If the test name is in the task list, it is not part
      // of a group, so no need to check if other members are done
      displayResult(test);
    } else {
      // Test if part of a test group, so check if all members
      // of the group have finished before producing output
      const string& group = test_functions[test_name].test_group;
      const vector<string>& vtests_group = test_groups[group];
      const uint ntests_group = vtests_group.size();
      uint i = 0;
      for (; i < ntests_group; i++) {
        const string& test_name_group = vtests_group[i];
        if (!test_functions[test_name_group].finished) {
          break;
        }
      }
      if (i == ntests_group) {
        // All finished - print output
        uint nsuccess = 0;
        for (size_t t = 0; t < vtests_group.size(); t++) {
          if (taskSuccessful(vtests_group[t])) {
            nsuccess++;
          }
        }
        stringstream message;
        const long double duration = aurostd::get_delta_seconds(start);
        if (nsuccess == ntests_group) {
          message << "Unit tests of group " << group << " passed successfully (passing " << ntests_group << " test" << ((ntests_group == 1) ? "" : "s") << " in " << std::setprecision(3) << duration << " s).";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
        } else {
          message << "Some unit tests of group " << group << " failed (" << (ntests_group - nsuccess) << " of " << ntests_group << " failed).";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
        }
        for (size_t t = 0; t < vtests_group.size(); t++) {
          displayResult(test_functions[vtests_group[t]]);
        }
      }
    }
  }

  /// @brief Check if a task has finished successfully.
  ///
  /// @param task The task to check.
  ///
  /// @return Whether a task has finished without errors.
  ///
  /// Criteria for returning true:
  ///   1) Task is finished.
  ///   2) All unit tests passed.
  ///   3) There are no additional errors.
  bool UnitTest::taskSuccessful(const string& task) {
    const std::map<string, vector<string>>::iterator it = test_groups.find(task);
    if (it != test_groups.end()) {
      const vector<string> members = (*it).second;
      for (size_t m = 0; m < members.size(); m++) {
        const string& member = members[m];
        const xcheck& xchk = test_functions[member];
        if (!xchk.finished || (!xchk.errors.empty()) || (xchk.passed_checks != xchk.results.size())) {
          return false;
        }
      }
      return true;
    } else {
      const xcheck& xchk = test_functions[task];
      return (xchk.finished && (xchk.errors.empty()) && (xchk.passed_checks == xchk.results.size()));
    }
  }
} // namespace unittest

// Output formatters
namespace unittest {
  /// @brief Convert results vector into a formatted table.
  ///
  /// @param table Structured table data.
  ///
  /// @return Formatted table string.
  ///
  /// This function makes no assumption about the table dimensions, so tables
  /// rows can have different sizes.
  /// Empty columns will be skipped, but not empty rows.
  /// The first column will have spaces prepended to indent the table.
  string UnitTest::formatResultsTable(const vector<vector<string>>& table) {
    const size_t nrows = table.size();
    if (nrows == 0) {
      return "";  // Empty table
    }

    // Determine dimensions of the table
    vector<size_t> col_sizes;
    size_t str_length = 0;
    size_t maxcol = 0;
    for (size_t r = 0; r < nrows; r++) {
      maxcol = std::max(table[r].size(), maxcol);
      for (size_t c = 0; c < table[r].size(); c++) {
        str_length = table[r][c].length();
        if (c == col_sizes.size()) {
          col_sizes.push_back(str_length);
        } else if (str_length > col_sizes[c]) {
          col_sizes[c] = str_length;
        }
      }
    }
    if (maxcol == 0) {
      return "";  // Empty rows
    }

    vector<string> output(nrows);
    vector<string> row;
    string col;
    for (size_t r = 0; r < nrows; r++) {
      row.clear();
      // Last column does not need padding
      for (size_t c = 0; c < maxcol - 1; c++) {
        if (col_sizes[c] > 0) {
          col = (c == 0 ? "  " : "") + aurostd::PaddedPOST((c < table[r].size() ? table[r][c] : ""), col_sizes[c]);
          row.push_back(col);
        }
      }
      if (table[r].size() == maxcol) {
        row.push_back(table[r].back());
      }
      output[r] = aurostd::joinWDelimiter(row, " | ");
    }
    return aurostd::joinWDelimiter(output, "\n");
  }

  /// @brief Display results of a unit test.
  ///
  /// @param xchk Unit test object containing all results.
  void UnitTest::displayResult(const xcheck& xchk) {
    stringstream message;
    const size_t check_num = xchk.results.size();
    if (xchk.passed_checks == check_num) {
      if (!xchk.errors.empty()) {
        // All attempted checks passed, but there were errors.
        // This happens when a prerequisite for a test fails (e.g. file loading).
        message << "FAIL " << xchk.task_description << " due to runtime errors" << std::endl;
        pflow::logger(__AFLOW_FILE__, xchk.function_name, message, aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
      } else {
        message << "SUCCESS " << xchk.task_description << " (passing " << check_num << " check" << ((check_num == 1) ? "" : "s") << ")" << std::endl;
        pflow::logger(__AFLOW_FILE__, xchk.function_name, message, aflags, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
      }
    } else {
      message << "FAIL " << xchk.task_description << " (" << (check_num - xchk.passed_checks) << " of " << check_num << " checks failed)" << std::endl;
      pflow::logger(__AFLOW_FILE__, xchk.function_name, message, aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
    }
    pflow::logger(__AFLOW_FILE__, xchk.function_name, formatResultsTable(xchk.results), aflags, *p_FileMESSAGE, *p_oss, _LOGGER_RAW_);
    if (!xchk.errors.empty()) {
      message << "\nAdditional error messages:\n" << aurostd::joinWDelimiter(xchk.errors, "\n");
      pflow::logger(__AFLOW_FILE__, xchk.function_name, message, aflags, *p_FileMESSAGE, *p_oss, _LOGGER_RAW_);
    }
  }
} // namespace unittest

// Collection of generic check functions, to streamline testing.
namespace unittest {
  template <typename utype>
  void UnitTest::checkEqual(const vector<utype>& calculated, const vector<utype>& expected, const string& check_function, const string& check_description, uint& passed_checks, vector<vector<string>>& results) {
    bool passed = (calculated.size() == expected.size());
    for (size_t i = 0; i < calculated.size() && passed; i++) {
      passed = aurostd::isequal(calculated[i], expected[i]);
    }
    check(passed, calculated, expected, check_function, check_description, passed_checks, results);
  }

  // AZ20231030 added vector of strings variant of check equal
  void UnitTest::checkEqual(const vector<string>& calculated, const vector<string>& expected, const string& check_function, const string& check_description, uint& passed_checks, vector<vector<string>>& results) {
    bool passed = (calculated.size() == expected.size());
    for (size_t i = 0; i < calculated.size() && passed; i++) {
      passed = (calculated[i] == expected[i]);
    }
    check(passed, calculated, expected, check_function, check_description, passed_checks, results);
  }

  void UnitTest::checkEqual(const vector<bool>& calculated, const vector<bool>& expected, const string& check_function, const string& check_description, uint& passed_checks, vector<vector<string>>& results) {
    bool passed = (calculated.size() == expected.size());
    for (size_t i = 0; i < calculated.size() && passed; i++) {
      passed = (calculated[i] == expected[i]);
    }
    check(passed, calculated, expected, check_function, check_description, passed_checks, results);
  }

  template <typename utype>
  void UnitTest::checkEqual(const utype& calculated, const utype& expected, const string& check_function, const string& check_description, uint& passed_checks, vector<vector<string>>& results) {
    const bool passed = (aurostd::isequal(calculated, expected));
    check(passed, calculated, expected, check_function, check_description, passed_checks, results);
  }

  void UnitTest::checkEqual(const string& calculated, const string& expected, const string& check_function, const string& check_description, uint& passed_checks, vector<vector<string>>& results) {
    const bool passed = (calculated == expected);
    check(passed, calculated, expected, check_function, check_description, passed_checks, results);
  }

  template <typename utype>
  void UnitTest::check(const bool passed, const vector<utype>& calculated, const vector<utype>& expected, const string& check_function, const string& check_description, uint& passed_checks, vector<vector<string>>& results) {
    check(passed, aurostd::joinWDelimiter(calculated, ","), aurostd::joinWDelimiter(expected, ","), check_function, check_description, passed_checks, results);
  }
  void UnitTest::check(const bool passed, const vector<double>& calculated, const vector<double>& expected, const string& check_function, const string& check_description, uint& passed_checks, vector<vector<string>>& results) {
    check(passed, aurostd::vecDouble2String(calculated), aurostd::vecDouble2String(expected), check_function, check_description, passed_checks, results);
  }

  template <typename utype>
  void UnitTest::check(const bool passed, const xmatrix<utype>& calculated, const xmatrix<utype>& expected, const string& check_function, const string& check_description, uint& passed_checks, vector<vector<string>>& results) {
    check(passed, aurostd::xmat2String(calculated), aurostd::xmat2String(expected), check_function, check_description, passed_checks, results);
  }
  void UnitTest::check(const bool passed, const xmatrix<double>& calculated, const xmatrix<double>& expected, const string& check_function, const string& check_description, uint& passed_checks, vector<vector<string>>& results) {
    check(passed, aurostd::xmatDouble2String(calculated), aurostd::xmatDouble2String(expected), check_function, check_description, passed_checks, results);
  }

  /// Checks the equality of two files by their hash.
  void UnitTest::checkFiles(const string& calculated, const string& expected, const string& check_function, const string& check_description, uint& passed_checks, vector<vector<string>>& results) {
    checkEqual(aurostd::file2hash(calculated), aurostd::file2hash(expected), check_function, check_description, passed_checks, results);
  }

  /// @brief Base function to check and update results.
  ///
  /// @param passed            Whether the test has passed.
  /// @param calculated        Calculated value.
  /// @param expected          Expected value.
  /// @param check_function    Function called for the test.
  /// @param check_description Description of the performed test.
  /// @param passed_checks     Number of passed checks.
  /// @param results           Results data - doubles as number of performed checks.
  template <typename utype>
  void UnitTest::check(const bool passed, const utype& calculated, const utype& expected, const string& check_function, const string& check_description, uint& passed_checks, vector<vector<string>>& results) {
    vector<string> result;
    const uint check_num = results.size() + 1;
    result.push_back(aurostd::utype2string<uint>(check_num));
    if (passed) {
      passed_checks++;
      result.emplace_back("pass");
    } else {
      result.emplace_back("FAIL");
    }
    result.push_back(check_function);
    result.push_back(check_description);
    if (!passed) {
      stringstream failstring;
      failstring << " (result: " << calculated << " | expected: " << expected << ")";
      result.back() += failstring.str();
    }
    results.push_back(result);
  }

  /// @brief Alternate function to check and update results.
  /// Avoids cluttered output in cases where displaying expected and calculated results is not helpful
  ///
  /// @param passed            Whether the test has passed.
  /// @param result_note       Note displayed explaining result
  /// @param check_function    Function called for the test.
  /// @param check_description Description of the performed test.
  /// @param passed_checks     Number of passed checks.
  /// @param results           Results data - doubles as number of performed checks.
  void UnitTest::check(const bool passed, const std::string& result_note, const std::string& check_function, const std::string& check_description, uint& passed_checks, std::vector<std::vector<std::string>>& results) {
    vector<string> result;
    const uint check_num = results.size() + 1;
    result.push_back(aurostd::utype2string<uint>(check_num));
    if (passed) {
      passed_checks++;
      result.emplace_back("pass");
    } else {
      result.emplace_back("FAIL");
    }
    result.push_back(check_function);
    result.push_back(check_description);
    if (!result_note.empty()) {
      result.back() += " (";
      result.back() += result_note;
      result.back() += ")";
    }
    results.push_back(result);
  }

} // namespace unittest

namespace unittest {

  void UnitTest::xscalarTest(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    (void) errors;  // Suppress compiler warnings
                    // setup test environment
    string check_function;
    string check_description;
    double calculated_dbl = 0.0;
    double expected_dbl = 0.0;
    int calculated_int = 0;
    int expected_int = 0;
    vector<int> calculated_vint;
    vector<int> expected_vint;

    // ---------------------------------------------------------------------------
    // Check | double2fraction conversion //DX20210908
    // ---------------------------------------------------------------------------
    check_function = "aurostd::double2fraction()";
    check_description = "convert a double to a fraction";

    const double test_double = 1.625;
    int numerator = 1;
    int denominator = 1;
    const string answer = "13/8";
    aurostd::double2fraction(test_double, numerator, denominator);
    stringstream result_ss;
    result_ss << numerator << "/" << denominator;

    checkEqual(result_ss.str(), answer, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | mod_floored (int) //SD20220124
    // ---------------------------------------------------------------------------
    check_function = "aurostd::mod_floored()";
    check_description = "floored mod; numbers as int";
    expected_int = -1;

    calculated_int = aurostd::mod_floored(5, -3);
    checkEqual(calculated_int, expected_int, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | mod_floored (double) //SD20220124
    // ---------------------------------------------------------------------------
    check_function = "aurostd::mod_floored()";
    check_description = "floored mod; numbers as double";
    expected_dbl = 1.4;

    calculated_dbl = aurostd::mod_floored(-5.2, 3.3);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | mod_floored (divisor 0) //SD20220124
    // ---------------------------------------------------------------------------
    check_function = "aurostd::mod_floored()";
    check_description = "floored mod; divisor is 0";
    expected_dbl = 11.11;

    calculated_dbl = aurostd::mod_floored(11.11, 0.0);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | mod_floored (divisor inf) //SD20220124
    // ---------------------------------------------------------------------------
    check_function = "aurostd::mod_floored()";
    check_description = "floored mod; divisor is inf";
    expected_dbl = 11.11;

    calculated_dbl = aurostd::mod_floored(11.11, (double) INFINITY);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | gcd //CO20190520
    // ---------------------------------------------------------------------------
    check_function = "aurostd::GCD()";
    int a = 0;
    int b = 0;
    int x1 = 0;
    int y1 = 0;
    int gcd = 0;

    check_description = "gcd(25,15)";
    a = 25;
    b = 15;
    expected_vint = {5, -1, 2};
    aurostd::GCD(a, b, gcd, x1, y1);
    calculated_vint = {gcd, x1, y1};
    checkEqual(calculated_vint, expected_vint, check_function, check_description, passed_checks, results);

    check_description = "gcd(25,0)";
    a = 25;
    b = 0;
    expected_vint = {25, 1, 0};
    aurostd::GCD(a, b, gcd, x1, y1);
    calculated_vint = {gcd, x1, y1};
    checkEqual(calculated_vint, expected_vint, check_function, check_description, passed_checks, results);

    check_description = "gcd(0,15)";
    a = 0;
    b = 15;
    expected_vint = {15, 0, 1};
    aurostd::GCD(a, b, gcd, x1, y1);
    calculated_vint = {gcd, x1, y1};
    checkEqual(calculated_vint, expected_vint, check_function, check_description, passed_checks, results);

    check_description = "gcd(-5100,30450)";
    a = -5100;
    b = 30450;
    expected_vint = {150, -6, -1};
    aurostd::GCD(a, b, gcd, x1, y1);
    calculated_vint = {gcd, x1, y1};
    checkEqual(calculated_vint, expected_vint, check_function, check_description, passed_checks, results);
  }

  void UnitTest::xvectorTest(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    (void) errors;  // Suppress compiler warnings
                    // setup test environment
    string check_function;
    string check_description;

    // HE20210511
    double expected_dbl = 0.0;
    double calculated_dbl = 0.0;
    xvector<double> calculated_xvecdbl;
    xvector<double> expected_xvecdbl;
    int expected_int = 0;
    string expected_str;
    vector<xvector<double>> points;
    vector<xvector<int>> ipoints;
    vector<vector<uint>> facets;
    vector<uint> facet;

    // Define test data
    points.clear();
    ipoints.clear();
    facets.clear();
    points.push_back({0.0, 0.0, 0.0});
    points.push_back({1.0, 0.0, 0.0});
    points.push_back({1.0, 1.0, 0.0});
    points.push_back({0.0, 1.0, 0.0});
    points.push_back({0.0, 0.0, 2.0});
    points.push_back({1.0, 0.0, 2.0});
    points.push_back({1.0, 1.0, 3.0});
    points.push_back({0.0, 1.0, 3.0});

    // covert points to integer to test special implementation
    for (size_t i_point = 0; i_point < points.size(); i_point++) {
      ipoints.push_back({(int) points[i_point][1], (int) points[i_point][2], (int) points[i_point][3]});
    }

    // define the facets
    facets.push_back({0, 1, 2, 3});
    facets.push_back({4, 5, 6, 7});
    facets.push_back({1, 2, 6, 5});
    facets.push_back({0, 3, 7});
    facets.push_back({0, 1, 5, 4});
    facets.push_back({3, 2, 6, 7});

    // ---------------------------------------------------------------------------
    // Check | convex solid volume (double)
    // ---------------------------------------------------------------------------
    check_function = "aurostd::volume()";
    check_description = "convex solid, points as doubles";
    expected_dbl = 2.5;

    calculated_dbl = aurostd::volume(points, facets, true);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | convex solid volume (int)
    // ---------------------------------------------------------------------------
    check_function = "aurostd::volume()";
    check_description = "convex solid, points as int";

    calculated_dbl = aurostd::volume(ipoints, facets, true);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    // define non convex solid
    points.clear();
    ipoints.clear();
    facets.clear();
    points.push_back({0.0, 0.0, 0.0});
    points.push_back({0.0, 4.0, 0.0});
    points.push_back({2.0, 4.0, 0.0});
    points.push_back({1.0, 1.0, 0.0});
    points.push_back({4.0, 2.0, 0.0});
    points.push_back({4.0, 0.0, 0.0});
    points.push_back({0.0, 0.0, 4.0});
    points.push_back({0.0, 4.0, 4.0});
    points.push_back({2.0, 4.0, 4.0});
    points.push_back({1.0, 1.0, 4.0});
    points.push_back({4.0, 2.0, 4.0});
    points.push_back({4.0, 0.0, 4.0});

    // covert points to integer to test special implementation
    for (size_t i_point = 0; i_point < points.size(); i_point++) {
      ipoints.push_back({(int) points[i_point][1], (int) points[i_point][2], (int) points[i_point][3]});
    }

    // define the facets
    facets.push_back({5, 4, 3, 2, 1, 0});
    facets.push_back({6, 7, 8, 9, 10, 11});
    facets.push_back({0, 6, 11, 5});
    facets.push_back({4, 5, 11, 10});
    facets.push_back({3, 4, 10, 9});
    facets.push_back({3, 9, 8, 2});
    facets.push_back({1, 2, 8, 7});
    facets.push_back({0, 1, 7, 6});

    // ---------------------------------------------------------------------------
    // Check | non convex solid volume (double)
    // ---------------------------------------------------------------------------
    check_function = "aurostd::volume()";
    check_description = "non convex solid, points as doubles";
    expected_dbl = 40.0;

    calculated_dbl = aurostd::volume(points, facets);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | error facet/normals mismatch
    // ---------------------------------------------------------------------------
    check_function = "aurostd::volume()";
    check_description = "error: facet/normals mismatch";
    const vector<xvector<double>> normals;
    expected_str = "xerror code 30 (VALUE_ERROR)";
    expected_int = _VALUE_ERROR_;

    try {
      calculated_dbl = aurostd::volume(points, facets, normals);
      check(false, std::string("no error"), expected_str, check_function, check_description, passed_checks, results);
    } catch (aurostd::xerror e) {
      if (e.whatCode() == expected_int) {
        check(true, "", "", check_function, check_description, passed_checks, results);
      } else {
        check(false, aurostd::utype2string(e.whatCode()), expected_str, check_function, check_description, passed_checks, results);
      }
    } catch (...) {
      check(false, std::string("not an xerror"), expected_str, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | non convex solid volume (int)
    // ---------------------------------------------------------------------------
    check_function = "aurostd::volume()";
    check_description = "non convex solid, points as int";
    expected_dbl = 40.0;

    calculated_dbl = aurostd::volume(ipoints, facets);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | error facet size
    // ---------------------------------------------------------------------------
    check_function = "aurostd::volume()";
    check_description = "error: wrong facet size";
    expected_str = "xerror code 30 (VALUE_ERROR)";
    expected_int = _VALUE_ERROR_;

    facet.clear();
    facets.push_back({1, 2});
    try {
      calculated_dbl = aurostd::volume(points, facets);
      check(false, std::string("no error"), expected_str, check_function, check_description, passed_checks, results);
    } catch (aurostd::xerror e) {
      if (e.whatCode() == expected_int) {
        check(true, "", "", check_function, check_description, passed_checks, results);
      } else {
        check(false, aurostd::utype2string(e.whatCode()), expected_str, check_function, check_description, passed_checks, results);
      }
    } catch (...) {
      check(false, std::string("not an xerror"), expected_str, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | numerical integration of convex solid with triangle faces and density function
    // ---------------------------------------------------------------------------
    {
      check_function = "convex_volume_integral()";
      check_description = "numerical integration of a convex solid";
      const std::function<double(const xvector<double>)> density = [](const xvector<double> point) -> double { return point[1] + point[2] + point[3]; };

      const vector<xvector<double>> points = {
          {0.0, 0.0, 0.0},
          {0.0, 0.0, 1.0},
          {0.0, 1.0, 1.0},
          {0.0, 1.0, 0.0},
          {1.0, 0.0, 0.0},
          {1.0, 0.0, 1.0},
          {1.0, 1.0, 1.0},
          {1.0, 1.0, 0.0}
      }; // cube
      vector<vector<uint>> facets = {
          {0, 1, 2},
          {0, 2, 3},
          {0, 1, 5},
          {0, 5, 4},
          {0, 3, 7},
          {0, 7, 4},
          {4, 5, 6},
          {4, 6, 7},
          {5, 6, 2},
          {5, 2, 1},
          {2, 3, 7},
          {2, 7, 6}
      };
      vector<xvector<double>> normals;
      const double expected = 1.5;
      double calculated;
      sortPolyhedronFacets(points, facets, normals);
      calculated = aurostd::convex_volume_integral_triangle(points, facets, density);
      checkEqual(calculated, expected, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | non convex area (double)
    // ---------------------------------------------------------------------------
    check_function = "aurostd::areaPointsOnPlane()";
    check_description = "non convex area; points as double";
    expected_dbl = 10.0;

    // fill vectors with data
    points.clear();
    ipoints.clear();
    facets.clear();
    points.push_back({0.0, 0.0, 0.0});
    points.push_back({0.0, 4.0, 0.0});
    points.push_back({2.0, 4.0, 0.0});
    points.push_back({1.0, 1.0, 0.0});
    points.push_back({4.0, 2.0, 0.0});
    points.push_back({4.0, 0.0, 0.0});

    // convert points to integer to test special implementation
    for (size_t i_point = 0; i_point < points.size(); i_point++) {
      ipoints.push_back({(int) points[i_point][1], (int) points[i_point][2], (int) points[i_point][3]});
    }
    //    points.push_back(p0); points.push_back(p1); points.push_back(p2); points.push_back(p3); points.push_back(p4);
    //    points.push_back(p5);
    //    ipoints.push_back(p0i); ipoints.push_back(p1i); ipoints.push_back(p2i); ipoints.push_back(p3i); ipoints.push_back(p4i);
    //    ipoints.push_back(p5i);

    calculated_dbl = aurostd::areaPointsOnPlane(points);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | non convex area (int)
    // ---------------------------------------------------------------------------
    check_function = "aurostd::areaPointsOnPlane()";
    check_description = "non convex area; points as int";
    expected_dbl = 10.0;

    calculated_dbl = aurostd::areaPointsOnPlane(ipoints);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    // define triangle in 3D to better test int handling
    points.clear();
    ipoints.clear();
    facets.clear();
    points.push_back({0.0, 0.0, 0.0});
    points.push_back({1.0, 1.0, 1.0});
    points.push_back({5.0, 0.0, 5.0});

    // convert points to integer to test special implementation
    for (size_t i_point = 0; i_point < points.size(); i_point++) {
      ipoints.push_back({(int) points[i_point][1], (int) points[i_point][2], (int) points[i_point][3]});
    }

    // ---------------------------------------------------------------------------
    // Check | 3d triangle area (double)
    // ---------------------------------------------------------------------------
    check_function = "aurostd::areaPointsOnPlane()";
    check_description = "3d triangle; points as double";
    expected_dbl = 3.5355339059;

    calculated_dbl = aurostd::areaPointsOnPlane(points);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | 3d triangle area (int)
    // ---------------------------------------------------------------------------
    check_function = "aurostd::areaPointsOnPlane()";
    check_description = "3d triangle; points as int";
    expected_dbl = 3.5355339059;

    calculated_dbl = aurostd::areaPointsOnPlane(ipoints);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | linspace //SD20220324
    // ---------------------------------------------------------------------------
    check_function = "aurostd::linspace()";
    check_description = "generate n linearly spaced points";
    expected_xvecdbl = {1.0, 1.375, 1.75, 2.125, 2.5};
    calculated_xvecdbl = aurostd::linspace(1.0, 2.5, 5);
    checkEqual(calculated_xvecdbl, expected_xvecdbl, check_function, check_description, passed_checks, results);

    {
      // ---------------------------------------------------------------------------
      // Check | histogram edges (manual bins) //AZ20230213
      // ---------------------------------------------------------------------------
      check_function = "aurostd::histogram()";
      check_description = "edges of histogram bins (manual bins)";
      expected_xvecdbl = {0., 1.8, 3.6, 5.4, 7.2, 9.};
      vector<xvector<double>> vec_xvecdbls;
      xvector<double> in_xvecdbl = {1., 3., 9., 8., 4., 7., 0., 2., 1., 9., 8., 3., 7., 4., 1., 0., 2., 8., 9., 7., 3., 4., 0., 1., 7., 8., 7., 6., 4., 3., 8., 7., 1., 2., 6., 4., 8., 1., 2., 0.,
                                    9., 3., 4., 8., 8., 7., 8., 1., 5., 3., 8., 3., 7., 4., 6., 1., 2., 8., 3., 4., 1., 9., 8., 3., 6., 4., 8., 7., 1., 6., 2., 4., 3., 0., 9., 7., 1., 2., 3., 4.};
      vec_xvecdbls = aurostd::histogram(in_xvecdbl, 5);
      checkEqual(vec_xvecdbls[1], expected_xvecdbl, check_function, check_description, passed_checks, results);

      // ---------------------------------------------------------------------------
      // Check | histogram bins //AZ20230213
      // ---------------------------------------------------------------------------
      check_description = "counts of histogram bins (manual bins)";
      expected_xvecdbl = {16., 18., 12., 15., 19.};
      checkEqual(vec_xvecdbls[0], expected_xvecdbl, check_function, check_description, passed_checks, results);

      // ---------------------------------------------------------------------------
      // Check | auto histogram edges (full vector) //AZ20230213
      // ---------------------------------------------------------------------------
      vec_xvecdbls = aurostd::histogram(in_xvecdbl, 5, 1);
      check_description = "automatic histogram bins (large vector)";
      expected_xvecdbl = {5., 11., 7., 11., 11., 1., 5., 10., 19.};
      checkEqual(vec_xvecdbls[0], expected_xvecdbl, check_function, check_description, passed_checks, results);

      // ---------------------------------------------------------------------------
      // Check | auto histogram edges (full vector) //AZ20230213
      // ---------------------------------------------------------------------------
      vec_xvecdbls = aurostd::histogram(in_xvecdbl, 5, 1);
      check_function = "aurostd::histogram()";
      check_description = "automatic histogram edges (large vector)";
      expected_xvecdbl = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9.};
      checkEqual(vec_xvecdbls[1], expected_xvecdbl, check_function, check_description, passed_checks, results);

      // ---------------------------------------------------------------------------
      // Check | auto histogram bins //AZ20230213
      // ---------------------------------------------------------------------------
      in_xvecdbl = {1., 3., 9., 8., 4., 7., 0., 2., 1., 9., 8., 3., 7., 4., 1., 0., 2.};
      vec_xvecdbls = aurostd::histogram(in_xvecdbl, 5, 1);
      check_function = "aurostd::histogram()";
      check_description = "automatic histogram bins (small vector)";
      expected_xvecdbl = {5., 4., 2., 2., 4.};
      checkEqual(vec_xvecdbls[0], expected_xvecdbl, check_function, check_description, passed_checks, results);

      check_function = "aurostd::histogram()";
      check_description = "automatic histogram edges (small vector)";
      expected_xvecdbl = {0., 1.8, 3.6, 5.4, 7.2, 9.};
      checkEqual(vec_xvecdbls[1], expected_xvecdbl, check_function, check_description, passed_checks, results);
    }
  }

  void UnitTest::xmatrixTest(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    (void) errors;  // Suppress compiler warnings
                    // setup test environment
    string check_function;
    string check_description;

    {
      // ---------------------------------------------------------------------------
      // Check convert to xvector //AZ20220627
      // ---------------------------------------------------------------------------
      check_function = "aurostd::getxvec()";
      xmatrix<int> full_xmatint;
      const xmatrix<int> xmatint;
      // need this matrix to test slicing
      full_xmatint = xmatrix<int>(3, 4);
      full_xmatint = {
          {1,  2,  3,  4},
          {5,  6,  7,  8},
          {9, 10, 11, 12}
      };
      check_description = "getxvec() test for type conversion";
      xvector<int> expected_xvecint(3);
      xvector<int> calculated_xvecint(3);
      expected_xvecint = {1, 5, 9};
      calculated_xvecint = full_xmatint.getxmat(1, 3, 1, 1).getxvec();
      checkEqual(calculated_xvecint, expected_xvecint, check_function, check_description, passed_checks, results);

      // ---------------------------------------------------------------------------
      // Check | column xvector //AZ20220627
      // ---------------------------------------------------------------------------
      check_description = "get column xvector from xmatrix";
      calculated_xvecint = full_xmatint.getxvec(1, 3, 1, 1);
      checkEqual(calculated_xvecint, expected_xvecint, check_function, check_description, passed_checks, results);

      // ---------------------------------------------------------------------------
      // Check | row xvector //AZ20220627
      // ---------------------------------------------------------------------------
      check_description = "get row xvector from xmatrix";
      expected_xvecint = {1, 2, 3};
      calculated_xvecint = full_xmatint.getxvec(1, 1, 1, 3);
      checkEqual(calculated_xvecint, expected_xvecint, check_function, check_description, passed_checks, results);

      // ---------------------------------------------------------------------------
      // Check | 1x1 xvector //AZ20220627
      // ---------------------------------------------------------------------------
      check_description = "get a 1x1 vector from xmatrix";
      expected_xvecint = {12};
      calculated_xvecint = full_xmatint.getxvec(3, 3, 4, 4);
      checkEqual(calculated_xvecint, expected_xvecint, check_function, check_description, passed_checks, results);
    }
    // ---------------------------------------------------------------------------
    // Check | ehermite //CO20190520
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::getEHermite()";
      check_description = "calculate elementary Hermite transformation";
      const xmatrix<int> expected_xmatint = {
          {  5, -2},
          {-12,  5}
      };
      xmatrix<int> calculated_xmatint = xmatrix<int>(2, 2);
      aurostd::getEHermite(5, 12, calculated_xmatint);
      checkEqual(calculated_xmatint, expected_xmatint, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | equilibrateMatrix //SD20220505
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::equilibrateMatrix()";
      xmatrix<long double> icm(5, 5);
      // define a Hilbert matrix
      for (int i = 1; i <= icm.rows; i++) {
        for (int j = 1; j <= icm.cols; j++) {
          icm(i, j) = 1.0 / (i + j - 1.0);
        }
      }
      check_description = "pre-condition a Hilbert matrix, lowers the condition number";
      xmatrix<long double> em;
      xmatrix<long double> rm;
      xmatrix<long double> cm;
      aurostd::equilibrateMatrix(icm, em, rm, cm);
      const bool expected_bool = true;
      const bool calculated_bool = aurostd::isequal((long double) 1.0, aurostd::sign(aurostd::condition_number(icm) - aurostd::condition_number(em)));
      checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);
      check_description = "pre-condition a Hilbert matrix, finds the inverse";
      checkEqual(icm, aurostd::inverse(rm) * em * aurostd::inverse(cm), check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | solve a linear system //HE20220912
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::inverse(xmatrix) * xvector";
      check_description = "solve a simple linear system with shifted matrices ";
      const xvector<double> expected_xvecdouble = {5.0, 3.0, -2.0};
      xvector<double> calculated_xvecdouble;
      const xvector<double> b = {6.0, -4, 27};
      xmatrix<double> A = {
          {1.0, 1.0,  1.0},
          {0.0, 2.0,  5.0},
          {2.0, 5.0, -1.0}
      };
      bool shift_check = true;
      // inv(A)*b should be solvable even when the xmatrix and xvector have different index boundaries
      for (const int shift_row : {-2, -1, 0, 1, 2}) {
        for (const int shift_col : {-2, -1, 0, 1, 2}) {
          A.shift(shift_row, shift_col);
          calculated_xvecdouble = aurostd::inverse(A) * b;
          if (not aurostd::isequal(expected_xvecdouble, calculated_xvecdouble)) {
            shift_check = false;
            break;
          }
        }
      }
      check(shift_check, calculated_xvecdouble, expected_xvecdouble, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | shifting xmatrix
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::xmatrix.shift()";
      check_description = "shift the lower bounds of an xmatrix to new values";
      bool overall_check = true;
      xmatrix<double> A = {
          {1.0, 2.0, 3.0},
          {4.0, 5.0, 6.0},
          {7.0, 8.0, 9.0}
      };
      xmatrix<double> B = A;

      for (const std::pair<int, int> p : vector<std::pair<int, int>>({
               {  -2,   -4},
               {   0,    3},
               {5461, 8461}
      })) {
        A.shift(p.first, p.second);
        for (const std::pair<int, int> t : vector<std::pair<int, int>>({
                 {1, 1},
                 {2, 2},
                 {3, 2}
        })) {
          overall_check = overall_check && (A[p.first + t.first - 1][p.second + t.second - 1] == B[t.first][t.second]);
          if (!overall_check) {
            check(overall_check, A[p.first + t.first - 1][p.second + t.second - 1], B[t.first][t.second], check_function, check_description, passed_checks, results);
            break;
          }
          if (!overall_check) {
            break;
          }
        }
      }
      if (overall_check) {
        check(overall_check, "", "", check_function, check_description, passed_checks, results);
      }
    }

    // ---------------------------------------------------------------------------
    // Check | reshaping xmatrix
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::xmatrix.reshape()";
      check_description = "reshaping a xmatrix";
      bool overall_check = true;
      const vector<int> B = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
      xmatrix<int> A = {
          {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}
      };
      int c;
      int r;
      for (const std::pair<int, int>& p : vector<std::pair<int, int>>({
               {4, 3},
               {2, 6},
               {3, 4}
      })) {
        A.reshape(p.first, p.second);
        for (const int t : B) {
          r = (t - 1) / p.second + 1;
          c = t - (r - 1) * (p.second);
          overall_check = overall_check && (A[r][c] == t);
          if (!overall_check) {
            check(overall_check, A[r][c], t, check_function, check_description, passed_checks, results);
            break;
          }
        }
        if (!overall_check) {
          break;
        }
      }
      if (overall_check) {
        check(overall_check, "", "", check_function, check_description, passed_checks, results);
      }
    }

    // ---------------------------------------------------------------------------
    // Check | reshaping xvector
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::reshape()";
      check_description = "reshaping xvectors";
      bool overall_check = true;
      const vector<int> B = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
      const xvector<int> xv = aurostd::vector2xvector(B);
      aurostd::xmatrix<int> A;
      int c;
      int r;
      for (const std::pair<int, int>& p : vector<std::pair<int, int>>({
               {4, 3},
               {2, 6},
               {3, 4}
      })) {
        A = aurostd::reshape(xv, p.first, p.second);
        for (const int t : B) {
          r = (t - 1) / p.second + 1;
          c = t - (r - 1) * (p.second);
          overall_check = overall_check && (A[r][c] == t);
          if (!overall_check) {
            check(overall_check, A[r][c], t, check_function, check_description, passed_checks, results);
            break;
          }
        }
        if (!overall_check) {
          break;
        }
      }
      if (overall_check) {
        check(overall_check, "", "", check_function, check_description, passed_checks, results);
      }
    }

    // ---------------------------------------------------------------------------
    // Check | stacking xvector
    // ---------------------------------------------------------------------------
    {
      const xvector<int> xv = {1, 2, 3, 4};
      xmatrix<int> smat;
      const xmatrix<int> vsmat = {
          {1, 1, 1},
          {2, 2, 2},
          {3, 3, 3},
          {4, 4, 4}
      };
      const xmatrix<int> hsmat = {
          {1, 2, 3, 4},
          {1, 2, 3, 4},
          {1, 2, 3, 4}
      };
      check_function = "aurostd::vstack()";
      check_description = "stack xvectors vertically into a xmatrix";
      smat = aurostd::vstack(vector<xvector<int>>({xv, xv, xv}));
      checkEqual(smat, vsmat, check_function, check_description, passed_checks, results);
      check_function = "aurostd::hstack()";
      check_description = "stack xvectors horizontal into a xmatrix";
      smat = aurostd::hstack(vector<xvector<int>>({xv, xv, xv}));
      checkEqual(smat, hsmat, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | getSmithNormalForm
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::getSmithNormalForm()";
      check_description = "calculate the Smith form of the rank 3 inverse Hilbert matrix";
      xmatrix<int> U1;
      xmatrix<int> V1;
      xmatrix<int> S1;
      const xmatrix<int> A1 = {
          {  9,  -36,   30},
          {-36,  192, -180},
          { 30, -180,  180}
      };
      aurostd::getSmithNormalForm(A1, U1, V1, S1);
      checkEqual(S1, U1 * A1 * V1, check_function, check_description, passed_checks, results);
      check_description = "calculate the Smith form of the rank 5 inverse Hilbert matrix";
      xmatrix<long long int> U2;
      xmatrix<long long int> V2;
      xmatrix<long long int> S2;
      const xmatrix<long long int> A2 = {
          {   25,   -300,    1050,   -1400,    630},
          { -300,   4800,  -18900,   26880, -12600},
          { 1050, -18900,   79380, -117600,  56700},
          {-1400,  26880, -117600,  179200, -88200},
          {  630, -12600,   56700,  -88200,  44100}
      };
      aurostd::getSmithNormalForm(A2, U2, V2, S2);
      checkEqual(S2, U2 * A2 * V2, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | SVDecomposition_Jacobi
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::SVDecomposition_Jacobi()";
      check_description = "calculate the singular value decomposition of a 3x2 matrix";
      xmatrix<double> U;
      xmatrix<double> S;
      xmatrix<double> S_inv;
      xmatrix<double> V;
      const xmatrix<double> A1 = {
          {3, 2,  2},
          {2, 3, -2},
      };
      aurostd::SVDecomposition_Jacobi(A1, U, S, S_inv, V);
      checkEqual(U * S * aurostd::trasp(V), A1, check_function, check_description, passed_checks, results);
      check_description = "calculate the singular value decomposition of a 2x3 matrix";
      const xmatrix<double> A2 = {
          {3,  2},
          {2,  3},
          {2, -2},
      };
      aurostd::SVDecomposition_Jacobi(A2, U, S, S_inv, V);
      checkEqual(U * S * aurostd::trasp(V), A2, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | pseudoinverse
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::pseudoinverse()";
      check_description = "calculate the pseudoinverse of the rank 4 Hilbert matrix";
      const xmatrix<long double> A1 = {
          {1.0 / 1.0, 1.0 / 2.0, 1.0 / 3.0, 1.0 / 4.0},
          {1.0 / 2.0, 1.0 / 3.0, 1.0 / 4.0, 1.0 / 5.0},
          {1.0 / 3.0, 1.0 / 4.0, 1.0 / 5.0, 1.0 / 6.0},
          {1.0 / 4.0, 1.0 / 5.0, 1.0 / 6.0, 1.0 / 7.0},
      };
      checkEqual(A1 * aurostd::pseudoinverse(A1), aurostd::identity(A1 * aurostd::trasp(A1)), check_function, check_description, passed_checks, results);
    }
  }

  void UnitTest::aurostdMainTest(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    (void) errors;  // Suppress compiler warnings
    string check_function;
    string check_description;
    bool calculated_bool = false;
    bool expected_bool = false;
    int calculated_int = 0;
    int expected_int = 0;
    bool multi_check = true;
    uint64_t expected_uint64 = 0;
    uint64_t calculated_uint64 = 0;

    // ---------------------------------------------------------------------------
    // Check | substringlist2bool
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::substringlist2bool()";
      string str = "hawk owl";
      vector<string> strlist;

      strlist.emplace_back("falcon");
      check_description = aurostd::joinWDelimiter(strlist, " and ") + " in " + str + " (match all)";
      expected_bool = false;
      calculated_bool = aurostd::substringlist2bool(str, strlist, true);
      checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

      check_description = aurostd::joinWDelimiter(strlist, " or ") + " in " + str + " (match one)";
      expected_bool = false;
      calculated_bool = aurostd::substringlist2bool(str, strlist, false);
      checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

      strlist.emplace_back("owl");
      check_description = aurostd::joinWDelimiter(strlist, " and ") + " in " + str + " (match all)";
      expected_bool = false;
      calculated_bool = aurostd::substringlist2bool(str, strlist, true);
      checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

      check_description = aurostd::joinWDelimiter(strlist, " or ") + " in " + str + " (match one)";
      expected_bool = true;
      calculated_bool = aurostd::substringlist2bool(str, strlist, false);
      checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

      strlist[0] = "hawk";
      check_description = aurostd::joinWDelimiter(strlist, " and ") + " in " + str + " (match all)";
      expected_bool = true;
      calculated_bool = aurostd::substringlist2bool(str, strlist, true);
      checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

      check_description = aurostd::joinWDelimiter(strlist, " or ") + " in " + str + " (match one)";
      expected_bool = true;
      calculated_bool = aurostd::substringlist2bool(str, strlist, false);
      checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

      str = "falcon";
      check_description = "Find substring " + str + " in [" + aurostd::joinWDelimiter(strlist, ", ") + "]";
      expected_int = -1;
      aurostd::SubstringWithinList(strlist, str, calculated_int);
      checkEqual(calculated_int, expected_int, check_function, check_description, passed_checks, results);

      strlist.emplace_back("peregrine falcon");
      check_description = "Find substring " + str + " in [" + aurostd::joinWDelimiter(strlist, ", ") + "]";
      expected_int = (int) strlist.size() - 1;
      aurostd::SubstringWithinList(strlist, str, calculated_int);
      checkEqual(calculated_int, expected_int, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | substring2string //SD20220525
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::substring2string()";
      check_description = "return the third match of the substring";
      const string test_string = "_FILE_START_\nIALGO==48\nALGO==FAST\nIALGO==49\nALGO==MEDIUM\nIALGO==50\nALGO==SLOW\n_FILE_END_";
      const string calculated_string = aurostd::substring2string(test_string, "ALGO", 3);
      const string expected_string = "==49";
      checkEqual(calculated_string, expected_string, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | substring2string //SG20240401
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::substring2string()";
      check_description = "return the second match in between 2 substrings of a text";
      string test_string = "_FILE_START_\nIALGO==48\nALGO==FAST\nIALGO==49 //comment\nALGO==MEDIUM\nIALGO==50\nALGO==SLOW\n_FILE_END_";
      string calculated_string = aurostd::substring2string(test_string, "IALGO", "ALGO", 2, true, true);
      string expected_string = "==49";
      checkEqual(calculated_string, expected_string, check_function, check_description, passed_checks, results);

      check_function = "aurostd::substring2string()";
      check_description = "return the second match in between 2 substrings with a nested substring between them";
      test_string = "_FILE_START_\nX==48\nY==FAST\nX==49 //comment\nX==MEDIUM\nY==50\nX==SLOW\n_FILE_END_";
      calculated_string = aurostd::substring2string(test_string, "X", "Y", 2, true, true);
      expected_string = "==49\nX==MEDIUM";
      checkEqual(calculated_string, expected_string, check_function, check_description, passed_checks, results);

      check_description = "return the third-to-last match in between 2 substrings of the substring";
      test_string = "_FILE_START_\nIALGO==48\nALGO==FAST\nIALGO==49 //comment\nALGO==MEDIUM\nIALGO==50\nALGO==SLOW\n_FILE_END_";
      calculated_string = aurostd::substring2string(test_string, "IALGO", "ALGO", -3, true, true);
      expected_string = "==48";
      checkEqual(calculated_string, expected_string, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | kvpair2string //SD20220525
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::kvpair2string()";
      check_description = "return the second match of the kvpair";
      const string test_string = "_FILE_START_\nIALGO==48\nALGO==FAST\nIALGO==49\nALGO==MEDIUM\nIALGO==50\nALGO==SLOW\n_FILE_END_";
      const string calculated_string = aurostd::kvpair2string(test_string, "ALGO", "==", 2);
      const string expected_string = "MEDIUM";
      checkEqual(calculated_string, expected_string, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | substring2string //SD20220525
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::substring2string()";
      check_description = "return the last match of the substring";
      const string test_string = "_FILE_START_\nIALGO==48\nALGO==FAST\nIALGO==49\nALGO==MEDIUM\nIALGO==50\nALGO==SLOW\n_FILE_END_";
      const string calculated_string = aurostd::substring2string(test_string, "ALGO", -1);
      const string expected_string = "==SLOW";
      checkEqual(calculated_string, expected_string, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | substring2strings //SD20230322
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::substring2strings()";
      check_description = "return all the matches as vector of strings";
      const string test_string = "_FILE_START_\nIALGO==48\nALGO==FAST\nIALGO==49\nALGO==MEDIUM\nIALGO==50\nALGO==SLOW\n_FILE_END_";
      vector<string> calculated_vstring;
      aurostd::substring2strings(test_string, calculated_vstring, "IALGO");
      const vector<string> expected_vstring = {"==48", "==49", "==50"};
      checkEqual(calculated_vstring, expected_vstring, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | substring2strings //CO20230502
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::substring2strings()";
      check_description = "return all the matches as vector of strings";
      const string test_string = "[START]1234\n 4321[STOP] 456[START] 789 \n 987 [STOP]\n[START]4321 [STOP]654[START]987[STOP]";
      vector<string> calculated_vstring;
      aurostd::substring2strings(test_string, calculated_vstring, string("[START]"), string("[STOP]"), -3); // CO20230502 - without string(), compiler picks the wrong overload...
      const vector<string> expected_vstring = {"789 \n 987", "4321", "987"};
      checkEqual(calculated_vstring, expected_vstring, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | kvpair2string //SD20220525
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::kvpair2string()";
      check_description = "return the last match of the kvpair";
      const string test_string = "_FILE_START_\nIALGO==48\nALGO==FAST\nIALGO==49\nALGO==MEDIUM\nIALGO==50\nALGO==SLOW\n_FILE_END_";
      const string calculated_string = aurostd::kvpair2string(test_string, "ALGO", "==", -1);
      const string expected_string = "SLOW";
      checkEqual(calculated_string, expected_string, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | string2utype //HE20220324
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::string2utype()";
      check_description = "bool - bases 16, 10, 8, 5, 2";
      multi_check = true;
      multi_check = (multi_check && (aurostd::string2utype<bool>("-420") == true));
      multi_check = (multi_check && (aurostd::string2utype<bool>("-420.23") == true));
      multi_check = (multi_check && (aurostd::string2utype<bool>("922337203685477432") == true));
      multi_check = (multi_check && (aurostd::string2utype<bool>("T") == true));
      multi_check = (multi_check && (aurostd::string2utype<bool>("true", 2) == true));
      multi_check = (multi_check && (aurostd::string2utype<bool>(".true.", 16) == true));
      multi_check = (multi_check && (aurostd::string2utype<bool>("F", 8) == false));
      multi_check = (multi_check && (aurostd::string2utype<bool>("false", 5) == false));
      multi_check = (multi_check && (aurostd::string2utype<bool>(".false.") == false));
      multi_check = (multi_check && (aurostd::string2utype<bool>("aflow", 16)) == true);
      multi_check = (multi_check && (aurostd::string2utype<bool>("0.0", 16)) == false);
      multi_check = (multi_check && (aurostd::string2utype<bool>("0.0000")) == false);
      multi_check = (multi_check && (aurostd::string2utype<bool>("0", 2)) == false);
      checkEqual(multi_check, true, check_function, check_description, passed_checks, results);

      check_function = "aurostd::string2utype()";
      check_description = "int - bases 16, 10, 8, 5, 2";
      multi_check = true;
      multi_check = (multi_check && (aurostd::string2utype<int>("-420") == -420));
      multi_check = (multi_check && (aurostd::string2utype<int>("-420.23") == -420));
      multi_check = (multi_check && (aurostd::string2utype<long long int>("922337203685477432") == 922337203685477432));
      multi_check = (multi_check && (aurostd::string2utype<int>("T") == 1));
      multi_check = (multi_check && (aurostd::string2utype<int>("true") == 1));
      multi_check = (multi_check && (aurostd::string2utype<int>(".true.") == 1));
      multi_check = (multi_check && (aurostd::string2utype<int>("-420", 16)) == -1056);
      multi_check = (multi_check && (aurostd::string2utype<int>("-0x420", 16)) == -1056);
      multi_check = (multi_check && (aurostd::string2utype<int>("-0X420", 16)) == -1056);
      multi_check = (multi_check && (aurostd::string2utype<int>("-420", 16)) == -1056);
      multi_check = (multi_check && (aurostd::string2utype<int>("-420", 8)) == -272);
      multi_check = (multi_check && (aurostd::string2utype<int>("-0420", 8)) == -272);
      multi_check = (multi_check && (aurostd::string2utype<int>("-420", 5)) == -110);
      multi_check = (multi_check && (aurostd::string2utype<int>("-110100100", 2)) == -420);
      multi_check = (multi_check && (aurostd::string2utype<int>("aflow", 2)) == 0);
      multi_check = (multi_check && (aurostd::string2utype<int>("   45.114455.2154.2145", 10)) == 45);
      multi_check = (multi_check && (aurostd::string2utype<int>("   45.1syfgdhen", 10)) == 45);
      multi_check = (multi_check && (aurostd::string2utype<int>("   45.1syfgdhen", 8)) == 37);
      multi_check = (multi_check && (aurostd::string2utype<int>("gf hsdg", 10)) == 0);

      checkEqual(multi_check, true, check_function, check_description, passed_checks, results);
      check_description = "float - bases 16, 10, 8, 5, 2";
      multi_check = true;
      multi_check = (multi_check && (aurostd::string2utype<float>("-4.20") == -4.20F));
      multi_check = (multi_check && (aurostd::string2utype<float>("-420", 16)) == -1056.0F);
      multi_check = (multi_check && (aurostd::string2utype<float>("-420", 8)) == -272.0F);
      multi_check = (multi_check && (aurostd::string2utype<float>("-420.35", 8)) == -272.0F);
      multi_check = (multi_check && (aurostd::string2utype<float>("-420", 5)) == -110.0F);
      multi_check = (multi_check && (aurostd::string2utype<float>("-110100100", 2)) == -420.0F);
      multi_check = (multi_check && (aurostd::string2utype<float>("   45.1syfgdhen", 10)) == 45.1F);
      checkEqual(multi_check, true, check_function, check_description, passed_checks, results);

      check_description = "double - bases 16, 10, 8, 5, 2";
      multi_check = true;
      multi_check = (multi_check && (aurostd::string2utype<double>("-4.20") == -4.20));
      multi_check = (multi_check && (aurostd::string2utype<double>("-420", 16)) == -1056.0);
      multi_check = (multi_check && (aurostd::string2utype<double>("-420", 8)) == -272.0);
      multi_check = (multi_check && (aurostd::string2utype<double>("-420", 5)) == -110.0);
      multi_check = (multi_check && (aurostd::string2utype<double>("-110100100", 2)) == -420.0);
      checkEqual(multi_check, true, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | crc64 //HE20220404
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::crc64()";
      check_description = "runtime hashing";
      expected_uint64 = 15013402708409085989UL;
      calculated_uint64 = aurostd::crc64("aflowlib_date");
      checkEqual(calculated_uint64, expected_uint64, check_function, check_description, passed_checks, results);

      check_function = "aurostd::ctcrc64()";
      check_description = "compiler hashing (constexpr)";
      static constexpr uint64_t calculated_const_uint64 = aurostd::ctcrc64("aflowlib_date");
      checkEqual(calculated_const_uint64, expected_uint64, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | crc2human //HE20230221
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::crc2human()";
      check_description = "create human readable hash";
      std::string expected = "NEPPW041L8FJ2";
      std::string result = aurostd::crc2human("Hello World!");
      checkEqual(result, expected, check_function, check_description, passed_checks, results);

      check_description = "create truncated human readable hash";
      expected = "NEPPW0";
      result = aurostd::crc2human("Hello World!", 6);
      checkEqual(result, expected, check_function, check_description, passed_checks, results);

      check_description = "check padding";
      expected = "2ZP300";
      result = aurostd::crc2human(145624, 6);
      checkEqual(result, expected, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | file2hash //HE20240908
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::file2hash()";
      check_description = "hash a file using openSSL with ";
      const std::string tmp_path = aurostd::TmpStrCreate();
      string result;
      aurostd::string2file("Some content to be hashed", tmp_path);
      const std::map<std::string, std::string> expected = {
          {   "MD5",                                 "7f9da88e81fdddfca1c7ee06efb65138"},
          {  "SHA1",                         "ff6c23447b2652e6b583f8f930f1114ca681aa6c"},
          {"SHA256", "69e89d1ea0e2b67fd21a9651a092a7c572ad4fb27ed6cf5cd3f93f4fca37026b"}
      };
      for (const auto& hash_system : expected) {
        result = aurostd::file2hash(tmp_path, hash_system.first);
        checkEqual(result, hash_system.second, check_function, check_description + hash_system.first, passed_checks, results);
      }
      aurostd::RemoveFile(tmp_path);
    }

    // ---------------------------------------------------------------------------
    // Check | sting2hash ///HE20240908
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::sting2hash()";
      check_description = "hash a string using openSSL with ";
      string result;
      const std::map<std::string, std::string> expected = {
          {   "MD5",                                 "7f9da88e81fdddfca1c7ee06efb65138"},
          {  "SHA1",                         "ff6c23447b2652e6b583f8f930f1114ca681aa6c"},
          {"SHA256", "69e89d1ea0e2b67fd21a9651a092a7c572ad4fb27ed6cf5cd3f93f4fca37026b"}
      };
      for (const auto& hash_system : expected) {
        result = aurostd::string2hash("Some content to be hashed", hash_system.first);
        checkEqual(result, hash_system.second, check_function, check_description + hash_system.first, passed_checks, results);
      }
    }

    // ---------------------------------------------------------------------------
    // Check | aurostd::execute2OutErrPair
    // ---------------------------------------------------------------------------

    // execute command with large std_out and std_err
    {
      check_function = "aurostd::execute2OutErrPair()";
      check_description = "reading 256kBytes from /dev/zero";
      const std::string std_out_expected(262144, 'a'); // 256 kBytes
      const std::string std_err_expected(262144, 'b');
      const std::pair<std::string, std::string> output = aurostd::execute2OutErrPair(R"(head -c 262144 /dev/zero | tr "\000" "\141"; head -c 262144 /dev/zero | tr "\000" "\142" 1>&2)");
      checkEqual(output.first, std_out_expected, check_function, check_description + "std_out", passed_checks, results);
      checkEqual(output.second, std_err_expected, check_function, check_description + "std_err", passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | aurostd::substring_present_file
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::substring_present_file()";
      check_description = "check if substring is present in file";
      const string placeholder
          = "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi "
            "ut aliquip ex ea commodo consequat.";
      const string substring = "exercitation";
      const string filename = aurostd::TmpFileCreate();
      aurostd::string2file(placeholder, filename);
      bool expected = true;
      bool result = aurostd::substring_present_file(filename, substring);
      aurostd::RemoveFile(filename);
      checkEqual(result, expected, check_function, check_description, passed_checks, results);
      check_description = "check if substring is not present in file";
      aurostd::string2file(filename, filename);
      expected = false;
      result = aurostd::substring_present_file(filename, substring);
      aurostd::RemoveFile(filename);
      checkEqual(result, expected, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | aurostd::substrings_present_file
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::substrings_present_file()";
      check_description = "check if substrings are present in file (vector)";
      const string placeholder
          = "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi "
            "ut aliquip ex ea commodo consequat.";
      const vector<string> substrings = {"exercitation", "dolor"};
      const string filename = aurostd::TmpFileCreate();
      aurostd::string2file(placeholder, filename);
      vector<bool> expected = {true, true};
      vector<bool> result = aurostd::substrings_present_file(filename, substrings);
      aurostd::RemoveFile(filename);
      checkEqual(result, expected, check_function, check_description, passed_checks, results);
      check_description = "check if substrings are not present in file (vector)";
      aurostd::string2file(filename, filename);
      expected = {false, false};
      result = aurostd::substrings_present_file(filename, substrings);
      aurostd::RemoveFile(filename);
      checkEqual(result, expected, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | aurostd::substrings_map_present_file
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::substrings_map_present_file()";
      check_description = "check if substrings are present in file (unordered_map)";
      const string placeholder
          = "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi "
            "ut aliquip ex ea commodo consequat.";
      std::unordered_map<string, bool> substrings_map = {
          {"exercitation", false},
          {       "dolor", false}
      };
      const string filename = aurostd::TmpFileCreate();
      aurostd::string2file(placeholder, filename);
      vector<bool> expected = {true, true};
      aurostd::substrings_map_present_file(filename, substrings_map);
      vector<bool> result;
      for (const auto& entry : substrings_map) {
        result.emplace_back(entry.second);
      }
      aurostd::RemoveFile(filename);
      checkEqual(result, expected, check_function, check_description, passed_checks, results);
      check_description = "check if substrings are not present in file (unordered_map)";
      aurostd::string2file(filename, filename);
      substrings_map = {
          {"exercitation", false},
          {       "dolor", false}
      };
      expected = {false, false};
      aurostd::substrings_map_present_file(filename, substrings_map);
      result.clear();
      for (const auto& entry : substrings_map) {
        result.emplace_back(entry.second);
      }
      aurostd::RemoveFile(filename);
      checkEqual(result, expected, check_function, check_description, passed_checks, results);
    }
  }

  void UnitTest::xfileTest(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    (void) errors;  // Suppress compiler warnings
    // setup test environment
    string check_function;
    string check_description;

    // ---------------------------------------------------------------------------
    // Check | CopyDirectory and CopyFile //ST20251208 //HE20260129
    // ---------------------------------------------------------------------------
    {
      bool check_passed = false;
      std::string check_note;
      const string path_from = aurostd::TmpDirectoryCreate();
      const string path_to = aurostd::TmpDirectoryCreate();
      aurostd::DirectoryMake(path_from + "/extra");
      const std::vector<std::string> filenames = {"/A.txt", "/B.txt", "/C.txt", "/extra/D.txt"};
      const string filename_from = path_from + filenames[0];
      const string filename_to = path_to + filenames[0];
      constexpr size_t COLS = 80;
      constexpr size_t LINES = 1000;
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_int_distribution<> dist(33, 126);
      for (const std::string& filename : filenames) {
        std::ofstream out{path_from + filename, std::ios::binary};
        for (size_t i = 0; i < LINES; i++) {
          std::string line;
          line.reserve(COLS + 1);
          for (size_t j = 0; j < COLS; j++) {
            line.push_back(static_cast<char>(dist(gen)));
          }
          line.push_back('\n');
          out << line;
        }
        out.close();
      }

      auto read_lambda = [](const std::string& filename) -> std::string {
        std::ifstream ifs{filename};
        std::ostringstream oss;
        oss << ifs.rdbuf();
        return oss.str();
      };

      // gives True if they don't match to enable the use of std::any_of for dir_lambda
      auto compare_lambda = [&path_from, &path_to, &read_lambda](const std::string& filename) -> bool { return read_lambda(path_from + filename) != read_lambda(path_to + filename); };

      // gives True if they don't match to be the same as compare_lambda
      auto dir_lambda = [&compare_lambda, &filenames]() -> bool { return std::any_of(filenames.cbegin(), filenames.cend(), compare_lambda); };

      check_description = "copy a file";
      check_function = "aurostd::CopyFile()";
      check_note = "";
      check_passed = aurostd::CopyFile(filename_from, filename_to);
      if (!check_passed) {
        check_note = "copy function signaled a problem";
      } else {
        check_passed = !compare_lambda(filenames[0]);
        if (!check_passed) {
          check_note = "file content does not match";
        }
      }
      check(check_passed, check_note, check_function, check_description, passed_checks, results);

      aurostd::RemoveFile(filename_to);

      check_function = "aurostd::CopyFile_safer()";
      check_note = "";
      check_passed = aurostd::CopyFile_safer(filename_from, filename_to);
      if (!check_passed) {
        check_note = "copy function signaled a problem";
      } else {
        check_passed = !compare_lambda(filenames[0]);
        if (!check_passed) {
          check_note = "file content does not match";
        }
      }
      check(check_passed, check_note, check_function, check_description, passed_checks, results);
      aurostd::RemoveDirectory(path_to);

      check_description = "copy a folder";
      check_function = "aurostd::CopyDirectory()";
      check_note = "";
      check_passed = aurostd::CopyDirectory(path_from, path_to);
      if (!check_passed) {
        check_note = "copy function signaled a problem";
      } else {
        check_passed = !dir_lambda();
        if (!check_passed) {
          check_note = "file content does not match";
        }
      }
      check(check_passed, check_note, check_function, check_description, passed_checks, results);
      aurostd::RemoveDirectory(path_to);

      check_description = "copy a folder";
      check_function = "aurostd::CopyDirectory_safer()";
      check_note = "";
      check_passed = aurostd::CopyDirectory_safer(path_from, path_to);
      if (!check_passed) {
        check_note = "copy function signaled a problem";
      } else {
        check_passed = !dir_lambda();
        if (!check_passed) {
          check_note = "file content does not match";
        }
      }
      check(check_passed, check_note, check_function, check_description, passed_checks, results);
      aurostd::RemoveDirectory(path_to);
      aurostd::RemoveDirectory(path_from);
    }

    // ---------------------------------------------------------------------------
    // Check | FileEmpty //SD20240313
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::FileEmpty()";
      check_description = "check if file is empty";
      const string filename = aurostd::TmpFileCreate();
      const bool expected = true;
      const bool result = aurostd::FileEmpty(filename);
      aurostd::RemoveFile(filename);
      checkEqual(result, expected, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | FileNotEmpty //SD20240313
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::FileNotEmpty()";
      check_description = "check if file is not empty";
      const string filename = aurostd::TmpFileCreate();
      aurostd::string2file(filename, filename);
      const bool expected = true;
      const bool result = aurostd::FileNotEmpty(filename);
      aurostd::RemoveFile(filename);
      checkEqual(result, expected, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | string2file and file2string variants //HE20241217
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::string2file()";
      check_description = "write + read - simple file";
      string filename = aurostd::TmpFileCreate();
      aurostd::string2file(filename, filename);
      string result = aurostd::file2string(filename);
      aurostd::RemoveFile(filename);
      checkEqual(result, filename, check_function, check_description, passed_checks, results);

      check_description = "write + read - compressed XZ file";
      aurostd::string2file(filename, filename, aurostd::compression_type::XZ);
      result = aurostd::file2string(filename + aurostd::compression_suffix[aurostd::compression_type::XZ]);
      aurostd::RemoveFile(filename + aurostd::compression_suffix[aurostd::compression_type::XZ]);
      checkEqual(result, filename, check_function, check_description, passed_checks, results);

      check_description = "write + read - create parent folder";
      const string folder = aurostd::TmpStrCreate("", "", false, true);
      filename = folder + "/xfile.test";
      aurostd::string2file(filename, filename);
      result = aurostd::file2string(filename);
      std::filesystem::remove_all(folder);
      checkEqual(result, filename, check_function, check_description, passed_checks, results);
    }

    {
      // ---------------------------------------------------------------------------
      // Check | CompressFile/DecompressFile //SD20240326
      // ---------------------------------------------------------------------------
      check_function = "aurostd::CompressFile()";
      check_description = "compress a file using ";
      const vector<aurostd::compression_type> compressions = {
          aurostd::compression_type::BZ2, aurostd::compression_type::GZ, aurostd::compression_type::XZ, aurostd::compression_type::ZIP, aurostd::compression_type::ZSTD};
      const string expected = aurostd::TmpStrCreate();
      string result;
      aurostd::string2file(expected, expected);
      for (const auto& compression : compressions) {
        try {
          aurostd::CompressFile(expected, compression);
          aurostd::DecompressFile(expected + aurostd::compression_suffix[compression]);
          result = aurostd::file2string(expected);
        } catch (aurostd::xerror& err) {
          aurostd::string2file(expected, expected);
          result = err.buildMessageString();
        }
        checkEqual(result, expected, check_function, check_description + aurostd::compression_suffix[compression], passed_checks, results);
      }
      aurostd::RemoveFile(expected);

      // ---------------------------------------------------------------------------
      // Check | string2compressfile/compressfile2string //SD20240326
      // ---------------------------------------------------------------------------
      check_function = "aurostd::string2compressfile()";
      check_description = "compress a string using ";
      for (const auto& compression : compressions) {
        try {
          aurostd::string2compressfile(expected, expected, compression);
          aurostd::compressfile2string(expected + aurostd::compression_suffix[compression], result);
        } catch (aurostd::xerror& err) {
          result = err.buildMessageString();
        }
        aurostd::RemoveFile(expected + aurostd::compression_suffix[compression]);
        checkEqual(result, expected, check_function, check_description + aurostd::compression_suffix[compression], passed_checks, results);
      }
    }

    // ---------------------------------------------------------------------------
    // Check | compressfile2compressfile //SD20240329
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::compressfile2compressfile()";
      check_description = "convert compression from bz2 to xz";
      const string expected = aurostd::TmpStrCreate();
      string result;
      aurostd::string2file(expected, expected);
      try {
        aurostd::CompressFile(expected, aurostd::compression_type::BZ2);
        aurostd::compressfile2compressfile(expected + ".bz2", aurostd::compression_type::XZ);
        result = aurostd::compressfile2string(expected + ".xz");
      } catch (aurostd::xerror& err) {
        result = err.buildMessageString();
      }
      aurostd::RemoveFile(expected + ".xz");
      checkEqual(result, expected, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | CompressFiles/DecompressFile //SD20240331
    // changed to check relativ path compression //HE20240908
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::CompressFiles()";
      check_description = "compress multiple files";
      const string expected;
      string result;
      vector<fs::path> filenames = {aurostd::TmpStrCreate(), aurostd::TmpStrCreate(), aurostd::TmpStrCreate() + ".zip"};
      const fs::path directory(aurostd::TmpStrCreate());
      aurostd::string2file(filenames[0], filenames[0]);
      aurostd::string2file(filenames[1], filenames[1]);
      try {
        const fs::path archive_path = aurostd::CompressFiles({filenames[0], filenames[1]}, filenames[0].parent_path(), filenames[2], aurostd::compression_type::ZIP);
        aurostd::DecompressFiles(archive_path, directory);
        result = (fs::exists(directory / filenames[0].filename()) && fs::exists(directory / filenames[1].filename())) ? "" : "Not an xerror";
      } catch (aurostd::xerror& err) {
        result = err.buildMessageString();
      }
      aurostd::RemoveDirectory(directory);
      checkEqual(result, expected, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | Check if compressed //HE20240506
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::IsCompressed()";
      check_description = "check if compressed";
      const string expected;
      string result;
      const vector<aurostd::compression_type> compressions = {
          aurostd::compression_type::BZ2, aurostd::compression_type::GZ, aurostd::compression_type::XZ, aurostd::compression_type::ZIP, aurostd::compression_type::ZSTD};
      const string base = aurostd::TmpStrCreate() + ".txt";
      checkEqual(aurostd::IsCompressed(base), false, check_function, check_description + " (txt)", passed_checks, results);
      aurostd::string2file(base, base);
      for (const auto& compression : compressions) {
        try {
          aurostd::CompressFile(base, compression, true);
          result = aurostd::IsCompressed(base + aurostd::compression_suffix[compression]) ? "" : "Not an xerror";
        } catch (aurostd::xerror& err) {
          result = err.buildMessageString();
        }
        checkEqual(result, expected, check_function, check_description + " (" + aurostd::compression_suffix[compression] + ")", passed_checks, results);
        aurostd::RemoveFile(base + aurostd::compression_suffix[compression]);
      }
      aurostd::RemoveFile(base);
    }

    // ---------------------------------------------------------------------------
    // Check | file modification timestamps
    // ---------------------------------------------------------------------------
    {
      check_function = "aurostd::SecondsSinceFileModified()";
      check_description = "check duration";
      const string folder = aurostd::TmpStrCreate("", "", false, true);
      string base = folder + "/xfile.test";
      constexpr long int test_time = 3;
      aurostd::string2file(base, base);
      std::this_thread::sleep_for(std::chrono::seconds(test_time));
      long int eclipsed = aurostd::SecondsSinceFileModified(base);
      // allow off by one error
      if (eclipsed >= test_time and eclipsed <= test_time + 1) {
        checkEqual(true, true, check_function, check_description, passed_checks, results);
      } else {
        checkEqual(eclipsed, test_time, check_function, check_description, passed_checks, results);
      }
      check_function = "aurostd::SecondsSinceDirectoryContentModified()";
      check_description = "check duration";
      base = folder + "/extra.test";
      aurostd::string2file(base, base);
      std::this_thread::sleep_for(std::chrono::seconds(test_time));
      eclipsed = aurostd::SecondsSinceDirectoryContentModified(folder);
      if (eclipsed >= test_time and eclipsed <= test_time + 1) {
        checkEqual(true, true, check_function, check_description, passed_checks, results);
      } else {
        checkEqual(eclipsed, test_time, check_function, check_description, passed_checks, results);
      }
      aurostd::RemoveDirectory(folder);
    }
  }

  void UnitTest::xparserTest(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    (void) errors; // Suppress compiler warnings

    const string task_description = "Test xparsers";
    string check_description;
    string check_function;

    // ---------------------------------------------------------------------------
    // Check | aurostd::JSON
    // ---------------------------------------------------------------------------

    // read strings
    // tests based on https://seriot.ch/projects/parsing_json.html
    {
      check_function = "JSON";
      check_description = "test string parsing";
      const std::string string_json = string(utd["xparser"]["string_json"]);

      const std::map<string, string> string_results({
          {             "1_2_3_bytes_UTF-8_sequences",       "\u0060\u012a\u12AB"},
          {                         "allowed_escapes",          "\"\\/\b\f\n\r\t"},
          {                                      "pi",                        "π"},
          {                    "escaped_noncharacter",                   "\uFFFF"},
          {                 "last_surrogates_1_and_2",                     "􏿿"},
          {                 "accepted_surrogate_pair",                        "𐐷"},
          {                                 "unicode",                   "\uA66D"},
          {                 "unicode_U+1FFFE_nonchar",                     "🿾"},
          {            "backslash_and_u_escaped_zero",                  "\\u0000"},
          {                        "three-byte-utf-8",                   "\u0821"},
          {                  "backslash_doublequotes",                       "\""},
          {                        "uescaped_newline",            "new\u000Aline"},
          {            "unicode_escaped_double_quote",                   "\u0022"},
          {                         "double_escape_a",                      "\\a"},
          {                                "comments",            "a/*b*/c/*d//e"},
          {                                   "space",                        " "},
          {                                 "uEscape", "\u0061\u30af\u30EA\u30b9"},
          {         "unicode_U+200B_ZERO_WIDTH_SPACE",                   "\u200B"},
          {                         "u+2028_line_sep",                   "\u2028"},
          {                          "two-byte-utf-8",                   "\u0123"},
          {                 "unicodeEscapedBackslash",                   "\u005C"},
          {           "unicode_U+2064_invisible_plus",                   "\u2064"},
          {               "escaped_control_character",                   "\u0012"},
          {                            "simple_ascii",                     "asd "},
          {                "unicode_U+10FFFE_nonchar",                     "􏿾"},
          {                                    "utf8",                       "€𝄞"},
          {                   "unescaped_char_delete",                        ""},
          {"surrogates_U+1D11E_MUSICAL_SYMBOL_G_CLEF",                        "𝄞"},
          {                         "double_escape_n",                      "\\n"},
          {                      "with_del_character",                      "aa"},
          {              "nonCharacterInUTF-8_U+FFFF",                      "￿"},
          {                "accepted_surrogate_pairs",                     "😹💍"},
          {             "in_array_with_leading_space",                      "asd"},
          {            "nonCharacterInUTF-8_U+10FFFF",                     "􏿿"},
          {                          "one-byte-utf-8",                   "\u002c"},
          {                               "unicode_2",                     "⍂㈴⍂"},
          {                  "unicode_U+FDD0_nonchar",                   "\uFDD0"},
          {                           "nbsp_uescaped",            "new\u00A0line"},
          {                  "unicode_U+FFFE_nonchar",                   "\uFFFE"},
          {        "reservedCharacterInUTF-8_U+1BFFF",                     "𛿿"},
          {                          "u+2029_par_sep",                   "\u2029"}
      });

      const std::map<string, string> escaped_results({
          {             "1_2_3_bytes_UTF-8_sequences",             "`\\u012a\\u12ab"},
          {                 "accepted_surrogate_pair",              "\\ud801\\udc37"},
          {                "accepted_surrogate_pairs", R"(\ud83d\ude39\ud83d\udc8d)"},
          {                         "allowed_escapes",         R"(\"\\\/\b\f\n\r\t)"},
          {            "backslash_and_u_escaped_zero",                   "\\\\u0000"},
          {                  "backslash_doublequotes",                        "\\\""},
          {                                "comments",       R"(a\/*b*\/c\/*d\/\/e)"},
          {                         "double_escape_a",                       "\\\\a"},
          {                         "double_escape_n",                       "\\\\n"},
          {               "escaped_control_character",                           ""},
          {                    "escaped_noncharacter",                     "\\uffff"},
          {             "in_array_with_leading_space",                         "asd"},
          {                 "last_surrogates_1_and_2",              "\\udbff\\udfff"},
          {                           "nbsp_uescaped",              "new\\u00a0line"},
          {            "nonCharacterInUTF-8_U+10FFFF",              "\\udbff\\udfff"},
          {              "nonCharacterInUTF-8_U+FFFF",                     "\\uffff"},
          {                          "one-byte-utf-8",                           ","},
          {                                      "pi",                     "\\u03c0"},
          {        "reservedCharacterInUTF-8_U+1BFFF",              "\\ud82f\\udfff"},
          {                            "simple_ascii",                        "asd "},
          {                                   "space",                           " "},
          {"surrogates_U+1D11E_MUSICAL_SYMBOL_G_CLEF",              "\\ud834\\udd1e"},
          {                        "three-byte-utf-8",                     "\\u0821"},
          {                          "two-byte-utf-8",                     "\\u0123"},
          {                         "u+2028_line_sep",                     "\\u2028"},
          {                          "u+2029_par_sep",                     "\\u2029"},
          {                                 "uEscape",      R"(a\u30af\u30ea\u30b9)"},
          {                        "uescaped_newline",                  "new\\nline"},
          {                   "unescaped_char_delete",                           ""},
          {                                 "unicode",                     "\\ua66d"},
          {                 "unicodeEscapedBackslash",                        "\\\\"},
          {                               "unicode_2",       R"(\u2342\u3234\u2342)"},
          {                "unicode_U+10FFFE_nonchar",              "\\udbff\\udffe"},
          {                 "unicode_U+1FFFE_nonchar",              "\\ud83f\\udffe"},
          {         "unicode_U+200B_ZERO_WIDTH_SPACE",                     "\\u200b"},
          {           "unicode_U+2064_invisible_plus",                     "\\u2064"},
          {                  "unicode_U+FDD0_nonchar",                     "\\ufdd0"},
          {                  "unicode_U+FFFE_nonchar",                     "\\ufffe"},
          {            "unicode_escaped_double_quote",                        "\\\""},
          {                                    "utf8",       R"(\u20ac\ud834\udd1e)"},
          {                      "with_del_character",                         "aa"}
      });

      aurostd::JSON::object jo = aurostd::JSON::loadString(string_json);
      bool overall_test = true;
      for (const auto& [key, value] : string_results) {
        if (static_cast<string>(jo[key]) != value) {
          check_description += "(at " + key + " subtest)";
          check(false, static_cast<string>(jo[key]), value, check_function, check_description, passed_checks, results);
          overall_test = false;
          break;
        }
      }
      if (overall_test) {
        check(overall_test, 0, 0, check_function, check_description, passed_checks, results);
      }

      check_description = "test string conversion to JSON";
      overall_test = true;
      for (const auto& [key, value] : escaped_results) {
        if (jo[key].toString(true, true) != "\"" + value + "\"") {
          check_description += "(at " + key + " subtest)";
          check(false, static_cast<string>(jo[key]), "\"" + value + "\"", check_function, check_description, passed_checks, results);
          overall_test = false;
          break;
        }
      }
      if (overall_test) {
        check(overall_test, 0, 0, check_function, check_description, passed_checks, results);
      }

      // save to file and read again
      check_function = "JSON";
      check_description = "file save + load (";

      const aurostd::compression_type ct = aurostd::compression_type::XZ;
      const std::string normal_tmp_path = aurostd::TmpStrCreate() + ".json";
      const std::string compressed_tmp_path = normal_tmp_path + aurostd::compression_suffix[ct];
      aurostd::JSON::saveFile(jo, normal_tmp_path);
      jo.saveFile(compressed_tmp_path, ct);

      const std::map<std::string, aurostd::JSON::object> jo_list = {
          {    "direct",     aurostd::JSON::loadFile(normal_tmp_path)},
          {"compressed", aurostd::JSON::loadFile(compressed_tmp_path)}
      };
      std::map<std::string, bool> test_result = {
          {    "direct", true},
          {"compressed", true}
      };
      fs::remove(compressed_tmp_path);
      fs::remove(normal_tmp_path);

      for (const auto& [name, test_jo] : jo_list) {
        for (const auto& [key, value] : escaped_results) {
          if (test_jo[key].toString(true, true) != "\"" + value + "\"") {
            check_description += "(at " + key + " subtest)";
            check(false, static_cast<string>(test_jo[key]), "\"" + value + "\"", check_function, check_description + name + ")", passed_checks, results);
            test_result[name] = false;
            break;
          }
        }
      }

      for (const auto& [name, res] : test_result) {
        if (res) {
          check(res, 0, 0, check_function, check_description + name + ")", passed_checks, results);
        }
      }
    }

    // read single chars

    {
      check_function = "JSON";
      check_description = "test char parsing";

      std::vector<char> test_chars{'!', '#', '$', '%', '&', '(', ')', '*', '+', '-', '.', '/', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', '>', '?', '@', 'A',
                                   'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '[', ']', '^', '_', '`',
                                   'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '{', '|', '}', '~'};

      std::vector<aurostd::JSON::object> test_json(test_chars.size());
      bool overall_test = true;

      for (int i = 0; i < test_chars.size(); i++) {
        test_json[i] = aurostd::JSON::object(test_chars[i]);
      }

      for (int i = 0; i < test_chars.size(); i++) {
        if (static_cast<char>(test_chars[i]) != static_cast<char>(test_json[i])) {
          checkEqual(static_cast<char>(test_chars[i]), static_cast<char>(test_json[i]), check_function, check_description, passed_checks, results);
          overall_test = false;
          break;
        }
      }

      if (overall_test) {
        check(overall_test, 0, 0, check_function, check_description, passed_checks, results);
      }
    }

    // read single numbers
    // tests based on https://seriot.ch/projects/parsing_json.html
    {
      check_function = "JSON";
      check_description = "test number parsing (double, long long)";
      const std::string number_json = string(utd["xparser"]["number_json"]);

      aurostd::JSON::object jo = aurostd::JSON::loadString(number_json);

      const std::map<string, double> double_results({
          {                     "0e+1",          0},
          {                      "0e1",          0},
          {              "after_space",          4},
          {     "double_close_to_zero",     -1e-78},
          {"double__too_close_to_zero",        NAN},
          {             "int_with_exp",        200},
          {               "minus_zero",         -0},
          {             "negative_int",       -123},
          {             "negative_one",         -1},
          {           "real_capital_e",      1e+22},
          {   "real_capital_e_neg_exp",       0.01},
          {   "real_capital_e_pos_exp",        100},
          {            "real_exponent",   1.23e+47},
          {   "real_fraction_exponent", 123.456e78},
          {             "real_neg_exp",       0.01},
          {        "real_pos_exponent",        100},
          {               "simple_int",        123},
          {              "simple_real", 123.456789},
          {                   "number",   1.23e+67},
          {                     "true",        1.0},
          {                    "false",        0.0},
          {                     "null",        NAN}
      });
      const std::map<string, long long> ll_results({
          {                  "0e+1",                                     0},
          {                   "0e1",                                     0},
          {           "after_space",                                     4},
          {  "double_close_to_zero",                                     0},
          {          "int_with_exp",                                   200},
          {            "minus_zero",                                     0},
          {          "negative_int",                                  -123},
          {          "negative_one",                                    -1},
          {        "real_capital_e", std::numeric_limits<long long>::max()},
          {"real_capital_e_neg_exp",                                     0},
          {"real_capital_e_pos_exp",                                   100},
          {         "real_exponent", std::numeric_limits<long long>::max()},
          {"real_fraction_exponent", std::numeric_limits<long long>::max()},
          {          "real_neg_exp",                                     0},
          {     "real_pos_exponent",                                   100},
          {            "simple_int",                                   123},
          {           "simple_real",                                   123},
          {                "number", std::numeric_limits<long long>::max()},
          {                  "true",                                     1},
          {                 "false",                                     0}
      });

      bool overall_test = true;
      for (const auto& entry : double_results) {
        if ((double) jo[entry.first] != entry.second) {
          if (std::isnan((double) jo[entry.first]) && std::isnan(entry.second)) {
            continue;
          }
          check_description += "(at " + entry.first + " subtest)";
          check(false, (double) jo[entry.first], entry.second, check_function, check_description, passed_checks, results);
          overall_test = false;
          break;
        }
      }
      for (const auto& [key, value] : ll_results) {
        if (static_cast<long long>(jo[key]) != value) {
          check_description += "(at " + key + " subtest)";
          check(false, static_cast<long long>(jo[key]), value, check_function, check_description, passed_checks, results);
          overall_test = false;
          break;
        }
      }
      if (overall_test) {
        check(overall_test, 0, 0, check_function, check_description, passed_checks, results);
      }
    }

    // read xvectors
    {
      check_function = "JSON";
      check_description = "test parsing xvector<";
      const std::string number_json = string(utd["xparser"]["vector_number_json"]);

      aurostd::JSON::object jo = aurostd::JSON::loadString(number_json);

      xvector<double> xvd_res = jo["xvector_double"];
      xvector<double> xvd_exp = {923.49445786, -441.74004105, 465.49355057, 96.15610686, 557.6834903, 147.6777196, 871.81485459, 287.89958188, 863.66132302, 876.36635155};
      checkEqual(xvd_res, xvd_exp, check_function, check_description + "double>", passed_checks, results);

      xvd_res = jo["xvector_nan"];
      xvd_exp = {449644208, -441.74004105, NAN, 725767871, 1.0, 946128279, 65015635, 287.89958188, 863.66132302, 762013152};
      checkEqual(xvd_res, xvd_exp, check_function, check_description + "double> mixed + NAN", passed_checks, results);

      const xvector<float> xvf_res = jo["xvector_nan"];
      const xvector<float> xvf_exp = {449644224, -441.74004105, NAN, 725767872, 1.0, 946128256, 65015636, 287.89958188, 863.66132302, 762013184}; // values different from double -> narrowing
      checkEqual(xvf_res, xvf_exp, check_function, check_description + "float> mixed + NAN", passed_checks, results);

      xvector<long long> xvll_res = jo["xvector_ll"];
      xvector<long long> xvll_exp = {449644208, -252515403, 601496576, 725767871, 502088591, 946128279, 65015635, 352203056, 717938486, 762013152};
      checkEqual(xvll_res, xvll_exp, check_function, check_description + "long long>", passed_checks, results);

      const xvector<unsigned long long> xvull_res = jo["xvector_ull"];
      const xvector<unsigned long long> xvull_exp = {449644208, 252515403, 601496576, 725767871, 502088591, 946128279, 65015635, 352203056, 717938486, 762013152};
      checkEqual(xvull_res, xvull_exp, check_function, check_description + "unsigned long long>", passed_checks, results);

      xvll_res = jo["xvector_mixed"];
      xvll_exp = {449644208, -441, 465, 725767871, 1, 946128279, 65015635, 287, 863, 762013152};
      checkEqual(xvll_res, xvll_exp, check_function, check_description + "long long> mixed", passed_checks, results);

      const xvector<int> xvi_res = jo["xvector_mixed"];
      const xvector<int> xvi_exp = {449644208, -441, 465, 725767871, 1, 946128279, 65015635, 287, 863, 762013152};
      checkEqual(xvi_res, xvi_exp, check_function, check_description + "int> mixed", passed_checks, results);

      // <uint> would fail as operator << xvector<uint> does not exist yet
    }

    // read xmatrix
    {
      check_function = "JSON";
      check_description = "test parsing xmatrix<";

      const std::string matrix_json = string(utd["xparser"]["matrix_json"]);
      aurostd::JSON::object jo = aurostd::JSON::loadString(matrix_json);

      const xmatrix<double> xmd_exp = {
          { 2.85899214e+07, 2.70836286e+04,              82,  7.47207733e+08,  4.67299748e+04, -6.25361399e+09, 8.00388892e+08},
          {-2.50841390e+04, 1.86048071e+07,  7.18691607e+04, -7.18492121e+00,  5.11971140e+02,  8.97364584e+06, 1.68446225e+02},
          { 4.94590604e+08, 1.89067109e+09,  4.64803095e+00,  7.74359909e+00,            5465,  1.37944663e+02, 1.37565108e+08},
          {           8499, 1.73734230e+03,  3.28957479e+08, -5.41490015e+03,  5.52617678e+00,  2.64034084e+09, 4.76561709e+01},
          { 3.31904798e+04,     2288175298, -2.96213531e+07,  1.22631199e+02, -8.83941811e+00, -2.17327450e+04, 3.77389069e+03},
          {-2.98263326e+04, 5.39915781e+05,  1.28198687e+09,  2.69529424e+07,  7.47855186e+01,      1538912783, 2.55027330e+02}
      };
      const xmatrix<double> xmd_res = jo["xmatrix_double_mix"];
      checkEqual(xmd_res, xmd_exp, check_function, check_description + "double> mixed", passed_checks, results);

      const xmatrix<int> xmll_exp = {
          { 28589921,       27083,         82, 747207733, 46729, -1958646694, 800388892},
          {   -25084,    18604807,      71869,        -7,   511,     8973645,       168},
          {494590604,  1890671090,          4,         7,  5465,         137, 137565108},
          {     8499,        1737,  328957479,     -5414,     5, -1654626456,        47},
          {    33190, -2006791998,  -29621353,       122,    -8,      -21732,      3773},
          {   -29826,      539915, 1281986870,  26952942,    74,  1538912783,       255}
      };
      const xmatrix<int> xmll_res = jo["xmatrix_double_mix"];
      checkEqual(xmll_res, xmll_exp, check_function, check_description + "int> mixed", passed_checks, results);
    }

    // read std::map
    {
      check_function = "JSON";
      check_description = "test parsing std:map<string,";
      const std::string map_json = string(utd["xparser"]["map_json"]);

      const std::map<string, double> md_exp = {
          {"VNZPC",                   2},
          {"HQPZP", -2472285305.2312846},
          {"DUKPG",               24171},
          {"NUJYO",          -108940230},
          {"MDPSG",     15011078.748e-5},
          {"UKYVK",    1072982.12214322},
          {"UCZSX",       825.8070164E2},
          {"SQBER",               -2056},
          {"XGPXD",   853427.8245249396},
          {"UASBY",          4474520688}
      };
      aurostd::JSON::object jo = aurostd::JSON::loadString(map_json);
      std::map<string, double> md_res = jo["map"];
      bool overall_test = true;
      for (const std::pair<string, double> entry : md_exp) {
        if (entry.second != md_res[entry.first]) {
          check_description += "double> (at " + entry.first + " subtest key)";
          check(false, md_res[entry.first], entry.second, check_function, check_description, passed_checks, results);
          overall_test = false;
          break;
        }
      }
      if (overall_test) {
        check(overall_test, 0, 0, check_function, check_description + "double>)", passed_checks, results);
      }

      const std::map<string, long long> mll_exp = {
          {"VNZPC",           2},
          {"HQPZP", -2472285305},
          {"DUKPG",       24171},
          {"NUJYO",  -108940230},
          {"MDPSG",         150},
          {"UKYVK",     1072982},
          {"UCZSX",       82580},
          {"SQBER",       -2056},
          {"XGPXD",      853427},
          {"UASBY",  4474520688}
      };
      std::map<string, long long> mll_res = jo["map"];
      overall_test = true;
      for (const std::pair<string, long long> entry : mll_exp) {
        if (entry.second != mll_res[entry.first]) {
          check_description += "long long> (at subtest key " + entry.first + ")";
          check(false, mll_res[entry.first], entry.second, check_function, check_description, passed_checks, results);
          overall_test = false;
          break;
        }
      }
      if (overall_test) {
        check(overall_test, 0, 0, check_function, check_description + "long long>)", passed_checks, results);
      }
    }

    // read std::vector
    {
      check_function = "JSON";
      check_description = "test parsing vector<";

      const std::string vector_json = string(utd["xparser"]["vector_mixed_json"]);

      aurostd::JSON::object jo = aurostd::JSON::loadString(vector_json);
      std::vector<std::string> vs_res = jo["string_vector"];
      std::vector<std::string> vs_exp = {"First w\"ord", "Second word", "Last word"};
      checkEqual(vs_res, vs_exp, check_function, check_description + "string>", passed_checks, results);

      vs_res = jo["mixed_vector"];
      vs_exp = {R"({"DUKPG":24171,"HQPZP":-2.47229e+09,"VNZPC":2})", "5.68e-06", "152.324", "784", "A string", "true", "false", "null", R"(["A","nested","list"])"};
      checkEqual(vs_res, vs_exp, check_function, check_description + "string> mixed", passed_checks, results);

      const std::vector<double> vd_res = jo["number_vector"];
      const std::vector<double> vd_exp = {449644208, 441.74004105, 465.49355057, 725767871, 1.0, 946128279, 65015635, 287.89958188, 863.66132302, 762013152};
      checkEqual(vd_res, vd_exp, check_function, check_description + "double> mixed", passed_checks, results);

      const std::vector<int> vll_res = jo["number_vector"];
      const std::vector<int> vll_exp = {449644208, 441, 465, 725767871, 1, 946128279, 65015635, 287, 863, 762013152};
      checkEqual(vll_res, vll_exp, check_function, check_description + "int> mixed", passed_checks, results);

      for (const string key : {"bool>", "bool> mixed"}) {// aurostd::joinWDelimiter does not work for bool
        std::vector<bool> vb_res = jo[key];
        std::vector<bool> vb_exp = {true, false, false, false, true};
        bool overall_test = true;
        for (size_t idx = 0; idx < vb_exp.size(); idx++) {
          if (vb_res[idx] != vb_exp[idx]) {
            check_description += key;
            check(false, vb_res[idx], vb_exp[idx], check_function, check_description, passed_checks, results);
            overall_test = false;
            break;
          }
        }
        if (overall_test) {
          check(overall_test, 0, 0, check_function, check_description + "bool>", passed_checks, results);
        }
      }
    }

    // create JSON storage objects
    {
      check_function = "JSON";
      check_description = "test converting ";
      aurostd::JSON::object so;
      {
        const int exp = 4652;
        so = exp;
        checkEqual((int) so, exp, check_function, check_description + "int", passed_checks, results);
      }
      {
        so = "Hello World";
        checkEqual((std::string) so, "Hello World", check_function, check_description + "std::string", passed_checks, results);
      }
      {
        so = nullptr;
        checkEqual((std::string) so, "null", check_function, check_description + "nullptr", passed_checks, results);
      }
      {
        const xcomplex<double> exp(51.25, 0.3485);
        so = exp;
        checkEqual((xcomplex<double>) so, exp, check_function, check_description + "xcomplex<double>", passed_checks, results);
      }
      {
        const xcomplex<float> exp(51.25, 0.3485);
        so = exp;
        checkEqual((xcomplex<float>) so, exp, check_function, check_description + "xcomplex<float>", passed_checks, results);
      }
      {
        const vector<float> exp = {12.33, 20.7, 34.23454};
        so = exp;
        checkEqual((std::string) so, "[12.33,20.7,34.2345]", check_function, check_description + "vector<float>", passed_checks, results);
      }
      {
        const deque<float> exp = {12.33, 20.7, 34.23454};
        so = exp;
        checkEqual((std::string) so, "[12.33,20.7,34.2345]", check_function, check_description + "deque<float>", passed_checks, results);
      }
      {
        const xvector<uint> exp = {3, 7, 10};
        so = exp;
        checkEqual((std::string) so, "[3,7,10]", check_function, check_description + "xvector<uint>", passed_checks, results);
      }
      {
        const xmatrix<float> exp = {
            {1.0, 2.0, 3.0, 4.0},
            {  5,   6,   7,   8},
            {  9,  10,  11,  12},
            { 13,  14,  15,  16}
        };
        so = exp;
        checkEqual((std::string) so, "[[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]]", check_function, check_description + "xmatrix<float>", passed_checks, results);
      }
      {
        xmatrix<xcomplex<float>> exp(3, 3, 1, 1);
        exp(1, 1) = {1.1, 1.11};
        exp(1, 2) = {1.2, 1.22};
        exp(1, 3) = {1.3, 1.33};
        exp(2, 1) = {2.1, 2.11};
        exp(2, 2) = {2.2, 2.22};
        exp(2, 3) = {2.3, 2.33};
        exp(3, 1) = {3.1, 3.11};
        exp(3, 2) = {3.2, 3.22};
        exp(3, 3) = {3.3, 3.33};
        so = exp;
        checkEqual((string) so, R"([[{"imag":1.11,"real":1.1},{"imag":1.22,"real":1.2},{"imag":1.33,"real":1.3}],[{"imag":2.11,"real":2.1},{"imag":2.22,"real":2.2},{"imag":2.33,"real":2.3}],[{"imag":3.11,"real":3.1},{"imag":3.22,"real":3.2},{"imag":3.33,"real":3.3}]])",
                   check_function, check_description + "xmatrix<xcomplex<float>>", passed_checks, results);
      }
      {
        const std::vector<std::vector<std::string>> exp = {
            { "Hello",  "World"},
            {"Hello2", "World2"}
        };
        so = exp;
        checkEqual(static_cast<std::string>(so), R"([["Hello","World"],["Hello2","World2"]])", check_function, check_description + "vector<vector<string>>", passed_checks, results);
      }
      {
        const std::deque<std::deque<std::string>> exp = {
            { "Hello",  "World"},
            {"Hello2", "World2"}
        };
        so = exp;
        checkEqual(static_cast<std::string>(so), R"([["Hello","World"],["Hello2","World2"]])", check_function, check_description + "deque<deque<string>>", passed_checks, results);
      }
      {
        const vector<float> a = {12.33, 20.7, 34.23454};
        const vector<float> b = {89.6, 2.243, 76.432};
        std::map<std::string, std::vector<float>> exp;
        exp.insert({"Hello1", a});
        exp.insert({"Hello2", b});
        so = exp;
        checkEqual(static_cast<std::string>(so), R"({"Hello1":[12.33,20.7,34.2345],"Hello2":[89.6,2.243,76.432]})", check_function, check_description + "map<string, vector<float>>", passed_checks, results);
      }
      {
        const deque<float> a = {12.33, 20.7, 34.23454};
        const deque<float> b = {89.6, 2.243, 76.432};
        std::map<std::string, std::deque<float>> exp;
        exp.insert({"Hello1", a});
        exp.insert({"Hello2", b});
        so = exp;
        checkEqual(static_cast<std::string>(so), R"({"Hello1":[12.33,20.7,34.2345],"Hello2":[89.6,2.243,76.432]})", check_function, check_description + "map<string, deque<float>>", passed_checks, results);
      }
    }

    // check type specific functions
    {
      check_function = "JSON";
      check_description = "type specific function ";
      aurostd::JSON::object so;
      size_t expected = 11;
      so = "Hello World";
      checkEqual(so.size(), expected, check_function, check_description + "size for string", passed_checks, results);
      checkEqual(so.empty(), false, check_function, check_description + "empty for string", passed_checks, results);
      so = vector<int>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11});
      checkEqual(so.size(), expected, check_function, check_description + "size for list", passed_checks, results);
      checkEqual(so.empty(), false, check_function, check_description + "empty for list", passed_checks, results);
      so = std::map<string, int>({
          {"eins", 1},
          {"zwei", 2},
          {"drei", 3}
      });
      expected = 3;
      checkEqual(so.size(), expected, check_function, check_description + "size for dictionary", passed_checks, results);
      checkEqual(so.empty(), false, check_function, check_description + "empty for dictionary", passed_checks, results);
    }

    // check join for LIST
    {
      check_function = "JSON";
      check_description = "joining two lists";
      const aurostd::JSON::object list_one(aurostd::JSON::object_types::LIST);
      const aurostd::JSON::object list_two(aurostd::JSON::object_types::LIST);
      const aurostd::JSON::object list_joined(aurostd::JSON::object_types::LIST);
      list_one.push_back(1);
      list_one.push_back(1.235);
      list_one.push_back("Hello One");
      list_two.push_back(2);
      list_two.push_back(2.235);
      list_two.push_back("Hello Two");
      list_joined.push_back(1);
      list_joined.push_back(1.235);
      list_joined.push_back("Hello One");
      list_joined.push_back(2);
      list_joined.push_back(2.235);
      list_joined.push_back("Hello Two");
      list_one.join(list_two);
      checkEqual(static_cast<std::string>(list_one), static_cast<std::string>(list_joined), check_function, check_description, passed_checks, results);
    }

    // check join for DICTIONARY
    {
      check_function = "JSON";
      check_description = "joining two dictionaries";
      aurostd::JSON::object dict_one(aurostd::JSON::object_types::DICTIONARY);
      aurostd::JSON::object dict_two(aurostd::JSON::object_types::DICTIONARY);
      aurostd::JSON::object dict_joined(aurostd::JSON::object_types::DICTIONARY);
      dict_one["one"] = 1;
      dict_one["both"] = 1;
      dict_two["two"] = 2;
      dict_two["both"] = 2;
      dict_joined["one"] = 1;
      dict_joined["two"] = 2;
      dict_joined["both"] = 2;
      dict_one.join(dict_two);
      checkEqual(static_cast<std::string>(dict_one), static_cast<std::string>(dict_joined), check_function, check_description, passed_checks, results);
    }
  }

  // -----------------------------------------------------------------------------
  // Check | JSON Serialization / Deserialization
  // -----------------------------------------------------------------------------

  /// @authors
  /// @mod{ST,20250108,created}
  void UnitTest::jsonSerialization(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    this->jsonMockSerialization(passed_checks, results, errors);
    this->xKpointJsonSerialization(passed_checks, results, errors);
    this->xPOTCARJsonSerialization(passed_checks, results, errors);
    this->xVASPRUNXMLJsonSerialization(passed_checks, results, errors);
    this->xDOSCARJsonSerialization(passed_checks, results, errors);
    this->xIBZKPTJsonSerialization(passed_checks, results, errors);
    this->xQMVASPJsonSerialization(passed_checks, results, errors);
    this->xEIGENVALJsonSerialization(passed_checks, results, errors);
    this->xOUTCARJsonSerialization(passed_checks, results, errors);
  }

  /// @brief set the datetime key to empty string
  void strip_datetime(aurostd::JSON::object& jo) {
    jo["_aflow_serialization"]["datetime"] = "";
  }

  /// @brief Test the serializable class T by running a (de)serialization cycle and checking json content.
  /// @tparam T The serializable class
  /// @tparam F An initializer function that returns a T and takes no args
  /// @authors
  /// @mod{ST,20250316,created}
  template <class T, aurostd::JSON::enable_serializable<T>, class F> void UnitTest::json_cycle_test(F initializer, uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    const string check_function = string("JSON serialization cycle for class ").append(typeid(T).name());
    const string check_description = string("Check JSON string equality for cycled ").append(typeid(T).name());
    try {
      const T t = initializer();
      aurostd::JSON::object json1 = t.dumpToJson();
      const T t_cycled = T::loadFromJson(json1);
      aurostd::JSON::object json2 = t_cycled.dumpToJson();
      strip_datetime(json1);
      strip_datetime(json2);
      checkEqual(json2.toString(), json1.toString(), check_function, check_description, passed_checks, results);
    } catch (const std::exception& e) {
      errors.emplace_back(e.what());
    } catch (aurostd::xerror& e) {
      errors.emplace_back(e.what());
    }
  }

  /// @authors
  /// @mod{NA,20250304,created}
  void UnitTest::xPOTCARJsonSerialization(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    json_cycle_test<xPOTCAR>(
        []() {
          xPOTCAR ret{};
          ret.GetProperties(aurostd::EmbData::get_test_file("POTCAR.txt"));
          return ret;
        },
        passed_checks, results, errors);
  }

  /// @authors
  /// @mod{NA,20250304,created}
  void UnitTest::xVASPRUNXMLJsonSerialization(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    json_cycle_test<xVASPRUNXML>(
        []() {
          xVASPRUNXML ret{};
          ret.GetProperties(aurostd::EmbData::get_test_file("VASPRUNXML.txt"));
          return ret;
        },
        passed_checks, results, errors);
  }

  /// @authors
  /// @mod{NA,20250304,created}
  void UnitTest::xDOSCARJsonSerialization(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    json_cycle_test<xDOSCAR>(
        []() {
          xDOSCAR ret{};
          ret.GetProperties(aurostd::EmbData::get_test_file("DOSCAR.txt"));
          return ret;
        },
        passed_checks, results, errors);
  }

  /// @authors
  /// @mod{NA,20250304,created}
  void UnitTest::xIBZKPTJsonSerialization(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    json_cycle_test<xIBZKPT>(
        []() {
          xIBZKPT ret{};
          ret.GetProperties(aurostd::EmbData::get_test_file("IBZKPT.txt"));
          return ret;
        },
        passed_checks, results, errors);
  }

  /// @authors
  /// @mod{NA,20250304,created}
  void UnitTest::xQMVASPJsonSerialization(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    json_cycle_test<xQMVASP>(
        []() {
          xQMVASP ret{};
          ret.GetProperties(aurostd::EmbData::get_test_file("QMVASP.txt"));
          return ret;
        },
        passed_checks, results, errors);
  }

  /// @authors
  /// @mod{NA,20250304,created}
  void UnitTest::xEIGENVALJsonSerialization(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    json_cycle_test<xEIGENVAL>(
        []() {
          xEIGENVAL ret{};
          ret.GetProperties(aurostd::EmbData::get_test_file("EIGENVAL.txt"));
          return ret;
        },
        passed_checks, results, errors);
  }

  /// @authors
  /// @mod{NA,20250320,created}
  void UnitTest::xOUTCARJsonSerialization(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    json_cycle_test<xOUTCAR>(
        []() {
          xOUTCAR ret{};
          ret.GetProperties(aurostd::EmbData::get_test_file("OUTCAR.txt"));
          return ret;
        },
        passed_checks, results, errors);
  }

  /// @authors
  /// @mod{ST,20250316,created}
  void UnitTest::jsonMockSerialization(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    json_cycle_test<MockSerializableMain>([]() { return MockSerializableMain{}; }, passed_checks, results, errors);
  }

  /// @authors
  /// @mod{ST,20250108,created}
  /// @mod{ST,20250316,use the json cycle test method}
  void UnitTest::xKpointJsonSerialization(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    json_cycle_test<xKPOINTS>(
        []() {
          xKPOINTS ret{};
          ret.GetProperties(string(aurostd::EmbData::get_test_file("KPOINTS.txt")));
          return ret;
        },
        passed_checks, results, errors);
  }

  void UnitTest::xfitTest(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    (void) errors;  // Suppress compiler warnings
    // setup test environment
    string check_function;
    string check_description;

    bool calculated_bool;
    bool expected_bool;
    double calculated_dbl;
    double expected_dbl;
    xvector<double> calculated_xvecdbl;
    xvector<double> expected_xvecdbl;
    xvector<double> calculated_xvecdbl_r;
    xvector<double> expected_xvecdbl_r;
    xvector<double> calculated_xvecdbl_i;
    xvector<double> expected_xvecdbl_i;
    xmatrix<double> calculated_xmatdbl;
    xmatrix<double> expected_xmatdbl;

    // ---------------------------------------------------------------------------
    // Check | companion matrix //SD20220318
    // ---------------------------------------------------------------------------
    check_function = "aurostd::companion_matrix()";
    check_description = "calculate the companion matrix of a univariate polynomial";
    const xvector<double> pc = {6.0, -5.0, -2.0, 3.0};
    expected_xmatdbl = {
        {       0.0,       1.0,       0.0},
        {       0.0,       0.0,       1.0},
        {-6.0 / 3.0, 5.0 / 3.0, 2.0 / 3.0}
    };
    calculated_xmatdbl = aurostd::companion_matrix(pc);
    checkEqual(calculated_xmatdbl, expected_xmatdbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | polynomialFindRoots //SD20220318
    // ---------------------------------------------------------------------------
    check_function = "aurostd::polynomialFindRoots()";
    expected_xvecdbl_r = {1.05576592536838, 1.05576592536838, -1.44486518407010};
    expected_xvecdbl_i = {-0.519201791550296, 0.519201791550296, 0.0};
    calculated_xvecdbl_r = xvector<double>(3), calculated_xvecdbl_i = xvector<double>(3);
    aurostd::polynomialFindRoots(pc, calculated_xvecdbl_r, calculated_xvecdbl_i);
    check_description = "calculate the roots of a univariate polynomial (real part)";
    checkEqual(calculated_xvecdbl_r, expected_xvecdbl_r, check_function, check_description, passed_checks, results);
    check_description = "calculate the roots of a univariate polynomial (imag part)";
    checkEqual(calculated_xvecdbl_i, expected_xvecdbl_i, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | polynomialCurveFit //SD20220422
    // ---------------------------------------------------------------------------
    check_function = "aurostd::polynomialCurveFit()";
    check_description = "calculate the coefficients for a quintic polynomial that fits the data";
    const xvector<double> xdata = {
        0.551878738140095, 0.815194385592879, -1.436650500869055, 0.837122604741089, 0.024588378708606, -1.521792147897101, -0.998568921765830, -0.222562772922286, 0.964723358008907, 0.986066878262693};
    const xvector<double> ydata = {
        -1.895746346279865, -2.426516802630079, -0.329326356733414, -2.476265605752963, -1.148143476025739, -0.304622388440250, -0.471352985639471, -0.911928012994989, -2.783990433247946, -2.838593165314427};
    const xvector<double> wdata = {
        0.757740130578333, 0.743132468124916, 0.392227019534168, 0.655477890177557, 0.171186687811562, 0.706046088019609, 0.031832846377421, 0.276922984960890, 0.046171390631154, 0.097131781235848};
    expected_xvecdbl = {-1.121837858162905, -1.056905122203599, -0.543168130435282, -0.145858244844311, -0.007682867520147, 0.000748981621075};
    calculated_xvecdbl = aurostd::polynomialCurveFit(xdata, ydata, 5, wdata);
    checkEqual(expected_xvecdbl, calculated_xvecdbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | findZeroBrent //SD20220517
    // ---------------------------------------------------------------------------
    check_function = "aurostd::findZeroBrent()";
    check_description = "find the zero of a univariate function";
    const std::function<double(double)> func = [](double x) { return std::pow(x, 2.0) - 3.0; };
    expected_dbl = 1.732050807568877;
    aurostd::findZeroBrent(0.0, 10.0, func, calculated_dbl);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | checkDerivatives //SD20220622
    // ---------------------------------------------------------------------------
    check_function = "aurostd::checkDerivatives()";
    check_description = "check that the analytical derivatives of a 2D function are correct";
    xvector<double> x0 = {aurostd::ran2(), aurostd::ran2()};
    const std::function<double(xvector<double>)> testf = [](xvector<double> x) { return std::sin(x(1)) + std::cos(x(2)); };
    vector<std::function<double(xvector<double>)>> testdf;
    testdf.emplace_back([](xvector<double> x) { return std::cos(x(1)); });
    testdf.emplace_back([](xvector<double> x) { return -1.0 * std::sin(x(2)); });
    expected_bool = true;
    calculated_bool = checkDerivatives(x0, testf, testdf);
    checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | calcNumericalJacobian //SD20220624
    // ---------------------------------------------------------------------------
    check_function = "aurostd::calcNumericalJacobian()";
    check_description = "calculate the numerical Jacobian of a 3D rectangular system";
    const xvector<double> dx = 10.0 * _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_ * aurostd::ones_xv<double>(3);
    x0 = {aurostd::ran2(), aurostd::ran2(), aurostd::ran2()};
    vector<std::function<double(xvector<double>)>> vfunc;
    vector<std::function<double(xvector<double>)>> vdfunc;
    vector<vector<std::function<double(xvector<double>)>>> jac;
    vfunc.emplace_back([](xvector<double> x) { return std::exp(1.0 * x(1) + 2.0 * x(2) + 3.0 * x(3)); });
    vfunc.emplace_back([](xvector<double> x) { return std::cos(1.0 * x(1) + 2.0 * x(2) + 3.0 * x(3)); });
    for (int i = 1; i <= 3; i++) {
      vdfunc.emplace_back([i](xvector<double> x) { return (double) i * std::exp(1.0 * x(1) + 2.0 * x(2) + 3.0 * x(3)); });
    }
    jac.push_back(vdfunc);
    vdfunc.clear();
    for (int i = 1; i <= 3; i++) {
      vdfunc.emplace_back([i](xvector<double> x) { return -1.0 * (double) i * std::sin(1.0 * x(1) + 2.0 * x(2) + 3.0 * x(3)); });
    }
    jac.push_back(vdfunc);
    expected_bool = true;
    calculated_bool = true;
    vector<vector<std::function<double(xvector<double>)>>> testjac = aurostd::calcNumericalJacobian(vfunc, dx);
    for (size_t i = 0; i < jac.size(); i++) {
      for (size_t j = 0; j < jac[0].size(); j++) {
        calculated_bool = calculated_bool && aurostd::isequal(jac[i][j](x0), testjac[i][j](x0));
      }
    }
    checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | findZeroNewtonRaphson //SD20220616
    // ---------------------------------------------------------------------------
    check_function = "aurostd::findZeroNewtonRaphson()";
    check_description = "find the zeros of a 2D nonlinear square system";
    x0 = {aurostd::ran2(), aurostd::ran2()};
    vfunc.clear();
    vdfunc.clear();
    jac.clear();
    vfunc.emplace_back([](xvector<double> x) { return std::exp(-std::exp(-(x(1) + x(2)))) - x(2) * (1.0 + std::pow(x(1), 2.0)); });
    vfunc.emplace_back([](xvector<double> x) { return x(1) * std::cos(x(2)) + x(2) * std::sin(x(1)) - 0.5; });
    vdfunc.emplace_back([](xvector<double> x) { return std::exp(-std::exp(-(x(1) - x(2))) - x(1) - x(2)) - 2 * x(1) * x(2); });
    vdfunc.emplace_back([](xvector<double> x) { return std::exp(-std::exp(-(x(1) - x(2))) - x(1) - x(2)) - (1.0 + std::pow(x(1), 2.0)); });
    jac.push_back(vdfunc);
    vdfunc.clear();
    vdfunc.emplace_back([](xvector<double> x) { return x(2) * std::cos(x(1)) + std::cos(x(2)); });
    vdfunc.emplace_back([](xvector<double> x) { return std::sin(x(1)) - x(1) * std::sin(x(2)); });
    jac.push_back(vdfunc);
    expected_xvecdbl = {0.353246561918931, 0.606082026502552};
    aurostd::findZeroNewtonRaphson(x0, vfunc, jac, calculated_xvecdbl);
    checkEqual(calculated_xvecdbl, expected_xvecdbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | findZeroDeflation //SD20220616
    // ---------------------------------------------------------------------------
    check_function = "aurostd::findZeroDeflation()";
    check_description = "find multiple zeros of a 2D nonlinear square system";
    vfunc.clear();
    vdfunc.clear();
    jac.clear();
    vfunc.emplace_back([](xvector<double> x) { return x(1) * x(2) - 1.0; });
    vfunc.emplace_back([](xvector<double> x) { return std::pow(x(1), 2.0) + std::pow(x(2), 2.0) - 4.0; });
    vdfunc.emplace_back([](xvector<double> x) { return x(2); });
    vdfunc.emplace_back([](xvector<double> x) { return x(1); });
    jac.push_back(vdfunc);
    vdfunc.clear();
    vdfunc.emplace_back([](xvector<double> x) { return 2.0 * x(1); });
    vdfunc.emplace_back([](xvector<double> x) { return 2.0 * x(2); });
    jac.push_back(vdfunc);
    expected_xmatdbl = {
        {0.517638090205041,  1.93185165257813,  -1.93185165257813},
        { 1.93185165257813, 0.517638090205041, -0.517638090205041}
    };
    aurostd::findZeroDeflation(aurostd::vector2xvector<double>({0.0, 1.0}), vfunc, jac, calculated_xmatdbl);
    checkEqual(calculated_xmatdbl, expected_xmatdbl, check_function, check_description, passed_checks, results);
  }

  void UnitTest::entryLoaderTest(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    (void) errors;  // Suppress compiler warnings

    // setup test environment
    const string task_description = "Testing EntryLoader";
    string check_description;
    stringstream check_description_helper;
    string check_function;

    const std::string test_alloy = "MnPdPt";
    bool recursive = false;
    aflowlib::EntryLoader el;
    aflowlib::_aflowlib_entry test_entry;
    xstructure test_structure;

    size_t expected_size_t = 0;
    const std::vector<std::string> test_AUIDs = {"aflow:2de63b1ebe0a1a83", "4d8cf7edb50d1901", "auid:6d47aa3f4f1286d0", "aflow:7dd846bc04c764e8", "9d84facf8161aa60", "broken"};
    const std::string test_AUID = "aflow:0d16c1946df2435c";

    const std::vector<std::string> test_AURLs = {"aflowlib.duke.edu:AFLOWDATA/LIB2_WEB/Ca_svCu_pv/84", "AFLOWDATA/LIB2_WEB/Ca_svCu_pv/546", "aurl:AFLOWDATA/LIB2_WEB/Ca_svCu_pv/724.BA",
                                                 "LIB2_WEB/Ca_svCu_pv/253", "aflowlib.duke.edu:AFLOWDATA/LIB2_WEB/Ca_svCu_pv/230"};
    const std::string test_AURL = "aflowlib.duke.edu:AFLOWDATA/LIB2_WEB/Ca_svCu_pv/539";

    std::map<std::string, aflowlib::EntryLoader::Source> test_sources = {
        {   "AUTO SELECT",           aflowlib::EntryLoader::Source::NONE},
        {        "SQLITE",         aflowlib::EntryLoader::Source::SQLITE},
        {         "AFLUX",          aflowlib::EntryLoader::Source::AFLUX},
        {    "FILESYSTEM",     aflowlib::EntryLoader::Source::FILESYSTEM},
        {       "RESTAPI",        aflowlib::EntryLoader::Source::RESTAPI},
        {"FILESYSTEM_RAW", aflowlib::EntryLoader::Source::FILESYSTEM_RAW},
        {   "RESTAPI_RAW",    aflowlib::EntryLoader::Source::RESTAPI_RAW}
    };

    std::map<std::string, aflowlib::EntryLoader::Source> short_test_sources = {
        {    "SQLITE",     aflowlib::EntryLoader::Source::SQLITE},
        {     "AFLUX",      aflowlib::EntryLoader::Source::AFLUX},
        {"FILESYSTEM", aflowlib::EntryLoader::Source::FILESYSTEM},
        {   "RESTAPI",    aflowlib::EntryLoader::Source::RESTAPI},
    };

    // ---------------------------------------------------------------------------
    // Check | load alloys
    for (std::map<std::string, aflowlib::EntryLoader::Source>::iterator source = test_sources.begin(); source != test_sources.end(); source++) {
      check_function = "EntryLoader::loadAlloy()";
      if (source->first == "RESTAPI" || source->first == "RESTAPI_RAW") {
        recursive = false;
      } else {
        recursive = true;
      }
      check_description = source->first + " - " + test_alloy;
      if (recursive) {
        check_description += " - recursive";
        expected_size_t = 3400;
      } else {
        expected_size_t = 90;
      }
      el.clear();
      el.m_out_silent = true;
      {
        long double start = aurostd::get_seconds();
        if (el.setSource(source->second)) { // don't test if basic requirements are not met for a source
          if (source->first == "FILESYSTEM" || source->first == "FILESYSTEM_RAW") {
            if (el.m_sqlite_alloy_db_ptr == nullptr) {
              expected_size_t = 1300;
            }
          }
          el.loadAlloy(test_alloy, recursive);
          const long double duration = aurostd::get_delta_seconds(start);
          aurostd::StringstreamClean(check_description_helper);
          check_description_helper << " | speed " << el.m_entries_flat->size() / duration << " entries/s; " << el.m_entries_flat->size() << " entries";
          check_description += check_description_helper.str();
          check((expected_size_t < el.m_entries_flat->size()), el.m_entries_flat->size(), expected_size_t, check_function, check_description, passed_checks, results);
        }
      }
    }

    // ---------------------------------------------------------------------------
    // Check | load AUID + Xstructure
    for (std::map<std::string, aflowlib::EntryLoader::Source>::iterator source = short_test_sources.begin(); source != short_test_sources.end(); source++) {
      check_function = "EntryLoader::loadAUID()";
      check_description = source->first + " + xstructure";
      expected_size_t = 6;
      el.clear();
      el.m_out_silent = true;
      el.m_xstructure_original = true;
      el.m_xstructure_relaxed = true;
      if (source->first == "RESTAPI" || source->first == "RESTAPI_RAW") {
        el.m_filesystem_path = "/fake/"; // force xstructure test to use REST API
      }
      long double start = aurostd::get_seconds();
      if (el.setSource(source->second)) { // don't test if basic requirements are not met for a source
        el.loadAUID(test_AUID);
        if (source->first == "AFLUX") {
          test_entry = *el.m_entries_flat->back();
        }
        el.loadAUID(test_AUIDs);

        const long double duration = aurostd::get_delta_seconds(start);
        aurostd::StringstreamClean(check_description_helper);
        check_description_helper << " | speed " << el.m_entries_flat->size() / duration << " entries/s; " << el.m_entries_flat->size() << " entries";
        check_description += check_description_helper.str();
        checkEqual(el.m_entries_flat->size(), expected_size_t, check_function, check_description, passed_checks, results);
      }
    }

    // ---------------------------------------------------------------------------
    // Check | load xstructure from file
    check_function = "EntryLoader::loadXstructureFile()";
    check_description = "load xstructure from CONTCAR.relax";
    if (!test_entry.auid.empty()) {
      el.loadXstructureFile(test_entry, test_structure, {"CONTCAR.relax"});
      checkEqual(test_structure.atoms.size(), (size_t) 6, check_function, check_description, passed_checks, results);
    } else {
      check_description += " | failed to load example structure form AFLUX in previous test";
      check(false, 0, 0, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | load AURL
    for (std::map<std::string, aflowlib::EntryLoader::Source>::iterator source = short_test_sources.begin(); source != short_test_sources.end(); source++) {
      check_function = "EntryLoader::loadAURL()";
      check_description = source->first;
      expected_size_t = 6;
      el.clear();
      el.m_out_silent = true;
      {
        long double start = aurostd::get_seconds();
        if (el.setSource(source->second)) { // don't test if basic requirements are not met for a source
          el.loadAURL(test_AURLs);
          el.loadAURL(test_AURL);
          const long double duration = aurostd::get_delta_seconds(start);
          aurostd::StringstreamClean(check_description_helper);
          check_description_helper << " | speed " << el.m_entries_flat->size() / duration << " entries/s; " << el.m_entries_flat->size() << " entries";
          check_description += check_description_helper.str();
          checkEqual(el.m_entries_flat->size(), expected_size_t, check_function, check_description, passed_checks, results);
        }
      }
    }
  }

  void UnitTest::ovaspTest(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    const string tmpdir = aurostd::TmpDirectoryCreate("unit_test_ovasp");

    const string qmvasp_tfile = tmpdir + "/QMVASP.txt";
    aurostd::string2file(aurostd::EmbData::get_test_file("QMVASP.txt"), qmvasp_tfile);
    const xQMVASP qmvasp(qmvasp_tfile);
    checkEqual(qmvasp.H_atom_relax, -3.742564, "xQMVASP()", "xQMVASP construction from file", passed_checks, results);
    checkEqual(qmvasp.H_atom_static, AUROSTD_NAN, "xQMVASP()", "xQMVASP construction from file", passed_checks, results);
    checkEqual(qmvasp.vforces.size(), 1UL, "xQMVASP()", "xQMVASP construction from file", passed_checks, results);

    const string doscar_tfile = tmpdir + "/DOSCAR.txt";
    aurostd::string2file(aurostd::EmbData::get_test_file("DOSCAR.txt"), doscar_tfile);
    const xDOSCAR doscar(doscar_tfile);
    checkEqual(doscar.Efermi, 8.00489507, "xDOSCAR()", "xDOSCAR construction from file", passed_checks, results);
    checkEqual(doscar.energy_max, 24.77735449, "xDOSCAR()", "xDOSCAR construction from file", passed_checks, results);
    checkEqual(doscar.energy_min, -5.13478138, "xDOSCAR()", "xDOSCAR construction from file", passed_checks, results);
    checkEqual(doscar.denergy, 0.1, "xDOSCAR()", "xDOSCAR construction from file", passed_checks, results);

    const string eig_tfile = tmpdir + "/EIGENVAL.txt";
    aurostd::string2file(aurostd::EmbData::get_test_file("EIGENVAL.txt"), eig_tfile);
    const xEIGENVAL eig(eig_tfile);
    checkEqual(eig.number_bands, 36U, "xEIGENVAL()", "xEIGENVAL construction from file", passed_checks, results);
    checkEqual(eig.number_kpoints, 256U, "xEIGENVAL()", "xEIGENVAL construction from file", passed_checks, results);
    checkEqual(eig.energy_max, 87.728, "xEIGENVAL()", "xEIGENVAL construction from file", passed_checks, results);
    checkEqual(eig.energy_min, -3.1278, "xEIGENVAL()", "xEIGENVAL construction from file", passed_checks, results);

    const string ibz_tfile = tmpdir + "/IBZKPT.txt";
    aurostd::string2file(aurostd::EmbData::get_test_file("IBZKPT.txt"), ibz_tfile);
    const xIBZKPT ibz(ibz_tfile);
    checkEqual(ibz.nweights, 8000U, "xIBZKPT()", "xIBZKPT construction from file", passed_checks, results);
    checkEqual(ibz.nkpoints_irreducible, 256U, "xIBZKPT()", "xIBZKPT construction from file", passed_checks, results);

    const string kpt_tfile = tmpdir + "/KPOINTS.txt";
    aurostd::string2file(aurostd::EmbData::get_test_file("KPOINTS.txt"), kpt_tfile);
    const xKPOINTS kpt(kpt_tfile);
    checkEqual(kpt.mode, 0, "xKPOINTS()", "xKPOINTS construction from file", passed_checks, results);
    checkEqual(kpt.nkpoints, 8000, "xKPOINTS()", "xKPOINTS construction from file", passed_checks, results);
    checkEqual(kpt.grid_type, "Gamma", "xKPOINTS()", "xKPOINTS construction from file", passed_checks, results);

    const string pot_tfile = tmpdir + "/POTCAR.txt";
    aurostd::string2file(aurostd::EmbData::get_test_file("POTCAR.txt"), pot_tfile);
    const xPOTCAR pot(pot_tfile);
    checkEqual(pot.ENMAX, 240.3, "xPOTCAR()", "xPOTCAR construction from file", passed_checks, results);
    checkEqual(pot.ENMIN, 180.225, "xPOTCAR()", "xPOTCAR construction from file", passed_checks, results);
    checkEqual(pot.RMAX_min, 2.974, "xPOTCAR()", "xPOTCAR construction from file", passed_checks, results);
    checkEqual(pot.RMAX_max, 2.974, "xPOTCAR()", "xPOTCAR construction from file", passed_checks, results);

    const string vrx_tfile = tmpdir + "/VASPRUNXML.txt";
    aurostd::string2file(aurostd::EmbData::get_test_file("VASPRUNXML.txt"), vrx_tfile);
    const xVASPRUNXML vrx(vrx_tfile);
    checkEqual(vrx.natoms, 1.0, "xVASPRUNXML()", "xVASPRUNXML construction from file", passed_checks, results);

    const string out_tfile = tmpdir + "/OUTCAR.txt";
    aurostd::string2file(aurostd::EmbData::get_test_file("OUTCAR.txt"), out_tfile);
    xOUTCAR out(out_tfile);
    checkEqual(out.NELM, 60, "xOUTCAR()", "xOUTCAR construction from file", passed_checks, results);
    checkEqual(out.NIONS, 1, "xOUTCAR()", "xOUTCAR construction from file", passed_checks, results);
    checkEqual(out.IALGO, 38, "xOUTCAR()", "xOUTCAR construction from file", passed_checks, results);
    checkEqual(out.ISMEAR, 1, "xOUTCAR()", "xOUTCAR construction from file", passed_checks, results);

    out.GetBandGap(out.Efermi);
    checkEqual(out.Egap.front(), -AUROSTD_NAN, "xOUTCAR()", "xOUTCAR construction from file", passed_checks, results);
    checkEqual(out.Egap_type.front(), "metal", "xOUTCAR()", "xOUTCAR construction from file", passed_checks, results);

    aurostd::RemoveDirectory(tmpdir);
  }

  /// @brief verify mathematical accuracy of XOUTCAR dielectric functions
  /// @note A dielectric calculation was performed on Si to produce frequency dependent dielectric matrix
  /// @see
  /// @doi{10.1103/PhysRevB.27.985}
  /// @authors
  /// @mod{NHA,20251116,created function}
  void UnitTest::dielectricTest(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    const string check_func = "XOUTCAR_dielectric";
    const vector<double> intra_real_verified = {
        8.9962900000000001, 8.9963619999999995, 8.9965759999999992, 8.9969330000000003, 8.9974329999999991, 8.9980759999999993, 8.9988620000000008, 8.9997910000000001, 9.000864,           9.0020790000000002,
        9.0034379999999992, 9.0049410000000005, 9.0065880000000007, 9.0083789999999997, 9.0103139999999993, 9.0123929999999994, 9.0146169999999994, 9.0169859999999993, 9.0195000000000007, 9.0221590000000003,
        9.0249649999999999, 9.0279159999999994, 9.0310140000000008, 9.0342590000000005, 9.0376499999999993, 9.0411900000000003, 9.0448769999999996, 9.0487129999999993, 9.0526970000000002, 9.0568310000000007};
    const vector<double> intra_im_verified = {
        0.000000000000000000, 0.001433000000000000, 0.0028670000000000002, 0.0043010000000000001, 0.0057349999999999996, 0.0071710000000000003, 0.0086060000000000008, 0.0100430000000,
        0.011481000000000000, 0.012921000000000000, 0.014361000000000002,  0.015803999999999999,  0.017247999999999999,  0.018693999999999999,  0.020143000000000001,  0.021593000000000001,
        0.023046000000000001, 0.024502000000000003, 0.025961000000000001,  0.027422000000000002,  0.028887000000000003,  0.030355000000000000,  0.031827000000000001,  0.033301999999999998,
        0.034780999999999999, 0.036263999999999998, 0.037752000000000001,  0.039244000000000001,  0.040739999999999998,  0.042241000000000001};
    const vector<double> EELS_verified = {
        0.00000000000000000000, 1.7705668732315631e-05, 3.542200521713607e-05,  5.3134961554959499e-05, 7.0842846514947559e-05, 8.8568671949350677e-05, 0.00010627369014402978, 0.00012399325865473552,
        0.00014171329986747774, 0.00015944451800580192, 0.00017716046876869128, 0.0001948964550911989,  0.00021262609215140624, 0.00023035999960744786, 0.00024810877836849019, 0.00026584607247050283,
        0.00028359473222434563, 0.00030135301302048597, 0.00031911916973442555, 0.00033687917238259675, 0.00035465569610725484, 0.00037243484627397109, 0.0003902269677087916,  0.00040801801740359528,
        0.00042581856825337188, 0.00044362654738172526, 0.00046145256872511561, 0.00047928238828890985, 0.00049711430441148874, 0.00051495848229323599};
    const vector<double> reflectivity_iso_verified = {
        0.24992269439112227, 0.24992419879231761, 0.24992867026239729, 0.2499361294539881,  0.24994657605311624, 0.24996000967406046, 0.24997642968671099, 0.24999583549726434,
        0.25001824710999604, 0.25004362198225621, 0.25007200063128066, 0.25010338203364119, 0.25013776478000782, 0.25017514750175229, 0.25021552872538766, 0.25025890662975137,
        0.25030530037352822, 0.25035470803116283, 0.25040712754970673, 0.25046255664640488, 0.25052103477922538, 0.25058251782831176, 0.25064704489142259, 0.25071461313334953,
        0.25078519907067992, 0.2508588620049399,  0.25093555745104196, 0.25101532352166839, 0.25109813608794512, 0.25118403323459781};

    xOUTCAR xout;
    xout.GetProperties(aurostd::EmbData::get_test_file("OUTCAR_dielectric.txt"));
    xout.GetOptical();

    // ---------------------------------------------------------------------------
    // Check 1 | dielectric_full_iso_real
    // ---------------------------------------------------------------------------
    vector<double> intra_re_sample(xout.dielectric_full_iso_real.begin(), xout.dielectric_full_iso_real.begin() + 30);
    checkEqual(intra_re_sample, intra_real_verified, check_func, "Intraband_re_accuracy", passed_checks, results);
    // ---------------------------------------------------------------------------
    // Check 2 | dielectric_full_iso_imag
    // ---------------------------------------------------------------------------
    vector<double> intra_im_sample(xout.dielectric_full_iso_imag.begin(), xout.dielectric_full_iso_imag.begin() + 30);
    checkEqual(intra_im_sample, intra_im_verified, check_func, "Intraband_im_accuracy", passed_checks, results);
    // ---------------------------------------------------------------------------
    // Check 3 | energy_loss_function_iso
    // ---------------------------------------------------------------------------
    vector<double> EELS_sample(xout.energy_loss_function_iso.begin(), xout.energy_loss_function_iso.begin() + 30);
    checkEqual(EELS_sample, EELS_verified, check_func, "EELS_accuracy", passed_checks, results);
    // ---------------------------------------------------------------------------
    // Check 4 | reflectivity_iso
    // ---------------------------------------------------------------------------
    vector<double> reflectivity_sample(xout.reflectivity_iso.begin(), xout.reflectivity_iso.begin() + 30);
    checkEqual(reflectivity_sample, reflectivity_iso_verified, check_func, "Reflectivity_accuracy", passed_checks, results);
  }

  // database
  void UnitTest::schemaTest(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    (void) errors;  // Suppress compiler warnings

    // Set up test environment
    string check_function;
    string check_description;

    // ---------------------------------------------------------------------------
    // Check | internal consistency
    // ---------------------------------------------------------------------------
    check_function = "XHOST.vschema";
    check_description = "Internal consistency of vschema";
    vector<string> vschema_keys;
    vector<string> vschema_types = {"UNIT", "TYPE"};
    string key;
    uint ninconsistent = 0;
    for (size_t i = 0; i < XHOST.vschema.vxsghost.size(); i += 2) {
      if (XHOST.vschema.vxsghost[i].find("::NAME:") != string::npos) {
        key = aurostd::RemoveSubString(XHOST.vschema.vxsghost[i], "SCHEMA::NAME:");
        vschema_keys.push_back(XHOST.vschema.getattachedscheme("SCHEMA::NAME:" + key));
        for (size_t j = 0; j < vschema_types.size(); j++) {
          if (!XHOST.vschema.isdefined("SCHEMA::" + vschema_types[j] + ":" + key)) {
            ninconsistent++;
            if (LDEBUG) {
              std::cerr << __AFLOW_FUNC__ << " SCHEMA::" << vschema_types[j] << ":" << key << " not found." << std::endl;
            }
          }
        }
      }
    }
    checkEqual(ninconsistent, (uint) 0, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | consistency between _aflowlib_entry and schema
    // ---------------------------------------------------------------------------
    check_function = "_aflowlib_entry";
    check_description = "Consistency between _aflowlib_entry json and schema";
    aflowlib::_aflowlib_entry aentry;
    const string aflowlib_json = aentry.aflowlib2string("JSON", true);
    vector<string> json_keys = aurostd::extractJsonKeysAflow(aflowlib_json);

    const vector<string> vkeys_ignore = {"data_language", "error_status", "natoms_orig", "density_orig", "volume_cell_orig", "volume_atom_orig", "spinD_magmom_orig", "freq_plasma", "dielectric_static"};
    for (size_t k = 0; k < json_keys.size(); k++) {
      const string& key = json_keys[k];
      if (!aurostd::WithinList(vkeys_ignore, key) && !aurostd::WithinList(vschema_keys, key)) {
        ninconsistent++;
        if (LDEBUG) {
          std::cerr << __AFLOW_FUNC__ << " " << key << " not found in schema." << std::endl;
        }
      }
    }
    checkEqual(ninconsistent, (uint) 0, check_function, check_description, passed_checks, results);
  }

  void UnitTest::qhullTest(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    (void) errors; // Suppress compiler warnings
    string check_function;
    string check_description;

    //Set-up test environments: create a convex hull of MnPdPt
    string test_alloy = "MnPdPt";
    string sc_auid = "aflow:0309029b2b3f8940";
    aurostd::xoption vpflow;

    vpflow.flag(nhull::convex_flags.m_nhull_print, false);
    vpflow.flag(nhull::convex_flags.m_interquartile_range, true);
    vpflow.flag(nhull::convex_flags.m_nminus1_flag, true);
    vpflow.flag(nhull::convex_flags.m_stability_criterion_flag, true);
    vpflow.addattachedscheme(nhull::convex_flags.m_alloy, test_alloy, true);
    vpflow.addattachedscheme(nhull::convex_flags.m_stability_criterion_flag, sc_auid, true);

    nhull::ConvexHull test_run = nhull::convexHull(vpflow);

    //set all tests to fail by default if we fail to initialize hull:
    //test1:
    bool hull_initialized_fail = true;
    string hull_initialized_reason = "hull_initialization_failed";
    //test2:
    bool phase_decomp_fail = true;
    bool dhull_fail = true;
    string phase_decomp_reason = "hull_initialization_failed";
    string dhull_reason = "hull_initialization_failed";
    //test3:
    bool nminus1_fail = true;
    bool stability_criterion_fail = true;
    string nminus1_reason = "hull_initialization_failed";
    string sc_reason = "hull_initialization_failed";
    //test4:
    bool below_facet_fail1 = true;
    bool below_facet_fail2 = true;
    bool below_facet_fail3 = true;
    bool below_facet_fail4 = true;
    bool below_facet_fail5 = true;
    bool below_facet_fail6 = true;
    string check_description1 = "hull_initialization_failed";
    string check_description2 = "hull_initialization_failed";
    string check_description3 = "hull_initialization_failed";
    string check_description4 = "hull_initialization_failed";
    string check_description5 = "hull_initialization_failed";
    string check_description6 = "hull_initialization_failed";
    //test5
    bool distance_fail = true;
    string distance_fail_reason = "hull_initialization_failed";

    // ---------------------------------------------------------------------------
    // Test 1: check to make sure hull was properly initialized
    // ---------------------------------------------------------------------------
    vector<nhull::Entry> calculated_points = test_run.getPoints();
    if (!calculated_points.empty()) {
      hull_initialized_fail = false;
      hull_initialized_reason = "hull_initialization_passed";
    }

    if (!hull_initialized_fail) {
     // ---------------------------------------------------------------------------
     // Test 2: verify hull equivalence with test data
     // ---------------------------------------------------------------------------
      phase_decomp_fail = false;
      dhull_fail = false;
      phase_decomp_reason = "";
      dhull_reason = "";
      aurostd::JSON::object jo = aurostd::JSON::parse(aurostd::EmbData::get_test_file("MnPdPt_hull.json"));
      vector<nhull::Entry> reference_points = static_cast<vector<nhull::Entry>>(jo);

      for (const nhull::Entry& entry : calculated_points) {
        for (const nhull::Entry& ref_entry : reference_points) {
          if (entry.auid == ref_entry.auid && !entry.flagged_entry && entry.nspecies != 1 && entry.is_hull_point == false) {
            auto phase_decomp = entry.nhull_phase_decomp;
            auto ref_phase_decomp = ref_entry.nhull_phase_decomp;
            double energy = entry.distance_hull_enthalpy_formation_atom;
            double ref_energy = ref_entry.distance_hull_enthalpy_formation_atom;

            //phase decomp test
            for (auto phase : ref_phase_decomp) {
              if (fabs(phase_decomp[phase.first] - phase.second) > NHULL_UNIT_TEST_TOL) {
                phase_decomp_fail = true;
                phase_decomp_reason += "Auid: " + entry.auid + " failed phase decomp test!\n";
              }
            }

            //distance to hull test
            //ref_energy is loaded from JSON which saves energy is meVs. Hull runs in eVs, so a conversion is done here:
            if (fabs(ref_energy / 1000 - energy) > NHULL_UNIT_TEST_TOL) {
              dhull_fail = true;
              dhull_reason += "Auid: " + entry.auid + " failed dhull test!\n";
            }
          }
        }
      }

      if (!phase_decomp_fail) {
        phase_decomp_reason = "passed phase decomp test";
      }
      if (!dhull_fail) {
        dhull_reason = "passed distance-to-hull test";
      }

      // ---------------------------------------------------------------------------
      // Test 3: verify stability criterion and nminus1
      // ---------------------------------------------------------------------------
      nminus1_fail = false;
      stability_criterion_fail = false;
      nminus1_reason = "";
      sc_reason = "";

      std::pair<string, double> verified_data_MnPdPt_auid_nminus1_1 = {"aflow:0309029b2b3f8940", 72.4623};
      std::pair<string, double> verified_data_MnPdPt_auid_nminus1_2 = {"aflow:fb9eaa58604ce774", 32.9375};
      std::pair<string, double> verified_data_MnPdPt_auid_sc_1 = {sc_auid, 72.4623};
      std::pair<string, double> verified_data_MnPdPt_auid_sc_2 = {"aflow:11ca5563e6e32342", 0.575667};

      for (const auto& entry : calculated_points) {
        //nminus1 and stability criterion is loaded from JSON which saves energy is meVs. Hull runs in eVs, so a conversion is done here:
        if (entry.auid == verified_data_MnPdPt_auid_nminus1_1.first) {
          if (fabs(verified_data_MnPdPt_auid_nminus1_1.second / 1000 - entry.m_nminus1_enthalpy_gain) > NHULL_UNIT_TEST_TOL) {
            nminus1_fail = true;
            nminus1_reason += "Auid: " + entry.auid + " failed nminus1 test!\n";
          }
        }
        if (entry.auid == verified_data_MnPdPt_auid_nminus1_2.first) {
          if (fabs(verified_data_MnPdPt_auid_nminus1_2.second / 1000 - entry.m_nminus1_enthalpy_gain) > NHULL_UNIT_TEST_TOL) {
            nminus1_fail = true;
            nminus1_reason += "Auid: " + entry.auid + " failed nminus1 test!\n";
          }
        }
        if (entry.auid == verified_data_MnPdPt_auid_sc_1.first) {
          if (fabs(verified_data_MnPdPt_auid_sc_1.second / 1000 - entry.m_stability_criterion) > NHULL_UNIT_TEST_TOL) {
            stability_criterion_fail = true;
            sc_reason += "Auid: " + entry.auid + " failed sc test!\n";
          }
        }
        if (entry.auid == verified_data_MnPdPt_auid_sc_2.first) {
          if (fabs(verified_data_MnPdPt_auid_sc_2.second / 1000 - entry.m_stability_criterion) > NHULL_UNIT_TEST_TOL) {
            stability_criterion_fail = true;
            sc_reason += "Auid: " + entry.auid + " failed sc test!\n";
          }
        }
      }

      if (!nminus1_fail) {
        nminus1_reason = "passed nminus1 test";
      }
      if (!stability_criterion_fail) {
        sc_reason = "passed sc test";
      }

      // ---------------------------------------------------------------------------
      // Test 4: facet identification
      // ---------------------------------------------------------------------------

      std::vector<aurostd::xvector<double>> v;
      aurostd::xvector<double> qhull_point;
      bool expected;

      //test1 2D
      check_description1 = "check below facet test1 passed";
      expected = true;
      below_facet_fail1 = false;
      v = {
          {1, 0, 5},
          {0, 1, 2},
          {0, 0, 8}
      };
      qhull_point = {.4, .5, 17};
      nhull::NhullFacet facet1(v, 0, 0, 3);
      if (nhull::checkPointBelowHull(facet1, qhull_point) != expected) {
        below_facet_fail1 = true;
        check_description1 = "check below facettest1 failed";
      }

      //test2 2D
      check_description2 = "check below facet test2 passed";
      expected = false;
      below_facet_fail2 = false;
      v = {
          {1, 0, 7},
          {0, 1, 2},
          {0, 0, 8}
      };
      qhull_point = {.6, .5, 17};
      nhull::NhullFacet facet2(v, 0, 0, 3);
      if (nhull::checkPointBelowHull(facet2, qhull_point) != expected) {
        below_facet_fail2 = true;
        check_description2 = "check below facet test2 failed";
      }

      //test3 2D
      check_description3 = "check below facet test3 passed";
      expected = false;
      below_facet_fail3 = false;
      v = {
          {1, 0, 2},
          {0, 1, 9},
          {0, 0, 7}
      };
      qhull_point = {1.1, 0, 17};
      nhull::NhullFacet facet3(v, 0, 0, 3);
      if (nhull::checkPointBelowHull(facet3, qhull_point) != expected) {
        below_facet_fail3 = true;
        check_description3 = "check below facet test3 failed";
      }

      //test4 2D
      check_description4 = "check below facet test4 passed";
      expected = false;
      below_facet_fail4 = false;
      v = {
          {1, 0, 4},
          {0, 1, 7},
          {0, 0, 5}
      };
      qhull_point = {-.1, 0, 17};
      nhull::NhullFacet facet4(v, 0, 0, 3);
      if (nhull::checkPointBelowHull(facet4, qhull_point) != expected) {
        below_facet_fail4 = true;
        check_description4 = "check below facet test4 failed";
      }

      //test5 3D
      check_description5 = "check below facet test5 passed";
      expected = true;
      below_facet_fail5 = false;
      v = {
          {1, 0, 0, 4},
          {0, 1, 0, 7},
          {0, 0, 0, 5},
          {0, 0, 1, 5}
      };
      qhull_point = {0.25, 0.25, 0.25, 17};
      nhull::NhullFacet facet5(v, 0, 0, 4);
      if (nhull::checkPointBelowHull(facet5, qhull_point) != expected) {
        below_facet_fail5 = true;
        check_description5 = "check below facet test5 failed";
      }

      //test6 3D
      check_description6 = "check below facet test6 passed";
      expected = false;
      below_facet_fail6 = false;
      v = {
          {1, 0, 0, 4},
          {0, 1, 0, 7},
          {0, 0, 0, 5},
          {0, 0, 1, 5}
      };
      qhull_point = {1, 0.25, 0.25, 17};
      nhull::NhullFacet facet6(v, 0, 0, 4);
      if (nhull::checkPointBelowHull(facet6, qhull_point) != expected) {
        below_facet_fail6 = true;
        check_description6 = "check below facet test6 failed";
      }

      // ---------------------------------------------------------------------------
      // Test 5: general distance formula
      // ---------------------------------------------------------------------------
      distance_fail = false;
      distance_fail_reason = "passed distance formula test";
      const aurostd::xvector<double> test_point = {0, 0.5862, -0.1137};
      std::vector<aurostd::xvector<double>> dummy_vertices;
      aurostd::xvector<double> facet_normal = {0.379, 0.417, 0.826};
      double facet_offset = 0.036734999999999962;
      const double expected_distance = 0.1873;

      nhull::NhullFacet test_facet(dummy_vertices, facet_normal, facet_offset, 3);

      if (std::fabs(test_facet.distance(test_point) - expected_distance) > 1e-4) {
        distance_fail = true;
        distance_fail_reason = "failed vertical distance to facet";
      }
    }

    //All checks done here:
    //test1:
    check_function = "check hull initialization";
    checkEqual(hull_initialized_fail, false, check_function, hull_initialized_reason, passed_checks, results);
    //test2:
    check_function = "calc_vertical_distance_and_phase_decomposition()";
    checkEqual(dhull_fail, false, check_function, dhull_reason, passed_checks, results);
    checkEqual(phase_decomp_fail, false, check_function, phase_decomp_reason, passed_checks, results);
    //test3:
    check_function = "calculateNminus1()";
    checkEqual(nminus1_fail, false, check_function, nminus1_reason, passed_checks, results);
    check_function = "calculateStabilityCriterion()";
    checkEqual(stability_criterion_fail, false, check_function, sc_reason, passed_checks, results);

    check_function = "nhull::checkPointBelowHull()";
    checkEqual(below_facet_fail1, false, check_function, check_description1, passed_checks, results);
    checkEqual(below_facet_fail2, false, check_function, check_description2, passed_checks, results);
    checkEqual(below_facet_fail3, false, check_function, check_description3, passed_checks, results);
    checkEqual(below_facet_fail4, false, check_function, check_description4, passed_checks, results);
    checkEqual(below_facet_fail5, false, check_function, check_description5, passed_checks, results);
    checkEqual(below_facet_fail6, false, check_function, check_description6, passed_checks, results);

    check_function = "NhullFacet::distance()";
    checkEqual(distance_fail, false, check_function, distance_fail_reason, passed_checks, results);
  }

  // structure tests
  void UnitTest::atomicEnvironmentTest(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    (void) errors;  // Suppress compiler warnings

    // setup test environment
    string check_function;
    string check_description;

    // ---------------------------------------------------------------------------
    // Test 1: create AE - mode 1
    // ---------------------------------------------------------------------------

    {
      // load test system
      aflowlib::_aflowlib_entry entry;
      xstructure str;
      aflowlib::EntryLoader el;
      el.m_xstructure_relaxed = true;
      el.loadAUID("aflow:d912e209c81aeb94");

      entry = *el.m_entries_flat->back();
      str = entry.vstr.back();
      vector<AtomEnvironment> AE = getAtomEnvironments(str, 1);
      check_function = "getAtomEnvironments()";

      // ---------------------------------------------------------------------------
      // Check | number of created AEs
      check_description = "number of created AEs";
      checkEqual(uint(AE.size()), uint(6), check_function, check_description, passed_checks, results);

      // ---------------------------------------------------------------------------
      // Check | point index mapping
      check_description = "point index mapping";
      checkEqual(AE[1].index2Point(10), AE[1].coordinates_neighbor[1][2], check_function, check_description, passed_checks, results);

      // ---------------------------------------------------------------------------
      // Check | coordinate matching
      check_description = "coordinate matching";
      const xvector<double> compare_point = {-2.9522760602661933, 1.020477106008512, 2.4272163095277017};
      checkEqual(AE[1].index2Point(2), compare_point, check_function, check_description, passed_checks, results);

      // ---------------------------------------------------------------------------
      // Check | center element
      check_description = "center element";
      checkEqual(std::string(AE[2].element_center), std::string("Ca"), check_function, check_description, passed_checks, results);

      // ---------------------------------------------------------------------------
      // Check | center element id
      check_description = "center element id";
      checkEqual(AE[4].type_center, uint(1), check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Test 2: create AE convex hull
    // ---------------------------------------------------------------------------

    // setup test environment
      const uint test_AE = 4;

    // create hull
      AE[test_AE].constructAtomEnvironmentHull();
      check_function = "constructAtomEnvironmentHull()";

    // ---------------------------------------------------------------------------
    // Check | hull bit set
      check_description = "hull bit set";
      checkEqual(AE[test_AE].has_hull, true, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | hull volume
      check_description = "hull volume";
      checkEqual(AE[test_AE].volume, 33.035348492423054, check_function, check_description, passed_checks, results);
    // ---------------------------------------------------------------------------
    // Check | hull area
      check_description = "hull area";
      checkEqual(AE[test_AE].area, 62.498084711616272, check_function, check_description, passed_checks, results);
    // ---------------------------------------------------------------------------
    // Check | triangle count
      check_description = "triangle count";
      checkEqual(AE[test_AE].facet_order[0], uint(6), check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | tetragon count
      check_description = "tetragon count";
      checkEqual(AE[test_AE].facet_order[1], uint(2), check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | pentagon count
      check_description = "pentagon count";
      checkEqual(AE[test_AE].facet_order[2], uint(0), check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Test 3: Soliquidy
    // ---------------------------------------------------------------------------

    {
      check_function = "soliquidy::Calculate()";
      std::vector<std::string> auid = {"aflow:d912e209c81aeb94", "aflow:03cb24fa8c754e9d"};
      aurostd::JSON::object result = soliquidy::CalculateAUID(auid[0]);
      check_description = "single AUID (" + auid[0] + ")";
      checkEqual(aurostd::round(static_cast<double>(result["aflow:d912e209c81aeb94"]["value"]), 2), aurostd::round(4.06346, 2), check_function, check_description, passed_checks, results);
      result = soliquidy::CalculateAUID(auid[1]);
      check_description = "single AUID (" + auid[1] + ")";
      checkEqual(aurostd::round(static_cast<double>(result["aflow:03cb24fa8c754e9d"]["value"]), 2), aurostd::round(24.2974, 2), check_function, check_description, passed_checks, results);

      result = soliquidy::CalculateAUID(auid);
      check_description = "multiple AUIDs (" + auid[0] + ")";
      checkEqual(aurostd::round(static_cast<double>(result["aflow:d912e209c81aeb94"]["value"]), 2), aurostd::round(4.06346, 2), check_function, check_description, passed_checks, results);
      check_description = "multiple AUIDs (" + auid[1] + ")";
      checkEqual(aurostd::round(static_cast<double>(result["aflow:03cb24fa8c754e9d"]["value"]), 2), aurostd::round(24.2974, 2), check_function, check_description, passed_checks, results);
    }
  }

  void UnitTest::xstructureParserTest(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    (void) errors;  // Suppress compiler warnings

    if (LDEBUG) {
      std::cerr << "Running " << __AFLOW_FUNC__ << std::endl;
    }

    // Set up test environment
    string check_function;
    string check_description;
    bool calculated_bool = false;
    bool expected_bool = false;

    // Set up structure variables
    _aflags aflags;
    aflags.Directory = aurostd::getPWD();
    bool check_passed = true;
    stringstream xstrss;
    xstructure xstr_cif;
    xstructure xstr_poscar;
    string str_cif;
    string str_poscar;

    const bool same_species = true;
    const bool scale_volume = false;
    const bool optimize_match = false;
    double misfit = 0.0;

    // ---------------------------------------------------------------------------
    // Check | CIF parser
    // ---------------------------------------------------------------------------

    // ---------------------------------------------------------------------------
    // Test 1: parse structure with known settings
    // ---------------------------------------------------------------------------
    if (LDEBUG) {
      std::cerr << __AFLOW_FUNC__ << " Test 1: parse structure with known settings" << std::endl;
    }

    // CrO3 was a problematic structure in the past
    str_cif = string(utd["xstructure_parser"]["cif_CrO3"]);
    str_poscar = string(utd["xstructure_parser"]["poscar_CrO3"]);

    check_function = "xstructure::operator<<";
    check_description = "Parsing CIF file with recognized setting (CrO3)";

    if (LDEBUG) {
      std::cerr << __AFLOW_FUNC__ << " [1] Parsing CIF file with recognized setting (CrO3)" << std::endl;
    }
    aurostd::StringstreamClean(xstrss);
    xstrss << str_cif;
    try {
      xstr_cif = xstructure(xstrss);
    } catch (aurostd::xerror& excpt) {
      errors.emplace_back("Could not generate CrO3");
    }
    if (LDEBUG) {
      std::cerr << __AFLOW_FUNC__ << " [2] Parsing CIF file with recognized setting (CrO3)" << std::endl;
    }

    aurostd::StringstreamClean(xstrss);
    xstrss << str_poscar;
    xstr_poscar = xstructure(xstrss);
    if (LDEBUG) {
      std::cerr << __AFLOW_FUNC__ << " [3] Parsing CIF file with recognized setting (CrO3)" << std::endl;
    }

    // ---------------------------------------------------------------------------
    // test: parse structure
    if (LDEBUG) {
      std::cerr << __AFLOW_FUNC__ << " test: parse structure" << std::endl;
    }
    expected_bool = true;
    calculated_bool = compare::aflowCompareStructure(xstr_cif, xstr_poscar, same_species, scale_volume, optimize_match, misfit);
    checkEqual(expected_bool, calculated_bool, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // test: compare Wyckoff positions
    if (LDEBUG) {
      std::cerr << __AFLOW_FUNC__ << " test: compare Wyckoff positions" << std::endl;
    }
    check_description = "Compare parsed Wyckoff positions of CrO3";
    vector<wyckoffsite_ITC> vwyckoff(4);
    const xvector<double> coords;
    vwyckoff[0].type = "Cr";
    vwyckoff[0].letter = "b";
    vwyckoff[0].site_symmetry = "m..";
    vwyckoff[0].multiplicity = 4;
    vwyckoff[0].coord[1] = 0.25;
    vwyckoff[0].coord[2] = 0.09676;
    vwyckoff[0].coord[3] = 0.5;
    vwyckoff[1].type = "O";
    vwyckoff[1].letter = "a";
    vwyckoff[1].site_symmetry = "..2";
    vwyckoff[1].multiplicity = 4;
    vwyckoff[1].coord[1] = 0.0;
    vwyckoff[1].coord[2] = 0.0;
    vwyckoff[1].coord[3] = 0.3841;
    vwyckoff[2].type = "O";
    vwyckoff[2].letter = "b";
    vwyckoff[2].site_symmetry = "m..";
    vwyckoff[2].multiplicity = 4;
    vwyckoff[2].coord[1] = 0.25;
    vwyckoff[2].coord[2] = 0.2677;
    vwyckoff[2].coord[3] = 0.3755;
    vwyckoff[3].type = "O";
    vwyckoff[3].letter = "b";
    vwyckoff[3].site_symmetry = "m..";
    vwyckoff[3].multiplicity = 4;
    vwyckoff[3].coord[1] = 0.25;
    vwyckoff[3].coord[2] = 0.6078;
    vwyckoff[3].coord[3] = 0.3284;
    check_passed = (vwyckoff.size() == xstr_cif.wyckoff_sites_ITC.size());
    for (size_t i = 0; i < vwyckoff.size() && check_passed; i++) {
      check_passed = (check_passed && (xstr_cif.wyckoff_sites_ITC[i].type == vwyckoff[i].type));
      check_passed = (check_passed && (xstr_cif.wyckoff_sites_ITC[i].letter == vwyckoff[i].letter));
      check_passed = (check_passed && (xstr_cif.wyckoff_sites_ITC[i].multiplicity == vwyckoff[i].multiplicity));
      check_passed = (check_passed && (xstr_cif.wyckoff_sites_ITC[i].coord == vwyckoff[i].coord));
      if (LDEBUG && !check_passed) {
        std::cerr << "Failed site:" << std::endl;
        std::cerr << "type = " << xstr_cif.wyckoff_sites_ITC[i].type << ", ";
        std::cerr << "letter = " << xstr_cif.wyckoff_sites_ITC[i].letter << ", ";
        std::cerr << "multiplicity = " << xstr_cif.wyckoff_sites_ITC[i].multiplicity << ", ";
        std::cerr << "coord = " << xstr_cif.wyckoff_sites_ITC[i].coord << std::endl << std::endl;
        std::cerr << "Should be:" << std::endl;
        std::cerr << "type = " << vwyckoff[i].type << ", ";
        std::cerr << "letter = " << vwyckoff[i].letter << ", ";
        std::cerr << "multiplicity = " << vwyckoff[i].multiplicity << ", ";
        std::cerr << "coord = " << vwyckoff[i].coord << std::endl;
      }
    }

    expected_bool = true;
    calculated_bool = check_passed;
    checkEqual(expected_bool, calculated_bool, check_function, check_description, passed_checks, results);

    // May need a better test case where the labels actually change
    check_description = "Compare calculated Wyckoff positions of CrO3";
    xstr_cif.SpaceGroup_ITC();
    vwyckoff[0].coord[1] = 0.25;
    vwyckoff[0].coord[2] = 0.59676;
    vwyckoff[0].coord[3] = 0.0000;
    vwyckoff[1].coord[1] = 0.00;
    vwyckoff[1].coord[2] = 0.00000;
    vwyckoff[1].coord[3] = 0.3841;
    vwyckoff[2].coord[1] = 0.25;
    vwyckoff[2].coord[2] = 0.76770;
    vwyckoff[2].coord[3] = 0.8755;
    vwyckoff[3].coord[1] = 0.25;
    vwyckoff[3].coord[2] = 0.60780;
    vwyckoff[3].coord[3] = 0.3284;
    check_passed = (vwyckoff.size() == xstr_cif.wyckoff_sites_ITC.size());
    for (size_t i = 0; i < vwyckoff.size() && check_passed; i++) {
      check_passed = (check_passed && (xstr_cif.wyckoff_sites_ITC[i].type == vwyckoff[i].type));
      check_passed = (check_passed && (xstr_cif.wyckoff_sites_ITC[i].letter == vwyckoff[i].letter));
      check_passed = (check_passed && (xstr_cif.wyckoff_sites_ITC[i].multiplicity == vwyckoff[i].multiplicity));
      check_passed = (check_passed && aurostd::isequal(xstr_cif.wyckoff_sites_ITC[i].coord, vwyckoff[i].coord));
      if (LDEBUG && !check_passed) {
        std::cerr << "Failed site:" << std::endl;
        std::cerr << "type = " << xstr_cif.wyckoff_sites_ITC[i].type << ", ";
        std::cerr << "letter = " << xstr_cif.wyckoff_sites_ITC[i].letter << ", ";
        std::cerr << "multiplicity = " << xstr_cif.wyckoff_sites_ITC[i].multiplicity << ", ";
        std::cerr << "coord = " << xstr_cif.wyckoff_sites_ITC[i].coord << std::endl << std::endl;
        std::cerr << "Should be:" << std::endl;
        std::cerr << "type = " << vwyckoff[i].type << ", ";
        std::cerr << "letter = " << vwyckoff[i].letter << ", ";
        std::cerr << "multiplicity = " << vwyckoff[i].multiplicity << ", ";
        std::cerr << "coord = " << vwyckoff[i].coord << std::endl;
      }
    }

    expected_bool = true;
    calculated_bool = check_passed;
    checkEqual(expected_bool, calculated_bool, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Test 2: parse structure with unrecognized (old) settings
    // ---------------------------------------------------------------------------
    if (LDEBUG) {
      std::cerr << __AFLOW_FUNC__ << " Test 2: parse structure with unrecognized (old) settings" << std::endl;
    }
    check_description = "Parsing CIF file with unrecognized setting (GePt3)";
    aurostd::StringstreamClean(xstrss);
    str_cif = string(utd["xstructure_parser"]["cif_Ge8Pt24"]);
    str_poscar = string(utd["xstructure_parser"]["poscar_Ge8Pt24"]);

    const bool quiet_tmp = XHOST.QUIET;
    XHOST.QUIET = !LDEBUG;  // Suppress warnings
    xstrss << str_cif;
    xstr_cif = xstructure(xstrss);
    XHOST.QUIET = quiet_tmp;

    // ---------------------------------------------------------------------------
    // test: parse structure
    aurostd::StringstreamClean(xstrss);
    xstrss << str_poscar;
    xstr_poscar = xstructure(xstrss);
    expected_bool = true;
    calculated_bool = compare::aflowCompareStructure(xstr_cif, xstr_poscar, same_species, scale_volume, optimize_match, misfit);
    checkEqual(expected_bool, calculated_bool, check_function, check_description, passed_checks, results);
  }

  void UnitTest::kpathTest(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    (void) errors;  // Suppress compiler warnings
    // Set up test environment
    string check_function;
    string check_description;
    // initializing variables for kpath
    const double grid = 20.0;
//    string stmp = "", line = "";
//      lattice = ""  , kpath_string = "";
    bool foundBZ = false;
    xstructure str_sp;
    xstructure str_sc;
    // list of structure types to test is in LATTICE::lattice_names

    // ---------------------------------------------------------------------------
    // Test 1: vasp kpath test
    // ---------------------------------------------------------------------------
    for (const std::pair<string, string> it : LATTICE::lattice_names) {
      foundBZ = false;
      const string lattice_name = it.first;
      string lat_xstr_str_out;
      const xstructure lat_xstr_out;
      lat_xstr_str_out = string(utd["kpath"][lattice_name]["kpath_vasp"]); // access vasp kpath from unit test data
      stringstream lat_xstrss((string(utd["kpath"][lattice_name]["poscar"]))); // access poscar from unit test data json
      xstructure lat_xstr_in;
      lat_xstrss >> lat_xstr_in;
      LATTICE::Standard_Lattice_StructureDefault(lat_xstr_in, str_sp, str_sc); // This functions format the str_in, and stores data in str_sp and str_sc
      const std::string lattice = str_sp.bravais_lattice_variation_type;
      const std::string kpath_string = LATTICE::KPOINTS_Directions(lattice, str_sp.lattice, grid, str_sp.iomode, foundBZ);
      check_function = "KPOINTS_Directions";
      check_description = "checks that correct vasp kpath is generated for " + lattice_name;
      std::vector<std::string> vkpath_string;
      std::vector<std::string> vlat_xstr_str_out;
      aurostd::string2vectorstring(kpath_string, vkpath_string);
      aurostd::string2vectorstring(lat_xstr_str_out, vlat_xstr_str_out);
      checkEqual(vkpath_string, vlat_xstr_str_out, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Test 2: quantum espresso kpath test
    // ---------------------------------------------------------------------------
    for (const std::pair<string, string> it : LATTICE::lattice_names) {
      foundBZ = false;
      const string lattice_name = it.first;
      const string lat_xstr_str_in;
      string lat_xstr_str_out;
      const xstructure lat_xstr_out;
      lat_xstr_str_out = string(utd["kpath"][lattice_name]["kpath_qe"]);
      stringstream lat_xstrss;
      lat_xstrss.str(string(utd["kpath"][lattice_name]["poscar"]));
      lat_xstrss << lat_xstr_str_in;
      xstructure lat_xstr_in;
      try {
        lat_xstrss >> lat_xstr_in;          // CO20200404 - this WILL throw an error because det(lattice)<0.0, leave alone
      } catch (aurostd::xerror& excpt) {} // CO20200404 - this WILL throw an error because det(lattice)<0.0, leave alone
      LATTICE::Standard_Lattice_StructureDefault(lat_xstr_in, str_sp, str_sc); // This functions format the str_in, and stores data in str_sp and str_sc
      const std::string lattice = str_sp.bravais_lattice_variation_type;
      str_sp.iomode = IOQE_AUTO; // QE IO necessary
      const std::string kpath_string = LATTICE::KPOINTS_Directions(lattice, str_sp.lattice, grid, str_sp.iomode, foundBZ);
      check_function = "KPOINTS_Directions";
      check_description = "checks that correct quantum espresso kpath is generated for " + lattice_name;
      std::vector<std::string> vkpath_string;
      std::vector<std::string> vlat_xstr_str_out;
      aurostd::string2vectorstring(kpath_string, vkpath_string);
      aurostd::string2vectorstring(lat_xstr_str_out, vlat_xstr_str_out);
      checkEqual(vkpath_string, vlat_xstr_str_out, check_function, check_description, passed_checks, results);
    }
  }

  void UnitTest::xstructureTest(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    (void) errors;  // Suppress compiler warnings
    const bool LDEBUG = (false || XHOST.DEBUG);
    // Set up test environment
    string check_function;
    string check_description;
    bool expected_bool = false;
    bool calculated_bool = false;
    uint expected_uint = 0;
    uint calculated_uint = 0;

    // ---------------------------------------------------------------------------
    // Check | getCoordinations() //CO20190520
    // ---------------------------------------------------------------------------
    xstructure xstr("aflowlib.duke.edu:AFLOWDATA/ICSD_WEB/FCC/Cl1Na1_ICSD_240599", "CONTCAR.relax.vasp", IOAFLOW_AUTO);
    deque<deque<uint>> coordinations;
    xstr.GetCoordinations(coordinations);

    check_function = "xstructure::getCoordinations()";
    check_description = "Number of iatoms";
    calculated_uint = coordinations.size();
    expected_uint = 2;
    checkEqual(calculated_uint, expected_uint, check_function, check_description, passed_checks, results);

    check_description = "Number of coordination environments atom 1 > 2";
    expected_bool = true;
    calculated_bool = (coordinations[0].size() > 2);
    checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

    check_description = "Number of coordination environments atom 2 > 2";
    expected_bool = true;
    calculated_bool = (coordinations[1].size() > 2);
    checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

    // first iatom
    // first shell
    check_description = "First shell atom 1";
    calculated_uint = coordinations[0][0];
    expected_uint = 6;
    checkEqual(calculated_uint, expected_uint, check_function, check_description, passed_checks, results);

    // second shell
    check_description = "Second shell atom 1";
    calculated_uint = coordinations[0][1];
    expected_uint = 12;
    checkEqual(calculated_uint, expected_uint, check_function, check_description, passed_checks, results);

    // second iatom
    // first shell
    check_description = "First shell atom 2";
    calculated_uint = coordinations[1][0];
    expected_uint = 6;
    checkEqual(calculated_uint, expected_uint, check_function, check_description, passed_checks, results);

    // second shell
    check_description = "Second shell atom 2";
    calculated_uint = coordinations[1][1];
    expected_uint = 12;
    checkEqual(calculated_uint, expected_uint, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | foldAtomdInCell //DX20210129
    // ---------------------------------------------------------------------------

    // set fold atoms in cell variables
    const bool skew = false;
    const double tol = 0.01;
    const bool check_min_dists = false;

    // ---------------------------------------------------------------------------
    // test 1: expand cell
    // create 3x1x1 supercell expansion matrix
    check_function = "xstructure::foldAtomsInCell()";
    check_description = "expand cell (3x1x1)";
    xmatrix<double> supercell_matrix = aurostd::eye<double>(3, 3);
    supercell_matrix(1, 1) = 3.0;
    const xmatrix<double> lattice_new = supercell_matrix * xstr.lattice; // transform lattice

    xstructure xstr_supercell = xstr;
    xstr_supercell.foldAtomsInCell(lattice_new, skew, tol, check_min_dists);

    const bool same_species = true;
    expected_bool = true;
    calculated_bool = compare::structuresMatch(xstr, xstr_supercell, same_species);
    checkEqual(expected_bool, calculated_bool, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // test 2: reduce cell
    // convert supercell back to original lattice
    check_function = "xstructure::foldAtomsInCell()";
    check_description = "reduce cell into primitive";
    xstructure xstr_reduced = xstr_supercell;
    xstr_reduced.foldAtomsInCell(xstr.lattice, skew, tol, check_min_dists);

    expected_bool = true;
    calculated_bool = compare::structuresMatch(xstr, xstr_reduced, same_species);
    checkEqual(expected_bool, calculated_bool, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | slab test //CO20190520
    // ---------------------------------------------------------------------------
    // See W. Sun and G. Ceder, Surface Science 617 (2013) 53-59
    string xstr_str;
    stringstream xstrss;
    double min_dist = 0.0;
    double min_dist_orig = 0.0;
    const xvector<int> hkl({1, 0, 4});
    //    hkl[1] = 1; hkl[2] = 0; hkl[3] = 4;

    // create input structure
    xstr_str = string(utd["xstructure"]["poscar_FeO"]);
    aurostd::StringstreamClean(xstrss);
    xstrss << xstr_str;
    xstructure xstr_in;
    try {
      xstrss >> xstr_in;          // CO20200404 - this WILL throw an error because det(lattice)<0.0, leave alone
    } catch (aurostd::xerror& excpt) {} // CO20200404 - this WILL throw an error because det(lattice)<0.0, leave alone
    min_dist = min_dist_orig = xstr_in.MinDist();
    if (LDEBUG) {
      std::cerr << __AFLOW_FUNC__ << " xstr_in=\n" << xstr_in << std::endl;
      std::cerr << __AFLOW_FUNC__ << " xstr_in.MinDist()=" << min_dist << std::endl;
    }

    // create xstr_slab (correct answer)
    xstr_str = string(utd["xstructure"]["poscar_FeO_slab"]);
    xstructure xstr_slab_correct;
    aurostd::StringstreamClean(xstrss);
    xstrss << xstr_str;
    try {
      xstrss >> xstr_slab_correct; // CO20200404 - this WILL throw an error because det(lattice)<0.0, leave alone
    } catch (aurostd::xerror& excpt) {} // CO20200404 - this WILL throw an error because det(lattice)<0.0, leave alone
                                        //  ---------------------------------------------------------------------------
                                        //  test 1: compare min distance of bulk and correct slab
    min_dist = xstr_slab_correct.MinDist();
    if (LDEBUG) {
      std::cerr << __AFLOW_FUNC__ << " xstr_slab_correct=\n" << xstr_slab_correct << std::endl;
      min_dist = xstr_slab_correct.MinDist();
      std::cerr << __AFLOW_FUNC__ << " xstr_slab_correct.MinDist()=" << min_dist << std::endl;
    }
    check_function = "xstructure::CreateSlab_SurfaceLattice()";
    check_description = "Compare minimum distance bulk vs. correct slab";
    checkEqual(min_dist, min_dist_orig, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // test 2: compare min distance of generated slab and correct slab
    const xmatrix<double> lattice_slab_origbasis;
    xstructure xstr_slab_test = slab::CreateSlab_SurfaceLattice(xstr_in, hkl, 1, 0, 5.0); // the v3len_max_strict is very important here, as the test from Sun et al. takes a shortcut here
    min_dist = xstr_slab_test.MinDist();
    check_function = "xstructure::CreateSlab_SurfaceLattice(hkl = 104)";
    check_description = "Compare minimum distance slab vs. correct slab";
    if (LDEBUG) {
      std::cerr << __AFLOW_FUNC__ << " xstr_slab_test=\n" << xstr_slab_test << std::endl;
      min_dist = xstr_slab_test.MinDist();
      std::cerr << __AFLOW_FUNC__ << " xstr_slab_test.MinDist()=" << min_dist << std::endl;
    }
    checkEqual(min_dist, min_dist_orig, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // test 3: compare structures of generated slab and correct slab
    check_description = "Match structures of generated slab and correct slab";
    expected_bool = true;
    calculated_bool = compare::structuresMatch(xstr_slab_correct, xstr_slab_test, true, false, false);
    checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | GetStandardPrimitive //CO20230220
    // ---------------------------------------------------------------------------

    // create xstr_str
    xstr_str
        = "O Co Sr\n"
          "1.00000000000000\n"
          "  5.4205419737648093    0.0000000000000000    0.0000000000000000\n"
          "  0.0000000000000000    5.4205419737648093    0.0000000000000000\n"
          "  0.0000000000000000    0.0000000000000000    7.6623790337915390\n"
          "O    Co   Sr\n"
          "12    4    4\n"
          "Direct\n"
          "  0.7499997640097504  0.2500002359902496  0.5000000000000000\n"
          "  0.2500002359902496  0.7499997640097504  0.5000000000000000\n"
          "  0.2499997640097504  0.7500002359902496  0.0000000000000000\n"
          "  0.7500002359902496  0.2499997640097504  0.0000000000000000\n"
          "  0.7499997640097504  0.7499997640097504  0.5000000000000000\n"
          "  0.2500002359902496  0.2500002359902496  0.5000000000000000\n"
          "  0.2499997640097504  0.2499997640097504  0.0000000000000000\n"
          "  0.7500002359902496  0.7500002359902496  0.0000000000000000\n"
          "  0.5000000000000000  0.5000000000000000  0.7499999723900999\n"
          "  0.0000000000000000  0.0000000000000000  0.7500000276099001\n"
          "  0.0000000000000000  0.0000000000000000  0.2499999723900999\n"
          "  0.5000000000000000  0.5000000000000000  0.2500000276099001\n"
          "  0.0000000000000000  0.0000000000000000  0.5000000000000000\n"
          "  0.5000000000000000  0.5000000000000000  0.0000000000000000\n"
          "  0.5000000000000000  0.5000000000000000  0.5000000000000000\n"
          "  0.0000000000000000  0.0000000000000000  0.0000000000000000\n"
          "  0.5000000000000000  0.0000000000000000  0.7500000000000000\n"
          "  0.0000000000000000  0.5000000000000000  0.7500000000000000\n"
          "  0.5000000000000000  0.0000000000000000  0.2500000000000000\n"
          "  0.0000000000000000  0.5000000000000000  0.2500000000000000\n";
    aurostd::StringstreamClean(xstrss);
    xstrss << xstr_str;
    try {
      xstrss >> xstr_in;          // CO20200404 - this WILL throw an error because det(lattice)<0.0, leave alone
    } catch (aurostd::xerror& excpt) {} // CO20200404 - this WILL throw an error because det(lattice)<0.0, leave alone

    // ---------------------------------------------------------------------------
    // test 1: expand cell
    // create 3x1x1 supercell expansion matrix
    check_function = "xstructure::GetStandardPrimitive()";
    check_description = "resolve primitive cell and find Bravais lattice type";
    xstr_in.GetStandardPrimitive();

    expected_bool = true;
    calculated_bool = (xstr_in.bravais_lattice_type == "CUB");
    checkEqual(expected_bool, calculated_bool, check_function, check_description, passed_checks, results);
  }

  // structure generation
  void UnitTest::ceramgenTest(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    (void) errors;  // Suppress compiler warnings

    // setup test environment
    string check_function;
    string check_description;
    bool calculated_bool = false;
    bool expected_bool = false;
    uint calculated_uint = 0;
    uint expected_uint = 0;

    //./aflow --generate_ceramics --nm=N,C --m=Co,Mo,Fe,Ru,Ni,Rh,Pt,Cu,Cr,V --N=5
    const vector<string> vnonmetals{"N", "C"};
    const vector<string> vmetals{"Co", "Mo", "Fe", "Ru", "Ni", "Rh", "Pt", "Cu", "Cr", "V"};

    check_function = "pflow::GENERATE_CERAMICS()";
    const vector<string> commands = pflow::GENERATE_CERAMICS(vnonmetals, vmetals, 5);

    check_description = "number of commands";
    calculated_uint = commands.size();
    expected_uint = 6;
    checkEqual(calculated_uint, expected_uint, check_function, check_description, passed_checks, results);

    // C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N
    // C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V
    // C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V
    // C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V
    // C:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V
    // C:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V

    string search;

    search = "C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N";
    check_description = "find " + search;
    expected_bool = true;
    calculated_bool = aurostd::WithinList(commands, search);
    checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

    search = "C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V";
    check_description = "find " + search;
    calculated_bool = aurostd::WithinList(commands, search);
    checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

    search = "C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V";
    check_description = "find " + search;
    calculated_bool = aurostd::WithinList(commands, search);
    checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

    search = "C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V";
    check_description = "find " + search;
    calculated_bool = aurostd::WithinList(commands, search);
    checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

    search = "C:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V";
    check_description = "find " + search;
    calculated_bool = aurostd::WithinList(commands, search);
    checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

    search = "C:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V";
    check_description = "find " + search;
    calculated_bool = aurostd::WithinList(commands, search);
    checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);
  }

  // DX20200925
  // ME20220324 - refactored to run in parallel
  void _testPrototype(uint i,
                      const vector<string>& vuid,
                      vector<uint>& nprotos,
                      vector<string>& errors
#ifdef AFLOW_MULTITHREADS_ENABLE
                      ,
                      std::mutex& m
#endif
  ) {
    double tolerance_sym = 0.0;
    string label_input;
    bool generated = false;
    bool sym = false;
    bool unique = false;
    xstructure xstr;
    ofstream ofs("/dev/null");
    const anrl::ProtoData pd = anrl::ProtoData::get();
    const std::string& uid = vuid[i];

    string error;
    try {
      // TODO use real parameters not string
      xstr = aflowlib::PrototypeLibraries(ofs, static_cast<string>(pd.content[uid]["label"]), static_cast<string>(pd.content[uid]["parameter_values"]), 1);
      generated = true;
    } catch (aurostd::xerror& excpt) {
      error = "Could not generate prototype=" + static_cast<string>(pd.content[uid]["label"]) + ", params=" + static_cast<string>(pd.content[uid]["parameter_values"]);
    }

    if (error.empty()) {
      stringstream label_input_ss;
      label_input_ss << static_cast<string>(pd.content[uid]["label"]);
      label_input = label_input_ss.str();
      tolerance_sym = anrl::specialCaseSymmetryTolerances(label_input);
      if (tolerance_sym != AUROSTD_MAX_DOUBLE) {
        xstr.sym_eps = tolerance_sym;
        xstr.sym_eps_calculated = true;
      }
    }

    if (error.empty()) {
      string updated_label_and_params;
      if (!anrl::structureAndLabelConsistent(xstr, static_cast<string>(pd.content[uid]["label"]), updated_label_and_params, tolerance_sym)) { // DX20201105 - added symmetry tolerance
        error = "The structure has a higher symmetry than indicated by the label (orig: proto="
              + static_cast<string>(pd.content[uid]["label"])
              + ", params="
              + static_cast<string>(pd.content[uid]["parameter_values"])
              + ")."
              + " The correct label and parameters for this structure are:\n"
              + updated_label_and_params;
      } else {
        sym = true;
      }
    }
    if (error.empty()) {
      const aurostd::xoption vpflow;
      // check if the prototype matches to more than one prototype
      // (i.e., a prototype should match with itself, but no others)
      const string catalog = "anrl";
      const vector<string> protos_matching = compare::getMatchingPrototypes(xstr, vpflow, catalog);
      // if it matches to more than one
      if (protos_matching.size() > 1 && !anrl::isSpecialCaseEquivalentPrototypes(protos_matching)) {
        error = static_cast<string>(pd.content[uid]["label"])
              + ", params="
              + static_cast<string>(pd.content[uid]["parameter_values"])
              + " matches multiple prototypes (and not a documented special case): "
              + aurostd::joinWDelimiter(protos_matching, ",")
              + "."
              + " If the prototype was newly added, ONLY include it in the encyclopedia"
              + " for a valid reason (e.g., historical, special designation, etc.)"
              + " and document this in anrl::isSpecialCaseEquivalentPrototypes().";
        // if it doesn't match with ITSELF
      } else if (protos_matching.empty()) {
        error = static_cast<string>(pd.content[uid]["label"])
              + ", params="
              + static_cast<string>(pd.content[uid]["parameter_values"])
              + " does not match to any prototypes"
              + " (requires special symmetry tolerance or there is a bug with XtalFinder).";
      } else {
        unique = true;
      }
    }
#ifdef AFLOW_MULTITHREADS_ENABLE
    std::lock_guard<std::mutex> const lk(m);
#endif
    if (generated) {
      nprotos[0]++;
    }
    if (sym) {
      nprotos[1]++;
    }
    if (unique) {
      nprotos[2]++;
    }
    if (!error.empty()) {
      errors.push_back(error);
    }
  }

  void UnitTest::prototypeGeneratorTest(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    const bool LDEBUG = (false || XHOST.DEBUG);

    // Set up test environment
    string check_function;
    string check_description;

    const anrl::ProtoData pd = anrl::ProtoData::get();
    vector<string> prototype_labels = pd.content.keys();
    const uint num_protos = prototype_labels.size();

    if (LDEBUG) {
      std::cerr << __AFLOW_FUNC__ << "Number of prototype labels = " << num_protos << std::endl;
    }

    // Test
    // 1: if the prototype can be generated,
    // 2: if symmetry and label are consistent
    // 3: if it is a unique prototye
    // Keep results in vector to simplify function input
    vector<uint> nprotos(3, 0);
#ifdef AFLOW_MULTITHREADS_ENABLE
    xthread::xThread xt(KBIN::get_NCPUS()); // Okay to be greedy - xThread will manage number of threads
    std::mutex m;
    xt.run(num_protos, _testPrototype, prototype_labels, nprotos, errors, m);
#else
    for (uint i = 0; i < num_protos; i++) _testPrototype(i, prototype_labels, nprotos, errors);
#endif

    check_function = "aflowlib::PrototypeLibraries()";
    check_description = "generate prototypes";
    checkEqual(nprotos[0], num_protos, check_function, check_description, passed_checks, results);

    check_function = "anrl::structureAndLabelConsistent()";
    check_description = "symmetry consistent with prototype label";
    checkEqual(nprotos[1], num_protos, check_function, check_description, passed_checks, results);

    check_function = "compare::getMatchingPrototypes()";
    check_description = "protoypes are unique";
    checkEqual(nprotos[2], num_protos, check_function, check_description, passed_checks, results);
  }

  /// Tests functionality of the --lib2raw option with the --local option.
  /// Note that test is intended to only be run where the database is accessible.
  /// Also note that LIB2RAW requires the following debian packages:
  /// gnuplot texlive-binaries imagemagick texlive-full
  void UnitTest::lib2rawTest(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    // ST20241107

    const string check_function = "aflowlib::LIB2RAW()";

    const string system = "/FCC/Au1_ICSD_53763";
    const string system_LIB_path = string("/common/ICSD/").append("LIB").append(system);
    const string system_RAW_path = string("/common/ICSD/").append("RAW").append(system);
    const string system_WEB_path = string("/common/ICSD/").append("WEB").append(system);
    const string tmpdir = aurostd::TmpDirectoryCreate("lib2raw_test");
    const string test_path = tmpdir + system_LIB_path;

    if (!std::filesystem::exists(system_LIB_path)) {
      errors.emplace_back("The system to test does not exist! Ensure the /common/ database is mounted and that " + system_LIB_path + " is available.");
      return;
    }
    if (!std::filesystem::exists(system_RAW_path) || !std::filesystem::exists(system_WEB_path)) {
      errors.emplace_back("The RAW and WEB locations corresponding to " + system_LIB_path + " do not exist!");
      return;
    }
    std::filesystem::create_directories(test_path);
    std::filesystem::copy(system_LIB_path, test_path);

    // -----------------------------------------------------------------------
    // Check usage of `--lib2raw dir` when the database already has the entry
    // -----------------------------------------------------------------------
    if (XHOST.hostname == XHOST.AFLOW_MATERIALS_SERVER || XHOST.hostname == XHOST.AFLOW_WEB_SERVER) {
      // Check for this because LIB2LIB requires it, otherwise we can skip this set of checks because we aren't on the server

      const string output = RedirectStream::run_with_capture_cout([&]() { aflowlib::LIB2RAW(test_path, false, false); }).str();
      std::cout << output << std::endl;
      vector<string> vlines;
      aurostd::string2vectorstring(output, vlines);
      string found_LIB;
      string found_RAW;
      string found_WEB;
      bool already_calculated = false;
      for (const auto& line : vlines) {
        for (const auto& [found, key] : vector<std::pair<string*, string>>{
                 {&found_LIB, "directory_LIB"},
                 {&found_RAW, "directory_RAW"},
                 {&found_WEB, "directory_WEB"}
        }) {
          if (aurostd::substring2bool(line, key)) {
            const size_t pos = line.find('=');
            *found = line.substr(pos + 1);
          }
        }
        if (aurostd::substring2bool(line, "ALREADY CALCULATED")) {
          already_calculated = true;
        }
      }
      checkEqual(found_LIB, system_LIB_path, check_function, "find directory_LIB", passed_checks, results);
      checkEqual(found_RAW, system_RAW_path, check_function, "find directory_RAW", passed_checks, results);
      checkEqual(found_WEB, system_WEB_path, check_function, "find directory_WEB", passed_checks, results);
      check(already_calculated, already_calculated, true, check_function, "The entry should already be calculated", passed_checks, results);
    }

    // ------------------------------------------------------------------------------------------
    // Check some consistency of files produced by `--lib2raw dir --local` with existing database
    // ------------------------------------------------------------------------------------------
    try {
      aflowlib::LIB2RAW(test_path, true, true);
    } catch (aurostd::xerror& e) {
      errors.push_back(e.what());
      std::filesystem::remove_all(tmpdir);
      return;
    }

    const string test_RAW = test_path + "/RAW";
    const vector<string> outs{
        DEFAULT_AFLOW_PGROUP_OUT,  DEFAULT_AFLOW_PGROUPK_OUT,  DEFAULT_AFLOW_FGROUP_OUT,  DEFAULT_AFLOW_PGROUP_XTAL_OUT,  DEFAULT_AFLOW_PGROUPK_XTAL_OUT,  DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT,
        DEFAULT_AFLOW_IATOMS_OUT,  DEFAULT_AFLOW_AGROUP_OUT,

        DEFAULT_AFLOW_PGROUP_JSON, DEFAULT_AFLOW_PGROUPK_JSON, DEFAULT_AFLOW_FGROUP_JSON, DEFAULT_AFLOW_PGROUP_XTAL_JSON, DEFAULT_AFLOW_PGROUPK_XTAL_JSON, DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON,
        DEFAULT_AFLOW_IATOMS_JSON, DEFAULT_AFLOW_AGROUP_JSON};
    const vector<string> pre_exts{"bands", "relax", "orig"};
    const string extension = ".xz";
    for (const string& base : outs) {
      for (const string& pre : pre_exts) {
        string file;
        aflowlib::AddFileNameBeforeExtension(base, pre, file);
        file += extension;
        checkFiles(string(test_RAW).append("/").append(file), string(system_RAW_path).append("/").append(file), check_function, "checking equality for file: " + file, passed_checks, results);
      }
    }

    const vector<string> other_checks{
        "Au1_ICSD_53763_bandsdata.json.xz",
        "Au1_ICSD_53763_dosdata.json.xz",
        "OSZICAR.static",
        "INCAR.bands",
        "LOCK",

        "CONTCAR.relax",
        "CONTCAR.relax1",
        "CONTCAR.relax.abinit",
        "CONTCAR.relax.aims",
        "CONTCAR.relax.vasp",

        "KPOINTS.relax",
        "KPOINTS.static",

        "OUTCAR.bands",
        "OUTCAR.relax.xz",

        "POSCAR.bands",
        "POSCAR.orig"};
    for (const string& file : other_checks) {
      checkFiles(string(test_RAW).append("/").append(file), string(system_RAW_path).append("/").append(file), check_function, "checking equality for file: " + file, passed_checks, results);
    }

    std::filesystem::remove_all(tmpdir);
  }

  void UnitTest::surface_TriangleArea(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    const string check_func = "surface::TriangleArea";

    // basic in the xy plane
    xvector<double> v1 = {0.0, 0.0, 0.0};
    xvector<double> v2 = {2.0, 0.0, 0.0};
    xvector<double> v3 = {1.0, 2.0, 0.0};
    double area = surface::TriangleArea(v1, v2, v3);
    checkEqual(area, 2.0, check_func, "Basic area of triangle in xy plane.", passed_checks, results);

    // zero area
    v3 = {0.0, 0.0, 0.0};
    area = surface::TriangleArea(v1, v2, v3);
    checkEqual(area, 0.0, check_func, "Zero-area triangle in xy plane.", passed_checks, results);

    // skinny triangle
    v3 = {1.0e2, 2.0, 0.0};
    area = surface::TriangleArea(v1, v2, v3);
    checkEqual(area, 2.0, check_func, "Skinny triangle in xy plane.", passed_checks, results);

    // very skinny triangle
    v3 = {1.0e6, 2.0, 0.0};
    area = surface::TriangleArea(v1, v2, v3);
    checkEqual(area, 2.0, check_func, "Very skinny triangle in xy plane.", passed_checks, results);

    // equilateral triangle through cube
    v1 = {0.0, 1.0, 1.0};
    v2 = {1.0, 0.0, 1.0};
    v3 = {1.0, 1.0, 0.0};
    area = surface::TriangleArea(v1, v2, v3);
    // for cube with edge of 1, triangle with edge of sqrt(2), A = (a^2 * sqrt(3)) / 4
    checkEqual(area, 0.5 * sqrt(3), check_func, "Equilateral triangle through cube.", passed_checks, results);
  }

  void UnitTest::surface_PlaneGetABCD(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    const string check_func = "surface::PlaneGetABCD";

    double a;
    double b;
    double c;
    double d;
    xvector<double> v1 = {0.0, 0.0, 0.0};
    xvector<double> v2 = {2.0, 0.0, 0.0};
    xvector<double> v3 = {1.0, 2.0, 0.0};
    surface::PlaneGetABCD(a, b, c, d, v1, v2, v3);
    checkEqual(xvector<double>{a, b, c, d}, xvector<double>{0, 0, 4, 0}, check_func, "check the coefficients for basic triangle", passed_checks, results);
    v3 = {1.0e2, 2.0, 0.0};
    surface::PlaneGetABCD(a, b, c, d, v1, v2, v3);
    checkEqual(xvector<double>{a, b, c, d}, xvector<double>{0, 0, 4, 0}, check_func, "check the coefficients for skinny triangle", passed_checks, results);
    v3 = {1.0e6, 2.0, 0.0};
    surface::PlaneGetABCD(a, b, c, d, v1, v2, v3);
    checkEqual(xvector<double>{a, b, c, d}, xvector<double>{0, 0, 4, 0}, check_func, "check the coefficients for very skinny triangle", passed_checks, results);
    v3 = {0.0, 0.0, 0.0};
    surface::PlaneGetABCD(a, b, c, d, v1, v2, v3);
    checkEqual(xvector<double>{a, b, c, d}, xvector<double>{0, 0, 0, 0}, check_func, "check the coefficients for zero area", passed_checks, results);
    v1 = {0.0, 1.0, 1.0};
    v2 = {1.0, 0.0, 1.0};
    v3 = {1.0, 1.0, 0.0};
    surface::PlaneGetABCD(a, b, c, d, v1, v2, v3);
    checkEqual(xvector<double>{a, b, c, d}, xvector<double>{1, 1, 1, -2}, check_func, "check the coefficients for equilateral through cube", passed_checks, results);
  }

  void UnitTest::surface_PlaneDistance(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    const string check_func = "surface::PlaneDistance";
    const xvector<double> r = {0, 0, 0};

    double dist = surface::PlaneDistance(r, 1.0, 1.0, 1.0, -2.0);
    // not a typo, the area of this triangle is also the same number as half the cube diagonal
    checkEqual(dist, 2.0 / 3.0 * sqrt(3.0), check_func, "check the distance for equilateral through cube", passed_checks, results);

    dist = surface::PlaneDistance(r, 0.0, 0.0, 0.0, 0.0);
    check(std::isnan(dist), dist, -std::numeric_limits<double>::signaling_NaN(), check_func, "check the distance for zero area", passed_checks, results);

    dist = surface::PlaneDistance(r, 1.0, 0.0, 0.0, -1.0);
    checkEqual(dist, 1.0, check_func, "check distance for simple point in x-direction", passed_checks, results);
  }

  void UnitTest::surface_PlaneGetProjection(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    const string check_func = "surface::PlaneGetProjection";
    const xvector<double> r = {0, 0, 0};

    xvector<double> proj = surface::PlaneGetProjection(r, 1, 1, 1, 0);
    checkEqual(proj, r, check_func, "check for point in plane", passed_checks, results);

    proj = surface::PlaneGetProjection(r, 1, 1, 1, -2);
    checkEqual(proj, {-2.0 / 3.0, -2.0 / 3.0, -2.0 / 3.0}, check_func, "check for point projection through cube's 111", passed_checks, results);

    proj = surface::PlaneGetProjection(r, 1, 0, 0, -1);
    checkEqual(proj, {-1.0, 0.0, 0.0}, check_func, "check for simple point projection in x-direction", passed_checks, results);
  }

  void UnitTest::surface_PlaneGetHKL(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    const string check_func = "surface::PlaneGetHKL";

    xvector<double> hkl = surface::PlaneGetHKL({1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1});
    checkEqual(hkl, {1.0, 1.0, 1.0}, check_func, "check first 111 hkl", passed_checks, results);
    hkl = surface::PlaneGetHKL({0, 1, 1}, {1, 0, 1}, {1, 1, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1});
    checkEqual(hkl, {0.5, 0.5, 0.5}, check_func, "check second 111 hkl", passed_checks, results);
  }

  void UnitTest::surface_PlaneGetVVV(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    const string check_func = "surface::surface_PlaneGetVVV";

    double area;
    xvector<double> v1;
    xvector<double> v2;
    xvector<double> v3;
    xvector<double> v4;
    bool rhombus = surface::PlaneGetVVV({1, 1, 1}, area, v1, v2, v3, v4, {1, 0, 0}, {0, 1, 0}, {0, 0, 1});
    checkEqual(rhombus, false, check_func, "check rhombus for 111", passed_checks, results);
    checkEqual(area, 2 * sqrt(3), check_func, "check area for 111", passed_checks, results);
    checkEqual(v1, xvector<double>{-1, 1, 1}, check_func, "v1 of 111", passed_checks, results);
    checkEqual(v2, xvector<double>{1, -1, 1}, check_func, "v2 of 111", passed_checks, results);
    checkEqual(v3, xvector<double>{1, 1, -1}, check_func, "v3 of 111", passed_checks, results);
    checkEqual(v4, xvector<double>{0, 0, 0}, check_func, "v4 of 111", passed_checks, results);

    rhombus = surface::PlaneGetVVV({1, 0, 0}, area, v1, v2, v3, v4, {1, 0, 0}, {0, 1, 0}, {0, 0, 1});
    checkEqual(rhombus, true, check_func, "check rhombus for 100", passed_checks, results);
    checkEqual(area, 1.0, check_func, "check area for 100", passed_checks, results);
    checkEqual(v1, xvector<double>{1, 0, 0}, check_func, "v1 of 100", passed_checks, results);
    checkEqual(v2, xvector<double>{1, 1, 0}, check_func, "v2 of 100", passed_checks, results);
    checkEqual(v3, xvector<double>{1, 0, 1}, check_func, "v3 of 100", passed_checks, results);
    checkEqual(v4, xvector<double>{1, 1, 1}, check_func, "v4 of 100", passed_checks, results);
  }

  void UnitTest::surface_PlaneGetVVV_V2(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    const string check_func = "surface::surface_PlaneGetVVV_V2";

    double area;
    double area_v2;
    xvector<double> v1;
    xvector<double> v2;
    xvector<double> v3;
    xvector<double> v4;
    xvector<double> v1_v2;
    xvector<double> v2_v2;
    xvector<double> v3_v2;
    xvector<double> v4_v2;
    // hkl, a1, a2, a3
    const vector<std::tuple<xvector<double>, xvector<double>, xvector<double>, xvector<double>>> combos{
        {   {1, 1, 1},   {1, 0, 0},   {0, 1, 0},   {0, 0, 1}},
        {   {1, 1, 1},   {0, 1, 0},   {1, 0, 0},   {0, 0, 1}},
        {   {1, 0, 0},   {1, 0, 0},   {0, 1, 0},   {0, 0, 1}},
        {   {0, 1, 0},   {1, 0, 0},   {0, 1, 0},   {0, 0, 1}},
        {   {0, 0, 1},   {1, 0, 0},   {0, 1, 0},   {0, 0, 1}},
        {   {1, 1, 0},   {1, 0, 0},   {0, 1, 0},   {0, 0, 1}},
        {   {2, 0, 0}, {1.2, 0, 0}, {0, 1.1, 0}, {0, 0, 1.3}},
        {   {1, 0, 2}, {1.2, 0, 0}, {0, 1.1, 0}, {0, 0, 1.3}},
        {   {0, 1, 2}, {1.2, 0, 0}, {0, 1.1, 0}, {0, 0, 1.3}},
        {  {0, -1, 2}, {1.2, 0, 0}, {0, 1.1, 0}, {0, 0, 1.3}},
        {   {1, 0, 0},   {0, 1, 1},   {1, 0, 1},   {1, 1, 0}},
        {{.5, .5, .5},   {0, 1, 1},   {1, 0, 1},   {1, 1, 0}}
    };
    for (const auto& [hkl, a1, a2, a3] : combos) {
      const string desc = string("combo of hkl: ")
                              .append("(")
                              .append(aurostd::xvecDouble2String(hkl))
                              .append(")")
                              .append(" in lattice a1 a2 a3: ")
                              .append(aurostd::xvecDouble2String(a1, 3))
                              .append(" ")
                              .append(aurostd::xvecDouble2String(a2, 3))
                              .append(" ")
                              .append(aurostd::xvecDouble2String(a3, 3));
      const bool rhombus = surface::PlaneGetVVV(hkl, area, v1, v2, v3, v4, a1, a2, a3);
      const bool rhombus_v2 = surface::PlaneGetVVV_V2(hkl, area_v2, v1_v2, v2_v2, v3_v2, v4_v2, a1, a2, a3);
      vector<aurostd::xvector<double>> vertices{v1, v2, v3, v4};
      vector<aurostd::xvector<double>> vertices_v2{v1_v2, v2_v2, v3_v2, v4_v2};
      // the vertices may be in a different order, so sort lexically to make the in the same order so we can compare them
      // this comparator is dirty, TODO xvec needs a dedicated comparator or a proper iterator
      auto comp_xvecdouble = [&](const aurostd::xvector<double>& a, const aurostd::xvector<double>& b) -> bool { return aurostd::xvecDouble2String(a) < aurostd::xvecDouble2String(b); };
      std::sort(vertices.begin(), vertices.end(), comp_xvecdouble);
      std::sort(vertices_v2.begin(), vertices_v2.end(), comp_xvecdouble);
      checkEqual(rhombus_v2, rhombus, check_func, "checking rhombus for " + desc, passed_checks, results);
      checkEqual(area_v2, area, check_func, "checking area for " + desc, passed_checks, results);
      for (size_t i = 0; i < vertices.size() && i < vertices_v2.size(); i++) {
        checkEqual(vertices_v2[i], vertices[i], check_func, "checking vertex " + std::to_string(i) + " for " + desc, passed_checks, results);
      }
    }
  }

  void UnitTest::surface_GetPlaneDensityAtoms(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {}

  void UnitTest::cumulantsTest(uint& passed_checks, vector<vector<string>>& results, vector<string>& errors) {
    const string check_func = "cumulants::CumulantsCalculator";
    const string tmpdir = aurostd::TmpDirectoryCreate("unit_test_cumulants");
    const vector<int> sol_chem = {283, 334, 382, 426, 464, 498, 524, 543, 555, 558, 552, 537, 514, 481, 440, 390, 331, 261, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    const vector<int> sol_chem_dvc = {218, 250, 279, 304, 327, 345, 359, 369, 373, 371, 363, 349, 328, 300, 265, 218, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    const vector<int> sol_coh = {276, 328, 374, 413, 445, 473, 494, 509, 517, 517, 506, 486, 456, 417, 368, 311, 244, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    const vector<int> sol_coh_dvc = {217, 251, 279, 301, 319, 333, 342, 348, 349, 343, 330, 309, 279, 238, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    aurostd::EmbData::save_to_file("cumulants_test.tar.xz", "TEST", tmpdir + "/cumulants_test.tar.xz");
    aurostd::DecompressFiles(tmpdir + "/cumulants_test.tar.xz", tmpdir);
    aurostd::xoption cumulants_opts;
    cumulants_opts.push_attached("CUMULANTS::DIRECTORY", tmpdir + "/AlZn");
    cumulants_opts.push_attached("CUMULANTS::NUM_CUMULANTS", "2");
    cumulants_opts.push_attached("CUMULANTS::POLY_ORDERS", "2,2");
    cumulants_opts.push_attached("CUMULANTS::CONC_START", "0.9,0.1");
    cumulants_opts.push_attached("CUMULANTS::CONC_STOP", "0.1,0.9");
    cumulants_opts.push_attached("CUMULANTS::CONC_STEP", "0.025");
    cumulants_opts.push_attached("CUMULANTS::POCC_TEMP", "2000");
    cumulants_opts.push_attached("CUMULANTS::MIN_TEMP", "200");
    {
      cumulants::CumulantsCalculator cc = cumulants::CumulantsCalculator::fromOptions(cumulants_opts);
      cc.calculateSpinodal();
      aurostd::JSON::object cumulants_jo = cc.getOutput();
      checkEqual(vector<int>(cumulants_jo["crit_temps"]), sol_chem, check_func, "chemical spinodal without DVC", passed_checks, results);
    }
    {
      cumulants_opts.flag("CUMULANTS::DVC", true);
      cumulants::CumulantsCalculator cc = cumulants::CumulantsCalculator::fromOptions(cumulants_opts);
      cc.calculateSpinodal();
      aurostd::JSON::object cumulants_jo = cc.getOutput();
      checkEqual(vector<int>(cumulants_jo["crit_temps"]), sol_chem_dvc, check_func, "chemical spinodal with DVC", passed_checks, results);
    }
    {
      cumulants_opts.flag("CUMULANTS::DVC", false);
      cumulants_opts.push_attached("CUMULANTS::KVEC", "1.0,0.0,0.0");
      cumulants_opts.push_attached("CUMULANTS::POLY_ORDERS_ELAS", "2,4");
      cumulants::CumulantsCalculator cc = cumulants::CumulantsCalculator::fromOptions(cumulants_opts);
      cc.calculateSpinodal();
      aurostd::JSON::object cumulants_jo = cc.getOutput();
      checkEqual(vector<int>(cumulants_jo["crit_temps"]), sol_coh, check_func, "coherent spinodal without DVC", passed_checks, results);
    }
    {
      cumulants_opts.flag("CUMULANTS::DVC", true);
      cumulants::CumulantsCalculator cc = cumulants::CumulantsCalculator::fromOptions(cumulants_opts);
      cc.calculateSpinodal();
      aurostd::JSON::object cumulants_jo = cc.getOutput();
      checkEqual(vector<int>(cumulants_jo["crit_temps"]), sol_coh_dvc, check_func, "coherent spinodal with DVC", passed_checks, results);
    }
    aurostd::RemoveDirectory(tmpdir);
  }

  void UnitTest::xpseudopotentialDataTest(uint& passed_checks, std::vector<std::vector<std::string>>& results, std::vector<std::string>& errors) {
    const std::vector<xPOTCAR> vxpseudopotential = get_pseudopotential_data();

    // collisions for potcars with same values but different AUIDs
    unsigned int collisions = 0;
    for (size_t i = 0; i < vxpseudopotential.size(); i++) {
      for (size_t j = i + 1; j < vxpseudopotential.size(); j++) {
        const bool collision = vxpseudopotential[i].AUID != vxpseudopotential[j].AUID
                            && vxpseudopotential[i].vTITEL.at(0) == vxpseudopotential[j].vTITEL.at(0)
                            && vxpseudopotential[i].vLEXCH.at(0) == vxpseudopotential[j].vLEXCH.at(0)
                            && std::abs(vxpseudopotential[i].vEATOM.at(0) - vxpseudopotential[j].vEATOM.at(0)) < 0.00001
                            && std::abs(vxpseudopotential[i].vRMAX.at(0) - vxpseudopotential[j].vRMAX.at(0)) < 0.0001;
        if (collision) {
          std::cerr
              << "xPOTCAR COLLISION: "
              << "first FILENAME="
              << vxpseudopotential[i].filename
              << " "
              << "second FILENAME="
              << vxpseudopotential[j].filename
              << " "
              << "AUID="
              << vxpseudopotential[i].AUID
              << " "
              << "TITEL="
              << vxpseudopotential[i].vTITEL.at(0)
              << " "
              << "LEXCH="
              << vxpseudopotential[i].vLEXCH.at(0)
              << " "
              << "EATOM="
              << vxpseudopotential[i].vEATOM.at(0)
              << " "
              << "RMAX="
              << vxpseudopotential[i].vRMAX.at(0)
              << " "
              << std::endl;
          collisions++;
        }
      }
    }
    checkEqual(collisions, 0U, "get_pseudopotential_data()", "check value, not auid, collisions in get_pseudopotential_data()", passed_checks, results);

    // collisions for potcars with same AUIDS
    collisions = 0;
    for (size_t i = 0; i < vxpseudopotential.size(); i++) {
      for (size_t j = i + 1; j < vxpseudopotential.size(); j++) {
        const bool collision = vxpseudopotential[i].AUID == vxpseudopotential[j].AUID;
        const bool skip_it = aurostd::substring2bool(vxpseudopotential[i].filename, "pot_GGA/potcar.Apr00/Xe")
                          || aurostd::substring2bool(vxpseudopotential[i].filename, "pot_GGA/potcar.Apr00/Kr")
                          || aurostd::substring2bool(vxpseudopotential[i].filename, "pot_GGA/potcar.Apr00/Li_pv")
                          || aurostd::substring2bool(vxpseudopotential[i].filename, "pot_LDA/potcar.Apr00/Li_pv");
        if (collision && !skip_it) {
          std::cerr
              << "xPOTCAR COLLISION: "
              << "first FILENAME="
              << vxpseudopotential[i].filename
              << " "
              << "second FILENAME="
              << vxpseudopotential[j].filename
              << " "
              << "AUID="
              << vxpseudopotential[i].AUID
              << std::endl;
          collisions++;
        }
      }
    }
    checkEqual(collisions, 0U, "get_pseudopotential_data()", "check only auid collisions in get_pseudopotential_data()", passed_checks, results);

    // collisions for potcars with same values
    collisions = 0;
    for (size_t i = 0; i < vxpseudopotential.size(); i++) {
      for (size_t j = i + 1; j < vxpseudopotential.size(); j++) {
        const bool collision = vxpseudopotential[i].vTITEL.at(0) == vxpseudopotential[j].vTITEL.at(0)
                            && vxpseudopotential[i].vLEXCH.at(0) == vxpseudopotential[j].vLEXCH.at(0)
                            && std::abs(vxpseudopotential[i].vEATOM.at(0) - vxpseudopotential[j].vEATOM.at(0)) < 0.00001
                            && std::abs(vxpseudopotential[i].vRMAX.at(0) - vxpseudopotential[j].vRMAX.at(0)) < 0.0001;
        const bool skip_it = aurostd::substring2bool(vxpseudopotential[i].filename, "pot_GGA/potcar.Apr00/Xe")
                          || aurostd::substring2bool(vxpseudopotential[i].filename, "pot_GGA/potcar.Apr00/Kr")
                          || aurostd::substring2bool(vxpseudopotential[i].filename, "pot_GGA/potcar.Apr00/Li_pv")
                          || aurostd::substring2bool(vxpseudopotential[i].filename, "pot_LDA/potcar.Apr00/Li_pv");
        if (collision && !skip_it) {
          std::cerr
              << "xPOTCAR COLLISION: "
              << "first FILENAME="
              << vxpseudopotential[i].filename
              << " "
              << "second FILENAME="
              << vxpseudopotential[j].filename
              << " "
              << "TITEL="
              << vxpseudopotential[i].vTITEL.at(0)
              << " "
              << "LEXCH="
              << vxpseudopotential[i].vLEXCH.at(0)
              << " "
              << "EATOM="
              << vxpseudopotential[i].vEATOM.at(0)
              << " "
              << "RMAX="
              << vxpseudopotential[i].vRMAX.at(0)
              << " "
              << std::endl;
          collisions++;
        }
      }
    }
    checkEqual(collisions, 0U, "get_pseudopotential_data()", "check only value collisions in get_pseudopotential_data()", passed_checks, results);
  }

} // namespace unittest

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
