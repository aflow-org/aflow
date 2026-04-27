
#ifndef AFLOW_UNIT_TEST_H
#define AFLOW_UNIT_TEST_H

#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xparser_json.h"

#include "flow/aflow_xclasses.h"

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Unit Tests - ME20220127

namespace unittest {

  typedef std::function<void(uint&, std::vector<std::vector<std::string>>&, std::vector<std::string>&)> unitTestFunction;

  struct xcheck {
    std::vector<std::string> errors;
    unitTestFunction func;
    std::string function_name;
    std::string task_description;
    std::string test_group;
    uint passed_checks;
    std::vector<std::vector<std::string>> results;
    bool finished;
  };

  class UnitTest : public xStream {
  public:
    UnitTest(std::ostream& oss = std::cout);
    UnitTest(std::ofstream& mf, std::ostream& oss = std::cout);
    UnitTest(const UnitTest& ut);
    const UnitTest& operator=(const UnitTest& ut);
    ~UnitTest();

    _aflags aflags;

    void clear();

    void resetUnitTest(const std::string& test_name);
    bool runTestSuites(const std::string& unit_test);
    bool runTestSuites(const std::vector<std::string>& unit_tests_in);

  private:
    aurostd::JSON::object utd;
    std::map<std::string, xcheck> test_functions;
    std::map<std::string, std::vector<std::string>> test_groups;
    std::map<std::string, std::string> test2group;
#ifdef AFLOW_MULTITHREADS_ENABLE
    std::mutex mtx;
#endif

    void free();
    void copy(const UnitTest& ut);

    void initialize();
    void initializeTestFunctions();
    void initializeTestGroups();

    template <class UTFunc> void register_test(const UTFunc&& test_function, const std::string& key, const std::string& name, const std::string& desc, const std::string& group);

    void resetUnitTest(xcheck&);
    xcheck initializeXCheck();

    void runUnitTest(std::vector<std::string>::iterator& it, const std::vector<std::string>& tasks);
    bool taskSuccessful(const std::string& task);

    std::string formatResultsTable(const std::vector<std::vector<std::string>>& table);
    void displayResult(const xcheck& xchk);

    template <typename utype>
    void check(bool passed,
               const std::vector<utype>& calculated,
               const std::vector<utype>& expected,
               const std::string& check_function,
               const std::string& checkDescription,
               uint& passed_checks,
               std::vector<std::vector<std::string>>& results);
    void check(bool passed,
               const std::vector<double>& calculated,
               const std::vector<double>& expected,
               const std::string& check_function,
               const std::string& checkDescription,
               uint& passed_checks,
               std::vector<std::vector<std::string>>& results);
    template <typename utype>
    void check(bool passed,
               const aurostd::xmatrix<utype>& calculated,
               const aurostd::xmatrix<utype>& expected,
               const std::string& check_function,
               const std::string& checkDescription,
               uint& passed_checks,
               std::vector<std::vector<std::string>>& results);
    void check(bool passed,
               const aurostd::xmatrix<double>& calculated,
               const aurostd::xmatrix<double>& expected,
               const std::string& check_function,
               const std::string& checkDescription,
               uint& passed_checks,
               std::vector<std::vector<std::string>>& results);

    void checkFiles(const std::string& calculated, const std::string& expected, const std::string& check_function, const std::string& check_description, uint& passed_checks, std::vector<std::vector<std::string>>& results);

    template <typename utype>
    void check(bool passed, const utype& calculated, const utype& expected, const std::string& check_function, const std::string& checkDescription, uint& passed_checks, std::vector<std::vector<std::string>>& results);
    void check(bool passed, const std::string& result_note, const std::string& check_function, const std::string& check_description, uint& passed_checks, std::vector<std::vector<std::string>>& results);

    template <typename utype>
    void checkEqual(const std::vector<utype>& calculated, const std::vector<utype>& expected, const std::string& check_function, const std::string& check_description, uint& passed_checks, std::vector<std::vector<std::string>>& results);
    void checkEqual(const std::vector<std::string>& calculated,
                    const std::vector<std::string>& expected,
                    const std::string& check_function, // AZ20231030 added this function to do std::string comparisions for kpath unit test
                    const std::string& check_description,
                    uint& passed_checks,
                    std::vector<std::vector<std::string>>& results);
    void checkEqual(const std::vector<bool>& calculated, const std::vector<bool>& expected, const std::string& check_function, const std::string& check_description, uint& passed_checks, std::vector<std::vector<std::string>>& results);
    template <typename utype>
    void checkEqual(const utype& calculated, const utype& expected, const std::string& check_function, const std::string& check_description, uint& passed_checks, std::vector<std::vector<std::string>>& results);
    void checkEqual(const std::string& calculated, const std::string& expected, const std::string& check_function, const std::string& check_description, uint& passed_checks, std::vector<std::vector<std::string>>& results);

    // Test functions ---------------------------------

      // aurostd
    void xscalarTest(uint&, std::vector<std::vector<std::string>>&, std::vector<std::string>&);
    void xvectorTest(uint&, std::vector<std::vector<std::string>>&, std::vector<std::string>&);
    void xmatrixTest(uint&, std::vector<std::vector<std::string>>&, std::vector<std::string>&);
    void xfitTest(uint&, std::vector<std::vector<std::string>>&, std::vector<std::string>&);
    void aurostdMainTest(uint&, std::vector<std::vector<std::string>>&, std::vector<std::string>&);
    void xfileTest(uint&, std::vector<std::vector<std::string>>&, std::vector<std::string>&);
    void xparserTest(uint&, std::vector<std::vector<std::string>>&, std::vector<std::string>&);

      // json
    template <class T, aurostd::JSON::enable_serializable<T> = true, class F>
    void json_cycle_test(F initializer, uint& passed_checks, std::vector<std::vector<std::string>>& results, std::vector<std::string>& errors);
    void jsonSerialization(uint& passed_checks, std::vector<std::vector<std::string>>& results, std::vector<std::string>& errors);
    void jsonMockSerialization(uint& passed_checks, std::vector<std::vector<std::string>>& results, std::vector<std::string>& errors);
    void xKpointJsonSerialization(uint& passed_checks, std::vector<std::vector<std::string>>& results, std::vector<std::string>& errors);
    void xPOTCARJsonSerialization(uint& passed_checks, std::vector<std::vector<std::string>>& results, std::vector<std::string>& errors);
    void xVASPRUNXMLJsonSerialization(uint& passed_checks, std::vector<std::vector<std::string>>& results, std::vector<std::string>& errors);
    void xDOSCARJsonSerialization(uint& passed_checks, std::vector<std::vector<std::string>>& results, std::vector<std::string>& errors);
    void xIBZKPTJsonSerialization(uint& passed_checks, std::vector<std::vector<std::string>>& results, std::vector<std::string>& errors);
    void xQMVASPJsonSerialization(uint& passed_checks, std::vector<std::vector<std::string>>& results, std::vector<std::string>& errors);
    void xPLASMONICSJsonSerialization(uint& passed_checks, std::vector<std::vector<std::string>>& results, std::vector<std::string>& errors);
    void xEIGENVALJsonSerialization(uint& passed_checks, std::vector<std::vector<std::string>>& results, std::vector<std::string>& errors);
    void xOUTCARJsonSerialization(uint& passed_checks, std::vector<std::vector<std::string>>& results, std::vector<std::string>& errors);

      // ovasp
    void ovaspTest(uint& passed_checks, std::vector<std::vector<std::string>>& results, std::vector<std::string>& errors);

    //dielectric testing
    void dielectricTest(uint& passed_checks, std::vector<std::vector<std::string>>& results, std::vector<std::string>& errors);

      // database
    void schemaTest(uint&, std::vector<std::vector<std::string>>&, std::vector<std::string>&);

      // xstructure
    void atomicEnvironmentTest(uint&, std::vector<std::vector<std::string>>&, std::vector<std::string>&);
    void xstructureParserTest(uint&, std::vector<std::vector<std::string>>&, std::vector<std::string>&);
    void xstructureTest(uint&, std::vector<std::vector<std::string>>&, std::vector<std::string>&);

      // kpath test, tests that all kpaths are consistent //AZ20231027
      /// @brief a unit test to make sure the kpaths are consistent
    void kpathTest(uint&, std::vector<std::vector<std::string>>&, std::vector<std::string>&); // AZ20231027

      // structure generation
    void ceramgenTest(uint&, std::vector<std::vector<std::string>>&, std::vector<std::string>&);
    void prototypeGeneratorTest(uint&, std::vector<std::vector<std::string>>&, std::vector<std::string>&);

      // ovasp
    void xoutcarTest(uint&, std::vector<std::vector<std::string>>&, std::vector<std::string>&);

    void lib2rawTest(uint&, std::vector<std::vector<std::string>>&, std::vector<std::string>&);

      // surface
    void surface_TriangleArea(uint& passed_checks, std::vector<std::vector<std::string>>& results, std::vector<std::string>& errors);
    void surface_PlaneGetABCD(uint& passed_checks, std::vector<std::vector<std::string>>& results, std::vector<std::string>& errors);
    void surface_PlaneDistance(uint& passed_checks, std::vector<std::vector<std::string>>& results, std::vector<std::string>& errors);
    void surface_PlaneGetProjection(uint& passed_checks, std::vector<std::vector<std::string>>& results, std::vector<std::string>& errors);
    void surface_PlaneGetHKL(uint& passed_checks, std::vector<std::vector<std::string>>& results, std::vector<std::string>& errors);
    void surface_PlaneGetVVV(uint& passed_checks, std::vector<std::vector<std::string>>& results, std::vector<std::string>& errors);
    void surface_PlaneGetVVV_V2(uint& passed_checks, std::vector<std::vector<std::string>>& results, std::vector<std::string>& errors);
    void surface_GetPlaneDensityAtoms(uint& passed_checks, std::vector<std::vector<std::string>>& results, std::vector<std::string>& errors);

      // entryLoader
    void entryLoaderTest(uint&, std::vector<std::vector<std::string>>&, std::vector<std::string>&);

    //qhull
    void qhullTest(uint&, std::vector<std::vector<std::string>>&, std::vector<std::string>&);
      // cumulants
    void cumulantsTest(uint&, std::vector<std::vector<std::string>>&, std::vector<std::string>&);

    // data
    void xpseudopotentialDataTest(uint&, std::vector<std::vector<std::string>>&, std::vector<std::string>&);
  };

  class RedirectStream {
  public:
    explicit RedirectStream(std::streambuf*, std::ostream* = &std::cout);
    ~RedirectStream();

    void capture(std::streambuf* buffer, std::ostream* captured_stream);
    void release();

    template <typename Function> static std::stringstream run_with_capture_cout(Function&& func);

  private:
    std::ostream* original_stream;
    std::streambuf* original_buffer;
    std::streambuf* capture_buffer;
    bool capturing;
  };

} // namespace unittest

#endif // AFLOW_UNIT_TEST_H
