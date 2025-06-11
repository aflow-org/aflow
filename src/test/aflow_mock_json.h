
#ifndef AFLOW_TEST_MOCK_JSON_H
#define AFLOW_TEST_MOCK_JSON_H

#include <deque>
#include <list>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "AUROSTD/aurostd_automatic_template.h"
#include "AUROSTD/aurostd_xserialization.h"

namespace unittest {

/// Mock for simple class with simple types.
  struct MockSerializableMember final : JsonSerializable<MockSerializableMember> {
    int member_x = 5;
    int member_y = 7;

#define JSON_MockSerializableMember_MEMBERS member_x, member_y
  protected:
    MockSerializableMember deserialize(const aurostd::JSON::object &jo) override {
      AST_JSON_SETTER(JSON_MockSerializableMember_MEMBERS)
      return *this;
    }
    [[nodiscard]] aurostd::JSON::object serialize() const override { return aurostd::JSON::object({AST_JSON_GETTER(JSON_MockSerializableMember_MEMBERS)}); }
    [[nodiscard]] std::string getJsonID() const override { return "MockSerializableMember"; }

  public:
    ~MockSerializableMember() override = default;
  };

/// Mock for class with list-like and dict-like members
  struct MockSerializableListDict final : JsonSerializable<MockSerializableListDict> {
    std::list<int> member_list{1, 2, 3};
    std::vector<int> member_vec{4, 5, 6};
    std::deque<int> member_deque{7, 8, 9};

    std::map<std::string, int> member_map{
        {"A", 1},
        {"B", 2},
        {"C", 3}
    };
    std::unordered_map<std::string, int> member_umap{
        {"a", 1},
        {"b", 2},
        {"c", 3}
    };
    std::multimap<std::string, int> member_multimap{
        {"1", 1},
        {"2", 2}
    };

#define JSON_MockSerializableListDict_MEMBERS member_list, member_vec, member_deque, member_map, member_umap, member_multimap
  protected:
    MockSerializableListDict deserialize(const aurostd::JSON::object &jo) override {
      AST_JSON_SETTER(JSON_MockSerializableListDict_MEMBERS)
      return *this;
    }
    [[nodiscard]] aurostd::JSON::object serialize() const override { return aurostd::JSON::object({AST_JSON_GETTER(JSON_MockSerializableListDict_MEMBERS)}); }
    [[nodiscard]] std::string getJsonID() const override { return "MockSerializableListDict"; }

  public:
    ~MockSerializableListDict() override = default;
  };

/// Mock for class with a JsonSerializable member.
  struct MockSerializableMain final : JsonSerializable<MockSerializableMain> {
    std::vector<std::vector<std::vector<int>>> member_3Dvec = {{{3}}};
    MockSerializableMember member_mock = {};
    MockSerializableListDict member_mock_listdict = {};

#define JSON_MockSerializableMain_MEMBERS member_3Dvec, member_mock, member_mock_listdict
  protected:
    MockSerializableMain deserialize(const aurostd::JSON::object &jo) override {
      AST_JSON_SETTER(JSON_MockSerializableMain_MEMBERS)
      return *this;
    }
    [[nodiscard]] aurostd::JSON::object serialize() const override { return aurostd::JSON::object({AST_JSON_GETTER(JSON_MockSerializableMain_MEMBERS)}); }
    [[nodiscard]] std::string getJsonID() const override { return "MockSerializableMain"; }

  public:
    ~MockSerializableMain() override = default;
  };

} // namespace unittest

#endif // AFLOW_TEST_MOCK_JSON_H
