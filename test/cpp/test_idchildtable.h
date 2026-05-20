#ifndef TEST_IDCHILDTABLE_H
#define TEST_IDCHILDTABLE_H

#include <gtest/gtest.h>
#include <fstream>
#include <nlohmann/json.hpp>
#include <unistd.h>
#include "json_module.h"
#include "idchildtable_m.h"

static std::string create_temp_json(const std::string& content) {
    std::string path = "/tmp/test_idchildtable_XXXXXX.json";
    char buf[256];
    strncpy(buf, path.c_str(), sizeof(buf));
    int fd = mkstemp(buf);
    close(fd);
    std::ofstream f(buf);
    f << content;
    f.close();
    return std::string(buf);
}

static void cleanup_temp(const std::string& path) {
    std::remove(path.c_str());
}

// Test IdChildTable_t construction and basic operations
TEST(idchildtable, basic_construction) {
    std::string json = R"({
        "materials": [
            {"id": 1, "type": "wire"},
            {"id": 2, "type": "terminal"}
        ]
    })";
    std::string path = create_temp_json(json);

    json_module::json_file jsonfile;
    jsonfile.initialize();
    jsonfile.load(path);
    json_module::json_core core;
    jsonfile.get_core(core);

    json_module::json_value root;
    jsonfile.get(".", root);

    idchildtable_m::IdChildTable_t table =
        idchildtable_m::IdChildTable_t::ctor(core, root, "materials");

    EXPECT_EQ(table.totalSize(), 2);
    EXPECT_EQ(table.checkId(1), 0);
    EXPECT_EQ(table.checkId(2), 0);
    EXPECT_EQ(table.checkId(999), 1);

    auto val1 = table.getId(1);
    EXPECT_NE(val1, nullptr);
    auto val2 = table.getId(2);
    EXPECT_NE(val2, nullptr);

    cleanup_temp(path);
}

// Test IdChildTable_t with non-existent path
TEST(idchildtable, empty_path) {
    std::string json = R"({
        "other": [1, 2, 3]
    })";
    std::string path = create_temp_json(json);

    json_module::json_file jsonfile;
    jsonfile.initialize();
    jsonfile.load(path);
    json_module::json_core core;
    jsonfile.get_core(core);

    json_module::json_value root;
    jsonfile.get(".", root);

    idchildtable_m::IdChildTable_t table =
        idchildtable_m::IdChildTable_t::ctor(core, root, "nonexistent");

    EXPECT_EQ(table.totalSize(), 0);

    cleanup_temp(path);
}

// Test IdChildTable_t with single entry
TEST(idchildtable, single_entry) {
    std::string json = R"({
        "items": [
            {"id": 42, "name": "test"}
        ]
    })";
    std::string path = create_temp_json(json);

    json_module::json_file jsonfile;
    jsonfile.initialize();
    jsonfile.load(path);
    json_module::json_core core;
    jsonfile.get_core(core);

    json_module::json_value root;
    jsonfile.get(".", root);

    idchildtable_m::IdChildTable_t table =
        idchildtable_m::IdChildTable_t::ctor(core, root, "items");

    EXPECT_EQ(table.totalSize(), 1);
    EXPECT_EQ(table.checkId(42), 0);
    EXPECT_EQ(table.checkId(41), 1);
    EXPECT_EQ(table.checkId(43), 1);

    auto val = table.getId(42);
    EXPECT_NE(val, nullptr);

    cleanup_temp(path);
}

// Test IdChildTable_t with many entries
TEST(idchildtable, many_entries) {
    std::string json = R"({
        "data": [
            {"id": 1, "val": "a"},
            {"id": 2, "val": "b"},
            {"id": 3, "val": "c"},
            {"id": 4, "val": "d"},
            {"id": 5, "val": "e"}
        ]
    })";
    std::string path = create_temp_json(json);

    json_module::json_file jsonfile;
    jsonfile.initialize();
    jsonfile.load(path);
    json_module::json_core core;
    jsonfile.get_core(core);

    json_module::json_value root;
    jsonfile.get(".", root);

    idchildtable_m::IdChildTable_t table =
        idchildtable_m::IdChildTable_t::ctor(core, root, "data");

    EXPECT_EQ(table.totalSize(), 5);
    for (int i = 1; i <= 5; i++) {
        EXPECT_EQ(table.checkId(i), 0);
        auto val = table.getId(i);
        EXPECT_NE(val, nullptr);
    }
    EXPECT_EQ(table.checkId(0), 1);
    EXPECT_EQ(table.checkId(6), 1);

    cleanup_temp(path);
}

// Test getId returns non-null for existing IDs and null for non-existing
TEST(idchildtable, get_id_null) {
    std::string json = R"({
        "items": [
            {"id": 10, "type": "foo"}
        ]
    })";
    std::string path = create_temp_json(json);

    json_module::json_file jsonfile;
    jsonfile.initialize();
    jsonfile.load(path);
    json_module::json_core core;
    jsonfile.get_core(core);

    json_module::json_value root;
    jsonfile.get(".", root);

    idchildtable_m::IdChildTable_t table =
        idchildtable_m::IdChildTable_t::ctor(core, root, "items");

    auto val10 = table.getId(10);
    EXPECT_NE(val10, nullptr);
    auto val99 = table.getId(99);
    EXPECT_EQ(val99, nullptr);

    cleanup_temp(path);
}

#endif // TEST_IDCHILDTABLE_H
