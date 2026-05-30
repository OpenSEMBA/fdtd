#ifndef TEST_IDCHILDTABLE_H
#define TEST_IDCHILDTABLE_H

#include <gtest/gtest.h>
#include <cstring>
#include <fstream>
#include <nlohmann/json.hpp>
#include <unistd.h>
#include "id_map_m.h"

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

static nlohmann::json load_json_file(const std::string& path) {
    std::ifstream input(path);
    nlohmann::json root;
    input >> root;
    return root;
}

TEST(id_map, basic_construction) {
    std::string json = R"({
        "materials": [
            {"id": 1, "type": "wire"},
            {"id": 2, "type": "terminal"}
        ]
    })";
    std::string path = create_temp_json(json);

    nlohmann::json root = load_json_file(path);

    id_map_m::id_map_t table = id_map_m::buildIdMap(root, "materials");

    EXPECT_EQ(table.size(), 2U);
    EXPECT_TRUE(id_map_m::containsId(table, 1));
    EXPECT_TRUE(id_map_m::containsId(table, 2));
    EXPECT_FALSE(id_map_m::containsId(table, 999));

    auto val1 = id_map_m::findById(table, 1);
    EXPECT_NE(val1, nullptr);
    auto val2 = id_map_m::findById(table, 2);
    EXPECT_NE(val2, nullptr);

    cleanup_temp(path);
}

TEST(id_map, empty_path) {
    std::string json = R"({
        "other": [1, 2, 3]
    })";
    std::string path = create_temp_json(json);

    nlohmann::json root = load_json_file(path);

    id_map_m::id_map_t table = id_map_m::buildIdMap(root, "nonexistent");

    EXPECT_EQ(table.size(), 0U);

    cleanup_temp(path);
}

TEST(id_map, single_entry) {
    std::string json = R"({
        "items": [
            {"id": 42, "name": "test"}
        ]
    })";
    std::string path = create_temp_json(json);

    nlohmann::json root = load_json_file(path);

    id_map_m::id_map_t table = id_map_m::buildIdMap(root, "items");

    EXPECT_EQ(table.size(), 1U);
    EXPECT_TRUE(id_map_m::containsId(table, 42));
    EXPECT_FALSE(id_map_m::containsId(table, 41));
    EXPECT_FALSE(id_map_m::containsId(table, 43));

    auto val = id_map_m::findById(table, 42);
    EXPECT_NE(val, nullptr);

    cleanup_temp(path);
}

TEST(id_map, many_entries) {
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

    nlohmann::json root = load_json_file(path);

    id_map_m::id_map_t table = id_map_m::buildIdMap(root, "data");

    EXPECT_EQ(table.size(), 5U);
    for (int i = 1; i <= 5; i++) {
        EXPECT_TRUE(id_map_m::containsId(table, i));
        auto val = id_map_m::findById(table, i);
        EXPECT_NE(val, nullptr);
    }
    EXPECT_FALSE(id_map_m::containsId(table, 0));
    EXPECT_FALSE(id_map_m::containsId(table, 6));

    cleanup_temp(path);
}

TEST(id_map, find_by_id_null_on_missing_key) {
    std::string json = R"({
        "items": [
            {"id": 10, "type": "foo"}
        ]
    })";
    std::string path = create_temp_json(json);

    nlohmann::json root = load_json_file(path);

    id_map_m::id_map_t table = id_map_m::buildIdMap(root, "items");

    auto val10 = id_map_m::findById(table, 10);
    EXPECT_NE(val10, nullptr);
    auto val99 = id_map_m::findById(table, 99);
    EXPECT_EQ(val99, nullptr);

    cleanup_temp(path);
}

TEST(id_map, duplicate_ids_keep_last_entry) {
    std::string json = R"({
        "items": [
            {"id": 1, "type": "wire"},
            {"id": 1, "type": "terminal"}
        ]
    })";
    std::string path = create_temp_json(json);

    nlohmann::json root = load_json_file(path);
    id_map_m::id_map_t table = id_map_m::buildIdMap(root, "items");

    EXPECT_EQ(table.size(), 1U);
    auto val = id_map_m::findById(table, 1);
    ASSERT_NE(val, nullptr);
    ASSERT_TRUE(val->contains("type"));
    EXPECT_EQ((*val)["type"].get<std::string>(), "terminal");

    cleanup_temp(path);
}

#endif // TEST_IDCHILDTABLE_H
