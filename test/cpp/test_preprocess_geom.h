#ifndef TEST_PREPROCESS_GEOM_H
#define TEST_PREPROCESS_GEOM_H

#include <gtest/gtest.h>

#include "preprocess_tags.h"
#include "nfde_types.h"

using namespace Preprocess_m;
using namespace NFDETypes_m;

namespace {

tagtype_t make_tagtype_with_capacity(int capacity) {
    tagtype_t tagtype;
    tagtype.numertags = 0;
    tagtype.tag.resize(static_cast<size_t>(capacity));
    return tagtype;
}

} // namespace

TEST(preprocess, searchtag) {
    tagtype_t test_tags;
    test_tags.numertags = 5;
    test_tags.tag = {
        "tag_alpha",
        "tag_beta",
        "tag_gamma",
        "tag_delta",
        "tag_epsilon",
    };

    EXPECT_EQ(searchtag(test_tags, "nonexistent"), -1);
    EXPECT_EQ(searchtag(test_tags, "tag_alpha"), 1);
    EXPECT_EQ(searchtag(test_tags, "tag_gamma"), 3);
    EXPECT_EQ(searchtag(test_tags, "tag_epsilon"), 5);
    EXPECT_EQ(searchtag(test_tags, "  tag_beta  "), 2);
    EXPECT_EQ(searchtag(test_tags, ""), -1);
    EXPECT_EQ(searchtag(test_tags, "TAG_ALPHA"), -1);
}

TEST(preprocess, searchtag_empty) {
    tagtype_t empty_tags;
    empty_tags.numertags = 0;
    EXPECT_EQ(searchtag(empty_tags, "any_tag"), -1);
}

TEST(preprocess, searchtag_single) {
    tagtype_t single_tags;
    single_tags.numertags = 1;
    single_tags.tag = {"only_tag"};

    EXPECT_EQ(searchtag(single_tags, "only_tag"), 1);
    EXPECT_EQ(searchtag(single_tags, "other_tag"), -1);
}

TEST(preprocess, searchtag_special_chars) {
    tagtype_t special_tags;
    special_tags.numertags = 4;
    special_tags.tag = {
        "tag_with_underscores",
        "tag-with-dashes",
        "tag.with.dots",
        "Tag123",
    };

    EXPECT_EQ(searchtag(special_tags, "tag_with_underscores"), 1);
    EXPECT_EQ(searchtag(special_tags, "tag-with-dashes"), 2);
    EXPECT_EQ(searchtag(special_tags, "tag.with.dots"), 3);
    EXPECT_EQ(searchtag(special_tags, "Tag123"), 4);
}

TEST(preprocess, checkDielectricTag_no_dup) {
    Dielectric_t diel_comp;
    diel_comp.n_C1P = 2;
    diel_comp.c1P.resize(2);
    diel_comp.c1P[0].tag = "dielectric_tag1";
    diel_comp.c1P[1].tag = "dielectric_tag2";
    diel_comp.n_C2P = 0;

    std::vector<Dielectric_t> prev_diel(1);
    prev_diel[0].n_C1P = 0;
    prev_diel[0].n_C2P = 0;

    tagtype_t tagtype = make_tagtype_with_capacity(10);
    const std::string error_msg = "Error in dielectric tag check";

    int32_t numertag = 0;
    numertag += 1;
    checkDielectricTagForDuplicate(diel_comp, prev_diel, 0, 1, "c1P", numertag, tagtype, 1, error_msg);
    EXPECT_EQ(numertag, 1);

    numertag += 1;
    checkDielectricTagForDuplicate(diel_comp, prev_diel, 0, 2, "c1P", numertag, tagtype, 1, error_msg);
    EXPECT_EQ(numertag, 2);
}

TEST(preprocess, checkLossyTag_basic) {
    LossyThinSurface_t lossy_comp;
    lossy_comp.nc = 2;
    lossy_comp.c.resize(2);
    lossy_comp.c[0].tag = "lossy_tag1";
    lossy_comp.c[1].tag = "lossy_tag2";

    std::vector<LossyThinSurface_t> prev_lossy(1);
    prev_lossy[0].nc = 0;

    tagtype_t tagtype = make_tagtype_with_capacity(10);

    int32_t numertag = 0;
    numertag += 1;
    checkLossyTagForDuplicate(lossy_comp, prev_lossy, 0, 1, numertag, tagtype, 1);
    EXPECT_EQ(numertag, 1);

    numertag += 1;
    checkLossyTagForDuplicate(lossy_comp, prev_lossy, 0, 2, numertag, tagtype, 1);
    EXPECT_EQ(numertag, 2);
}

#endif
