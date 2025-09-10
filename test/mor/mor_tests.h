#include <gtest/gtest.h>

extern "C" int test_evolution_operator_dimension_Field_basis();
extern "C" int test_evolution_operator_poisition_E_basis();
extern "C" int test_evolution_operator_position_H_basis();
extern "C" int test_evolution_operator_E_indices_map();
extern "C" int test_evolution_operator_H_indices_map();
extern "C" int test_evolution_operator_indices_map_all_fields();
extern "C" int test_evolution_operator_comparison_with_solver();

TEST(mor, evolutionOperator_BasisDimension)        { EXPECT_EQ(0, test_evolution_operator_dimension_Field_basis()); }
TEST(mor, evolutionOperator_PositionEBasis)        { EXPECT_EQ(0, test_evolution_operator_poisition_E_basis()); }
TEST(mor, evolutionOperator_PositionHBasis)        { EXPECT_EQ(0, test_evolution_operator_position_H_basis()); }
TEST(mor, evolutionOperator_EIndicesMap)           { EXPECT_EQ(0, test_evolution_operator_E_indices_map()); }
TEST(mor, evolutionOperator_HIndicesMap)           { EXPECT_EQ(0, test_evolution_operator_H_indices_map()); }
TEST(mor, evolutionOperator_IndicesMapAllFields)   { EXPECT_EQ(0, test_evolution_operator_indices_map_all_fields()); }
TEST(mor, evolutionOperator_ComparisonWithSolver)  { EXPECT_EQ(0, test_evolution_operator_comparison_with_solver()); }
