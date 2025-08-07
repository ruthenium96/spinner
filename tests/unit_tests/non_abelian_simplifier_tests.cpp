#include <cmath>

#include "gtest/gtest.h"
#include "src/common/Logger.h"
#include "src/space/optimization/OptimizedSpaceConstructor.h"


const auto factories = quantum::linear_algebra::FactoriesList();

void compare_two_spaces(const space::Space& one, const space::Space& two) {
   ASSERT_EQ(one.getBlocks().size(), two.getBlocks().size())
       << "Numbers of subspaces do not equal";
   for (size_t i = 0; i < one.getBlocks().size(); ++i) {
       uint32_t one_dimensionality = one.getBlocks()[i].properties.dimensionality;
       uint32_t two_dimensionality = two.getBlocks()[i].properties.dimensionality;
       uint32_t one_degeneracy = one.getBlocks()[i].properties.degeneracy;
       uint32_t two_degeneracy = two.getBlocks()[i].properties.degeneracy;

       ASSERT_EQ(
           one.getBlocks()[i].decomposition->size_cols() * one_degeneracy,
           two.getBlocks()[i].decomposition->size_cols() * two_degeneracy)
           << "Size of subspace " << i << " does not equal";

       EXPECT_EQ(one_dimensionality * one_degeneracy, two_dimensionality * two_degeneracy)
           << "First space :" << one.getBlocks()[i].properties
           << "Second space :" << two.getBlocks()[i].properties;

       std::unordered_map<uint32_t, uint32_t> one_met;
       std::unordered_map<uint32_t, uint32_t> two_met;

       // counting:
       for (uint32_t index_of_vector = 0;
            index_of_vector < one.getBlocks()[i].decomposition->size_cols();
            ++index_of_vector) {
           auto iterator = one.getBlocks()[i].decomposition->GetNewIterator(index_of_vector);
           while (iterator->hasNext()) {
               auto item = iterator->getNext();
               one_met[item.index] += 1;
           }
       }
       for (uint32_t index_of_vector = 0;
            index_of_vector < two.getBlocks()[i].decomposition->size_cols();
            ++index_of_vector) {
           auto iterator = two.getBlocks()[i].decomposition->GetNewIterator(index_of_vector);
           while (iterator->hasNext()) {
               auto item = iterator->getNext();
               two_met[item.index] += 1;
           }
       }

       // ceiling to closest n * dimensionality, because some projectors has some coefficients equal to 0
       for (auto& p : one_met) {
           p.second = (int)ceil(((double)p.second) / one_dimensionality) * one_dimensionality;
       }
       for (auto& p : two_met) {
           p.second = (int)ceil(((double)p.second) / two_dimensionality) * two_dimensionality;
       }

       // comparing spaces
       for (const auto& p : one_met) {
           EXPECT_EQ(p.second * one_degeneracy, two_met[p.first] * two_degeneracy)
               << "Lex index [" << p.first
               << "] is met different times in spaces, representation: "
               << one.getBlocks()[i].properties.get_representation_name();
       }
   }
}

TEST(nonAbelianSimplifier, 333_Dih3) {
    std::vector<spin_algebra::Multiplicity> mults = {3, 3, 3};
    model::ModelInput model(mults);

    common::physical_optimization::OptimizationList optimizationListFull;
    optimizationListFull.Symmetrize({group::Group::Dihedral, 3}, {{1, 2, 0}, {0, 2, 1}});
    common::physical_optimization::OptimizationList optimizationListSimplified;
    optimizationListSimplified.Symmetrize({group::Group::Dihedral, 3}, {{1, 2, 0}, {0, 2, 1}});
    optimizationListSimplified.NonAbelianSimplify();

    space::Space space_full = space::optimization::OptimizedSpaceConstructor::construct(
        {model, optimizationListFull},
        factories);
    space::Space space_simplified = space::optimization::OptimizedSpaceConstructor::construct(
        {model, optimizationListFull},
        factories);

    compare_two_spaces(space_full, space_simplified);
}

TEST(nonAbelianSimplifier, 333_doubleDih3) {
   std::vector<spin_algebra::Multiplicity> mults = {3, 3, 3};
   model::ModelInput model(mults);

   common::physical_optimization::OptimizationList optimizationListFull;
   optimizationListFull.Symmetrize({group::Group::Dihedral, 3}, {{1, 2, 0}, {0, 2, 1}});
   optimizationListFull.Symmetrize({group::Group::Dihedral, 3}, {{1, 2, 0}, {0, 2, 1}});
   common::physical_optimization::OptimizationList optimizationListSimplified;
   optimizationListSimplified.Symmetrize({group::Group::Dihedral, 3}, {{1, 2, 0}, {0, 2, 1}});
   optimizationListSimplified.Symmetrize({group::Group::Dihedral, 3}, {{1, 2, 0}, {0, 2, 1}});
   optimizationListSimplified.NonAbelianSimplify();

   space::Space space_full = space::optimization::OptimizedSpaceConstructor::construct(
       {model, optimizationListFull},
       factories);
   space::Space space_simplified = space::optimization::OptimizedSpaceConstructor::construct(
       {model, optimizationListFull},
       factories);

   compare_two_spaces(space_full, space_simplified);
}

TEST(nonAbelianSimplifier, 333_doubleDih3_tricky) {
   std::vector<spin_algebra::Multiplicity> mults = {3, 3, 3};
   model::ModelInput model(mults);

   common::physical_optimization::OptimizationList optimizationListFull;
   optimizationListFull.Symmetrize({group::Group::Dihedral, 3}, {{1, 2, 0}, {0, 2, 1}});
   optimizationListFull.Symmetrize({group::Group::Dihedral, 3}, {{2, 0, 1}, {1, 0, 2}});
   common::physical_optimization::OptimizationList optimizationListSimplified;
   optimizationListSimplified.Symmetrize({group::Group::Dihedral, 3}, {{1, 2, 0}, {0, 2, 1}});
   optimizationListSimplified.Symmetrize({group::Group::Dihedral, 3}, {{2, 0, 1}, {1, 0, 2}});
   optimizationListSimplified.NonAbelianSimplify();

   space::Space space_full = space::optimization::OptimizedSpaceConstructor::construct(
       {model, optimizationListFull},
       factories);
   space::Space space_simplified = space::optimization::OptimizedSpaceConstructor::construct(
       {model, optimizationListFull},
       factories);

   compare_two_spaces(space_full, space_simplified);
}

TEST(nonAbelianSimplifier, 333333_Dih3) {
   std::vector<spin_algebra::Multiplicity> mults = {3, 3, 3, 3, 3, 3};
   model::ModelInput model(mults);

   common::physical_optimization::OptimizationList optimizationListFull;
   optimizationListFull.Symmetrize({group::Group::Dihedral, 3}, {{1, 2, 0, 4, 5, 3}, {0, 2, 1, 3, 5, 4}});
   common::physical_optimization::OptimizationList optimizationListSimplified;
   optimizationListSimplified.Symmetrize({group::Group::Dihedral, 3}, {{1, 2, 0, 4, 5, 3}, {0, 2, 1, 3, 5, 4}});
   optimizationListSimplified.NonAbelianSimplify();

   space::Space space_full = space::optimization::OptimizedSpaceConstructor::construct(
       {model, optimizationListFull},
       factories);
   space::Space space_simplified = space::optimization::OptimizedSpaceConstructor::construct(
       {model, optimizationListFull},
       factories);

   compare_two_spaces(space_full, space_simplified);
}

TEST(nonAbelianSimplifier, 333333_Dih3xS2) {
   std::vector<spin_algebra::Multiplicity> mults = {3, 3, 3, 3, 3, 3};
   model::ModelInput model(mults);

   common::physical_optimization::OptimizationList optimizationListFull;
   optimizationListFull.Symmetrize({group::Group::Dihedral, 3}, {{1, 2, 0, 4, 5, 3}, {0, 2, 1, 3, 5, 4}});
   optimizationListFull.Symmetrize(group::Group::S2, {{3, 4, 5, 0, 1, 2}});
   common::physical_optimization::OptimizationList optimizationListSimplified;
   optimizationListSimplified.Symmetrize({group::Group::Dihedral, 3}, {{1, 2, 0, 4, 5, 3}, {0, 2, 1, 3, 5, 4}});
   optimizationListSimplified.Symmetrize(group::Group::S2, {{3, 4, 5, 0, 1, 2}});
   optimizationListSimplified.NonAbelianSimplify();

   space::Space space_full = space::optimization::OptimizedSpaceConstructor::construct(
       {model, optimizationListFull},
       factories);
   space::Space space_simplified = space::optimization::OptimizedSpaceConstructor::construct(
       {model, optimizationListFull},
       factories);

   compare_two_spaces(space_full, space_simplified);
}

TEST(nonAbelianSimplifier, 222222222) {
   std::vector<spin_algebra::Multiplicity> mults = {2, 2, 2, 2, 2, 2, 2, 2, 2};
   model::ModelInput model(mults);

   common::physical_optimization::OptimizationList optimizationListFull;
   optimizationListFull.Symmetrize({group::Group::Dihedral, 3}, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
   optimizationListFull.Symmetrize({group::Group::Dihedral, 3}, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});
   common::physical_optimization::OptimizationList optimizationListSimplified;
   optimizationListSimplified.Symmetrize({group::Group::Dihedral, 3}, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
   optimizationListSimplified.Symmetrize({group::Group::Dihedral, 3}, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});
   EXPECT_THROW(optimizationListSimplified.NonAbelianSimplify(), std::invalid_argument);
}

TEST(nonAbelianSimplifier, 333333333) {
   std::vector<spin_algebra::Multiplicity> mults = {3, 3, 3, 3, 3, 3, 3, 3, 3};
   model::ModelInput model(mults);

   common::physical_optimization::OptimizationList optimizationListFull;
   optimizationListFull.Symmetrize({group::Group::Dihedral, 3}, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
   optimizationListFull.Symmetrize({group::Group::Dihedral, 3}, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});
   common::physical_optimization::OptimizationList optimizationListSimplified;
   optimizationListSimplified.Symmetrize({group::Group::Dihedral, 3}, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
   optimizationListSimplified.Symmetrize({group::Group::Dihedral, 3}, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});
   EXPECT_THROW(optimizationListSimplified.NonAbelianSimplify(), std::invalid_argument);
}