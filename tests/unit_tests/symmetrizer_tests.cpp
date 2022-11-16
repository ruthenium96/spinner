#include "gtest/gtest.h"
#include "src/common/Logger.h"
#include "src/space/optimization/OptimizedSpaceConstructor.h"

size_t number_of_vectors(const space::Space& space) {
    size_t acc = 0;
    for (const auto& subspace : space.getBlocks()) {
        acc += subspace.decomposition->size();
    }
    return acc;
}

bool orthogonality_of_basis(const space::Space& space) {
    std::unique_ptr<quantum::linear_algebra::AbstractSparseMatrix> unitary_matrix =
        quantum::linear_algebra::AbstractSparseMatrix::defaultSparseMatrix();
    for (const auto& subspace : space.getBlocks()) {
        unitary_matrix->copy_all_from(subspace.decomposition);
    }
    bool answer = true;
#pragma omp parallel for shared(space, unitary_matrix, answer) default(none)
    for (size_t index_of_vector_i = 0; index_of_vector_i < unitary_matrix->size();
         ++index_of_vector_i) {
        for (size_t index_of_vector_j = index_of_vector_i + 1;
             index_of_vector_j < unitary_matrix->size();
             ++index_of_vector_j) {
            double accumulator = 0;
            auto iterator = unitary_matrix->GetNewIterator(index_of_vector_i);
            while (iterator->hasNext()) {
                auto item = iterator->getNext();
                if (!unitary_matrix->is_zero(index_of_vector_j, item.index)) {
                    accumulator += item.value * unitary_matrix->at(index_of_vector_j, item.index);
                }
            }
            // TODO: epsilon
            if (std::abs(accumulator) > 1e-9) {
                answer = false;
                break;
            }
        }
    }
    return answer;
}

size_t calculateTotalSpaceSize(const std::vector<int>& mults) {
    size_t acc = 1;
    for (const auto& mult : mults) {
        acc *= mult;
    }
    return acc;
}

TEST(symmetrizer, 4444) {
    std::vector<int> mults = {4, 4, 4, 4};
    model::ModelInput model(mults);
    uint32_t totalSpaceSize = calculateTotalSpaceSize(mults);

    // S2 group
    {
        common::physical_optimization::OptimizationList optimizationList;
        optimizationList.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});
        space::Space space =
            space::optimization::OptimizedSpaceConstructor::construct({model, optimizationList});
        EXPECT_EQ(totalSpaceSize, number_of_vectors(space));
        EXPECT_TRUE(orthogonality_of_basis(space)) << "Vectors are not orthogonal";
    }
    // S2 * the same S2 group
    {
        common::physical_optimization::OptimizationList optimizationList;
        optimizationList.Symmetrize(group::Group::S2, {{1, 0, 3, 2}})
            .Symmetrize(group::Group::S2, {{1, 0, 3, 2}});
        space::Space space =
            space::optimization::OptimizedSpaceConstructor::construct({model, optimizationList});
        EXPECT_EQ(totalSpaceSize, number_of_vectors(space));
        EXPECT_TRUE(orthogonality_of_basis(space)) << "Vectors are not orthogonal";
    }
}

TEST(symmetrizer, 333) {
    std::vector<int> mults = {3, 3, 3};
    model::ModelInput model(mults);
    uint32_t totalSpaceSize = calculateTotalSpaceSize(mults);

    // S3
    {
        common::physical_optimization::OptimizationList optimizationList;
        optimizationList.Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}});
        space::Space space =
            space::optimization::OptimizedSpaceConstructor::construct({model, optimizationList});

        EXPECT_EQ(totalSpaceSize, number_of_vectors(space));
        EXPECT_TRUE(orthogonality_of_basis(space)) << "Vectors are not orthogonal";
    }
    // TODO: move these testes to another place!
    // S3 * the same S3
    {
        common::physical_optimization::OptimizationList optimizationList;
        optimizationList.Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}})
            .Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}});
        space::Space space =
            space::optimization::OptimizedSpaceConstructor::construct({model, optimizationList});
        EXPECT_EQ(totalSpaceSize, number_of_vectors(space));
        EXPECT_TRUE(orthogonality_of_basis(space)) << "Vectors are not orthogonal";
    }
    // S3 * the same S3 (but different generators)
    {
        common::physical_optimization::OptimizationList optimizationList;
        optimizationList.Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}})
            .Symmetrize(group::Group::S3, {{2, 0, 1}, {1, 0, 2}});
        space::Space space =
            space::optimization::OptimizedSpaceConstructor::construct({model, optimizationList});

        EXPECT_EQ(totalSpaceSize, number_of_vectors(space));
        EXPECT_TRUE(orthogonality_of_basis(space)) << "Vectors are not orthogonal";
    }
}

TEST(symmetrizer, 333333) {
    std::vector<int> mults = {3, 3, 3, 3, 3, 3};
    model::ModelInput model(mults);
    uint32_t totalSpaceSize = calculateTotalSpaceSize(mults);

    // S3
    {
        common::physical_optimization::OptimizationList optimizationList;
        optimizationList.Symmetrize(group::Group::S3, {{1, 2, 0, 4, 5, 3}, {0, 2, 1, 3, 5, 4}});
        space::Space space =
            space::optimization::OptimizedSpaceConstructor::construct({model, optimizationList});

        EXPECT_EQ(totalSpaceSize, number_of_vectors(space));
        EXPECT_TRUE(orthogonality_of_basis(space)) << "Vectors are not orthogonal";
    }
    // S3 * S2, S2 * S3, their commutativity
    {
        // S3 * S2
        common::physical_optimization::OptimizationList optimizationList_first;
        optimizationList_first
            .Symmetrize(group::Group::S3, {{1, 2, 0, 4, 5, 3}, {0, 2, 1, 3, 5, 4}})
            .Symmetrize(group::Group::S2, {{3, 4, 5, 0, 1, 2}});
        space::Space space_first = space::optimization::OptimizedSpaceConstructor::construct(
            {model, optimizationList_first});

        EXPECT_EQ(totalSpaceSize, number_of_vectors(space_first));
        EXPECT_TRUE(orthogonality_of_basis(space_first)) << "Vectors are not orthogonal";

        // S2 * S3
        common::physical_optimization::OptimizationList optimizationList_second;
        optimizationList_second.Symmetrize(group::Group::S2, {{3, 4, 5, 0, 1, 2}})
            .Symmetrize(group::Group::S3, {{1, 2, 0, 4, 5, 3}, {0, 2, 1, 3, 5, 4}});
        space::Space space_second = space::optimization::OptimizedSpaceConstructor::construct(
            {model, optimizationList_second});

        EXPECT_EQ(totalSpaceSize, number_of_vectors(space_second));
        EXPECT_TRUE(orthogonality_of_basis(space_second)) << "Vectors are not orthogonal";

        // check equivalence of S2*S3 and S3*S2:
        for (const auto& subspace_first : space_first.getBlocks()) {
            for (const auto& subspace_second : space_second.getBlocks()) {
                if (subspace_first.properties.representation[0]
                        == subspace_second.properties.representation[1]
                    && subspace_first.properties.representation[1]
                        == subspace_second.properties.representation[0]) {
                    EXPECT_TRUE(subspace_first.decomposition->is_equal_up_to_vector_order(
                        subspace_second.decomposition));
                }
            }
        }
    }
}

TEST(symmetrizer, 222222222_S3xS3) {
    std::vector<int> mults = {2, 2, 2, 2, 2, 2, 2, 2, 2};
    model::ModelInput model(mults);
    uint32_t totalSpaceSize = calculateTotalSpaceSize(mults);

    common::physical_optimization::OptimizationList optimizationList;
    optimizationList
        .Symmetrize(group::Group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}})
        .Symmetrize(group::Group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});
    space::Space space =
        space::optimization::OptimizedSpaceConstructor::construct({model, optimizationList});

    EXPECT_EQ(totalSpaceSize, number_of_vectors(space));
    EXPECT_TRUE(orthogonality_of_basis(space)) << "Vectors are not orthogonal";
}

TEST(symmetrizer, 333333333_S3xS3) {
    std::vector<int> mults = {3, 3, 3, 3, 3, 3, 3, 3, 3};
    model::ModelInput model(mults);
    uint32_t totalSpaceSize = calculateTotalSpaceSize(mults);

    common::physical_optimization::OptimizationList optimizationList;
    optimizationList
        .Symmetrize(group::Group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}})
        .Symmetrize(group::Group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});
    space::Space space =
        space::optimization::OptimizedSpaceConstructor::construct({model, optimizationList});

    EXPECT_EQ(totalSpaceSize, number_of_vectors(space));
    EXPECT_TRUE(orthogonality_of_basis(space)) << "Vectors are not orthogonal";
}
