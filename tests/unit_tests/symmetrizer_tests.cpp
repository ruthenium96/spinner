#include "gtest/gtest.h"
#include "src/common/Logger.h"
#include "src/space/optimization/OptimizedSpaceConstructor.h"

const auto factories = quantum::linear_algebra::FactoriesList();

size_t number_of_vectors(const space::Space& space) {
    size_t acc = 0;
    for (const auto& subspace : space.getBlocks()) {
        acc += subspace.size() * subspace.properties.degeneracy;
    }
    return acc;
}

void copySparseUnitaryMatrixToSparseUnitaryMatrix(
    std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>& matrix_to,
    const std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>& matrix_from,
    size_t shift) {
    size_t matrix_in_space_basis_size = matrix_from->size_cols();

    for (uint32_t index_of_space_vector_i = 0; index_of_space_vector_i < matrix_in_space_basis_size;
         ++index_of_space_vector_i) {
        auto iterator = matrix_from->GetNewIterator(index_of_space_vector_i);
        while (iterator->hasNext()) {
            auto item = iterator->getNext();
            uint32_t index_of_lexicographic_vector_k = item.index;
            double value = item.value;
            matrix_to->add_to_position(
                value,
                index_of_space_vector_i + shift,
                index_of_lexicographic_vector_k);
        }
    }
}

std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>
concatenateSpace(const space::Space& space) {
    EXPECT_FALSE(space.getBlocks().empty());
    auto unitary_matrix = factories.createSparseSemiunitaryMatrix(
        number_of_vectors(space),
        space.getBlocks()[0].decomposition->size_rows());

    size_t vector_number = 0;
    for (const auto& subspace : space.getBlocks()) {
        copySparseUnitaryMatrixToSparseUnitaryMatrix(
            unitary_matrix,
            subspace.decomposition,
            vector_number);
        vector_number += subspace.decomposition->size_cols();
    }
    return unitary_matrix;
}

bool orthogonality_of_basis(const space::Space& space) {
    auto unitary_matrix = concatenateSpace(space);
    bool answer = true;
#pragma omp parallel for shared(space, unitary_matrix, answer) default(none) collapse(2)
    for (size_t index_of_vector_i = 0; index_of_vector_i < unitary_matrix->size_cols();
         ++index_of_vector_i) {
        for (size_t index_of_vector_j = index_of_vector_i + 1;
             index_of_vector_j < unitary_matrix->size_cols();
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
            }
        }
    }
    return answer;
}

size_t calculateTotalSpaceSize(const std::vector<spin_algebra::Multiplicity>& mults) {
    size_t acc = 1;
    for (const auto& mult : mults) {
        acc *= mult;
    }
    return acc;
}

bool isEqualUpToVectorOrder(
    const std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>& lhs,
    const std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>& rhs) {
    std::vector<std::map<size_t, double>> lhs_(lhs->size_cols());
    std::vector<std::map<size_t, double>> rhs_(rhs->size_cols());

    for (uint32_t index_of_space_vector_i = 0; index_of_space_vector_i < lhs->size_cols();
         ++index_of_space_vector_i) {
        auto iterator = lhs->GetNewIterator(index_of_space_vector_i);
        while (iterator->hasNext()) {
            auto item = iterator->getNext();
            lhs_[index_of_space_vector_i][item.index] = item.value;
        }
    }
    for (uint32_t index_of_space_vector_i = 0; index_of_space_vector_i < rhs->size_cols();
         ++index_of_space_vector_i) {
        auto iterator = rhs->GetNewIterator(index_of_space_vector_i);
        while (iterator->hasNext()) {
            auto item = iterator->getNext();
            rhs_[index_of_space_vector_i][item.index] = item.value;
        }
    }

    std::sort(lhs_.begin(), lhs_.end());
    std::sort(rhs_.begin(), rhs_.end());

    return (lhs_ == rhs_);
}

TEST(symmetrizer, 2222_S2_broke_unitary_matrices) {
    std::vector<spin_algebra::Multiplicity> mults = {2, 2, 2, 2};
    model::ModelInput model(mults);
    uint32_t totalSpaceSize = calculateTotalSpaceSize(mults);

    {
        common::physical_optimization::OptimizationList optimizationList;
        optimizationList.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});
        space::Space space = space::optimization::OptimizedSpaceConstructor::construct(
            {model, optimizationList},
            factories);
        for (auto& subspace : space.getBlocks()) {
            subspace.decomposition->add_to_position(1, 0, 0);
        }
        EXPECT_FALSE(orthogonality_of_basis(space)) << "These tests does not work "
                                                       "and, probably, always returns true";
    }
}

TEST(symmetrizer, 4444) {
    std::vector<spin_algebra::Multiplicity> mults = {4, 4, 4, 4};
    model::ModelInput model(mults);
    uint32_t totalSpaceSize = calculateTotalSpaceSize(mults);

    // S2 group
    {
        common::physical_optimization::OptimizationList optimizationList;
        optimizationList.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});
        space::Space space = space::optimization::OptimizedSpaceConstructor::construct(
            {model, optimizationList},
            factories);
        EXPECT_EQ(totalSpaceSize, number_of_vectors(space));
        EXPECT_TRUE(orthogonality_of_basis(space)) << "Vectors are not orthogonal";
    }
    // S2 * the same S2 group
    {
        common::physical_optimization::OptimizationList optimizationList;
        optimizationList.Symmetrize(group::Group::S2, {{1, 0, 3, 2}})
            .Symmetrize(group::Group::S2, {{1, 0, 3, 2}});
        space::Space space = space::optimization::OptimizedSpaceConstructor::construct(
            {model, optimizationList},
            factories);
        EXPECT_EQ(totalSpaceSize, number_of_vectors(space));
        EXPECT_TRUE(orthogonality_of_basis(space)) << "Vectors are not orthogonal";
    }
    // S2 * other S2
    {
        common::physical_optimization::OptimizationList optimizationList;
        optimizationList.Symmetrize(group::Group::S2, {{1, 0, 3, 2}})
            .Symmetrize(group::Group::S2, {{3, 2, 1, 0}});
        space::Space space = space::optimization::OptimizedSpaceConstructor::construct(
            {model, optimizationList},
            factories);
        EXPECT_EQ(totalSpaceSize, number_of_vectors(space));
        EXPECT_TRUE(orthogonality_of_basis(space)) << "Vectors are not orthogonal";
    }
}

TEST(symmetrizer, 333) {
    std::vector<spin_algebra::Multiplicity> mults = {3, 3, 3};
    model::ModelInput model(mults);
    uint32_t totalSpaceSize = calculateTotalSpaceSize(mults);

    // S3
    {
        common::physical_optimization::OptimizationList optimizationList;
        optimizationList.Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}});
        space::Space space = space::optimization::OptimizedSpaceConstructor::construct(
            {model, optimizationList},
            factories);

        EXPECT_EQ(totalSpaceSize, number_of_vectors(space));
        EXPECT_TRUE(orthogonality_of_basis(space)) << "Vectors are not orthogonal";
    }
    // TODO: move these testes to another place!
    // S3 * the same S3
    {
        common::physical_optimization::OptimizationList optimizationList;
        optimizationList.Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}})
            .Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}});
        space::Space space = space::optimization::OptimizedSpaceConstructor::construct(
            {model, optimizationList},
            factories);
        EXPECT_EQ(totalSpaceSize, number_of_vectors(space));
        EXPECT_TRUE(orthogonality_of_basis(space)) << "Vectors are not orthogonal";
    }
    // S3 * the same S3 (but different generators)
    {
        common::physical_optimization::OptimizationList optimizationList;
        optimizationList.Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}})
            .Symmetrize(group::Group::S3, {{2, 0, 1}, {1, 0, 2}});
        space::Space space = space::optimization::OptimizedSpaceConstructor::construct(
            {model, optimizationList},
            factories);

        EXPECT_EQ(totalSpaceSize, number_of_vectors(space));
        EXPECT_TRUE(orthogonality_of_basis(space)) << "Vectors are not orthogonal";
    }
}

TEST(symmetrizer, 222222) {
    std::vector<spin_algebra::Multiplicity> mults = {2, 2, 2, 2, 2, 2};
    model::ModelInput model(mults);
    uint32_t totalSpaceSize = calculateTotalSpaceSize(mults);

    // S3
    {
        common::physical_optimization::OptimizationList optimizationList;
        optimizationList.Symmetrize(group::Group::S3, {{1, 2, 0, 4, 5, 3}, {0, 2, 1, 3, 5, 4}});
        space::Space space = space::optimization::OptimizedSpaceConstructor::construct(
            {model, optimizationList},
            factories);

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
            {model, optimizationList_first},
            quantum::linear_algebra::FactoriesList());

        EXPECT_EQ(totalSpaceSize, number_of_vectors(space_first));
        EXPECT_TRUE(orthogonality_of_basis(space_first)) << "Vectors are not orthogonal";

        // S2 * S3
        common::physical_optimization::OptimizationList optimizationList_second;
        optimizationList_second.Symmetrize(group::Group::S2, {{3, 4, 5, 0, 1, 2}})
            .Symmetrize(group::Group::S3, {{1, 2, 0, 4, 5, 3}, {0, 2, 1, 3, 5, 4}});
        space::Space space_second = space::optimization::OptimizedSpaceConstructor::construct(
            {model, optimizationList_second},
            factories);

        EXPECT_EQ(totalSpaceSize, number_of_vectors(space_second));
        EXPECT_TRUE(orthogonality_of_basis(space_second)) << "Vectors are not orthogonal";

        // check equivalence of S2*S3 and S3*S2:
        for (const auto& subspace_first : space_first.getBlocks()) {
            for (const auto& subspace_second : space_second.getBlocks()) {
                if (subspace_first.properties.representation[0]
                        == subspace_second.properties.representation[1]
                    && subspace_first.properties.representation[1]
                        == subspace_second.properties.representation[0]) {
                    EXPECT_TRUE(isEqualUpToVectorOrder(
                        subspace_first.decomposition,
                        subspace_second.decomposition));
                }
            }
        }
    }
}

TEST(symmetrizer, 333333) {
    std::vector<spin_algebra::Multiplicity> mults = {3, 3, 3, 3, 3, 3};
    model::ModelInput model(mults);
    uint32_t totalSpaceSize = calculateTotalSpaceSize(mults);

    // S3
    {
        common::physical_optimization::OptimizationList optimizationList;
        optimizationList.Symmetrize(group::Group::S3, {{1, 2, 0, 4, 5, 3}, {0, 2, 1, 3, 5, 4}});
        space::Space space = space::optimization::OptimizedSpaceConstructor::construct(
            {model, optimizationList},
            factories);

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
            {model, optimizationList_first},
            quantum::linear_algebra::FactoriesList());

        EXPECT_EQ(totalSpaceSize, number_of_vectors(space_first));
        EXPECT_TRUE(orthogonality_of_basis(space_first)) << "Vectors are not orthogonal";

        // S2 * S3
        common::physical_optimization::OptimizationList optimizationList_second;
        optimizationList_second.Symmetrize(group::Group::S2, {{3, 4, 5, 0, 1, 2}})
            .Symmetrize(group::Group::S3, {{1, 2, 0, 4, 5, 3}, {0, 2, 1, 3, 5, 4}});
        space::Space space_second = space::optimization::OptimizedSpaceConstructor::construct(
            {model, optimizationList_second},
            factories);

        EXPECT_EQ(totalSpaceSize, number_of_vectors(space_second));
        EXPECT_TRUE(orthogonality_of_basis(space_second)) << "Vectors are not orthogonal";

        // check equivalence of S2*S3 and S3*S2:
        for (const auto& subspace_first : space_first.getBlocks()) {
            for (const auto& subspace_second : space_second.getBlocks()) {
                if (subspace_first.properties.representation[0]
                        == subspace_second.properties.representation[1]
                    && subspace_first.properties.representation[1]
                        == subspace_second.properties.representation[0]) {
                    EXPECT_TRUE(isEqualUpToVectorOrder(
                        subspace_first.decomposition,
                        subspace_second.decomposition));
                }
            }
        }
    }
}

TEST(symmetrizer, 3333_D4) {
    std::vector<spin_algebra::Multiplicity> mults = {3, 3, 3, 3};
    model::ModelInput model(mults);
    uint32_t totalSpaceSize = calculateTotalSpaceSize(mults);

    common::physical_optimization::OptimizationList optimizationList;
    optimizationList.Symmetrize(group::Group::D4, {{1, 2, 3, 0}, {1, 0, 3, 2}});
    space::Space space = space::optimization::OptimizedSpaceConstructor::construct(
        {model, optimizationList},
        factories);

    EXPECT_EQ(totalSpaceSize, number_of_vectors(space));
    EXPECT_TRUE(orthogonality_of_basis(space)) << "Vectors are not orthogonal";

}

TEST(symmetrizer, 222222222_S3xS3) {
    std::vector<spin_algebra::Multiplicity> mults = {2, 2, 2, 2, 2, 2, 2, 2, 2};
    model::ModelInput model(mults);
    uint32_t totalSpaceSize = calculateTotalSpaceSize(mults);

    common::physical_optimization::OptimizationList optimizationList;
    optimizationList
        .Symmetrize(group::Group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}})
        .Symmetrize(group::Group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});
    space::Space space = space::optimization::OptimizedSpaceConstructor::construct(
        {model, optimizationList},
        factories);

    EXPECT_EQ(totalSpaceSize, number_of_vectors(space));
    EXPECT_TRUE(orthogonality_of_basis(space)) << "Vectors are not orthogonal";
}

TEST(symmetrizer, 222222222_S3xS3_different_one) {
    std::vector<spin_algebra::Multiplicity> mults = {2, 2, 2, 2, 2, 2, 2, 2, 2};
    model::ModelInput model(mults);
    uint32_t totalSpaceSize = calculateTotalSpaceSize(mults);

    common::physical_optimization::OptimizationList optimizationList;
    optimizationList
        .Symmetrize(group::Group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}})
        .Symmetrize(group::Group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {6, 7, 8, 3, 4, 5, 0, 1, 2}});
    space::Space space = space::optimization::OptimizedSpaceConstructor::construct(
        {model, optimizationList},
        factories);

    EXPECT_EQ(totalSpaceSize, number_of_vectors(space));
    EXPECT_TRUE(orthogonality_of_basis(space)) << "Vectors are not orthogonal";
}

TEST(symmetrizer, 222222222_S3xS3_different_two) {
    std::vector<spin_algebra::Multiplicity> mults = {2, 2, 2, 2, 2, 2, 2, 2, 2};
    model::ModelInput model(mults);
    uint32_t totalSpaceSize = calculateTotalSpaceSize(mults);

    common::physical_optimization::OptimizationList optimizationList;
    optimizationList
        .Symmetrize(group::Group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}})
        .Symmetrize(group::Group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {3, 4, 5, 0, 1, 2, 6, 7, 8}});
    space::Space space = space::optimization::OptimizedSpaceConstructor::construct(
        {model, optimizationList},
        factories);

    EXPECT_EQ(totalSpaceSize, number_of_vectors(space));
    EXPECT_TRUE(orthogonality_of_basis(space)) << "Vectors are not orthogonal";
}

TEST(symmetrizer, 333333333_S3xS3) {
    std::vector<spin_algebra::Multiplicity> mults = {3, 3, 3, 3, 3, 3, 3, 3, 3};
    model::ModelInput model(mults);
    uint32_t totalSpaceSize = calculateTotalSpaceSize(mults);

    common::physical_optimization::OptimizationList optimizationList;
    optimizationList
        .Symmetrize(group::Group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}})
        .Symmetrize(group::Group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});
    space::Space space = space::optimization::OptimizedSpaceConstructor::construct(
        {model, optimizationList},
        factories);

    EXPECT_EQ(totalSpaceSize, number_of_vectors(space));
    EXPECT_TRUE(orthogonality_of_basis(space)) << "Vectors are not orthogonal";
}

TEST(symmetrizer, 333333333_S3xS3_different_one) {
    std::vector<spin_algebra::Multiplicity> mults = {3, 3, 3, 3, 3, 3, 3, 3, 3};
    model::ModelInput model(mults);
    uint32_t totalSpaceSize = calculateTotalSpaceSize(mults);

    common::physical_optimization::OptimizationList optimizationList;
    optimizationList
        .Symmetrize(group::Group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}})
        .Symmetrize(group::Group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {6, 7, 8, 3, 4, 5, 0, 1, 2}});
    space::Space space = space::optimization::OptimizedSpaceConstructor::construct(
        {model, optimizationList},
        factories);

    EXPECT_EQ(totalSpaceSize, number_of_vectors(space));
    EXPECT_TRUE(orthogonality_of_basis(space)) << "Vectors are not orthogonal";
}

TEST(symmetrizer, 333333333_S3xS3_different_two) {
    std::vector<spin_algebra::Multiplicity> mults = {3, 3, 3, 3, 3, 3, 3, 3, 3};
    model::ModelInput model(mults);
    uint32_t totalSpaceSize = calculateTotalSpaceSize(mults);

    common::physical_optimization::OptimizationList optimizationList;
    optimizationList
        .Symmetrize(group::Group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}})
        .Symmetrize(group::Group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {3, 4, 5, 0, 1, 2, 6, 7, 8}});
    space::Space space = space::optimization::OptimizedSpaceConstructor::construct(
        {model, optimizationList},
        factories);

    EXPECT_EQ(totalSpaceSize, number_of_vectors(space));
    EXPECT_TRUE(orthogonality_of_basis(space)) << "Vectors are not orthogonal";
}