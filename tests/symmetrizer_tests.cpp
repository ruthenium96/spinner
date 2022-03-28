#include "gtest/gtest.h"
#include "src/common/Logger.h"
#include "src/common/runner/Runner.h"

size_t number_of_vectors(const space::Space& space) {
    size_t acc = 0;
    for (const auto& subspace : space.getBlocks()) {
        acc += subspace.decomposition.size();
    }
    return acc;
}

bool orthogonality_of_basis(const space::Space& space) {
    UnitarySparseMatrix unitary_matrix;
    for (const auto& subspace : space.getBlocks()) {
        unitary_matrix.copy_all_from(subspace.decomposition);
    }
    bool answer = true;
#pragma omp parallel for shared(space, unitary_matrix, answer) default(none)
    for (size_t index_of_vector_i = 0; index_of_vector_i < unitary_matrix.size();
         ++index_of_vector_i) {
        for (size_t index_of_vector_j = index_of_vector_i + 1;
             index_of_vector_j < unitary_matrix.size();
             ++index_of_vector_j) {
            double accumulator = 0;
            auto iterator = unitary_matrix.GetNewIterator(index_of_vector_i);
            while (iterator->hasNext()) {
                auto item = iterator->getNext();
                if (!unitary_matrix.is_zero(index_of_vector_j, item.index)) {
                    accumulator += item.value * unitary_matrix(index_of_vector_j, item.index);
                }
            }
            if (accumulator != 0) {
                answer = false;
                break;
            }
        }
    }
    return answer;
}

TEST(symmetrizer, throw_wrong_size_of_pemutation) {
    std::vector<int> mults = {4, 4, 4};
    model::Model model(mults);

    runner::Runner runner(model);
    EXPECT_THROW(runner.Symmetrize(group::Group::S2, {{1, 0, 3, 2}}), std::length_error);
}

TEST(symmetrizer, throw_permutes_different_multiplicities) {
    std::vector<int> mults = {4, 4, 4, 3};
    model::Model model(mults);

    runner::Runner runner(model);
    EXPECT_THROW(runner.Symmetrize(group::Group::S2, {{1, 0, 3, 2}}), std::invalid_argument);
}

TEST(symmetrizer, 4444_S2) {
    std::vector<int> mults = {4, 4, 4, 4};
    model::Model model(mults);

    runner::Runner runner(model);
    runner.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});
    EXPECT_EQ(
        runner.getIndexConverter().get_total_space_size(),
        number_of_vectors(runner.getSpace()));
    EXPECT_TRUE(orthogonality_of_basis(runner.getSpace())) << "Vectors are not orthogonal";
}

TEST(symmetrizer, 4444_doubleS2) {
    std::vector<int> mults = {4, 4, 4, 4};
    model::Model model(mults);

    runner::Runner runner(model);
    runner.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});
    runner.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});
    EXPECT_EQ(
        runner.getIndexConverter().get_total_space_size(),
        number_of_vectors(runner.getSpace()));
    EXPECT_TRUE(orthogonality_of_basis(runner.getSpace())) << "Vectors are not orthogonal";
}

TEST(symmetrizer, 333_S3) {
    std::vector<int> mults = {3, 3, 3};
    model::Model model(mults);

    runner::Runner runner(model);
    runner.Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}});
    EXPECT_EQ(
        runner.getIndexConverter().get_total_space_size(),
        number_of_vectors(runner.getSpace()));
    EXPECT_TRUE(orthogonality_of_basis(runner.getSpace())) << "Vectors are not orthogonal";
}

TEST(symmetrizer, 333_doubleS3) {
    std::vector<int> mults = {3, 3, 3};
    model::Model model(mults);

    runner::Runner runner(model);
    runner.Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}});
    runner.Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}});
    EXPECT_EQ(
        runner.getIndexConverter().get_total_space_size(),
        number_of_vectors(runner.getSpace()));
    EXPECT_TRUE(orthogonality_of_basis(runner.getSpace())) << "Vectors are not orthogonal";
}

TEST(symmetrizer, 333_doubleS3_tricky) {
    std::vector<int> mults = {3, 3, 3};
    model::Model model(mults);

    runner::Runner runner(model);
    runner.Symmetrize(group::Group::S3, {{1, 2, 0}, {0, 2, 1}});
    runner.Symmetrize(group::Group::S3, {{2, 0, 1}, {1, 0, 2}});
    EXPECT_EQ(
        runner.getIndexConverter().get_total_space_size(),
        number_of_vectors(runner.getSpace()));
    EXPECT_TRUE(orthogonality_of_basis(runner.getSpace())) << "Vectors are not orthogonal";
}

TEST(symmetrizer, 333333_S3) {
    std::vector<int> mults = {3, 3, 3, 3, 3, 3};
    model::Model model(mults);

    runner::Runner runner(model);
    runner.Symmetrize(group::Group::S3, {{1, 2, 0, 4, 5, 3}, {0, 2, 1, 3, 5, 4}});
    EXPECT_EQ(
        runner.getIndexConverter().get_total_space_size(),
        number_of_vectors(runner.getSpace()));
    EXPECT_TRUE(orthogonality_of_basis(runner.getSpace())) << "Vectors are not orthogonal";
}

TEST(symmetrizer, 333333_S3xS2) {
    std::vector<int> mults = {3, 3, 3, 3, 3, 3};
    model::Model model(mults);

    runner::Runner runner(model);
    runner.Symmetrize(group::Group::S3, {{1, 2, 0, 4, 5, 3}, {0, 2, 1, 3, 5, 4}});
    runner.Symmetrize(group::Group::S2, {{3, 4, 5, 0, 1, 2}});

    EXPECT_EQ(
        runner.getIndexConverter().get_total_space_size(),
        number_of_vectors(runner.getSpace()));
    EXPECT_TRUE(orthogonality_of_basis(runner.getSpace())) << "Vectors are not orthogonal";
}

TEST(symmetrizer, 333333_S2xS3) {
    std::vector<int> mults = {3, 3, 3, 3, 3, 3};
    model::Model model(mults);

    runner::Runner runner(model);
    runner.Symmetrize(group::Group::S2, {{3, 4, 5, 0, 1, 2}});
    runner.Symmetrize(group::Group::S3, {{1, 2, 0, 4, 5, 3}, {0, 2, 1, 3, 5, 4}});

    EXPECT_EQ(
        runner.getIndexConverter().get_total_space_size(),
        number_of_vectors(runner.getSpace()));
    EXPECT_TRUE(orthogonality_of_basis(runner.getSpace())) << "Vectors are not orthogonal";
}

TEST(symmetrizer, 333333_S2xS3_S3xS2_direct_product_commutativity) {
    std::vector<int> mults = {3, 3, 3, 3, 3, 3};
    model::Model model(mults);

    runner::Runner runner_first(model);
    runner_first.Symmetrize(group::Group::S2, {{3, 4, 5, 0, 1, 2}});
    runner_first.Symmetrize(group::Group::S3, {{1, 2, 0, 4, 5, 3}, {0, 2, 1, 3, 5, 4}});

    runner::Runner runner_second(model);
    runner_second.Symmetrize(group::Group::S3, {{1, 2, 0, 4, 5, 3}, {0, 2, 1, 3, 5, 4}});
    runner_second.Symmetrize(group::Group::S2, {{3, 4, 5, 0, 1, 2}});

    for (const auto& subspace_first : runner_first.getSpace().getBlocks()) {
        for (const auto& subspace_second : runner_second.getSpace().getBlocks()) {
            if (subspace_first.properties.representation[0]
                    == subspace_second.properties.representation[1]
                && subspace_first.properties.representation[1]
                    == subspace_second.properties.representation[0]) {
                EXPECT_TRUE(subspace_first.decomposition.is_equal_up_to_vector_order(
                    subspace_second.decomposition));
            }
        }
    }
}

TEST(symmetrizer, 222222222_S3xS3) {
    std::vector<int> mults = {2, 2, 2, 2, 2, 2, 2, 2, 2};
    model::Model model(mults);

    runner::Runner runner(model);
    runner.Symmetrize(group::Group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
    runner.Symmetrize(group::Group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});
    EXPECT_EQ(
        runner.getIndexConverter().get_total_space_size(),
        number_of_vectors(runner.getSpace()));
    EXPECT_TRUE(orthogonality_of_basis(runner.getSpace())) << "Vectors are not orthogonal";
}

TEST(symmetrizer, 333333333_S3xS3) {
    std::vector<int> mults = {3, 3, 3, 3, 3, 3, 3, 3, 3};
    model::Model model(mults);

    runner::Runner runner(model);
    runner.Symmetrize(group::Group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
    runner.Symmetrize(group::Group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});
    EXPECT_EQ(
        runner.getIndexConverter().get_total_space_size(),
        number_of_vectors(runner.getSpace()));
    EXPECT_TRUE(orthogonality_of_basis(runner.getSpace())) << "Vectors are not orthogonal";
}
