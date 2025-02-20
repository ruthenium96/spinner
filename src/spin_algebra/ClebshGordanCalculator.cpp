#include "ClebshGordanCalculator.h"

#include <omp.h>
#include <wignerSymbols.h>
#include <string>

namespace spin_algebra {

ClebshGordanCalculator::ClebshGordanCalculator() {
    hashed_CGs_for_all_threads = std::vector<Map>(omp_get_max_threads());
    hashed_W6Js_for_all_threads = std::vector<Map>(omp_get_max_threads());
}

double ClebshGordanCalculator::clebsh_gordan_coefficient(
    double l1,
    double l2,
    double l3,
    double m1,
    double m2) const {
    if (l1 + l2 < l3 || std::abs(l1 - l2) > l3) {
        return 0;
    }
    auto key = to_key(l1, l2, l3, m1, m2);
    // todo: calculate omp_get_thread_num() outside of this function?
    Map& hashed_CGs_per_thread = hashed_CGs_for_all_threads.at(omp_get_thread_num());

    if (hashed_CGs_per_thread.contains(key)) {
        return hashed_CGs_per_thread.at(key);
    } else {
        double value = WignerSymbols::clebschGordan(l1, l2, l3, m1, m2, m1 + m2);
        hashed_CGs_per_thread.emplace(key, value);
        return value;
    }
}

double ClebshGordanCalculator::wigner6j_coefficient(
    double l1,
    double l2,
    double l3,
    double l4,
    double l5,
    double l6) const {
    auto key = to_key(l1, l2, l3, l4, l5, l6);
    Map& hashed_W6J_per_thread = hashed_W6Js_for_all_threads.at(omp_get_thread_num());

    if (hashed_W6J_per_thread.contains(key)) {
        return hashed_W6J_per_thread.at(key);
    } else {
        double value = WignerSymbols::wigner6j(l1, l2, l3, l4, l5, l6);
        hashed_W6J_per_thread.emplace(key, value);
        return value;
    }
}

uint64_t
ClebshGordanCalculator::to_key(double l1, double l2, double l3, double m1, double m2) noexcept {
    // NB: each value cannot be bigger than BASE
    const uint64_t BASE = 4096;
    auto ul1 = (uint64_t)(2 * l1);
    auto ul2 = (uint64_t)(2 * l2);
    auto ul3 = (uint64_t)(2 * l3);
    auto um1 = (uint64_t)(m1 + l1);
    auto um2 = (uint64_t)(m2 + l2);
    uint64_t answer = 0;
    answer |= ul1;
    answer |= ul2 * BASE;
    answer |= ul3 * BASE * BASE;
    answer |= um1 * BASE * BASE * BASE;
    answer |= um2 * BASE * BASE * BASE * BASE;
    return answer;
}

uint64_t ClebshGordanCalculator::to_key(double l1, double l2, double l3,
                                        double l4, double l5, double l6) noexcept {
    // NB: each value cannot be bigger than BASE
    const uint64_t BASE = 1024;
    auto ul1 = (uint64_t)(2 * l1);
    auto ul2 = (uint64_t)(2 * l2);
    auto ul3 = (uint64_t)(2 * l3);
    auto ul4 = (uint64_t)(2 * l4);
    auto ul5 = (uint64_t)(2 * l5);
    auto ul6 = (uint64_t)(2 * l6);
    uint64_t answer = 0;
    answer |= ul1;
    answer |= ul2 * BASE;
    answer |= ul3 * BASE * BASE;
    answer |= ul4 * BASE * BASE * BASE;
    answer |= ul5 * BASE * BASE * BASE * BASE;
    answer |= ul6 * BASE * BASE * BASE * BASE * BASE;
    return answer;
}

double ClebshGordanCalculator::ninej_element(
    double left_one,
    double left_two,
    double left_fin,
    double right_one,
    double right_two,
    double right_fin,
    uint8_t rank_one,
    uint8_t rank_two,
    uint8_t rank_fin) const {
    if (rank_one == 0 && rank_two == 0 && rank_fin == 0) {
        if (left_one != right_one || left_two != right_two || left_fin != right_fin) {
            return 0;
        }
        return 1 / sqrt((2 * left_one + 1) * (2 * left_two + 1) * (2 * left_fin + 1));
    } else if (rank_one == 0 && rank_two == rank_fin) {
        if (left_one != right_one) {
            return 0;
        }
        double answer = 1;
        if ((int) (left_one + right_two + left_fin + rank_two) % 2 == 1) {
            answer *= -1;
        }
        answer /= sqrt((2 * rank_two + 1) * (2 * left_one + 1));
        answer *= wigner6j_coefficient(right_fin, left_fin, rank_two, left_two, right_two, left_one);
        return answer;
    } else if (rank_one == rank_fin && rank_two == 0) {
        if (left_two != right_two) {
            return 0;
        }
        double answer = 1;
        if ((int) (left_one + right_fin + left_two + rank_one) % 2 == 1) {
            answer *= -1;
        }
        answer /= sqrt((2 * rank_one + 1) * (2 * left_two + 1));
        answer *= wigner6j_coefficient(right_one, left_one, rank_one, left_fin, right_fin, left_two);
        return answer;
    } else if (rank_one == rank_two && rank_fin == 0) {
        if (left_fin != right_fin) {
            return 0;
        }
        double answer = 1;
        if ((int) (right_one + left_two + left_fin + rank_one) % 2 == 1) {
            answer *= -1;
        }
        answer /= sqrt((2 * rank_one + 1) * (2 * left_fin + 1));
        answer *= wigner6j_coefficient(left_one, right_one, rank_one, right_two, left_two, left_fin);
        return answer;
    } else if (rank_one == 1 && rank_two == 1 && rank_fin == 2) {
        // x can be right_one - 1, right_one and right_one + 1 due triangle rule
        double answer = 0;
        for (double x = right_one - 1; x <= right_one + 1; ++x) {
            if (x < 0) {
                continue;
            }
            double temp_answer = 2 * x + 1;
            temp_answer *= wigner6j_coefficient(left_one, right_one, 1, 1, 2, x);
            temp_answer *= wigner6j_coefficient(left_two, right_two, 1, right_one, x, right_fin);
            temp_answer *= wigner6j_coefficient(left_fin, right_fin, 2, x, left_one, left_two);
            answer += temp_answer;
        }
        if ((int)(2 * right_one) % 2 == 1) {
            answer *= -1;
        }
        return answer;
    } else {
        std::string ranks = std::to_string(rank_one) + " " + std::to_string(rank_two) + " " + std::to_string(rank_fin); 
        throw std::invalid_argument("Invalid value of rank in 9j-symbol: " + ranks);
    }
}
}  // namespace spin_algebra