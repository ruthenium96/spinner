#ifndef SPINNER_ONEORMANY_H
#define SPINNER_ONEORMANY_H

#include <cassert>
#include <cstddef>
#include <functional>
#include <optional>
#include <stdexcept>
#include <variant>
#include <vector>

template <typename T>
using OneOrMany = std::variant<T, std::vector<T>>;

template <typename T> 
bool holdsOne(const OneOrMany<T>& t) {
    return std::holds_alternative<T>(t);
}

template <typename T> 
bool holdsMany(const OneOrMany<T>& t) {
    return std::holds_alternative<std::vector<T>>(t);
}

template <typename T> 
const T& getOneRef(const OneOrMany<T>& t) {
    return std::get<T>(t);
}

template <typename T> 
const std::vector<T>& getManyRef(const OneOrMany<T>& t) {
    return std::get<std::vector<T>>(t);
}

template <typename T, typename U>
OneOrMany<U> copyRef(const OneOrMany<T>& t) {
    if (holdsOne(t)) {
        return getOneRef(t);
    } else {
        const auto& many = getManyRef(t);
        std::vector<U> answer;
        for (const auto& el : many) {
            answer.push_back(el);
        }
        return answer;
    }
}

template <typename T, typename U>
OneOrMany<U> transform_one_or_many(std::function<U(T)> f, const OneOrMany<T>& one_or_many) {
    if (holdsOne(one_or_many)) {
        const T& one = getOneRef(one_or_many);
        return f(one);
    } else {
        const std::vector<T>& many = getManyRef(one_or_many);
        std::vector<U> answer;
        for (const auto& el : many) {
            answer.emplace_back(f(el));
        }
        return answer;
    }
}

template <typename T>
const T& get_element(const OneOrMany<T>& one_or_many, size_t i) {
    if (holdsOne(one_or_many)) {
        return getOneRef(one_or_many);
    } else {
        return getManyRef(one_or_many).at(i);
    }
}

template <typename T>
std::optional<size_t> get_size(const OneOrMany<T>& one_or_many) {
    if (holdsOne(one_or_many)) {
        return std::nullopt;
    } else {
        return getManyRef(one_or_many).size();
    }
}

template <typename T>
inline void update_size(const OneOrMany<T>& t, std::optional<size_t>& size) {
    auto mb_t_size = get_size(t);
    if (mb_t_size.has_value()) {
        if (size.has_value() && mb_t_size.value() != size.value()) {
            throw std::invalid_argument("Manys with different sizes were passed to transform");
        } else {
            size = mb_t_size.value();
        }
    }
}

template <typename T, typename U, typename W>
OneOrMany<W> transform_one_or_many(
    std::function<W(T, U)> f, 
    const OneOrMany<T>& one_or_many_first,
    const OneOrMany<U>& one_or_many_second) {
    std::optional<size_t> mb_size_of_many;
    update_size(one_or_many_first, mb_size_of_many);
    update_size(one_or_many_second, mb_size_of_many);

    if (!mb_size_of_many.has_value()) {
        const T& one_first = getOneRef(one_or_many_first);
        const U& one_second = getOneRef(one_or_many_second);
        return f(one_first, one_second);
    } else {
        std::vector<W> answer;
        for (int i = 0; i < mb_size_of_many.value(); ++i) {
            const auto& el_first = get_element(one_or_many_first, i);
            const auto& el_second = get_element(one_or_many_second, i);
            answer.emplace_back(f(el_first, el_second));
        }
        return answer;
    }
}

template <typename T, typename U, typename W>
OneOrMany<W> transform_one_or_many(
    std::function<W(const T&, const U&)> f, 
    const OneOrMany<T>& one_or_many_first,
    const OneOrMany<U>& one_or_many_second) {
    std::optional<size_t> mb_size_of_many;
    update_size(one_or_many_first, mb_size_of_many);
    update_size(one_or_many_second, mb_size_of_many);

    if (!mb_size_of_many.has_value()) {
        const T& one_first = getOneRef(one_or_many_first);
        const U& one_second = getOneRef(one_or_many_second);
        return f(one_first, one_second);
    } else {
        std::vector<W> answer;
        for (int i = 0; i < mb_size_of_many.value(); ++i) {
            const auto& el_first = get_element(one_or_many_first, i);
            const auto& el_second = get_element(one_or_many_second, i);
            answer.emplace_back(f(el_first, el_second));
        }
        return answer;
    }
}

template <typename T, typename U, typename W, typename Z>
OneOrMany<Z> transform_one_or_many(
    std::function<Z(T, U, W)> f, 
    const OneOrMany<T>& one_or_many_first,
    const OneOrMany<U>& one_or_many_second,
    const OneOrMany<W>& one_or_many_third) {
    std::optional<size_t> mb_size_of_many;
    update_size(one_or_many_first, mb_size_of_many);
    update_size(one_or_many_second, mb_size_of_many);
    update_size(one_or_many_third, mb_size_of_many);

    if (!mb_size_of_many.has_value()) {
        const T& one_first = getOneRef(one_or_many_first);
        const U& one_second = getOneRef(one_or_many_second);
        const W& one_third = getOneRef(one_or_many_third);
        return f(one_first, one_second, one_third);
    } else {
        std::vector<Z> answer;
        for (int i = 0; i < mb_size_of_many.value(); ++i) {
            const auto& el_first = get_element(one_or_many_first, i);
            const auto& el_second = get_element(one_or_many_second, i);
            const auto& el_third = get_element(one_or_many_third, i);
            answer.emplace_back(f(el_first, el_second, el_third));
        }
        return answer;
    }
}

template <typename T>
void apply_to_one_or_many(std::function<void(const T&)> f, const OneOrMany<T>& one_or_many) {
    if (holdsOne(one_or_many)) {
        const T& one = getOneRef(one_or_many);
        f(one);
        return;
    } else {
        const std::vector<T>& many = getManyRef(one_or_many);
        for (const auto& el : many) {
            f(el);
        }
        return;
    }
}

#endif // SPINNER_ONEORMANY_H