#ifndef SPINNER_ONEORMANY_H
#define SPINNER_ONEORMANY_H

#include <cassert>
#include <functional>
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

template <typename T, typename U, typename W>
OneOrMany<W> transform_one_or_many(
    std::function<W(T, U)> f, 
    const OneOrMany<T>& one_or_many_first,
    const OneOrMany<U>& one_or_many_second) {
    if (holdsOne(one_or_many_first) && holdsOne(one_or_many_second)) {
        const T& one_first = getOneRef(one_or_many_first);
        const U& one_second = getOneRef(one_or_many_second);
        return f(one_first, one_second);
    } else if (holdsMany(one_or_many_first) && holdsMany(one_or_many_second)) {
        const std::vector<T>& many_first = getManyRef(one_or_many_first);
        const std::vector<U>& many_second = getManyRef(one_or_many_second);
        assert(many_first.size() == many_second.size());

        std::vector<W> answer;
        for (int i = 0; i < many_first.size(); ++i) {
            const auto& el_first = many_first[i];
            const auto& el_second = many_second[i];
            answer.emplace_back(f(el_first, el_second));
        }
        return answer;
    } else if (holdsMany(one_or_many_first) && holdsOne(one_or_many_second)) {
        const std::vector<T>& many_first = getManyRef(one_or_many_first);
        const U& one_second = getOneRef(one_or_many_second);

        std::vector<W> answer;
        for (int i = 0; i < many_first.size(); ++i) {
            const T& el_first = many_first[i];
            answer.emplace_back(f(el_first, one_second));
        }
        return answer;
    } else { // holdsOne(one_or_many_first) && holdsMany(one_or_many_second)
        const T& one_first = getOneRef(one_or_many_first);
        const std::vector<U>& many_second = getManyRef(one_or_many_second);

        std::vector<W> answer;
        for (int i = 0; i < many_second.size(); ++i) {
            const U& el_second = many_second[i];
            answer.emplace_back(f(one_first, el_second));
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