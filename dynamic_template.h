#pragma once
#include <memory>

// FOR CLANG_FORMAT constexpr-if FORMATING BUG!!!
#define static_if                                                                                  \
    if                                                                                             \
    constexpr

template <typename T, T t>
struct static_template {
    using type = T;
    static constexpr T value = t;
    constexpr operator T() { return t; }
};

template <typename T, T low, T high>
struct dynamic_template {
    template <typename F>
    dynamic_template(T target, F f)
    {
        eval(target, f);
    }
    template <typename F>
    static void eval(T target, F f)
    {
        eval(target, f, [] {});
    }
    template <typename F, typename R, typename result_t = void>
    static result_t eval(T target, F f, R r)
    {
        static_if(low < high)
        {
            if (target < low || target >= high) return r();
            if (target == low) return f(static_template<T, low>{});
            return dynamic_template<T, low + 1, high>::template eval<F, R, result_t>(
                std::forward<T>(target), std::forward<F>(f), std::forward<R>(r));
        }
        else return r();
    }
};

template <typename T, T low, T high, bool reverse = false>
struct static_for {
    template <typename F>
    static_for(F f)
    {
        static_if(low < high)
        {
            f(static_template<T, (reverse ? (high - low - 1) : (low))>{});
            static_for<T, low + 1, high, reverse>(std::forward<F>(f));
        }
    }
};


